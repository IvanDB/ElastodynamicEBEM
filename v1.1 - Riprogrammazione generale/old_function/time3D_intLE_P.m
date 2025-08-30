function [ris1,ris2,ris3,ris4] = time3D_intLE_P(zeta,R_P,pR_P,children,c,sign_prod)

%La function time3D_intLE_P() prende in INPUT:
%- la variabile zeta contenente la distanza tra il punto sorgente e il
%  piano individuato dal triangolo di campo
%
%- la variabile R_P=c_P*Delta contenente il raggio della superifice 
%  sferica individuata dal fronte d'onda relativo alle ONDE P
%
%- la variabile pR_P contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda P nel piano del triangolo di campo
%  pR_P = sqrt([c_P*Delta]^2-zeta^2);
%
%- l'array 3D 3x3x2 children formssato da 3 matrici 3x2 children(j,:,:) 
%  per j=1,2,3 che su ogni riga contengono le coordinate di uno dei tre 
%  vertici del j-esimo triangolo figlio nel sistema di riferimento 
%  bidimensionale
%
%- il vettore 1x3 c contenente i coefficienti associati a
%  ciascuno dei tre triangoli figli che permetteranno di determinare se il
%  triangolo figlio è degenere oppure se ha i vertici ordinati in senso
%  orario o antiorario
%
%- la variabile sign_prod contenente il segno del prodotto scalare tra il 
%  versore normale al piano del triangolo di campo e il vettore distanza
%  tra il punto sorgente sp e la sua proiezione sp_plane nel piano del 
%  triangolo di campo
%
%e restituisce in OUTPUT:
%- le matrici 3x3 ris1, ris2, ris3 e ris4 contenenti, rispettivamente,
%  il risultato dell'integrale di 1/r, r_i*r_j/r^3, 1/r^3 e r_i*r_j/r^5 
%  per i,j=1,2,3 nel caso in cui è presente solamente l'onda P
%
%-------------------------------------------------------------------------

%INIZIALIZZIAMO il risultato degli integrali
ris1 = zeros(3,3); %integrale di 1/r
ris2 = zeros(3,3); %integrale di r_i*r_j/r^3
ris3 = zeros(3,3); %integrale di 1/r^3
ris4 = zeros(3,3); %integrale di r_i*r_j/r^5 

%CICLO sui tre TRIANGOLI FIGLI del TRIANGOLO di CAMPO
for ind_child = 1:3
    
    %Estraiamo le COORDINATE dei VERTICI del TRIANGOLO FIGLIO corrente
    child = squeeze(children(ind_child,:,:)); 
    %La function squeeze() riduce l'array passato come parametro di una
    %dimensione

    %COEFFICIENTE assegnato al TRIANGOLO FIGLIO corrente
    C = c(ind_child);

%         %Rapporto tra l'area del triangolo figlio corrente e l'area del 
%         %triangolo padre
%         eta_child = eta_children(ind_child);
%         %Area del triangolo figlio corrente
%         area_child = area*eta_child;

    %Controlliamo che il TRIANGOLO FIGLIO sia NON DEGENERE
    if(abs(C)>1.0e-06)   

        %Se il coefficiente C assegnato al triangolo figlio è negativo, 
        %SCAMBIO il SECONDO e il TERZO VERTICE, in modo tale che i tre
        %vertici siano disposti in SENSO ANTIORARIO
        if(C < 0)
            temp = child(2,:);
            child(2,:) = child(3,:);
            child(3,:) = temp;
        end
%             if(area_child<0)
%                 temp = child(2,:);
%                 child(2,:) = child(3,:);
%                 child(3,:) = temp;
%             end

%-----------------------------------------------------------------------        

        %STEP 1: RUOTIAMO il TRIANGOLO FIGLIO CORRENTE in modo tale da
        %ALLINEARE il suo PRIMO LATO al PRIMO ASSE del nuovo sistema di
        %riferimento bidimensionale

        %Determiniamo la COORDINATA POLARE theta0 che permette di
        %individuare il SECONDO VERTICE child(2,:) del triangolo figlio 
        %corrente nel piano (angolo rispetto al semiasse positivo 
        %delle ascisse) 
        theta0 = atan2(child(2,2),child(2,1));
        %N.B. L'angolo thetat0 appartiene all'intervallo (-pi,pi]
        
        %Definiamo l'ANGOLO rispetto al quale dobbiamo ruotare in SENSO
        %ANTIORARIO il TRIANGOLO FIGLIO CORRENTE in modo tale da portare 
        %il suo PRIMO LATO child(1,:)--child(2,:) lungo il semiasse
        %positivo delle ascisse nel nuovo sistema di riferimento 
        %bidimensionale introdotto
        if (theta0<0)
            theta0=abs(theta0);
        else
            theta0=2*pi-theta0;
        end

        %Memorizziamo nella struct InfoTri2D le INFORMAZIONI GEOMETRICHE
        %relative al triangolo figlio corrente (NON DEGENERE!)
        InfoTri2D=time3D_InfoTri2D(child);

        %COORDINATA POLARE r relativa al SECONDO VERTICE del TRIANGOLO
        %FIGLIO CORRENTE
        %(LUNGHEZZA del PRIMO LATO a del triangolo figlio corrente)
        L1=InfoTri2D.a;
        %L1 = sqrt(sum(child(2,:).^2));

        %COORDINATA POLARE r relativa al TERZO VERTICE del TRIANGOLO
        %FIGLIO CORRENTE
        %(LUNGHEZZA del TERZO LATO b del triangolo figlio corrente)
        L3=InfoTri2D.b;
        %L3 = sqrt(sum(child(3,:).^2));

        %SENO e COSENO dell' ANGOLO gamma compreso tra il PRIMO e il 
        %TERZO LATO del TRIANGOLO FIGLIO corrente
        sin_gamma=InfoTri2D.sin_gamma;
        cos_gamma=InfoTri2D.cos_gamma;

%             %angolo tra il primo e il terzo lato del triangolo figlio corrente
%             %(sfrutto la formula del modulo del prodotto vettoriale tra 
%             %child(2,:) e child(3,:):
%             %2A=|child(2,:)xchild(3,:)|=|child(2,:)|*|child(3,:)|*sin(theta)=
%             %L1*L3*sin(theta)
%             pvec = (child(2,1)*child(3,2)-child(2,2)*child(3,1))/L1;
%             theta = asin(pvec/L3);
%           
%             %faccio in modo che theta sia un angolo positivo 
%             if(theta<0) 
%                 theta = theta+pi; 
%             end

%-----------------------------------------------------------------------        

        %STEP 2: Individuiamo le NUOVE COORDINATE dei VERTICI del 
        %TRIANGOLO FIGLIO CORRENTE dopo la rotazione e memorizziamole
        %all'interno della matrice 3x2 tri2D

        %Il primo vertice del triangolo ruotato è sempre l'origine
        tri2D(1,:) = zeros(1,2);

        %Il secondo vertice del triangolo ruotato è sull'asse delle 
        %ascisse
        tri2D(2,:) = [L1 0];

        %Coordinate del terzo vertice del triangolo ruotato
        tri2D(3,:) = [L3*cos_gamma L3*sin_gamma];

%             %ascissa del terzo vertice del triangolo rotato
%             tri2D(3,1) = L3*cos(theta);
%             %ordinata del terzo vertice del triangolo rotato = altezza
%             %h=(2*A)/b=|child(2,:)xchild(3,:)|/L1
%             tri2D(3,2) = pvec;        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure(4)
%         fill([tri2D(1,1) tri2D(2,1) tri2D(3,1)],[tri2D(1,2) tri2D(2,2) tri2D(3,2)],'c')
%         hold on 
%         plot(tri2D(1,1),tri2D(1,2),'ks','markersize',10,'markerfacecolor','k')
%         text(tri2D(1,1)+0.01,tri2D(1,2),'$O$','interpreter','latex','fontsize',12)
% 
%         hold on
%         
%         plot(tri2D(2,1),tri2D(2,2),'ro','markerfacecolor','r','markersize',8)
%         text(tri2D(2,1)+0.01,tri2D(2,2)+0.01,'$N_1$','interpreter','latex','fontsize',12)
%         hold on
%         plot(tri2D(3,1),tri2D(3,2),'ro','markerfacecolor','r','markersize',8)
%         text(tri2D(3,1)+0.01,tri2D(3,2)+0.01,'$N_2$','interpreter','latex','fontsize',12)
%         hold on
%         t =linspace(0,pi);
%         x1 = pR_P*cos(t);
%         y1= pR_P*sin(t);
%         plot(x1,y1,'r','linewidth',1.5)
% 
%         axis('equal')
%         xlabel('$x$','interpreter','latex','fontsize' ,14)
%         ylabel('$y$','interpreter','latex','fontsize' ,14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure(ind_child+3)
%         fill([tri2D(1,1) tri2D(2,1) tri2D(3,1)],[tri2D(1,2) tri2D(2,2) tri2D(3,2)],'c')
%         hold on 
%         plot(tri2D(1,1),tri2D(1,2),'ks','markersize',10,'markerfacecolor','k')
%         %text(tri2D(1,1)+0.01,tri2D(1,2),'$O$','interpreter','latex','fontsize',12)
% 
%         hold on
%         
%         %plot(tri2D(2,1),tri2D(2,2),'ro','markerfacecolor','r','markersize',8)
%         %text(tri2D(2,1)+0.01,tri2D(2,2)+0.01,'$N_1$','interpreter','latex','fontsize',12)
%         %hold on
%         %plot(tri2D(3,1),tri2D(3,2),'ro','markerfacecolor','r','markersize',8)
%         %text(tri2D(3,1)+0.01,tri2D(3,2)+0.01,'$N_2$','interpreter','latex','fontsize',12)
%         hold on
%         t =linspace(0,pi);
%         x1_P = pR_P*cos(t);
%         y1_P= pR_P*sin(t);
%         plot(x1_P,y1_P,'r','linewidth',1.5)
%          
%         axis('equal')
%         xlabel('$x$','interpreter','latex','fontsize' ,14)
%         ylabel('$y$','interpreter','latex','fontsize' ,14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-----------------------------------------------------------------------        

        %STEP 3: Calcoliamo l'INTEGRALE sul TRIANGOLO FIGLIO CORRENTE
        %        RUOTATO in presenza delle SOLE ONDE P
        [ris1t,ris2t,ris3t,ris4t] = ...
            integr_vertexLE_P(tri2D,InfoTri2D,zeta,R_P,pR_P,sign_prod);
%       [ris1t,ris2t,ris3t,ris4t] = integr_vertexLE_P(tri2D,theta0,theta,L1,L3,zeta,pR_P);

%----------------------------------------------------------------      

        %STEP 4: TRASFORMIAMO gli integrali sul triangolo figlio ruotato 
        %        in integrali sul triangolo figlio non ruotato    
        ris2t=rotazione_inversa_attorno_asse_z(ris2t,theta0);
        ris4t=rotazione_inversa_attorno_asse_z(ris4t,theta0);
        
        %- la variabile theta0 contenente l'ampiezza dell'angolo rispetto al
        %  quale dobbiamo ruotare in senso antiorario il triangolo figlio
        %  corrente modo tale da portare il suo primo lato lungo il 
        %  semiasse positivo delle ascisse nel nuovo sistema di riferimento 

        
        %----------------------------------------------------------------     

        %STEP 5: Sommiamo il CONTRIBUTO ALGEBRICO dell'INTEGRALE sul
        %        TRIANGOLO FIGLIO corrente


        %Calcoliamo il SEGNO del CONTRIBUTO del TRIANGOLO FIGLIO 
        %CORRENTE nella SOMMA ALGEBRICA che ci da l'integrale
        coeff = sign(C);
        %N.B. C è sicuramente un valore non prossimo a zero 
        %in quanto siamo all'interno dell'if che ci garantisce il  
        %fatto che abs(C)> 1.0e-16 

        %coeff=dsign(1.d0,area_child)  
        ris1 = ris1+ris1t*coeff; %integrale di 1/r
        ris2 = ris2+ris2t*coeff; %integrale di r_i*r_j/r^3
        ris3 = ris3+ris3t*coeff; %integrale di 1/r^3
        ris4 = ris4+ris4t*coeff; %integrale di r_i*r_j/r^5     

    end %Fine if(abs(C)> 1.0e-16)
     
end %Fine ciclo sui ind_child (indice triangoli figli)

return

% %FORSE da TOGLIERE in QUESTO MOMENTO
% %Coefficienti della prima funzione di base (relativa al primo vertice 
% %del triangolo di campo, vedi ref.1) 
% x2a(1,:) = [parent(2,1)*parent(3,2)-parent(3,1)*parent(2,2) ...
%     parent(2,2)-parent(3,2)...
%     parent(3,1)-parent(2,1)];
% %Coefficienti della seconda funzione di base (relativa al secondo vertice 
% %del triangolo di campo, vedi ref.1)
% x2a(2,:) = [parent(3,1)*parent(1,2)-parent(1,1)*parent(3,2) ...
%     parent(3,2)-parent(1,2) ...
%     parent(1,1)-parent(3,1)];
% %Divido per 2 volte l'area (vedi ref.1) 
% x2a = x2a/(2.d0*area);