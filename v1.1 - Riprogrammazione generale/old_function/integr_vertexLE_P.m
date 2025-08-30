function [ris1,ris2,ris3,ris4] = integr_vertexLE_P(tri2D,Info_Tri2D,zeta,R_P,pR_P,sign_prod)

%La function integr_vertexLE_P() prende in INPUT:
%- la matrice 3x2 tri2D in cui l'i-esima riga contiene le coordinate 
%  dell'i-esimo vertice del triangolo figlio corrente (già opportunamente 
%  ruotato) nel sistema di riferimento bidimensionale 
%
%- la struct Info_Tri2D che contiene le informazioni geoemtriche sul
%  triangolo figlio corrente (lunghezza dei lati, ampiezza degli angoli,
%  valore del seno e del coseno degli angoli)
%
%- la variabile zeta contenente la distanza del punto sorgente dalla 
%  sua proiezione nel piano in cui giace il triangolo di campo
%
%- la variabile R_P=c_P*Delta contenente il raggio della superifice 
%  sferica individuata dal fronte d'onda relativo alle ONDE P
%
%- la variabile pR_P contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda P nel piano del triangolo di campo
%  pR_P = sqrt([c_P*Delta]^2-zeta^2);
%
%- la variabile sign_prod contenente il segno del prodotto scalare tra il 
%  versore normale al piano del triangolo di campo e il vettore distanza
%  tra il punto sorgente sp e la sua proiezione sp_plane nel piano del 
%  triangolo di campo
%
%e restituisce in OUTPUT:
%- le matrici 3x3 ris1, ris2, ris3 e ris4 contenenti, rispettivamente, il
%  valore degli integrali delle funzioni 1/r, r_ir_j/r^3, 1/r^3 e
%  r_ir_j/r^5 per i,j=1,2,3 calcolati in presenza della sola ONDA P
%
%-------------------------------------------------------------------------
     
%% STEP 1: CALCOLIAMO L'INTERSEZIONE tra LA RETTA SU CUI GIACE IL SECONDO 
%          LATO DEL TRIANGOLO e la CIRCONFERENZA DI RAGGIO pR_P

% %Valore assoluto di z (distanza del punto sorgente dalla sua proiezione nel
% %piano in cui giace il triangolo di campo)
% abszeta = abs(zeta);

%VETTORE DISTANZA tra i vertici 2 e 3 del TRIANGOLO FIGLIO corrente nel
%sistema di riferimento bidimensionale
vt = tri2D(3,:)-tri2D(2,:);

%LUNGHEZZA del SECONDO LATO c (che ha come estremi i vertici 2 e 3) del 
%TRIANGOLO FIGLIO corrente 
L2=Info_Tri2D.c;
%L2 = sqrt(sum(vt.^2));

%DIREZIONE del SEGMENTO che ha come ESTREMI i VERTICI 2 e 3 del TRIANGOLO
%FIGLIO corrente
vt = vt/L2;

%LUNGHEZZA del PRIMO LATO a (che ha come estremi i vertici 1 e 2) del 
%TRIANGOLO FIGLIO corrente 
L1=Info_Tri2D.a;

%SISTEMA tra l'EQUAZIONE 
% xi_1=L1+vt(1)*s; xi_2=vt(2)*s, con s parametro 
%della RETTA PASSANTE per il VERTICE 2 del TRIANGOLO FIGLIO CORRENTE e 
%PARALLELA al VERSORE vt con l'EQUAZIONE della CIRCONFERENZA 
%xi_1^2+xi_2^2=R^2

%Così facendo si ricava l'equazione s^2+2*L1*vt(1)*s+(L1^2-R^2)=0
b_mezzi = L1*vt(1);
c = L1^2-pR_P^2;
Delta = b_mezzi^2-c;

%% STEP 2: CALCOLIAMO le MATRICI B_i

%Richiamiamo la funzione time3D_matB() che permette di calcolare:
%- l'array 3D matB di formato 3x3x6 che contiene le 6 matrici B_i di 
%  fomato 3x3 in cui sono memeorizzati i coefficienti che permettono di 
%  esprimere il prodotto r_i*r_k 
[matB] = time3D_matB(zeta,Info_Tri2D,sign_prod);
%[matB,gamma] = time3D_matB(zeta,p1,p2);

%% STEP 3: RICAVIAMO le INFORMAZIONI GEOMETRICHE RELATIVE al TRIANGOLO

%LATI del triangolo
a=Info_Tri2D.a;
%b=Info_tri2D.b;
%c=Info_tri2D.c;

%SENO degli ANGOLI del triangolo
sin_beta=Info_Tri2D.sin_beta;
%sin_alpha=Info_Tri2D.sin_alpha;
%sin_gamma=Info_tri2D.sin_gamma;

%AMPIEZZA deli ANGOLI del triangolo
%alpha=Info_Tri2D.alpha;
beta=Info_Tri2D.beta;
gamma=Info_Tri2D.gamma;

% %Inizializzazione del vettore 1x2 p1 contenente le COORDINATE del 
% %SECONDO VERTICE del TRIANGOLO FIGLIO corrente 
% p1 = tri2D(2,:);
% 
% %Inizializzazione del vettore 1x2 p2 contenente le COORDINATE del 
% %TERZO VERTICE del TRIANGOLO FIGLIO corrente 
% p2 = tri2D(3,:);


%% STEP 4: INTEGRAZIONE

if(Delta<=1.0e-06) 
%if(Delta<1.0e-06*pR_P)

    %Siamo nel caso in cui NON ci sono INTERSEZIONI (il fronte d'onda è 
    %ancora piccolo) OPPURE ci sono due INTERSEZIONI COINCIDENTI

    %ANGOLI che individuano il SETTORE CIRCOLARE
    Theta1 = 0;
    Theta2 = gamma;
    
    %Calcolo gli INTEGRALI ANALITICI sul SETTORE CIRCOLARE (Regione E6)
    %in presenza della sola ONDA P
    [ris1,ris2,ris3,ris4] = time3D_intSCEX_P(R_P,pR_P,zeta,Theta1,Theta2,Info_Tri2D,matB);

%     %Angoli che individuano il settore circolare
%     theta1 = theta0;
%     theta2 = theta+theta0;
%     %Calcolo gli integrali sul settore circolare considerato
%     [ris1,ris2,ris3,ris4] = time3D_intSCEX_P(R,abszeta,angc,theta1,theta2,matB);
     
else
    
    %Siamo nel caso in cui Delta>1.0e-06  e quindi ci sono DUE 
    %INTERSEZIONI DISTINTE
    
    %Inizializziamo il VALORE degli INTEGRALI GLOBALI
    ris1 = zeros(3,3); %Integrale di 1/r (Onda P)
    ris2 = zeros(3,3); %Integrale di r_i*r_j/r^3 (Onda P)
    ris3 = zeros(3,3); %Integrale di 1/r^3 (Onda P)
    ris4 = zeros(3,3); %Integrale di r_i*r_j/r^5 (Onda P)
    
    %Definiamo una COSTANTE RICORRENTE
    r=a*sin_beta/pR_P;
    
    %Definiamo gli ANGOLI
    Omega1=asin(r)-beta;
    Omega2=pi-asin(r)-beta;
    
    %Definiamo gli angoli
    theta1=min(Omega1,Omega2);
    theta2=max(Omega1,Omega2);
    
    THETA1=0;
    THETA2=min(max(0,theta1),gamma);
    %THETA3=max(0,min(theta1,gamma));
    THETA4=max(0,min(theta2,gamma));
    THETA5=gamma;
    
    
    %Calcoliamo gli INTEGRALI sul SETTORE CIRCOLARE (Regione E6)
    %in presenza della sola ONDA P
    if (THETA1<THETA2)
        %Controlliamo che il SETTORE CIRCOLARE sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo
        
        [ris1t,ris2t,ris3t,ris4t] = time3D_intSCEX_P(R_P,pR_P,zeta,THETA1,THETA2,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;
    end
    
    %Calcoliamo gli INTEGRALI sul TRIANGOLO (Regioni E7+E8)
    %in presenza della sola ONDA P
    if (THETA2<THETA4)
        %Controlliamo che il TRIANGOLO sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo

        [ris1t,ris2t,ris3t,ris4t] = time3D_intT3EX_P(zeta,THETA2,THETA4,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;
    end
    
    %Calcoliamo gli INTEGRALI sul SETTORE CIRCOLARE (Regione E9)
    %in presenza della sola ONDA P
    if (THETA4<THETA5)
        %Controlliamo che il SETTORE CIRCOLARE sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo
        
        [ris1t,ris2t,ris3t,ris4t] = time3D_intSCEX_P(R_P,pR_P,zeta,THETA4,THETA5,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;
    end
       
%     %Calcolo gli integrali sul settore circolare considerato
%     [ris1,ris2,ris3,ris4] = time3D_intSCEX_P(R,abszeta,angc,theta1,theta2,matB);
%     
%     %Calcolo la radice quadrata del discriminate
%     Delta = sqrt(Delta);
%     %Calcolo i valori del parametro s che individuano le due intersezioni
%     s1 = -b_mezzi-Delta;
%     s2 = -b_mezzi+Delta;
%      
%     if ((s1>L2) || (s2<0)) %devo integrare solo su un settore circolare 
%         
%         %Angoli che individuano il settore circolare
%         theta1 = theta0;
%         theta2 = theta+theta0;
% % % % % % %         %Calcolo gli integrali sul settore circolare considerato
% % % % % % %         %[ris1,ris2,ris3,ris4] = time3D_intSCEX_P(R,abszeta,angc,theta1,theta2,matB);
%         
%     else
%         
%         %CALCOLO degli INTEGRALI su SETTORI CIRCOLARI
%         
%         if(s1>0) %integrale sul primo settore circolare
%             
%             %sposto il punto p1
%             p1 = tri2D(2,:)+s1*vt;
%             %calcolo gli angoli che descrivono il settore circolare
%             theta1 = theta0;
%             theta2 =  asin(p1(2)/pR_P);
%             %faccio in modo che theta2 sia un angolo positivo
%             if(theta2<0) theta2 = theta2+pi; end
%             theta2 = theta2+theta0;
% % % % % % %             %calcolo l'integrale sul settore circolare 
% % % % % % %             %[ris1,ris2,ris3,ris4] = time3D_intSCEX_P(R,abszeta,angc,theta1,theta2,matB);
%             
%         end
%         
%         if(s2<L2) %integrale sul secondo settore circolare
%             
%             %sposto il punto p2
%             p2 = tri2D(2,:)+s2*vt;
%             %calcolo il primo angolo che descrive il settore circolare
%             theta1 = asin(p2(2)/pR_P);
%             %faccio in modo che theta2 sia un angolo positivo
%             if(theta1<0) theta1 = theta1+pi; end
%             theta1 = theta1+theta0;
%             %calcolo il secondo angolo che descrive il settore circolare
%             theta2 = theta+theta0;
% % % % % % %             %calcolo degli integrali sul settore circolare 
% % % % % % %             %[ris1t,ris2t,ris3t,ris4t] = time3D_intSCEX_P(R,abszeta,angc,theta1,theta2,matB);
%             %aggiorno il valore degli integrali globali
%             ris1 = ris1+ris1t;
%             ris2 = ris2+ris2t;
%             ris3 = ris3+ris3t;
%             ris4 = ris4+ris4t;
%         end
%         
%         %CALCOLO degli INTEGRALI sul TRIANGOLO
%         
% % % % % % %         %calcolo degli integrali
% % % % % % %         %[ris1t,ris2t,ris3t,ris4t] = time3D_intT3EX_P(abszeta,p1,p2,matB);
%         
%         %aggiorno il valore degli integrali globali
%         ris1 = ris1+ris1t;
%         ris2 = ris2+ris2t;
%         ris3 = ris3+ris3t; 
%         ris4 = ris4+ris4t; 

%    end  
    
end %Fine if(Delta <=1.0e-6)

end