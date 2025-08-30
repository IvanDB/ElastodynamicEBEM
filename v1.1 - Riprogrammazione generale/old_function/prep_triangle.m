function [zeta,tri2D,children,eta_children,c,sign_prod] = prep_triangle(sp,T3,vnF,Sistema_rif) 

%La function prep_triangle() prende in INPUT:
%- il vettore 1x3 sp contenente le coordinate del punto sorgente
%
%- la matrice 3x3 T3 in cui l'i-esima riga contiene le 3 componenti
%  dell'i-esimo vertice del triangolo di campo
%
%- il vettore vnF contenente le coordinate del versore normale al 
%  triangolo di campo
%
%- la struct Sistema_rif contenente le informazioni utili a determinare 
%  gli assi del nuovo sistema di riferimento
%
%e restituisce in OUTPUT:
%- la variabile zeta contenente la distanza tra il punto sorgente e il
%  piano individuato dal triangolo di campo
%
%- la matrice 3x2 tri2D in cui l'i-esima riga contiene le coordinate del 
%  vertice i-esimo del triangolo di campo nel nuovo sistema di 
%  riferimento
%
%- l'array 3D 3x3x2 children formato da 3 matrici 3x2 children(j,:,:) 
%  per j=1,2,3 che su ogni riga contengono le coordinate di uno dei tre 
%  vertici del j-esimo triangolo figlio nel sistema di riferimento 
%  bidimensionale.
%
%- il vettore 1x3 eta_children contenente le coordinate baricentriche 
%  della proiezione del punto sorgente rispetto ai tre vertici del
%  triangolo di campo nel piano individuato da tale triangolo
%
%- il vettore 1x3 c contenente i coefficienti associati a ciascuno dei 
%  tre triangoli figli che permetteranno di determinare se il triangolo
%  figlio è degenere oppure se ha i vertici ordinati in senso orario
%  o antiorario
%
%- la variabile sign_prod contenente il segno del prodotto scalare tra il 
%  versore normale al piano del triangolo di campo e il vettore distanza
%  tra il punto sorgente sp e la sua proiezione sp_plane nel piano del 
%  triangolo di campo


%% STEP 1: PROIETTIAMO IL PUNTO SORGENTE NEL PIANO INDIVIDUATO DAL 
%          TRIANGOLO DI CAMPO

%Determiniamo l'INTERSEZIONE tra il PIANO individuato dal TRIANGOLO di
%CAMPO e la RETTA passante per il PUNTO SORGENTE che ha VETTORE 
%DIRETTORE PARALLELO al VERSORE NORMALE al TRIANGOLO di CAMPO

abscissa = vnF*(sp-T3(1,:))';

%Vettore 1x3 sp_plane contenente le coordinate della PROIEZIONE del
%PUNTO SORGENTE sp nel PIANO del TRIANGOLO di CAMPO
sp_plane = sp-abscissa*vnF;
%sp_plane = abscissa*vn;

%Calcoliamo il SEGNO del PRODOTTO SCALARE tra il VERSORE NORMALE al piano
%del TRIANGOLO DI CAMPO e il VETTORE DISTANZA tra il PUNTO SORGENTE sp e 
%la sua PROIEZIONE sp_plane nel piano del triangolo di campo
sign_prod=(sp-sp_plane)*vnF';
if abs(sign_prod)<=1.0e-6
    sign_prod=0;
end
sign_prod=sign(sign_prod);

%DISTANZA del PUNTO SORGENTE dal PIANO individuato dal TRIANGOLO di CAMPO
zeta = abs(abscissa);
%zeta = -abscissa;
        
%% STEP 2: INDIVIDUIAMO GLI ASSI DEL NUOVO SISTEMA DI RIFERIMENTO

% %VETTORE P_2-P_1 (PRIMO LATO P1P2 del TRIANGOLO di CAMPO) 
% vL1 = T3(2,:)-T3(1,:);

%LUNGHEZZA del VETTORE P_2-P_1
L1 = Sistema_rif.L1;
%VERSORE e_1 PARALLELO al VETTORE P_2-P_1
e1 = Sistema_rif.e1;

% %VETTORE P_3-P_1 (TERZO LATO P1P3 del TRIANGOLO di CAMPO)
% vL3 = T3(3,:)-T3(1,:); 

%LUNGHEZZA del VETTORE P_3-P_1
L3 = Sistema_rif.L3;

% %VERSORE e_3 PARALLELO al VETTORE P_3-P_1
% e3 = Sistema_rif.e3;

%COSENO dell'ANGOLO compreso TRA e_1 e e_3
proj = Sistema_rif.cos_e1_e3;
%ANGOLO compreso TRA i VERSORI e_1 e e_3
delta = Sistema_rif.delta;

%VERSORE e_2 (che è la direzione perpendicolare al versore e_1 nel piano 
%del triangolo di campo) è dato dal PRODOTTO VETTORIALE tra il VERSORE vn 
%NORMALE al TRIANGOLO di CAMPO e il VERSORE e_1 
e2 = Sistema_rif.e2;

% %VETTORE P_2-P_1 (PRIMO LATO P1P2 del TRIANGOLO di CAMPO) 
% vL1 = T3(2,:)-T3(1,:);  
% %LUNGHEZZA del VETTORE P_2-P_1
% L1 = sqrt(sum(vL1.^2));
% %VERSORE e_1 PARALLELO al VETTORE P_2-P_1
% e1 = vL1/L1;
% 
% %VETTORE P_3-P_1 (TERZO LATO P1P3 del TRIANGOLO di CAMPO)
% vL3 = T3(3,:)-T3(1,:); 
% %LUNGHEZZA del VETTORE P_3-P_1
% L3 = sqrt(sum(vL3.^2));
% %VERSORE e_3 PARALLELO al VETTORE P_3-P_1
% e3 = vL3/L3;
% 
% %COSENO dell'ANGOLO compreso TRA e_1 e e_3
% proj = e3*e1';
% %ANGOLO compreso TRA i VERSORI e_1 e e_3
% delta = acos(proj);
% 
% %VERSORE e_2 (che è la direzione perpendicolare al versore e_1 nel piano 
% %del triangolo di campo) è dato dal PRODOTTO VETTORIALE tra il VERSORE vn 
% %NORMALE al TRIANGOLO di CAMPO e il VERSORE e_1 
% e2 = cross(vnF,e1);

%% STEP 3: DETERMINIAMO LE INFORMAZIONI UTILI SUL TRIANGOLO DI CAMPO NEL 
%          NUOVO SISTEMA DI RIFERIMENTO

%VETTORE DISTANZA tra il PUNTO P_1 e l'ORIGINE sp_plane del NUOVO SISTEMA 
%di RIFERIMENTO
diff = T3(1,:)-sp_plane;
%diff = triangle(1,:)-sp_plane;

%Inizializzazione della matrice 3x2 tri2D contenente le coordinate dei 
%vertici del triangolo di campo nel nuovo sistema di riferimento.
%In particolare, l'i-esima riga di questa matrice contiene le coordinate 
%del vertice i-esimo del triangolo di campo
tri2D = zeros(3,2);

%COORDINATE dei VERTICI del TRIANGOLO di CAMPO nel NUOVO SISTEMA 
%di RIFERIMENTO
tri2D(1,1) = diff*e1';
tri2D(1,2) = diff*e2';
% tri2D(1,1) = sum(diff*e1');
% tri2D(1,2) = sum(diff*e2');
tri2D(2,:) = tri2D(1,:)+[L1 0];
tri2D(3,:) = tri2D(1,:)+L3*[proj sin(delta)];

% plot3(sp_plane(1), sp_plane(2), sp_plane(3),'k*')
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% fill3([T3(1,1) T3(2,1) T3(3,1)],[T3(1,2) T3(2,2) T3(3,2)],[T3(1,3) T3(2,3) T3(3,3)],'g')
% hold on
% 
% plot3(T3(1,1),T3(1,2),T3(1,3),'bo','markerfacecolor','r','markerfacecolor','r','markersize',8)
% text(T3(1,1)+0.01,T3(1,2)+0.01,T3(1,3)+0.01,'$N_1$','interpreter','latex','fontsize',12)
% hold on
% plot3(T3(2,1),T3(2,2),T3(2,3),'bo','markerfacecolor','r','markerfacecolor','r','markersize',8)
% text(T3(2,1)+0.01,T3(2,2)+0.01,T3(2,3)+0.01,'$N_2$','interpreter','latex','fontsize',12)
% hold on
% plot3(T3(3,1),T3(3,2),T3(3,3),'bo','markerfacecolor','r','markerfacecolor','r','markersize',8)
% text(T3(3,1)+0.01,T3(3,2)+0.01,T3(3,3)+0.01,'$N_3$','interpreter','latex','fontsize',12)
% hold on

% plot3(sp_plane(1),sp_plane(2),sp_plane(3),'k*','markersize',10)
% hold on
% plot3([sp(1) sp_plane(1)],[sp(2) sp_plane(2)],[sp(3) sp_plane(3)],'k--')
% hold on
% plot3(sp(1),sp(2),sp(3),'bo','markerfacecolor','b','markersize',8)
% 
% xlabel('$x$','interpreter','latex','fontsize' ,14)
% ylabel('$y$','interpreter','latex','fontsize' ,14)
% zlabel('$z$','interpreter','latex','fontsize' ,14)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% fill([tri2D(1,1) tri2D(2,1) tri2D(3,1)],[tri2D(1,2) tri2D(2,2) tri2D(3,2)],'g')
% hold on 
% 
% plot(tri2D(1,1),tri2D(1,2),'ro','markerfacecolor','r','markersize',8)
% text(tri2D(1,1)+0.01,tri2D(1,2)+0.01,'$N_1$','interpreter','latex','fontsize',12)
% hold on
% plot(tri2D(2,1),tri2D(2,2),'ro','markerfacecolor','r','markersize',8)
% text(tri2D(2,1)+0.01,tri2D(2,2)+0.01,'$N_2$','interpreter','latex','fontsize',12)
% hold on
% plot(tri2D(3,1),tri2D(3,2),'ro','markerfacecolor','r','markersize',8)
% text(tri2D(3,1)+0.01,tri2D(3,2)+0.01,'$N_3$','interpreter','latex','fontsize',12)
% hold on
% plot([0 tri2D(1,1)],[0 tri2D(1,2)],'b--','LineWidth',1.5)
% hold on
% plot([0 tri2D(2,1)],[0 tri2D(2,2)],'b--','LineWidth',1.5)
% hold on
% plot([0 tri2D(3,1)],[0 tri2D(3,2)],'b--','LineWidth',1.5)
% hold on
% 
% plot(0,0,'ks','markersize',10,'markerfacecolor','k')
% hold on
% 
% xlabel('$x$','interpreter','latex','fontsize' ,14)
% ylabel('$y$','interpreter','latex','fontsize' ,14)
% 
%% STEP 4: INDIVIDUIAMO I TRE TRIANGOLI FIGLI DEL TRIANGOLO DI CAMPO

%Inizializzazione dell'array 3D children di formato 3x3x2 contenente le 
%coordinate dei vertici dei triangoli figli nel piano del triangolo di 
%campo
children = zeros(3,3,2);

%Individuo le COORDINATE dei VERTICI  di ciascuno dei TRE TRIANGOLI FIGLI 
%del triangolo di campo nel nuovo sistema di riferimento e memorizziamole
%nell'array 3D children di formato 3x3x2.

%Primo triangolo figlio: (0,0) -- 1° vertice -- 2° vertice
%children(1,1,:) = zeros(1,2);
children(1,2,:) = tri2D(1,:);
children(1,3,:) = tri2D(2,:);


%Secondo triangolo figlio: (0,0) -- 2° vertice -- 3° vertice
%children(2,1,:) = zeros(1,2);
children(2,2,:) = tri2D(2,:);  
children(2,3,:) = tri2D(3,:);

%Terzo triangolo figlio: (0,0) -- 3° vertice -- 1° vertice
%children(3,1,:) = zeros(1,2);
children(3,2,:) = tri2D(3,:);
children(3,3,:) = tri2D(1,:);

% children(1,2,:) = tri2D(2,:);
% children(1,3,:) = tri2D(3,:);
% children(2,2,:) = tri2D(3,:);
% children(2,3,:) = tri2D(1,:);
% children(3,2,:) = tri2D(1,:);
% children(3,3,:) = tri2D(2,:);

%Quindi children(j,:,:) per j=1,2,3 è una matrice 3x2 che su ogni riga 
%contiene le coordinate di uno dei tre vertici del j-esimo triangolo 
%figlio. 


%% STEP 5: CALCOLIAMO I COEFFICIENTI DA ASSEGNARE A CIASCUNO DEI TRE 
%          TRIANGOLI FIGLI 

%COEFFICIENTE relativo al PRIMO TRIANGOLO FIGLIO
c(1)=tri2D(1,1)*tri2D(2,2)-tri2D(1,2)*tri2D(2,1);
%c(1) = det([tri2D(1,:);tri2D(2,:)]);

%COEFFICIENTE relativo al SECONDO TRIANGOLO FIGLIO
c(2)=tri2D(2,1)*tri2D(3,2)-tri2D(2,2)*tri2D(3,1);
%c(2) = det([tri2D(2,:);tri2D(3,:)]);

%COEFFICIENTE relativo al TERZO TRIANGOLO FIGLIO
c(3)=tri2D(3,1)*tri2D(1,2)-tri2D(3,2)*tri2D(1,1);
%c(3) = det([tri2D(3,:);tri2D(1,:)]);

% %SECONDA COORDINATA BARICENTRICA del PUNTO sp_pllane
% eta_children(2) = det([tri2D(3,:); tri2D(1,:)])/(A_doppia);
% %TERZA COORDINATA BARICENTRICA del PUNTO sp_pllane
% eta_children(3) = det([tri2D(1,:); tri2D(2,:)])/(A_doppia);
% %PRIMA COORDINATA BARICENTRICA del PUNTO sp_pllane
% eta_children(1) = 1-eta_children(2)-eta_children(3)/(A_doppia);

%% STEP 6: CALCOLIAMO LE COORDINATE BARICENTRICHE DELLA PRORIEZIONE DEL 
%          PUNTO SORGENTE sp_plane

%Valore del DOPPIO dell'AREA CON SEGNO del TRIANGOLO di CAMPO CORRENTE
A_doppia=(tri2D(2,1)-tri2D(1,1))*(tri2D(3,2)-tri2D(1,2))-(tri2D(2,2)-tri2D(1,2))*(tri2D(3,1)-tri2D(1,1));
%A_doppia=det([(tri2D(2,:)-tri2D(1,:))' (tri2D(3,:)-tri2D(1,:))']);

%SECONDA COORDINATA BARICENTRICA del PUNTO PROIEZIONE sp_plane
eta_children(2) = c(3)/(A_doppia);

%TERZA COORDINATA BARICENTRICA del PUNTO PROIEZIONE sp_plane
eta_children(3) = c(1)/(A_doppia);

%PRIMA COORDINATA BARICENTRICA del PUNTO PROIEZIONE sp_plane
eta_children(1) = 1-eta_children(2)-eta_children(3);

% %PRIMA COORDINATA BARICENTRICA del PUNTO sp_pllane
% eta_children(1) = dot(cross((T3(3,:)-T3(2,:)),sp_plane-T3(2,:)),vnF)/(norm(vnF,2)^2);
% %SECONDA COORDINATA BARICENTRICA del PUNTO sp_pllane
% eta_children(2) = dot(cross((T3(1,:)-T3(3,:)),sp_plane-T3(3,:)),vnF)/(norm(vnF,2)^2);
% %TERZA COORDINATA BARICENTRICA del PUNTO sp_pllane
% eta_children(3) = dot(cross((T3(2,:)-T3(1,:)),sp_plane-T3(1,:)),vnF)/(norm(vnF,2)^2);
% 
% %Rapporto tra l'area del primo triangolo figlio e l'area del triangolo padre
% eta_children(1) = (tri2D(2,1)*tri2D(3,2)-tri2D(2,2)*tri2D(3,1))/(2*area);
% %Rapporto tra l'area del secondo triangolo figlio e l'area del triangolo padre
% eta_children(2) = (tri2D(3,1)*tri2D(1,2)-tri2D(3,2)*tri2D(1,1))/(2*area);
% %Rapporto tra l'area del terzo triangolo figlio e l'area del triangolo padre
% eta_children(3) = 1-eta_children(1)-eta_children(2);

return