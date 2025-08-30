function [ris1,ris2,ris3,ris4,ris5,ris6] = integr_vertexLE_S_P_1(tri2D,Info_Tri2D,zeta,R_P,pR_P,R_S,pR_S,sign_prod)

%La function integr_vertexLE_S_P_1() prende in INPUT:
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
%  sferica individuata dal fronte d'onda relativo alle onde P
%
%- la variabile pR_P contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda P nel piano del triangolo di campo
%  pR_P = sqrt([c_P*Delta]^2-zeta^2)
%
%- la variabile R_S=c_S*Delta contenente il raggio della superifice 
%  sferica individuata dal fronte d'onda relativo alle onde S
%
%- la variabile pR_S contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda S nel piano del triangolo di campo
%  pR_S = sqrt([c_S*Delta]^2-zeta^2)
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
%- le matrici 3x3 ris5 e ris6 contenenti, rispettivamente, il
%  valore degli integrali delle funzioni 1/r e r_ir_j/r^3 per i,j=1,2,3 
%  calcolati in presenza delle onde S e P
%
%-------------------------------------------------------------------------

%% STEP 1: CALCOLIAMO L'INTERSEZIONE tra LA RETTA SU CUI GIACE IL SECONDO 
%          LATO DEL TRIANGOLO e le CIRCONFERENZE DI RAGGIO pR_P e pR_S

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
%xi_1^2+xi_2^2=pR_P^2

%Così facendo si ricava l'equazione s^2+2*L1*vt(1)*s+(L1^2-pR_P^2)=0
b_mezzi = L1*vt(1);
c = L1^2-pR_P^2;
Delta_P = b_mezzi^2-c;

%SISTEMA tra l'EQUAZIONE 
% xi_1=L1+vt(1)*s; xi_2=vt(2)*s, con s parametro 
%della RETTA PASSANTE per il VERTICE 2 del TRIANGOLO FIGLIO CORRENTE e 
%PARALLELA al VERSORE vt con l'EQUAZIONE della CIRCONFERENZA 
%xi_1^2+xi_2^2=pR_S^2

%Così facendo si ricava l'equazione s^2+2*L1*vt(1)*s+(L1^2-pR_S^2)=0
c = L1^2-pR_S^2;
Delta_S = b_mezzi^2-c;

%% STEP 2: CALCOLIAMO le MATRICI B_i

%Richiamiamo la funzione time3D_matB() che permette di calcolare:
%- l'array 3D matB di formato 3x3x6 che contiene le 6 matrici B_i di 
%  formato 3x3 in cui sono memeorizzati i coefficienti che permettono di 
%  esprimere il prodotto r_i*r_k 

[matB] = time3D_matB(zeta,Info_Tri2D,sign_prod);
% if (zeta>1.0e-6)
%     [matB] = time3D_matB(zeta,Info_Tri2D);
% else
%     [matB] = time3D_matB(0,Info_Tri2D);
% end

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

%% STEP 4: INTEGRAZIONE

%Inizializziamo il VALORE degli INTEGRALI GLOBALI
ris1 = zeros(3,3); %Integrale di 1/r (Onda P)
ris2 = zeros(3,3); %Integrale di r_i*r_j/r^3 (Onda P)
ris3 = zeros(3,3); %Integrale di 1/r^3 (Onda P)
ris4 = zeros(3,3); %Integrale di r_i*r_j/r^5 (Onda P)
ris5 = zeros(3,3); %Integrale di 1/r (Onde S e P)
ris6 = zeros(3,3); %Integrale di r_i*r_j/r^3 (Onde S e P)

if(Delta_P<=1.0e-06) 
    
    %Siamo nel caso in cui NON ci sono INTERSEZIONI tra la CIRCONFERENZA 
    %di RAGGIO pR_P e la retta passante per il LATO c OPPURE siamo nel 
    %caso in cui ci sono due INTERSEZIONI COINCIDENTI

    %ANGOLI che individuano il SETTORE e la CORONA CIRCOLARE
    theta1 = 0;
    theta2 = gamma;
    
%     hold on
%     plot(pR_P*cos(theta1),pR_P*sin(theta1),'b*');
%     hold on
%     plot(pR_P*cos(theta2),pR_P*sin(theta2),'b*');
%     
%     hold on
%     plot(pR_S*cos(theta1),pR_S*sin(theta1),'r*');
%     hold on
%     plot(pR_S*cos(theta2),pR_S*sin(theta2),'r*');

    
    %Calcolo gli INTEGRALI ANALITICI sul SETTORE CIRCOLARE (Regione E5) 
    %in presenza delle ONDE S e P
    [ris5,ris6] = time3D_intSCEX_S_P(R_S,pR_S,zeta,theta1,theta2,Info_Tri2D,matB);
    
    %Calcolo gli INTEGRALI ANALITICI sulla CORONA CIRCOLARE (Regione E6) 
    %in presenza dell'ONDA P
    [ris1,ris2,ris3,ris4] = time3D_intCCEX_S_P(R_S,pR_S,R_P,pR_P,zeta,theta1,theta2,Info_Tri2D,matB);

elseif (Delta_P>1.0e-06 && Delta_S<=1.0e-06) 
  
    %Siamo nel caso in cui la CIRCONFERENZA di RAGGIO pR_P INTERSECA la 
    %retta passante per il LATO c in due PUNTI DISTINTI, mentre la 
    %CIRCONFERENZA di RAGGIO pR_S NON ha INTERSEZIONI con la retta
    %passante per il LATO c oppure ha DUE INTERSEZIONI COINCIDENTI
    
    %Definiamo una COSTANTE RICORRENTE
    rP=a*sin_beta/pR_P;
    
    %Definiamo gli ANGOLI
    Omega1=asin(rP)-beta;
    Omega2=pi-asin(rP)-beta;
    
    %Definiamo gli angoli
    theta1=min(Omega1,Omega2);
    theta2=max(Omega1,Omega2);
    
    THETA1=0;
    THETA2=min(max(0,theta1),gamma);
    %THETA3=max(0,min(theta1,gamma));
    THETA4=max(0,min(theta2,gamma));
    THETA5=gamma;
    
    %Calcolo gli INTEGRALI ANALITICI sul SETTORE CIRCOLARE 
    %(Regioni E1+E4+E5) in presenza delle ONDE S e P
    [rist5,rist6] = time3D_intSCEX_S_P(R_S,pR_S,zeta,THETA1,THETA5,Info_Tri2D,matB);
    %Aggiorniamo il valore degli integrali globali
    ris5 = ris5+rist5;
    ris6 = ris6+rist6;    
    
    %Calcoliamo gli INTEGRALI sulla CORONA CIRCOLARE (Regione E6) 
    %in presenza dell'ONDA P
    if (THETA1<THETA2)
         %Controlliamo che la CORONA CIRCOLARE sia NON DEGENERE.
         %In caso contrario, il risultato degli integrali sarà nullo
    
        [rist1,rist2,rist3,rist4] = time3D_intCCEX_S_P(R_S,pR_S,R_P,pR_P,zeta,THETA1,THETA2,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+rist1;
        ris2 = ris2+rist2;
        ris3 = ris3+rist3;
        ris4 = ris4+rist4;
    end
    
    %Calcoliamo gli INTEGRALI sul TRAPEZIO/TRIANGOLO CIRCOLARE 
    %(Regioni E7+E8) in presenza dell'ONDA P
    if (THETA2<THETA4)
        %Controlliamo che il TRAPEZIO/TRIANGOLO CIRCOLARE sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo
     
        [ris1t,ris2t,ris3t,ris4t] = time3D_intTCEX_S_P(zeta,R_S,pR_S,THETA2,THETA4,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;
    end
    
    %Calcoliamo gli INTEGRALI sulla CORONA CIRCOLARE (Regione E9) 
    %in presenza dell'ONDA P
    
    if (THETA4<THETA5)
         %Controlliamo che la CORONA CIRCOLARE sia NON DEGENERE.
         %In caso contrario, il risultato degli integrali sarà nullo
 
        [ris1t,ris2t,ris3t,ris4t] = time3D_intCCEX_S_P(R_S,pR_S,R_P,pR_P,zeta,THETA4,THETA5,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;
    end
    
else
    
    %Siamo nel caso in cui le CIRCONFERENZE di RAGGIO pR_P e pR_S 
    %INTERSECANO la retta passante per il LATO c in due PUNTI DISTINTI
    
    %Definiamo una COSTANTE RICORRENTE
    rP=a*sin_beta/pR_P;
    rS=a*sin_beta/pR_S;
    
    %Definiamo gli ANGOLI
    Omega1=asin(rP)-beta;
    Omega2=pi-asin(rP)-beta;
    Omega3=asin(rS)-beta;
    Omega4=pi-asin(rS)-beta;
    
    %Definiamo gli angoli
    theta1=min(Omega1,Omega2);
    theta2=max(Omega1,Omega2);
    theta3=min(Omega3,Omega4);
    theta4=max(Omega3,Omega4);
    
    THETA1=0;
    THETA2=min(max(0,theta1),gamma);
    THETA3=min(max(0,theta3),gamma);
    THETA4=max(0,min(theta4,gamma));
    THETA5=max(0,min(theta2,gamma));
    THETA6=gamma;
    
%     hold on
%     plot(pR_S*cos(THETA1),pR_S*sin(THETA1),'b*');
%     hold on
%     plot(pR_P*cos(THETA1),pR_P*sin(THETA1),'b*');
%     hold on
%     plot(pR_S*cos(THETA2),pR_S*sin(THETA2),'b*');
%     hold on
%     plot(pR_S*cos(THETA3),pR_S*sin(THETA3),'b*');
%     hold on
%     plot(pR_S*cos(THETA4),pR_S*sin(THETA4),'b*');
%     hold on
%     plot([0 pR_S*cos(THETA4)],[0 pR_S*sin(THETA4)],'g-')
%     hold on
%     plot(pR_S*cos(THETA5),pR_S*sin(THETA5),'b*');
%     hold on
%     plot(pR_S*cos(THETA6),pR_S*sin(THETA6),'b*');

      
    %Calcolo gli INTEGRALI ANALITICI 
    %- sul SETTORE CIRCOLARE (Regione E1) in presenza delle ONDE S e P
    %- sulla CORONA CIRCOLARE (Regione E6) in presenza dell'ONDA P
    if (THETA1<THETA2)
        %Controlliamo che il SETTORE CIRCOLARE sia NON DEGENERE.
        %Controlliamo che la CORONA CIRCOLARE sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo

        [rist5,rist6] = time3D_intSCEX_S_P(R_S,pR_S,zeta,THETA1,THETA2,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris5 = ris5+rist5;
        ris6 = ris6+rist6;
        
        %--------------------------------------------------------------      
        
        [rist1,rist2,rist3,rist4] = time3D_intCCEX_S_P(R_S,pR_S,R_P,pR_P,zeta,THETA1,THETA2,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+rist1;
        ris2 = ris2+rist2;
        ris3 = ris3+rist3;
        ris4 = ris4+rist4;
    end
    
    %Calcolo gli INTEGRALI ANALITICI 
    %- sul SETTORE CIRCOLARE (Regione E2) in presenza delle ONDE S e P
    %- su un TRAPEZIO/TRIANGOLO CIRCOLARE (Regione E7) in presenza 
    %  dell'ONDA P
    if (THETA2<THETA3)
        %Controlliamo che il SETTORE CIRCOLARE sia NON DEGENERE.
        %Controlliamo che il TRAPEZIO/TRIANGOLO CIRCOLARE sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo

        [rist5,rist6] = time3D_intSCEX_S_P(R_S,pR_S,zeta,THETA2,THETA3,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris5 = ris5+rist5;
        ris6 = ris6+rist6;
        
        %--------------------------------------------------------------      
        
        [ris1t,ris2t,ris3t,ris4t] = time3D_intTCEX_S_P(zeta,R_S,pR_S,THETA2,THETA3,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;
        
    end
       
   
    %Calcolo gli INTEGRALI ANALITICI sul TRIANGOLO (Regione E3)
    %in presenza delle ONDE S e P
    if (THETA3<THETA4)
        %Controlliamo che il TRIANGOLO sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo
        
        [rist5,rist6] = time3D_intT3EX_S_P(zeta,THETA3,THETA4,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris5 = ris5+rist5;
        ris6 = ris6+rist6;
    end
    

    %Calcolo gli INTEGRALI ANALITICI 
    %- sul SETTORE CIRCOLARE (Regione E4) in presenza delle ONDE S e P
    %- sul TRAPEZIO/TRIANGOLO CIRCOLARE (Regione E8) in presenza 
    %  dell'ONDA P
    if (THETA4<THETA5)
        %Controlliamo che il SETTORE CIRCOLARE sia NON DEGENERE.
        %Controlliamo che il TRAPEZIO/TRIANGOLO CIRCOLARE sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo

        [rist5,rist6] = time3D_intSCEX_S_P(R_S,pR_S,zeta,THETA4,THETA5,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris5 = ris5+rist5;
        ris6 = ris6+rist6;
       
        %--------------------------------------------------------------            
        
        [ris1t,ris2t,ris3t,ris4t] = time3D_intTCEX_S_P(zeta,R_S,pR_S,THETA4,THETA5,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;
    end
    
     
    %Calcolo gli INTEGRALI ANALITICI 
    %- sul SETTORE CIRCOLARE (Regione E5) in presenza delle ONDE S e P
    %- sulla CORONA CIRCOLARE (Regione E9) in presenza dell'ONDA P
    if (THETA5<THETA6)
        %Controlliamo che il SETTORE CIRCOLARE sia NON DEGENERE.
        %Controlliamo che la CORONA CIRCOLARE sia NON DEGENERE.
        %In caso contrario, il risultato degli integrali sarà nullo
        
        [rist5,rist6] = time3D_intSCEX_S_P(R_S,pR_S,zeta,THETA5,THETA6,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris5 = ris5+rist5;
        ris6 = ris6+rist6;
        
        %--------------------------------------------------------------            
         
        [ris1t,ris2t,ris3t,ris4t] = time3D_intCCEX_S_P(R_S,pR_S,R_P,pR_P,zeta,THETA5,THETA6,Info_Tri2D,matB);
        %Aggiorniamo il valore degli integrali globali
        ris1 = ris1+ris1t;
        ris2 = ris2+ris2t;
        ris3 = ris3+ris3t;
        ris4 = ris4+ris4t;      
    end
       
end %Fine if(Delta <=1.0e-6)

end