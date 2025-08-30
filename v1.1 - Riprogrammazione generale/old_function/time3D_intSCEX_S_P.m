function [rist5,rist6] = time3D_intSCEX_S_P(R_S,pR_S,zeta,Theta1,Theta2,Info_Tri2D,B)

%La function time3D_intSCEX_S_P() prende in INPUT:
%- la variabile R_S=c_S*Delta contenente il raggio della superificie 
%  sferica individuata dal fronte d'onda relativo alle onde S
%
%- la variabile pR_S contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda S nel piano del triangolo di campo
%  pR_S = sqrt([c_S*Delta]^2-zeta^2);
%
%- la variabile zeta contenente la distanza del punto sorgente dalla 
%  sua proiezione nel piano in cui giace il triangolo di campo
%
%- la variabile Theta1 che contiene lampiezza del primo angolo che 
%  individua il settore circolare nel piano
%
%- la variabile Theta2 che contiene l'ampiezza del secondo angolo che 
%  individua il settore circolare nel piano
%
%- la struct Info_Tri2D che contiene le informazioni geometriche sul
%  triangolo figlio corrente (lunghezza dei lati, ampiezza degli angoli,
%  valore del seno e del coseno degli angoli)
%
%- l'array 3D B di formato 3x3x6, costituito da 6 matrici 3x3 ciascuna 
%  delle quali contiene una delle 6 tipologie di coefficienti che
%  permettono di esprimere il prodotto r_ir_k al variare di i,k=1,2,3
%
%e restituisce in OUTPUT:
%- le matrici 3x3 rist5 e rist6 contenenti, rispettivamente, 
%  il valore degli integrali delle funzioni 1/r e r_ir_j/r^3 per i,j=1,2,3 
%  calcolati sul SETTORE CIRCOLARE in presenza delle onde S e P 
%
%-------------------------------------------------------------------------
     
%STEP1: CALCOLO delle COSTANTI RICORRENTI NEGLI INTEGRALI e delle 
%       INFORMAZIONI GEOMETRICHE NECESSARIE

%Calcolo delle COSTANTI RICORRENTI negli INTEGRALI in THETA
Theta1_tilde = Theta2-Theta1;
Theta2_tilde = Theta2+Theta1;

%COSENO dell'ANGOLO GAMMA
cos_gamma=Info_Tri2D.cos_gamma;

%AMPIEZZA dell'ANGOLO GAMMA del triangolo
gamma=Info_Tri2D.gamma;

%-------------------------------------------------------------------------

% I(3) = 

if (zeta>1.0e-6) 
    %In questo caso la distanza zeta tra il punto sorgente e il piano del
    %triangolo di campo è maggiore di zero e quindi il triangolo sorgente 
    %e il triangolo di campo NON sono COMPLANARI
    
    %STEP 2: CALCOLO degli INTEGRALI in RHO nel caso zeta>0
    I(1) = R_S-zeta;
    I(2) = R_S+zeta^2/R_S-2*zeta;
    I(3) = log((pR_S+R_S)/zeta)-pR_S/R_S;
    I(4) = -1/R_S+1/zeta;  
    %--------------------------------------------------------------------

    %STEP 3: CAlCOLO degli INTEGRALI in THETA
    Phi(1) = 2*sin(Theta2_tilde/2)*sin(Theta1_tilde/2);
    Phi(2) = 2*sin(gamma-Theta2_tilde/2)*sin(Theta1_tilde/2);
    Phi(3) = 1/2*(Theta1_tilde-cos(Theta2_tilde)*sin(Theta1_tilde));
    Phi(4) = -1/2*(-sin(Theta1_tilde)*cos(gamma-Theta2_tilde)+Theta1_tilde*cos_gamma);
    Phi(5) = 1/2*(Theta1_tilde-cos(2*gamma-Theta2_tilde)*sin(Theta1_tilde));
    %--------------------------------------------------------------------

    %STEP 4: CALCOLO delle MATRICI B_1, B_2 e B_3 nel caso zeta > 0

    %Calcolo della matrice B_1
    B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+Phi(3)*B(:,:,6);
    %B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+2*Phi(3)*B(:,:,6);

    %Calcolo della matrice B_2
    B_2 = Phi(2)*B(:,:,2)+Phi(1)*B(:,:,3);

    %Calcolo della matrice B_3
    B_3 = B(:,:,1);  
    %--------------------------------------------------------------------

    %STEP 5: CALCOLO degli INTEGRALI sul SETTORE CIRCOLARE 
    %        nel caso zeta>0

    %INTEGRALE di 1/r sul SETTORE CIRCOLARE in presenza delle ONDE S e P
    %nel caso zeta maggiore di zero
    rist5 = Theta1_tilde*I(1)*eye(3,3);
%     f5=@(rho,theta) (rho)./(rho.^2+zeta^2).^(1/2);
%     RIS5=integral2(f5,0,pR_S,Theta1,Theta2,'Method','iterated','AbsTol',1e-15,'RelTol',1e-15)*eye(3,3);
%     diff5=abs(RIS5-rist5)
    
    %INTEGRALE di r_i*r_j/r^3 sul SETTORE CIRCOLARE in presenza delle 
    %ONDE S e P nel caso zeta maggiore di zero
    rist6 = I(2)*B_1+I(3)*B_2+Theta1_tilde*I(4)*B_3;      
%     f6=@(rho) (rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS6=integral(f6,0,pR_S,'AbsTol',1e-15,'RelTol',1e-15)*B_1;
%     f6=@(rho) (rho.^2)./(rho.^2+zeta^2).^(3/2);
%     RIS6=RIS6+integral(f6,0,pR_S,'AbsTol',1e-15,'RelTol',1e-15)*B_2;
%     f6=@(rho) (rho)./(rho.^2+zeta^2).^(3/2);
%     RIS6=RIS6+integral(f6,0,pR_S,'AbsTol',1e-15,'RelTol',1e-15)*B_3*Theta1_tilde;
%     diff6=abs(RIS6-rist6)
    %--------------------------------------------------------------------
    
else
    %In questo caso la distanza zeta tra il punto sorgente e il piano del
    %triangolo di campo è uguale a zero e quindi il triangolo sorgente 
    %e il triangolo di campo sono COMPLANARI
    
    %STEP 2: CALCOLO degli INTEGRALI in RHO nel caso zeta=0
    
    I(1) = R_S;
    I(2) = R_S;
    %--------------------------------------------------------------------
    
    %STEP 3: CAlCOLO degli INTEGRALI in THETA
    
    Phi(3) = 1/2*(Theta1_tilde-cos(Theta2_tilde)*sin(Theta1_tilde));
    Phi(4) = -1/2*(Theta1_tilde*cos_gamma-sin(Theta1_tilde)*cos(gamma-Theta2_tilde));
    Phi(5) = 1/2*(Theta1_tilde-cos(2*gamma-Theta2_tilde)*sin(Theta1_tilde));
    
    %--------------------------------------------------------------------

    %STEP 4: CALCOLO della MATRICE B_1 nel caso zeta=0

    %Calcolo della matrice B_1
    B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+Phi(3)*B(:,:,6);
    %--------------------------------------------------------------------

    %STEP 5: CALCOLO degli INTEGRALI sul SETTORE CIRCOLARE 
    %        nel caso zeta=0

    %INTEGRALE di 1/r sul SETTORE CIRCOLARE in presenza delle onde S e P 
    %nel caso in cui zeta=0
    rist5 = Theta1_tilde*I(1)*eye(3,3);
     
%     f5=@(rho,theta) (rho)./(rho.^2+zeta^2).^(1/2);
%     RIS5=integral2(f5,0,pR_S,Theta1,Theta2,'Method','iterated','AbsTol',1e-15,'RelTol',1e-15)*eye(3,3);
%     diff5=abs(RIS5-rist5)

    %INTEGRALE di r_i*r_j/r^3 sul SETTORE CIRCOLARE in presenza delle 
    %onde S e P nel caso in cui zeta=0
    rist6 = I(2)*B_1;
    
%     f6=@(rho) (rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS6=integral(f6,0,pR_S,'AbsTol',1e-15,'RelTol',1e-15)*B_1;
%     diff6=abs(RIS6-rist6)
    
    %--------------------------------------------------------------------

end %Fine if(zeta>1.0e-6)
   
end