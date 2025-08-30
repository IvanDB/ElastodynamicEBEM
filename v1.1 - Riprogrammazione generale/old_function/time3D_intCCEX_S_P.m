function [rist1,rist2,rist3,rist4] =time3D_intCCEX_S_P(R_S,pR_S,R_P,pR_P,zeta,Theta1,Theta2,Info_Tri2D,B)

%La function time3D_intCCEX_S_P() prende in INPUT:
%- la variabile R_S=c_S*Delta contenente il raggio della superificie 
%  sferica individuata dal fronte d'onda relativo alle ONDE S
%
%- la variabile pR_S contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda S nel piano del triangolo di campo
%  pR_S = sqrt([c_S*Delta]^2-zeta^2)
%
%- la variabile R_P=c_P*Delta contenente il raggio della superificie 
%  sferica individuata dal fronte d'onda relativo alle ONDE P
%
%- la variabile pR_P contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda P nel piano del triangolo di campo
%  pR_P = sqrt([c_S*Delta]^2-zeta^2)
%
%- la variabile zeta contenente la distanza del punto sorgente dalla 
%  sua proiezione nel piano in cui giace il triangolo di campo
%
%- la variabile Theta1 che contiene lampiezza del primo angolo che 
%  individua la corona circolare nel piano
%
%- la variabile Theta2 che contiene l'ampiezza del secondo angolo che 
%  individua la corona circolare nel piano
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
%- le matrici 3x3 rist1, rist2, rist3 e rist4 contenenti, rispettivamente, 
%  il valore degli integrali delle funzioni 1/r, r_ir_j/r^3, 1/r^3 e
%  r_ir_j/r^5 per i,j=1,2,3 calcolati sulla corona circolare in presenza 
%  dell'onda P
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

if (zeta>1.0e-6)
    %In questo caso la distanza zeta tra il punto sorgente e il piano del
    %triangolo di campo è maggiore di zero e quindi il triangolo sorgente 
    %e il triangolo di campo NON sono COMPLANARI
    
    %STEP 2: CALCOLO degli INTEGRALI in RHO nel caso zeta>0
    I(1) = R_P-R_S;
    I(2) = R_P-R_S+zeta^2*(1/R_P-1/R_S);
    I(3) = log((pR_P+R_P)/(pR_S+R_S))-pR_P/R_P+pR_S/R_S;
    I(4) = -1/R_P+1/R_S;
    I(5) = -1/R_P+1/R_S+(zeta^2/3)*(1/R_P^3-1/R_S^3);
    I(6) = (1/(3*zeta^2))*((pR_P/R_P)^3-(pR_S/R_S)^3);
    I(7) = -1/3*(1/R_P^3-1/R_S^3);
    %--------------------------------------------------------------------

    %STEP 3: CAlCOLO degli INTEGRALI in THETA
    Phi(1) = 2*sin(Theta2_tilde/2)*sin(Theta1_tilde/2);
    Phi(2) = 2*sin(gamma-Theta2_tilde/2)*sin(Theta1_tilde/2);
    Phi(3) = 1/2*(Theta1_tilde-cos(Theta2_tilde)*sin(Theta1_tilde));
    Phi(4) = -1/2*(-sin(Theta1_tilde)*cos(gamma-Theta2_tilde)+Theta1_tilde*cos_gamma);
    Phi(5) = 1/2*(Theta1_tilde-cos(2*gamma-Theta2_tilde)*sin(Theta1_tilde));
    %--------------------------------------------------------------------

    %STEP 4: CALCOLO delle MATRICI B_1, B_2 e B_3 nel caso zeta>0

    %Calcolo della matrice B_1
    B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+Phi(3)*B(:,:,6);
    %B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+2*Phi(3)*B(:,:,6);

    %Calcolo della matrice B_2
    B_2 = Phi(2)*B(:,:,2)+Phi(1)*B(:,:,3);

    %Calcolo della matrice B_3
    B_3 = B(:,:,1);  
    %--------------------------------------------------------------------

    %STEP 5: CALCOLO degli INTEGRALI nel caso zeta>0

    %INTEGRALE di 1/r sulla CORONA CIRCOLARE in presenza dell'ONDA P nel
    %caso zeta>0
    rist1 = Theta1_tilde*I(1)*eye(3,3);
%     f1=@(rho,theta) rho./(rho.^2+zeta^2).^(1/2);
%     RIS1=integral2(f1,pR_S,pR_P,Theta1,Theta2,'Method','iterated','AbsTol',1e-15,'RelTol',1e-15)*eye(3,3);
%     diff1=abs(RIS1-rist1)
    
    %INTEGRALE di r_i*r_j/r^3 sulla CORONA CIRCOLARE in presenza 
    %dell'ONDA P nel caso zeta>0
    rist2 = I(2)*B_1+I(3)*B_2+Theta1_tilde*I(4)*B_3;
%     f2=@(rho,theta) (rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS2=integral(f2,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_1;
%     f2=@(rho,theta) (rho.^2)./(rho.^2+zeta^2).^(3/2);
%     RIS2=RIS2+integral(f2,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_2;
%     f2=@(rho,theta) (rho)./(rho.^2+zeta^2).^(3/2);
%     RIS2=RIS2+integral(f2,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_3*Theta1_tilde;
%     diff2=abs(RIS2-rist2)
    
    %INTEGRALE di 1/r^3 sulla CORONA CIRCOLARE in presenza dell'ONDA P
    %nel caso zeta>0
    rist3 = Theta1_tilde*I(4)*eye(3,3);   
%     f3=@(rho,theta) rho./(rho.^2+zeta^2).^(3/2);
%     RIS3=integral2(f3,pR_S,pR_P,Theta1,Theta2,'Method','iterated','AbsTol',1e-15,'RelTol',1e-15)*eye(3,3);
%     diff3=abs(RIS3-rist3)

    %INTEGRALE di r_i*r_j/r^5 sulla CORONA CIRCOLARE in presenza 
    %dell'ONDA P nel caso zeta>0
    rist4 = I(5)*B_1+I(6)*B_2+Theta1_tilde*I(7)*B_3;     
%     f4=@(rho) (rho.^3)./(rho.^2+zeta^2).^(5/2);
%     RIS4=integral(f4,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_1;
%        f4=@(rho) (rho.^2)./(rho.^2+zeta^2).^(5/2);
%     RIS4=RIS4+integral(f4,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_2;
%        f4=@(rho) (rho)./(rho.^2+zeta^2).^(5/2);
%     RIS4=RIS4+integral(f4,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_3*Theta1_tilde;
%     diff4=abs(RIS4-rist4)
    %--------------------------------------------------------------------

else
    %In questo caso la distanza zeta tra il punto sorgente e il piano del
    %triangolo di campo è uguale a zero e quindi il triangolo sorgente 
    %e il triangolo di campo sono COMPLANARI
    
    %STEP 2: CALCOLO degli INTEGRALI in RHO nel caso zeta=0
    
    I(1) = R_P-R_S;
    I(2) = I(1);
    I(4) = -1/R_P+1/R_S;
    I(5) = I(4);
    %--------------------------------------------------------------------

   
    %STEP 3: CALCOLO degli INTEGRALI in THETA

    Phi(3) = 1/2*(Theta1_tilde-cos(Theta2_tilde)*sin(Theta1_tilde));
    Phi(4) = -1/2*(-sin(Theta1_tilde)*cos(gamma-Theta2_tilde)+Theta1_tilde*cos_gamma);
    Phi(5) = 1/2*(Theta1_tilde-cos(2*gamma-Theta2_tilde)*sin(Theta1_tilde));
    
    
%     Phi3=@(theta) (sin(theta)).^2;
%     RIS3=integral(Phi3,Theta1,Theta2,'AbsTol',1e-15,'RelTol',1e-15);
%     diff3=abs(RIS3-Phi(3))
%     
%     Phi4=@(theta) (sin(Info_Tri2D.gamma-theta)).*(sin(theta));
%     RIS4=integral(Phi4,Theta1,Theta2,'AbsTol',1e-15,'RelTol',1e-15);
%     diff4=abs(RIS4-Phi(4))
%     
%     Phi5=@(theta) (sin(Info_Tri2D.gamma-theta)).^2;
%     RIS5=integral(Phi5,Theta1,Theta2,'AbsTol',1e-15,'RelTol',1e-15);
%     diff5=abs(RIS5-Phi(5))

    %--------------------------------------------------------------------

    %STEP 4: CALCOLO della MATRICE B_1 nel caso zeta = 0

    %Calcolo della matrice B_1
    B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+Phi(3)*B(:,:,6);
    %--------------------------------------------------------------------

    % STEP 5: CALCOLO degli INTEGRALI 

    %INTEGRALE di 1/r sulla CORONA CIRCOLARE in presenza dell'ONDA P
    %nel caso zeta=0
    rist1 = Theta1_tilde*I(1)*eye(3,3);
    
%     f1=@(rho,theta) rho./(rho.^2+zeta^2).^(1/2);
%     RIS1=integral2(f1,pR_S,pR_P,Theta1,Theta2,'Method','iterated','AbsTol',1e-15,'RelTol',1e-15)*eye(3,3);
%     diff1=abs(RIS1-rist1)

    %INTEGRALE di r_i*r_j/r^3 sulla CORONA CIRCOLARE in presenza 
    %dell'ONDA P nel caso zeta=0
    rist2 = I(2)*B_1;
    
%     f2=@(rho,theta) (rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS2=integral(f2,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_1;
%     diff2=abs(RIS2-rist2)

    %INTEGRALE di 1/r^3 sulla CORONA CIRCOLARE in presenza dell'ONDA P
    %nel caso zeta=0
    rist3 = Theta1_tilde*I(4)*eye(3,3);
    
%     f3=@(rho,theta) rho./(rho.^2+zeta^2).^(3/2);
%     RIS3=integral2(f3,pR_S,pR_P,Theta1,Theta2,'Method','iterated','AbsTol',1e-15,'RelTol',1e-15)*eye(3,3);
%     diff3=abs(RIS3-rist3)
    
    %INTEGRALE di r_i*r_j/r^5 sulla CORONA CIRCOLARE in presenza 
    %dell'ONDA P nel caso zeta=0
    rist4 = I(5)*B_1;
    
%     f4=@(rho) (rho.^3)./(rho.^2+zeta^2).^(5/2);
%     RIS4=integral(f4,pR_S,pR_P,'AbsTol',1e-15,'RelTol',1e-15)*B_1;
%     diff4=abs(RIS4-rist4)
    
    %--------------------------------------------------------------------
   
end