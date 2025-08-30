function [rist1,rist2,rist3,rist4] = time3D_intSCEX_P(R_P,pR_P,zeta,Theta1,Theta2,Info_Tri2D,B)

%La function time3D_intSCEX_P() prende in INPUT:
%- la variabile R_P=c_P*Delta contenente il raggio della superifice 
%  sferica individuata dal fronte d'onda relativo alle ONDE P
%
%- la variabile pR_P contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda P nel piano del triangolo di campo
%  pR_P = sqrt([c_P*Delta]^2-zeta^2)
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
%- le matrici 3x3 rist1, rist2, rist3 e rist4 contenenti, rispettivamente, 
%  il valore degli integrali delle funzioni 1/r, r_ir_j/r^3, 1/r^3 e
%  r_ir_j/r^5 per i,j=1,2,3 calcolati sul SETTORE CIRCOLARE in presenza 
%  della sola ONDA P
%
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%NOTA BENE: in presenza di sole onde P il punto sorgente non si trova mai
%sullo stesso piano del triangolo di campo e quindi zeta > 0
%-------------------------------------------------------------------------
   
%% STEP1: CALCOLO delle COSTANTI RICORRENTI NEGLI INTEGRALI e delle 
%  INFORMAZIONI GEOMETRICHE NECESSARIE

%Calcolo delle COSTANTI RICORRENTI negli INTEGRALI in THETA
Theta1_tilde = Theta2-Theta1;
Theta2_tilde = Theta2+Theta1;

%COSENO dell'ANGOLO GAMMA
cos_gamma=Info_Tri2D.cos_gamma;

%AMPIEZZA dell'ANGOLO GAMMA del triangolo
gamma=Info_Tri2D.gamma;

%% STEP 2: CALCOLO degli INTEGRALI in RHO

I(1) = R_P-zeta;
I(2) = R_P+zeta^2/R_P-2*zeta;
I(3) = log((pR_P+R_P)/zeta)-pR_P/R_P;
I(4) = -1/R_P+1/zeta;
I(5) = -1/R_P+zeta^2/(3*R_P^3)+2/(3*zeta);
I(6) = 1/(3*zeta^2)*(pR_P/R_P)^3;
I(7) = -1/(3*R_P^3)+1/(3*zeta^3);

%% STEP 3: CALCOLO degli INTEGRALI in THETA

Phi(1) = 2*sin(Theta2_tilde/2)*sin(Theta1_tilde/2);
Phi(2) = 2*sin(gamma-Theta2_tilde/2)*sin(Theta1_tilde/2);
Phi(3) = 1/2*(Theta1_tilde-cos(Theta2_tilde)*sin(Theta1_tilde));
Phi(4) = -1/2*(-sin(Theta1_tilde)*cos(gamma-Theta2_tilde)+Theta1_tilde*cos_gamma);
Phi(5) = 1/2*(Theta1_tilde-cos(2*gamma-Theta2_tilde)*sin(Theta1_tilde));

%% STEP 4: CALCOLO delle MATRICI B_1, B_2 e B_3

%Calcolo della matrice B_1
B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+Phi(3)*B(:,:,6);
%B_1 = Phi(5)*B(:,:,4)+Phi(4)*B(:,:,5)+2*Phi(3)*B(:,:,6);

%Calcolo della matrice B_2
B_2 = Phi(2)*B(:,:,2)+Phi(1)*B(:,:,3);

%Calcolo della matrice B_3
B_3 = B(:,:,1);  

%% STEP 5: CALCOLO degli INTEGRALI

%INTEGRALE di 1/r sul SETTORE CIRCOLARE in presenza della sola ONDA P
rist1 = Theta1_tilde*I(1)*eye(3,3);

%INTEGRALE di r_i*r_j/r^3 sul SETTORE CIRCOLARE in presenza della sola ONDA P
rist2 = I(2)*B_1+I(3)*B_2+Theta1_tilde*I(4)*B_3;

%INTEGRALE di 1/r^3 sul SETTORE CIRCOLARE in presenza della sola ONDA P
rist3 = Theta1_tilde*I(4)*eye(3,3);

%INTEGRALE di r_i*r_j/r^5 sul SETTORE CIRCOLARE in presenza della sola ONDA P
rist4 = I(5)*B_1+I(6)*B_2+Theta1_tilde*I(7)*B_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f1=@(rho,theta)(rho./((rho.^2+zeta^2).^(1/2)));
% RIST1=integral2(f1,0,pR_P,Theta1,Theta2,'Method','iterated','AbsTol',0,'RelTol',1e-10)*eye(3,3);
% diff1=abs(RIST1-rist1);
% 
% f3=@(rho,theta)(rho./((rho.^2+zeta^2).^(3/2)));
% RIST3=integral2(f3,0,pR_P,Theta1,Theta2,'Method','iterated','AbsTol',0,'RelTol',1e-10)*eye(3,3);
% diff3=abs(RIST3-rist3);

  
end