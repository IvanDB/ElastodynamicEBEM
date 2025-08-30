function GG = time3D_coeffG_TC_zeta_uguale_0(R_S,Theta1,Theta2,Info_tri2D)

%La function time3D_coeffG_TC_zeta_uguale_0() prende in INPUT:
%- la variabile R_S=c_S*Delta contenente il raggio della superficie
%  sferica individuata dal fronte d'onda relativo alle onde S
%
%- la variabile Theta1 contenente il primo angolo di integrazione che
%  individua il trapezio/triangolo circolare
%
%- la variabile Theta2 contenente il secondo angolo di integrazione che
%  individua il trapezio/triangolo circolare
%
%- la struct Info_Tri2D che contiene le informazioni geometriche sul
%  triangolo figlio corrente (lunghezza dei lati, ampiezza degli angoli,
%  valore del seno e del coseno degli angoli)
%
%e restituisce in OUTPUT:
%- il vettore GG contenente il risultato degli integrali analitici G_i 
%  nel caso in cui zeta=0
%
%-------------------------------------------------------------------------

%% STEP 1: RICAVIAMO le INFORMAZIONI GEOMETRICHE RELATIVE al TRIANGOLO

%LATI del triangolo
a=Info_tri2D.a;
% b=Info_tri2D.b;
% c=Info_tri2D.c;

%SENO degli ANGOLI del triangolo
sin_beta=Info_tri2D.sin_beta;
sin_alpha=Info_tri2D.sin_alpha;
% sin_gamma=Info_tri2D.sin_gamma;

%COSENO degli ANGOLI del triangolo
cos_alpha=Info_tri2D.cos_alpha;
cos_beta=Info_tri2D.cos_beta;
% cos_gamma=Info_tri2D.cos_gamma;

%AMPIEZZA deli ANGOLI del triangolo
alpha=Info_tri2D.alpha;
beta=Info_tri2D.beta;
gamma=Info_tri2D.gamma;

%% STEP 2: CALCOLIAMO I COEFFICIENTI 

Theta1_tilde=Theta2-Theta1;
Theta2_tilde=Theta2+Theta1;

Phi(1) = 2*sin(Theta2_tilde/2)*sin(Theta1_tilde/2);
Phi(2) = 2*sin(gamma-Theta2_tilde/2)*sin(Theta1_tilde/2);
Phi(3) = 1/2*(Theta1_tilde-cos(Theta2_tilde)*sin(Theta1_tilde));
Phi(4) = -1/2*(-sin(Theta1_tilde)*cos(gamma-Theta2_tilde)+Theta1_tilde*cos(gamma));
Phi(5) = 1/2*(Theta1_tilde-cos(2*gamma-Theta2_tilde)*sin(Theta1_tilde));

F=a*sin_beta;

%Calcolo dei COEFFICIENTI g_i (vedi Milroy)
g3=1/2*(log((1+cos(Theta1+beta))*(1-cos(Theta2+beta)))-log((1-cos(Theta1+beta))*(1+cos(Theta2+beta))));

%% STEP 3: CALCOLIAMO IL RISULTATO DEGLI INTEGRALI ANALITICI G_i sul
%          TRAPEZIO CIRCOLARE NEL CASO zeta UGUALE A ZERO

%Calcolo degli INTEGRALI ANALITICI G_i nel caso zeta=0(vedi Milroy)
GG(1)=F*g3-R_S*Theta1_tilde;

GG(4)=1/F*(cos(Theta2+beta)-cos(Theta1+beta))+1/(R_S)*Theta1_tilde;

GG(7)=-R_S*Phi(5)+F*(cos(alpha-gamma+Theta2)-cos(alpha-gamma+Theta1)+sin_alpha^2*g3);

GG(8)=-R_S*Phi(4)+F*(cos(Theta1+alpha)-cos(Theta2+alpha)-sin_alpha*sin_beta*g3);

GG(9)=-R_S*Phi(3)+F*(cos(beta-Theta1)-cos(beta-Theta2)+sin_beta^2*g3);

GG(17)=1/(3*F)*(((1+cos_alpha^2)*cos(Theta2+beta)-((sin(Theta2+beta))^2)*cos(alpha-gamma+Theta2))...
      -((1+cos_alpha^2)*cos(Theta1+beta)-((sin(Theta1+beta))^2)*cos(alpha-gamma+Theta1)))+1/(R_S)*Phi(5);
  
GG(18)=1/(3*F)*((cos(Theta2+beta)*(cos_alpha*cos_beta+cos(alpha+beta))+((sin(Theta2+beta))^2)*cos(alpha+Theta2))...
      -(cos(Theta1+beta)*(cos_alpha*cos_beta+cos(alpha+beta))+((sin(Theta1+beta))^2)*cos(alpha+Theta1)))...
      +1/(R_S)*Phi(4);
     
GG(19)=1/(3*F)*(((1+cos_beta^2)*cos(Theta2+beta)+((sin(Theta2+beta))^2)*cos(Theta2-beta))...
      -((1+cos_beta^2)*cos(Theta1+beta)+((sin(Theta1+beta))^2)*cos(beta-Theta1)))...
      +1/(R_S)*Phi(3);

  
  
%  funz1=@(theta,rho)(rho)./((rho.^2).^(1/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_1=integral2(funz1,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_1=abs(G_hat_1-GG(1))   
%  
% funz4=@(theta,rho)(rho)./((rho.^2).^(3/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_4=integral2(funz4,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_4=abs(G_hat_4-GG(4)) 
%   
% funz7=@(theta,rho)((rho.^3).*(sin(gamma-theta)).^2)./((rho.^2).^(3/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_7=integral2(funz7,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_7=abs(G_hat_7-GG(7))   
%  
% funz8=@(theta,rho)(rho.^3).*(sin(theta)).*(sin(gamma-theta))./((rho.^2).^(3/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_8=integral2(funz8,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_8=abs(G_hat_8-GG(8))  
%   
% funz9=@(theta,rho)(rho.^3).*(sin(theta)).^2./((rho.^2).^(3/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_9=integral2(funz9,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_9=abs(G_hat_9-GG(9))  
%   
% funz17=@(theta,rho)((rho.^3).*(sin(gamma-theta)).^2)./((rho.^2).^(5/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_17=integral2(funz17,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_17=abs(G_hat_17-GG(17))   
%  
% funz18=@(theta,rho)(rho.^3).*(sin(theta)).*(sin(gamma-theta))./((rho.^2).^(5/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_18=integral2(funz18,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_18=abs(G_hat_18-GG(18))  
%   
% funz19=@(theta,rho)(rho.^3).*(sin(theta)).^2./((rho.^2).^(5/2));  
% R=@(theta)F./(sin(theta+beta));
% G_hat_19=integral2(funz19,Theta1,Theta2,R_S,R,'Method','iterated','AbsTol',1e-16,'RelTol',1e-16);
% E_hat_19=abs(G_hat_19-GG(19))

     
end