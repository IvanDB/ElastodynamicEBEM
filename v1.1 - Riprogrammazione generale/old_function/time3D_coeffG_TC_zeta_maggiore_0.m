function GG = time3D_coeffG_TC_zeta_maggiore_0(zeta,R_S,pR_S,Theta1,Theta2,Info_tri2D)

%La function time3D_coeffG_TC_zeta_maggiore_0() prende in INPUT:
%- la variabile zeta contenente la distanza tra il punto sorgente e il
%  piano individuato dal triangolo di campo
%
%- la variabile R_S=c_S*Delta contenente il raggio della superifice 
%  sferica individuata dal fronte d'onda relativo alle onde S
%
%- la variabile pR_S contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda S nel piano del triangolo di campo
%  pR_S = sqrt([c_S*Delta]^2-zeta^2)
%
%- la variabile RS_tilde = sqrt(pR_S^2+zeta^2)=c_S*Delta
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
%  sul trapezio/triangolo circolare nel caso in cui zeta>0
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
cos_gamma=Info_tri2D.cos_gamma;

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
Phi(4) = -1/2*(-sin(Theta1_tilde)*cos(gamma-Theta2_tilde)+Theta1_tilde*cos_gamma);
Phi(5) = 1/2*(Theta1_tilde-cos(2*gamma-Theta2_tilde)*sin(Theta1_tilde));

F=a*sin_beta;

R1=F/sin(Theta1+beta);
R2=F/sin(Theta2+beta);

%Calcolo delle COSTANTI (vedi Milroy)
zeta_frac_R1 = zeta/R1;
zeta_frac_R2 = zeta/R2;
z2_plus_F2 = zeta^2+F^2;

%Calcolo dei COEFFICIENTI g_i (vedi Milroy)
g(1)=sqrt(1+zeta_frac_R1^2);
g(2)=sqrt(1+zeta_frac_R2^2);
g(3)=1/2*(log((g(1)+cos(Theta1+beta))*(g(2)-cos(Theta2+beta)))-log((g(1)-cos(Theta1+beta))*(g(2)+cos(Theta2+beta))));
g(4)=log(1/zeta_frac_R1+sqrt(1+1/zeta_frac_R1^2));
g(5)=log(1/zeta_frac_R2+sqrt(1+1/zeta_frac_R2^2));
g(7)=pi-acos((-zeta*cos(Theta2+beta))/sqrt(z2_plus_F2))-acos((zeta*cos(Theta1+beta))/sqrt(z2_plus_F2));
g(8)=F^2/(z2_plus_F2);

%% STEP 3: CALCOLIAMO IL RISULTATO DEGLI INTEGRALI ANALITICI G_i sul
%          TRAPEZIO CIRCOLARE NEL CASO Z MAGGIORE DI ZERO

%Calcolo degli INTEGRALI ANALITICI G_i sul TRAPEZIO CIRCOLARE (vedi Milroy)
GG(1)=-R_S*Theta1_tilde+zeta*g(7)+F*g(3);

GG(4)=1/(R_S)*Theta1_tilde-1/zeta*g(7);

GG(5)=cos(gamma-Theta2)*g(5)-cos(gamma-Theta1)*g(4)-cos_alpha*g(3)...
    +(pR_S/R_S-log((pR_S+R_S)/zeta))*(cos(gamma-Theta2)-cos(gamma-Theta1));

GG(6)=-cos(Theta2)*g(5)+cos(Theta1)*g(4)-cos_beta*g(3)...
    -(pR_S/R_S-log((pR_S+R_S)/zeta))*(cos(Theta2)-cos(Theta1));

GG(7)=zeta*g(7)-((R_S^2+zeta^2)/R_S)*Phi(5)+F*(-g(1)*cos(alpha-gamma+Theta1)+g(2)*cos(alpha-gamma+Theta2)+sin_alpha^2*g(3));

GG(8)=-zeta*cos_gamma*g(7)-((R_S^2+zeta^2)/R_S)*Phi(4)...
    +F*(-sin_alpha*sin_beta*g(3)+g(1)*cos(Theta1+alpha)-g(2)*cos(Theta2+alpha));

GG(9)=zeta*g(7)-((R_S^2+zeta^2)/R_S)*Phi(3)...
    +F*(sin_beta^2*g(3)+g(1)*cos(Theta1-beta)-g(2)*cos(Theta2-beta));

GG(14)=1/(3*R_S^3)*Theta1_tilde...
    +1/(3*zeta^2)*((F/(z2_plus_F2)*(cos(Theta1+beta)/g(1)-cos(Theta2+beta)/g(2))-1/zeta*g(7)));

GG(15)=-pR_S^3/(3*zeta^2*R_S^3)*Phi(2)+1/(3*zeta^2)*((sin_alpha*sin(Theta2+beta)...
    -g(8)*cos_alpha*cos(Theta2+beta))/g(2)+(-sin_alpha*sin(Theta1+beta)+g(8)*cos_alpha*cos(Theta1+beta))/g(1));

GG(16)=-pR_S^3/(3*zeta^2*R_S^3)*Phi(1)+1/(3*zeta^2)*((sin_beta*sin(Theta1+beta)...
    +g(8)*cos_beta*cos(Theta1+beta))/g(1)+(-sin_beta*sin(Theta2+beta)-g(8)*cos_beta*cos(Theta2+beta))/g(2));

GG(17)=1/3*(-(sin(Theta2+beta))^2*cos(alpha-gamma+Theta2)+g(8)*(cos_alpha)^2*cos(Theta2+beta))/(F*g(2))...
         -1/3*(-(sin(Theta1+beta))^2*cos(alpha-gamma+Theta1)+g(8)*(cos_alpha)^2*cos(Theta1+beta))/(F*g(1))...
         -g(7)/(3*zeta)+(1/R_S-zeta^2/(3*R_S^3))*Phi(5);

GG(18)=1/3*((sin(Theta2+beta))^2*cos(alpha+Theta2)+g(8)*(cos_alpha*cos_beta*cos(Theta2+beta)))/(F*g(2))...
         -1/3*((sin(Theta1+beta))^2*cos(alpha+Theta1)+g(8)*(cos_alpha*cos_beta*cos(Theta1+beta)))/(F*g(1))...
         +g(7)/(3*zeta)*cos_gamma+(1/R_S-zeta^2/(3*R_S^3))*Phi(4);
     
GG(19)=1/3*((sin(Theta2+beta))^2*cos(beta-Theta2)+g(8)*(cos_beta^2*cos(Theta2+beta)))/(F*g(2))...
         -1/3*((sin(Theta1+beta))^2*cos(beta-Theta1)+g(8)*(cos_beta^2*cos(Theta1+beta)))/(F*g(1))...
         -g(7)/(3*zeta)+(1/R_S-zeta^2/(3*R_S^3))*Phi(3);  
     
end