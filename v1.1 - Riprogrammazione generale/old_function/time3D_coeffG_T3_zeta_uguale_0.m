function GG = time3D_coeffG_T3_zeta_uguale_0(Theta1,Theta2,Info_tri2D)

%La function time3D_coeffG_T3_zeta_uguale_0() prende in INPUT:
%- la variabile Theta1 contenente il primo angolo che individua il
%  triangolo di integrazione
%
%- la variabile Theta2 contenente il secondo angolo che individua il
%  triangolo di integrazione
%
%- la struct Info_Tri2D che contiene le informazioni geometriche sul
%  triangolo figlio corrente (lunghezza dei lati, ampiezza degli angoli,
%  valore del seno e del coseno degli angoli)
%
%e restituisce in OUTPUT:
%- il vettore GG contenente il risultato degli integrali analitici G_i 
%  in presenza di onde S e P sul triangolo nel caso zeta=0
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

%AMPIEZZA deli ANGOLI del triangolo
alpha=Info_tri2D.alpha;
beta=Info_tri2D.beta;
gamma=Info_tri2D.gamma;

%% STEP 2: CALCOLIAMO I COEFFICIENTI 

F=a*sin_beta;

%Calcolo dei COEFFICIENTI g_i (vedi Milroy)
% g(1)=sqrt(1+zeta_frac_R1^2);
% g(2)=sqrt(1+zeta_frac_R2^2);
g_3=1/2*(log((1+cos(Theta1+beta))*(1-cos(Theta2+beta)))-log((1-cos(Theta1+beta))*(1+cos(Theta2+beta))));
% g(4)=log(1/zeta_frac_R1+sqrt(1+1/zeta_frac_R1^2));
% g(5)=log(1/zeta_frac_R2+sqrt(1+1/zeta_frac_R2^2));
% g(7)=pi-acos((-z*cos(Theta2+beta))/sqrt(z2_plus_F2))-acos((z*cos(Theta1+beta))/sqrt(z2_plus_F2));
% g(8)=F^2/(z2_plus_F2);

%% STEP 3: CALCOLIAMO IL RISULTATO DEGLI INTEGRALI ANALITICI G_i 
%          NEL CASO zeta=0

%Calcolo degli INTEGRALI ANALITICI G_i NEL CASO zeta=0 (vedi Milroy)
GG(1)=F*g_3;

GG(7)=F*(-cos(alpha-gamma+Theta1)+cos(alpha-gamma+Theta2)+sin_alpha^2*g_3);

GG(8)=F*(-sin_alpha*sin_beta*g_3+cos(Theta1+alpha)-cos(Theta2+alpha));

GG(9)=F*(sin_beta^2*g_3+cos(Theta1-beta)-cos(Theta2-beta));
     
end