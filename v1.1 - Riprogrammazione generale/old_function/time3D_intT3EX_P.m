function [rist1,rist2,rist3,rist4] = time3D_intT3EX_P(zeta,Theta1,Theta2,Info_tri2D,B)
%function [ris1t,ris2t,ris3t,ris4t] = time3D_intT3EX_P(zeta,p1,p2,B)

%La function time3D_intT3EX_P() prende in INPUT:
%- la variabile zeta contenente la distanza del punto sorgente dalla 
%  sua proiezione nel piano in cui giace il triangolo di campo
%
%- la variabile Theta1 che contiene lampiezza del primo angolo che 
%  individua il triangolo nel piano
%
%- la variabile Theta2 che contiene l'ampiezza del secondo angolo che 
%  individua il triangolo nel piano
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
%  r_ir_j/r^5 per i,j=1,2,3 calcolati sul TRIANGOLO in presenza 
%  della sola ONDA P
%
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%NOTA BENE: in presenza di sole onde P il punto sorgente non si trova mai
%sullo stesso piano del triangolo di campo e quindi zeta > 0
%-------------------------------------------------------------------------


% %% STEP 1: RICAVIAMO le INFORMAZIONI GEOMETRICHE RELATIVE al TRIANGOLO
% 
% %coordinate del primo vertice del triangolo
% T3(1,:) = zeros(1,2);
% %coordinate del secondo vertice del triangolo
% T3(2,:) = p1;
% %coordinate del terzo vertice del triangolo
% T3(3,:) = p2;
% 
% %calcolo dell'area del triangolo
% area = (T3(2,1)*T3(3,2)-T3(2,2)*T3(3,1))/2;
% %Modulo del determinante della matrice jacobina della trasformazione che 
% %mappa il triangolo di riferimento nel triangolo considerato
% absarea2 = 2*abs(area);
% 
% %Lunghezza dei lati del triangolo di integrazione (NB: il primo vertice è (0,0))
% b =  sqrt(T3(2,1)^2+T3(2,2)^2);
% c =  sqrt((T3(3,1)-T3(2,1))^2+(T3(3,2)-T3(2,2))^2);
% a =  sqrt(T3(3,1)^2+T3(3,2)^2);
%   
% %Applico il THM di CARNOT per calcolare i coseni degli angoli del triangolo 
% %di integrazione
% cosa = (b^2+c^2-a^2)/(2*b*c);
% cosb = (a^2+c^2-b^2)/(2*a*c);
% cosc = (a^2+b^2-c^2)/(2*a*b);
%   
% %Applico la relazione che lega i seni degli angoli del triangolo di
% %integrazione all'area
% sina = absarea2/(c*b);
% sinb = absarea2/(a*c);
% sinc = absarea2/(b*a);
% 
% %Calcolo gli angoli del triangolo di integrazione
% anga = asin(sina);
% if(anga<0) anga = anga+pi; end 
% angc = asin(sinc);
% if(angc<0) angc = angc+pi; end
% angb = pi-anga-angc;
   

%% STEP 1: CALCOLIAMO I COEFFICIENTI 

%Calcoliamo i COEFFICIENTI G_i utili a calcolare gli integrali analitici
%in presenza della sola ONDA P (vedi Milroy)
G = time3D_coeff_G_P(zeta,Theta1,Theta2,Info_tri2D);

%% STEP 2: CALCOLIAMO GLI INTEGRALI

%INTEGRALE di 1/r sul TRIANGOLO in presenza della sola ONDA P
rist1 = G(1)*eye(3,3); 

%INTEGRALE di r_i*r_j/r^3 sul TRIANGOLO in presenza della sola ONDA P
rist2 = G(4)*B(:,:,1)+G(5)*B(:,:,2)+G(6)*B(:,:,3)+G(7)*B(:,:,4)+...
    +G(8)*B(:,:,5)+G(9)*B(:,:,6);

%INTEGRALE di 1/r^3 sul TRIANGOLO in presenza della sola ONDA P 
rist3 = G(4)*eye(3,3); 

%INTEGRALE di r_i*r_j/r^5 sul TRIANGOLO in presenza della sola ONDA P
rist4 = G(14)*B(:,:,1)+G(15)*B(:,:,2)+G(16)*B(:,:,3)+G(17)*B(:,:,4)+...
    +G(18)*B(:,:,5)+G(19)*B(:,:,6);
    
end