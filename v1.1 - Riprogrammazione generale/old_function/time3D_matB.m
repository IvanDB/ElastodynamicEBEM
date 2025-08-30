function B = time3D_matB(zeta,Info_Tri2D,sign_prod)
%function [B,gamma] = time3D_matB(zeta,p1,p2)

%La function time3D_matB() prende in INPUT:
%- la variabile zeta contenente la distanza tra il punto sorgente corrente
%  e il piano contenente il triangolo di campo corrente
%
%- la struct Info_Tri2D che contiene le informazioni geoemtriche sul
%  triangolo figlio corrente (lunghezza dei lati, ampiezza degli angoli,
%  valore del seno e del coseno degli angoli)
%
%- la variabile sign_prod contenente il segno del prodotto scalare tra il 
%  versore normale al piano del triangolo di campo e il vettore distanza
%  tra il punto sorgente sp e la sua proiezione sp_plane nel piano del 
%  triangolo di campo
%
%e restituisce in OUTPUT:
%- l'array 3D B di formato 3x3x6, costituito da 6 matrici 3x3 ciascuna 
%  delle quali contiene una delle 6 tipologie di coefficienti che
%  permettono di esprimere il prodotto r_ir_k al variare di i,k=1,2,3
%
%-----------------------------------------------------------------------

%% STEP 1: RICAVO le INFORMAZIONI GEOMETRICHE RELATIVE al TRIANGOLO

%COSENO dell'ANGOLO gamma OPPOSTO al LATO c 
cos_gamma=Info_Tri2D.cos_gamma;
%SENO dell'ANGOLO gamma OPPOSTO al LATO c 
sin_gamma=Info_Tri2D.sin_gamma;

% %COORDINATE del PRIMO VERTICE del TRIANGOLO 
% %NB: il primo vertice è (0,0)
% T3(1,:) = zeros(1,2);
% %COORDINATE del SECONDO VERTICE del TRIANGOLO
% T3(2,:) = p1;
% %COORDINATE del TERZO VERTICE del TRIANGOLO
% T3(3,:) = p2;
%LUNGHEZZA dei LATI del TRIANGOLO di integrazione 
%LATO a: Primo vertice -- Secondo vertice
%a =  sqrt(T3(2,1)^2 + T3(2,2)^2);
%LATO b: Primo vertice -- Terzo vertice
%b =  sqrt(T3(3,1)^2 + T3(3,2)^2);
%LATO c: Terzo vertice -- Secondo vertice
% c =  sqrt((T3(3,1) - T3(2,1))^2 + (T3(3,2)-T3(2,2))^2);
%COSENO dell'ANGOLO gamma OPPOSTO al LATO c 
%cos_gamma=(a^2+b^2-c^2)/(2*a*b);
%SENO dell'ANGOLO gamma OPPOSTO al LATO c 
%sin_gamma=sqrt(1-cos_gamma^2);
%ANGOLO gamma RELATIVO al PRIMO VERTICE del TRIANGOLO di integrazione
%gamma = acos(cos_gamma);
% %Calcolo dell'area del triangolo
% area = (T3(2,1)*T3(3,2)-T3(2,2)*T3(3,1))/2;
% %Doppio del modulo dell'area
% absarea2 = 2*abs(area);
%
% %Lunghezza dei lati del triangolo di integrazione (NB: il primo vertice è (0,0))
% b =  sqrt(T3(2,1)^2+T3(2,2)^2);
% c =  sqrt((T3(3,1)-T3(2,1))^2+(T3(3,2)-T3(2,2))^2);
% a =  sqrt(T3(3,1)^2+T3(3,2)^2);
%   
% %Applico la relazione che lega l'area e il seno dell'angolo relativo al 
% %primo vertice del triangolo di integrazione
% sinc = absarea2/(b*a);
% 
% %Calcolo l'angolo relativo al primo vertice del triangolo di integrazione
% angc = asin(sinc);
% if(angc<0) angc = angc+pi; end


%% STEP 2: CALCOLO della MATRICE A

%Inizializzazione della matrice 3x3 A
A = zeros(3,3);

%Calcolo della prima riga
A(1,1)=1;
%A(1,2:3)=zeros(1,2);

%Calcolo della seconda riga
A(2,1:2) = [cos_gamma sin_gamma];
%A(2,3)=0;

%Calcolo della terza riga
A(3,3) = -sign_prod;
%A(3,1:2)=zeros(1,2);

% %Inizializzazione della matrice A
% A = zeros(3,3);
% %Calcolo della prima riga
% A(1,1:2) = T3(2,:)/b;
% %Calcolo della seconda riga
% A(2,1:2) = T3(3,:)/a;
% %Calcolo della terza riga
% A(3,3) = 1;

%% STEP 3: CALCOLO delle MATRICI B_i

%Matrice B_1
B(:,:,1) = zeta^2*(A(3,:)'*A(3,:));
%Matrice B_2
B(:,:,2) = zeta*(A(1,:)'*A(3,:)+A(3,:)'*A(1,:))/sin_gamma;
%Matrice B_3
B(:,:,3) = zeta*(A(2,:)'*A(3,:)+A(3,:)'*A(2,:))/sin_gamma;
%Matrice B_4
B(:,:,4) = (A(1,:)'*A(1,:))/sin_gamma^2;
%Matrice B_5
B(:,:,5) = (A(1,:)'*A(2,:)+A(2,:)'*A(1,:))/sin_gamma^2;
%Matrice B_6
B(:,:,6) = (A(2,:)'*A(2,:))/sin_gamma^2;

% %Matrice B_5
% B(:,:,5) = [2*cos_gamma/sin_gamma^2 1/sin_gamma 0; 1/sin_gamma 0 0;0 0 0];
% %Matrice B_6
% B(:,:,6) = [(cos_gamma/sin_gamma)^2 cos_gamma/sin_gamma 0; cos_gamma/sin_gamma 1 0;0 0 0];

% %Matrice B_1
% B(:,:,1) = zeta^2*A(3,:)'*A(3,:);
% %Matrice B_2
% B(:,:,2) = zeta*(A(1,:)'*A(3,:)+A(3,:)'*A(1,:))/sinc;
% %Matrice B_3
% B(:,:,3) = zeta*(A(2,:)'*A(3,:)+A(3,:)'*A(2,:))/sinc;
% %Matrice B_4
% B(:,:,4) = A(1,:)'*A(1,:)/sinc^2;
% %Matrice B_5
% B(:,:,5) = (A(1,:)'*A(2,:)+A(2,:)'*A(1,:))/sinc^2;
% %Matrice B_6
% B(:,:,6) = A(2,:)'*A(2,:)/sinc^2;

return