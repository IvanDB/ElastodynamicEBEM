function Info_tri2D = time3D_InfoTri2D(T3)

%La function time3D_InfoTri2D() prende in INPUT:
%- la matrice 3x2 T3 sulla cui i-esima riga sono presenti le coordinate 
%  dell'i-esimo vertice del triangolo corrente nel sistema di riferimento
%  bidimensionale
%
%e restituisce in OUTPUT:
%- la struct Info_tri2D che contiene le informazioni sul triangolo 
%  corrente: 
%     * lati a = 1°--2° vertice, b = 1°--3° vertice, c = 3°--2° vertice
%     * coseno degli angoli alpha, beta, gamma
%     * seno degli angoli alpha, beta, gamma
%     * angoli alpha, beta, gamma
%
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%N.B. Si suppone che il triangolo considerato sia non degenere!
%Pertanto è necessario verificare prima di richiamare questa function che
%il triangolo sia non degenere
%-------------------------------------------------------------------------


%LUNGHEZZA dei LATI del triangolo 
Info_tri2D.a =  sqrt((T3(2,1)-T3(1,1))^2+(T3(2,2)-T3(1,2))^2);
Info_tri2D.b =  sqrt((T3(3,1)-T3(1,1))^2+(T3(3,2)-T3(1,2))^2);
Info_tri2D.c =  sqrt((T3(3,1)-T3(2,1))^2+(T3(3,2)-T3(2,2))^2);

%Applico il TEOREMA di CARNOT per calcolare i COSENI degli ANGOLI del 
%triangolo 
Info_tri2D.cos_alpha = (Info_tri2D.b^2+Info_tri2D.c^2-Info_tri2D.a^2)/(2*Info_tri2D.b*Info_tri2D.c);
Info_tri2D.cos_beta = (Info_tri2D.a^2+Info_tri2D.c^2-Info_tri2D.b^2)/(2*Info_tri2D.a*Info_tri2D.c);
Info_tri2D.cos_gamma = (Info_tri2D.a^2+Info_tri2D.b^2-Info_tri2D.c^2)/(2*Info_tri2D.a*Info_tri2D.b);

%Info_tri2D.cos_gamma=((T3(2,1)-T3(1,1))*(T3(3,1)-T3(1,1))+(T3(2,2)-T3(1,2))*(T3(3,2)-T3(1,2)))/(Info_tri2D.a*Info_tri2D.b);

%Applico la PRIMA RELAZIONE FONDAMENTALE della trigonometria che lega il 
%seno e il coseno di un angolo: (sin(x))^2+(cos(x))^2=1
%N.B. Il seno degli angoli di un triangolo è sempre positivo
Info_tri2D.sin_alpha = sqrt(1-Info_tri2D.cos_alpha^2);
Info_tri2D.sin_beta = sqrt(1-Info_tri2D.cos_beta^2);
Info_tri2D.sin_gamma = sqrt(1-Info_tri2D.cos_gamma^2);


%Calcoliamo l'AMPIEZZA degli ANGOLI del triangolo
Info_tri2D.alpha=acos(Info_tri2D.cos_alpha);
Info_tri2D.beta=acos(Info_tri2D.cos_beta);
Info_tri2D.gamma=acos(Info_tri2D.cos_gamma);


% p=(Info_tri2D.a+Info_tri2D.b+Info_tri2D.c)/2;
% A=sqrt(p*(p-Info_tri2D.a)*(p-Info_tri2D.b)*(p-Info_tri2D.c));
% Info_tri2D.sin_alpha = 2*A/(Info_tri2D.b*Info_tri2D.c);
% Info_tri2D.sin_beta = 2*A/(Info_tri2D.a*Info_tri2D.c);
% Info_tri2D.sin_gamma = 2*A/(Info_tri2D.a*Info_tri2D.b);
% 
% Info_tri2D.sin_gamma =max(Info_tri2D.c*Info_tri2D.sin_alpha/Info_tri2D.a,sqrt(1-Info_tri2D.cos_gamma^2));
% Info_tri2D.gamma=asin(Info_tri2D.sin_gamma);
% if (Info_tri2D.gamma<0)
%     Info_tri2D.gamma=Info_tri2D.gamma+pi;
% end

% Info_tri2D.gamma2 = Info_tri2D.cos_gamma^2+Info_tri2D.sin_gamma^2
% Info_tri2D.alpha2 = Info_tri2D.cos_alpha^2+Info_tri2D.sin_alpha^2
% Info_tri2D.beta2 = Info_tri2D.cos_beta^2+Info_tri2D.sin_beta^2

return