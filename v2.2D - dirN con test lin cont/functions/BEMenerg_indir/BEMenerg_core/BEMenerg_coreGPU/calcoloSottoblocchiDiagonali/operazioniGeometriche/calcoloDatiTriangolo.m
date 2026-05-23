function [tri2D, infoTri2D] = calcoloDatiTriangolo(verticiTriangolo)
%INPUT
% - 
%
%OUTPUT
% - 
% - 

%% CALCOLO INFORMAZIONI GEOMETRICHE TRIANGOLO
%Calcolo LUNGHEZZA LATI 
infoTri2D.a =  sqrt((verticiTriangolo(2, 1) - verticiTriangolo(1, 1))^2 + (verticiTriangolo(2, 2) - verticiTriangolo(1, 2))^2);
infoTri2D.b =  sqrt((verticiTriangolo(3, 1) - verticiTriangolo(1, 1))^2 + (verticiTriangolo(3, 2) - verticiTriangolo(1, 2))^2);
infoTri2D.c =  sqrt((verticiTriangolo(3, 1) - verticiTriangolo(2, 1))^2 + (verticiTriangolo(3, 2) - verticiTriangolo(2, 2))^2);

%Calcolo COSENO ANGOLI mediante TEOREMA di CARNOT
infoTri2D.cos_alpha = (infoTri2D.b^2 + infoTri2D.c^2 - infoTri2D.a^2)/(2*infoTri2D.b*infoTri2D.c);
infoTri2D.cos_beta  = (infoTri2D.a^2 + infoTri2D.c^2 - infoTri2D.b^2)/(2*infoTri2D.a*infoTri2D.c);
infoTri2D.cos_gamma = (infoTri2D.a^2 + infoTri2D.b^2 - infoTri2D.c^2)/(2*infoTri2D.a*infoTri2D.b);

%Calcolo SENO ANGOLI mediante PRIMA RELAZIONE FONDAMENTALE della trigonometria
infoTri2D.sin_alpha = sqrt(1 - infoTri2D.cos_alpha^2);
infoTri2D.sin_beta  = sqrt(1 - infoTri2D.cos_beta^2);
infoTri2D.sin_gamma = sqrt(1 - infoTri2D.cos_gamma^2);

%Calcolo AMPIEZZA ANGOLI
infoTri2D.alpha = acos(infoTri2D.cos_alpha);
infoTri2D.beta  = acos(infoTri2D.cos_beta);
infoTri2D.gamma = acos(infoTri2D.cos_gamma);

%% CALCOLO VERTICI POST-ROTAZIONE
%Calcolo coordinate primo vertice
tri2D(1, :) = zeros(1, 2);

%Calcolo coordinate secondo vertice
tri2D(2, :) = [infoTri2D.a 0];

%Calcolo coordinate terzo vertice
tri2D(3, :) = [infoTri2D.b * infoTri2D.cos_gamma infoTri2D.b * infoTri2D.sin_gamma];

return