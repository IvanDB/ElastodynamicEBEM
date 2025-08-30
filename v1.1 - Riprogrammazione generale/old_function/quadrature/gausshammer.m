function [gha, ghw] = gausshammer(mxghp)
% ****FUNZIONE NON REVISIONATA****
% INPUT 
%   - mxghp: intero contenente il massimo numero di nodi

% OUTPUT:
%   - gha: matrice mxghp x mxghp x 3 contenente le coordinate dei nodi di
%   Gauss-Hammer per il triangolo di default
%   - ghw: matrice mxghp x mxghp contenente i pesi Gauss-Hammer


%% Allocazione matrici
gha = zeros(mxghp, mxghp, 3);
ghw = zeros(mxghp, mxghp);

%Allocazione vettore ausiliario necessario per il calcolo
val = zeros(1, 3);

%% FORMULA a 1 NODO

val(1) = 1/3;

gha(1,1,1) = val(1);
gha(1,1,2) = val(1);
gha(1,1,3) = val(1);

ghw(1,1) = 1/2;

%% FORMULA di quadratura a 3 NODI

val(1) = 2/3;
val(2) = 1/6;
val(3) = 1/6;

%Numero di nodi - Indice 
gha(3,1,1) = val(1);
gha(3,1,2) = val(2);
gha(3,1,3) = val(3);
gha(3,2,1) = val(3);
gha(3,2,2) = val(1);
gha(3,2,3) = val(2);
gha(3,3,1) = val(2);
gha(3,3,2) = val(3);
gha(3,3,3) = val(1);

val(1) = 1/6;

ghw(3,1) = val(1);
ghw(3,2) = val(1);
ghw(3,3) = val(1);

%% FORMULA di quadratura a 7 NODI

gha(7,1,1) = 0.1012865073235;
gha(7,2,1) = 0.7974269853531;
gha(7,3,1) = 0.1012865073235;
gha(7,4,1) = 0.4701420641051;
gha(7,5,1) = 0.4701420641051;
gha(7,6,1) = 0.0597158717898;
gha(7,7,1) = 0.3333333333333;
gha(7,1,2) = 0.1012865073235;
gha(7,2,2) = 0.1012865073235;
gha(7,3,2) = 0.7974269853531;
gha(7,4,2) = 0.0597158717898;
gha(7,5,2) = 0.4701420641051;
gha(7,6,2) = 0.4701420641051;
gha(7,7,2) = 0.3333333333333;

for i1=1:7
    gha(7,i1,3) = 1-gha(7,i1,1)-gha(7,i1,2);
end

ghw(7,1) = 0.1259391805448/2;
ghw(7,2) = 0.1259391805448/2;
ghw(7,3) = 0.1259391805448/2;
ghw(7,4) = 0.1323941527885/2;
ghw(7,5) = 0.1323941527885/2;
ghw(7,6) = 0.1323941527885/2;
ghw(7,7) = 0.225/2;

%% FORMULA di quadratura a 12 NODI

val(1) = 0.873821971016996;
val(2) = (1-val(1))/2;
val(3) = val(2);

gha(12,1,1) = val(1);
gha(12,1,2) = val(2);
gha(12,1,3) = val(3);
gha(12,2,1) = val(3);
gha(12,2,2) = val(1);
gha(12,2,3) = val(2);
gha(12,3,1) = val(2);
gha(12,3,2) = val(3);
gha(12,3,3) = val(1);

val(1) = 0.501426509658179;
val(2) = (1-val(1))/2;
val(3) = val(2);

gha(12,4,1) = val(1);
gha(12,4,2) = val(2);
gha(12,4,3) = val(3);
gha(12,5,1) = val(3);
gha(12,5,2) = val(1);
gha(12,5,3) = val(2);
gha(12,6,1) = val(2);
gha(12,6,2) = val(3);
gha(12,6,3) = val(1);

val(1) = 0.636502499121399;
val(2) = 0.310352451033784;
val(3) = 1-val(1)-val(2);

gha(12,7,1)  = val(1);
gha(12,7,2)  = val(2);
gha(12,7,3)  = val(3);
gha(12,8,1)  = val(1);
gha(12,8,2)  = val(3);
gha(12,8,3)  = val(2);
gha(12,9,1)  = val(2);
gha(12,9,2)  = val(1);
gha(12,9,3)  = val(3);
gha(12,10,1) = val(2);
gha(12,10,2) = val(3);
gha(12,10,3) = val(1);
gha(12,11,1) = val(3);
gha(12,11,2) = val(1);
gha(12,11,3) = val(2);
gha(12,12,1) = val(3);
gha(12,12,2) = val(2);
gha(12,12,3) = val(1);

val(1) = 0.050844906370207/2;
val(2) = 0.116786275726379/2;
val(3) = (0.5-3*val(1)-3*val(2))/6;

ghw(12,1)  = val(1);
ghw(12,2)  = val(1);
ghw(12,3)  = val(1);
ghw(12,4)  = val(2);
ghw(12,5)  = val(2);
ghw(12,6)  = val(2);
ghw(12,7)  = val(3);
ghw(12,8)  = val(3);
ghw(12,9)  = val(3);
ghw(12,10) = val(3);
ghw(12,11) = val(3);
ghw(12,12) = val(3);

%% FORMULA di quadratura a 19 NODI

val(1) = 1/3;

gha(19,1,1:3) = val(1);

ghw(19,1) = 0.09713579628279610/2;

val(1) = 0.02063496160252593;
val(2) = 0.48968251919873700;

gha(19,2,1) = val(1);
gha(19,2,2) = val(2);
gha(19,2,3) = val(2);

gha(19,3,1) = val(2);
gha(19,3,2) = val(1);
gha(19,3,3) = val(2);

gha(19,4,1) = val(2);
gha(19,4,2) = val(2);
gha(19,4,3) = val(1);

ghw(19,2:4) = 0.09400410068141950/6;

val(1) = 0.1258208170141290;
val(2) = 0.4370895914929355;

gha(19,5,1) = val(1);
gha(19,5,2) = val(2);
gha(19,5,3) = val(2);

gha(19,6,1) = val(2);
gha(19,6,2) = val(1);
gha(19,6,3) = val(2);

gha(19,7,1) = val(2);
gha(19,7,2) = val(2);
gha(19,7,3) = val(1);

ghw(19,5:7) = 0.2334826230143263/6;

val(1) = 0.6235929287619356;
val(2) = 0.1882035356190322;

gha(19,8,1) = val(1);
gha(19,8,2) = val(2);
gha(19,8,3) = val(2);

gha(19,9,1) = val(2);
gha(19,9,2) = val(1);
gha(19,9,3) = val(2);

gha(19,10,1) = val(2);
gha(19,10,2) = val(2);
gha(19,10,3) = val(1);

ghw(19,8:10) = 0.2389432167816273/6;

val(1) = 0.91054097321109410;
val(2) = 0.04472951339445297;

gha(19,11,1) = val(1);
gha(19,11,2) = val(2);
gha(19,11,3) = val(2);

gha(19,12,1) = val(2);
gha(19,12,2) = val(1);
gha(19,12,3) = val(2);

gha(19,13,1) = val(2);
gha(19,13,2) = val(2);
gha(19,13,3) = val(1);

ghw(19,11:13) = 0.07673302697609430/6;

val(1) = 0.03683841205473626;
val(2) = 0.74119859878449800;
val(3) = 1-val(1)-val(2);

gha(19,14,1) = val(1);
gha(19,14,2) = val(2);
gha(19,14,3) = val(3);

gha(19,15,1) = val(1);
gha(19,15,2) = val(3);
gha(19,15,3) = val(2);

gha(19,16,1) = val(2);
gha(19,16,2) = val(1);
gha(19,16,3) = val(3);

gha(19,17,1) = val(2);
gha(19,17,2) = val(3);
gha(19,17,3) = val(1);

gha(19,18,1) = val(3);
gha(19,18,2) = val(1);
gha(19,18,3) = val(2);

gha(19,19,1) = val(3);
gha(19,19,2) = val(2);
gha(19,19,3) = val(1);

ghw(19,14:19) = 0.2597012362637364/12;

%% FORMULA di quadratura a 28 NODI

val(1) = 1/3;

gha(28,1,1:3) = val(1);

ghw(28,1) = 0.08797730116222190/2;

val(1) = 0.94802171814342330;
val(2) = 0.02598914092828833;

gha(28,2,1) = val(1);
gha(28,2,2) = val(2);
gha(28,2,3) = val(2);

gha(28,3,1) = val(2);
gha(28,3,2) = val(1);
gha(28,3,3) = val(2);

gha(28,4,1) = val(2);
gha(28,4,2) = val(2);
gha(28,4,3) = val(1);

ghw(28,2:4) = 0.02623293466120857/6;

val(1) = 0.81142499470415460;
val(2) = 0.09428750264792270;

gha(28,5,1) = val(1);
gha(28,5,2) = val(2);
gha(28,5,3) = val(2);

gha(28,6,1) = val(2);
gha(28,6,2) = val(1);
gha(28,6,3) = val(2);

gha(28,7,1) = val(2);
gha(28,7,2) = val(2);
gha(28,7,3) = val(1);

ghw(28,5:7) = 0.1142447159818060/6;

val(1) = 0.01072644996557060;
val(2) = 0.49463677501721470;

gha(28,8,1) = val(1);
gha(28,8,2) = val(2);
gha(28,8,3) = val(2);

gha(28,9,1) = val(2);
gha(28,9,2) = val(1);
gha(28,9,3) = val(2);

gha(28,10,1) = val(2);
gha(28,10,2) = val(2);
gha(28,10,3) = val(1);

ghw(28,8:10) = 0.05656634416839376/6;

val(1) = 0.5853132347709715;
val(2) = 0.2073433826145142;

gha(28,11,1) = val(1);
gha(28,11,2) = val(2);
gha(28,11,3) = val(2);

gha(28,12,1) = val(2);
gha(28,12,2) = val(1);
gha(28,12,3) = val(2);

gha(28,13,1) = val(2);
gha(28,13,2) = val(2);
gha(28,13,3) = val(1);

ghw(28,11:13) = 0.2164790926342230/6;

val(1) = 0.1221843885990187;
val(2) = 0.4389078057004907;

gha(28,14,1) = val(1);
gha(28,14,2) = val(2);
gha(28,14,3) = val(2);

gha(28,15,1) = val(2);
gha(28,15,2) = val(1);
gha(28,15,3) = val(2);

gha(28,16,1) = val(2);
gha(28,16,2) = val(2);
gha(28,16,3) = val(1);

ghw(28,14:16) = 0.2079874161166116/6;

val(1) = 0;
val(2) = 0.8588702812826364;
val(3) = 1-val(1)-val(2);

gha(28,17,1) = val(1);
gha(28,17,2) = val(2);
gha(28,17,3) = val(3);

gha(28,18,1) = val(1);
gha(28,18,2) = val(3);
gha(28,18,3) = val(2);

gha(28,19,1) = val(2);
gha(28,19,2) = val(1);
gha(28,19,3) = val(3);

gha(28,20,1) = val(2);
gha(28,20,2) = val(3);
gha(28,20,3) = val(1);

gha(28,21,1) = val(3);
gha(28,21,2) = val(1);
gha(28,21,3) = val(2);

gha(28,22,1) = val(3);
gha(28,22,2) = val(2);
gha(28,22,3) = val(1);

ghw(28,17:22) = 0.04417430269980344/12;

val(1) = 0.04484167758913055;
val(2) = 0.67793765488259020;
val(3) = 1-val(1)-val(2);

gha(28,23,1) = val(1);
gha(28,23,2) = val(2);
gha(28,23,3) = val(3);

gha(28,24,1) = val(1);
gha(28,24,2) = val(3);
gha(28,24,3) = val(2);

gha(28,25,1) = val(2);
gha(28,25,2) = val(1);
gha(28,25,3) = val(3);

gha(28,26,1) = val(2);
gha(28,26,2) = val(3);
gha(28,26,3) = val(1);

gha(28,27,1) = val(3);
gha(28,27,2) = val(1);
gha(28,27,3) = val(2);

gha(28,28,1) = val(3);
gha(28,28,2) = val(2);
gha(28,28,3) = val(1);

ghw(28,23:28) = 0.2463378925757316/12;

return