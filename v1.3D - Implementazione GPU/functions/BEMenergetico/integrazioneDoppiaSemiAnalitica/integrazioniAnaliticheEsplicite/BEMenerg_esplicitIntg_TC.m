function [ris1, ris2, ris3, ris4] = BEMenerg_esplicitIntg_TC(zeta, RS, rhoS, THETA1, THETA2, infoTri2D, B)
% INPUT
% - 
% - 
% OUTPUT
% - 

%Controllo in funzione di zeta

%Caso zeta = 0 (triangoli complanari)
G = BEMenerg_coeffG_TC(0, RS, rhoS, THETA1, THETA2, infoTri2D);

%Calcolo INTEGRALI

%INTEGRALE di 1/r
ris1 = G(1) * eye(3, 3);

%INTEGRALE di r_i*r_j/r^3
ris2 = G(7) * B(:, :, 4) + G(8) * B(:, :, 5) + G(9) * B(:, :, 6);

%INTEGRALE di 1/r^3
ris3 = G(4) * eye(3, 3);

%INTEGRALE di r_i*r_j/r^5
ris4 = G(17) * B(:, :, 4) + G(18) * B(:, :, 5) + G(19) * B(:, :, 6);

return
end