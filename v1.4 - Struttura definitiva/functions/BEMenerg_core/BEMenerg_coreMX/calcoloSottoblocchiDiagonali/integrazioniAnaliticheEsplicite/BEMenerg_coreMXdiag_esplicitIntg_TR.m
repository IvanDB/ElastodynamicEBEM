function [ris5, ris6] = BEMenerg_coreMXdiag_esplicitIntg_TR(THETA1, THETA2, infoTri2D, B)
% INPUT
% - 
% - 
% OUTPUT
% - 

%% CALCOLO VALORI SINGOLI NUCLEI INTEGRALI 
%Calcolo coeffincienti G_i
G = BEMenerg_coreMXdiag_coeffG_TR(THETA1, THETA2, infoTri2D);

%INTEGRALE di 1/r
ris5 = G(1) * eye(3, 3);

%INTEGRALE di r_i*r_j/r^3 
ris6 = G(7) * B(:, :, 4) + G(8) * B(:, :, 5) + G(9) * B(:, :, 6); 

return