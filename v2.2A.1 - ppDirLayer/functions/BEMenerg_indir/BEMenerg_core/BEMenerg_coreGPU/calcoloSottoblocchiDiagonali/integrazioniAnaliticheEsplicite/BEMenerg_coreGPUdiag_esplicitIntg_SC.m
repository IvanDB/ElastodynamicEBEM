function [ris5, ris6] = BEMenerg_coreGPUdiag_esplicitIntg_SC(rhoS, THETA1, THETA2, infoTri2D, B)
% INPUT
% - 
% - 
% OUTPUT
% - 

%% CALCOLO COSTANTI ANGOLARI RICORRENTI ed ESTRAZIONE INFORMAZIONI GEOMETRICHE TRIANGOLO

%Calcolo delle COSTANTI RICORRENTI negli INTEGRALI in THETA
THETA1tilde = THETA2 - THETA1;
THETA2tilde = THETA2 + THETA1;

%Estrazione COSENO e AMPIEZZA ANGOLO GAMMA
cos_gamma = infoTri2D.cos_gamma;
gamma = infoTri2D.gamma;

%% CALCOLO VALORI SINGOLI NUCLEI INTEGRALI 
%Calcolo integrali in rho
I(1) = rhoS;
I(2) = rhoS;

%Calcolo integrali in theta
Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
Phi(4) = -1/2 * ((THETA1tilde * cos_gamma) - (sin(THETA1tilde) * cos(gamma-THETA2tilde)));
Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));

%Calcolo MATRICE B_1 (B_2 e B_3 non utilizzate) 
B_1 = Phi(5) * B(:, :, 4) + Phi(4) * B(:, :, 5) + Phi(3) * B(:, :, 6);

%INTEGRALE di 1/r
ris5 = THETA1tilde * I(1) * eye(3, 3);

%INTEGRALE di r_i*r_j/r^3
ris6 = I(2) * B_1;

return