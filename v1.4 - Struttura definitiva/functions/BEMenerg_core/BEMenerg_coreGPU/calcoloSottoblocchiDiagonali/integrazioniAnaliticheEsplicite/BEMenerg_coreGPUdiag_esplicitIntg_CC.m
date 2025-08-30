function [ris1, ris2, ris3, ris4] = BEMenerg_coreGPUdiag_esplicitIntg_CC(rhoS, rhoP, THETA1, THETA2, infoTri2D, B)
% INPUT
% - 
% - 
% OUTPUT
% - 

%% CALCOLO COSTANTI ANGOLARI RICORRENTI ed ESTRAZIONE INFORMAZIONI GEOMETRICHE TRIANGOLO

%Calcolo COSTANTI RICORRENTI
THETA1tilde = THETA2 - THETA1;
THETA2tilde = THETA2 + THETA1;

%Estrazione COSENO e AMPIEZZA ANGOLO GAMMA
cos_gamma = infoTri2D.cos_gamma;
gamma = infoTri2D.gamma;

%% CALCOLO VALORI SINGOLI NUCLEI INTEGRALI 

%Calcolo integrali parziali in rho
I(1) = rhoP - rhoS;
I(2) = I(1);
I(4) = -(1/rhoP) + (1/rhoS);
I(5) = I(4);

%Calcolo integrali parziali in theta
Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
Phi(4) = -1/2 * (-(sin(THETA1tilde) * cos(gamma-THETA2tilde)) + (THETA1tilde * cos_gamma));
Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));

%Calcolo MATRICE B_1 (B_2 e B_3 non utilizzate) 
B_1 = (Phi(5) * B(:, :, 4)) + (Phi(4)*B(:, :, 5)) + (Phi(3)*B(:, :, 6));

%INTEGRALE di 1/r
ris1 = THETA1tilde * I(1) * eye(3, 3);

%INTEGRALE di r_i*r_j/r^3
ris2 = I(2) * B_1;

%INTEGRALE di 1/r^3
ris3 = THETA1tilde * I(4) * eye(3, 3);

%INTEGRALE di r_i*r_j/r^5
ris4 = I(5) * B_1;

return