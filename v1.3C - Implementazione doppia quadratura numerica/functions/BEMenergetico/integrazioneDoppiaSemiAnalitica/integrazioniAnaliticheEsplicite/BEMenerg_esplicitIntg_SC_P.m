function [ris1, ris2, ris3, ris4] = BEMenerg_esplicitIntg_SC_P(RP, rhoP, zeta, THETA1, THETA2, infoTri2D,B)
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
I(1) = RP - zeta;
I(2) = RP + (zeta^2/RP) - (2*zeta);
I(3) = log((rhoP+RP) / zeta) - (rhoP / RP);
I(4) = -(1/RP) + (1/zeta);
I(5) = -(1/RP) + (zeta^2 / (3*RP^3)) + (2 / (3*zeta));
I(6) = (1 / (3*zeta^2)) * (rhoP/RP)^3;
I(7) = -(1 / (3*RP^3)) + (1 / (3*zeta^3));

%Calcolo integrali in theta

Phi(1) = 2 * sin(THETA2tilde / 2) * sin(THETA1tilde / 2);
Phi(2) = 2 * sin(gamma - THETA2tilde / 2) * sin(THETA1tilde / 2);
Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
Phi(4) = -1/2 * (-(sin(THETA1tilde) * cos(gamma - THETA2tilde)) + (THETA1tilde*cos_gamma));
Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));

%Calcolo MATRICI B_1, B_2 e B_3 
B_1 = Phi(5) * B(:, :, 4) + Phi(4) * B(:, :, 5) + Phi(3)*B(:, :, 6);
B_2 = Phi(2) * B(:, :, 2) + Phi(1) * B(:, :, 3);
B_3 = B(:, :, 1);  

%Calcolo integrali

%Integrale di 1/r 
ris1 = THETA1tilde * I(1) * eye(3, 3);

%Integrale di r_i*r_j/r^3
ris2 = (I(2) * B_1) + (I(3) * B_2) + (THETA1tilde * I(4) * B_3);

%Integrale di 1/r^3
ris3 = THETA1tilde * I(4) * eye(3, 3);

%Integrale di r_i*r_j/r^5
ris4 = (I(5) * B_1) + (I(6) * B_2) + (THETA1tilde * I(7) * B_3);
end