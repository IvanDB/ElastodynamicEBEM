function [ris1, ris2, ris3, ris4] = BEMenerg_esplicitIntg_CC(RS, rhoS, RP, rhoP, zeta, THETA1, THETA2, infoTri2D, B)
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

%Controllo in funzione di zeta
if (zeta > 1.0e-6)
    %Caso zeta > 0 (triangoli non complanari)

    %Calcolo integrali parziali in rho
    I(1) = RP - RS;
    I(2) = RP - RS + (zeta^2 * (1/RP - 1/RS));
    I(3) = log((rhoP + RP) / (rhoS + RS)) - (rhoP/RP) + (rhoS/RS);
    I(4) = -(1 / RP) + (1 / RS);
    I(5) = -(1 / RP) + (1 / RS) + (zeta^2 / 3) * ((1/RP^3) - (1 / RS^3));
    I(6) = (1 / (3*zeta^2)) * ((rhoP / RP)^3 - (rhoS / RS)^3);
    I(7) = -1/3 * ((1/RP^3) - (1/RS^3));

    %Calcolo integrali parziali in theta
    Phi(1) = 2 * sin(THETA2tilde / 2) * sin(THETA1tilde / 2);
    Phi(2) = 2 * sin(gamma - (THETA2tilde/2)) * sin(THETA1tilde/2);
    Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
    Phi(4) = -1/2 * (-sin(THETA1tilde) * cos(gamma - THETA2tilde) + THETA1tilde * cos_gamma);
    Phi(5) = 1/2 * (THETA1tilde - cos(2*gamma - THETA2tilde) * sin(THETA1tilde));

    %Calcolo MATRICI B_1, B_2 e B_3 
    B_1 = Phi(5) * B(:, :, 4) + Phi(4) * B(:, :, 5) + Phi(3) * B(:, :, 6);
    B_2 = Phi(2) * B(:, :, 2) + Phi(1) * B(:, :, 3);
    B_3 = B(:, :, 1);  

    %Calcolo INTEGRALI

    %INTEGRALE di 1/r
    ris1 = THETA1tilde * I(1) * eye(3, 3);
    
    %INTEGRALE di r_i*r_j/r^3
    ris2 = (I(2) * B_1) + (I(3) * B_2) + (THETA1tilde * I(4) * B_3);

    %INTEGRALE di 1/r^3
    ris3 = THETA1tilde * I(4) * eye(3,3);

    %INTEGRALE di r_i*r_j/r^5
    ris4 = (I(5) * B_1) + (I(6) * B_2) + (THETA1tilde * I(7) * B_3); 

else
    %Caso zeta = 0 (triangoli complanari)
    
    %Calcolo integrali parziali in rho
    I(1) = RP - RS;
    I(2) = I(1);
    I(4) = -(1/RP) + (1/RS);
    I(5) = I(4);

    %Calcolo integrali parziali in theta
    Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
    Phi(4) = -1/2 * (-(sin(THETA1tilde) * cos(gamma-THETA2tilde)) + (THETA1tilde * cos_gamma));
    Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));
    
    

    %Calcolo MATRICE B_1 (B_2 e B_3 non utilizzate) 
    B_1 = (Phi(5) * B(:, :, 4)) + (Phi(4)*B(:, :, 5)) + (Phi(3)*B(:, :, 6));

    %Calcolo INTEGRALI

    %INTEGRALE di 1/r
    ris1 = THETA1tilde * I(1) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^3
    ris2 = I(2) * B_1;

    %INTEGRALE di 1/r^3
    ris3 = THETA1tilde * I(4) * eye(3, 3);
    
    %INTEGRALE di r_i*r_j/r^5
    ris4 = I(5) * B_1;
end

return
end