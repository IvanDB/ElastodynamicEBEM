function [ris5, ris6] = BEMenerg_analiticIntg_SC_SP(RS, rhoS, zeta, THETA1, THETA2, infoTri2D, B)
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

%Controllo in funzione di zeta
if (zeta > 1.0e-6)
    %Caso zeta > 0 (triangoli non complanari)

    %Calcolo integrali in rho
    I(1) = RS - zeta;
    I(2) = RS + (zeta^2 / RS) - (2*zeta);
    I(3) = log((rhoS+RS) / zeta) - (rhoS / RS);
    I(4) = -(1 / RS) + (1 / zeta);  

    %Calcolo integrali in theta
    Phi(1) = 2 * sin(THETA2tilde / 2) * sin(THETA1tilde / 2);
    Phi(2) = 2 * sin(gamma - (THETA2tilde / 2)) * sin(THETA1tilde / 2);
    Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
    Phi(4) = -1/2 * (-(sin(THETA1tilde) * cos(gamma - THETA2tilde)) + (THETA1tilde * cos_gamma));
    Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));
    
    %Calcolo MATRICI B_1, B_2 e B_3
    B_1 = Phi(5) * B(:, :, 4) + Phi(4) * B(:, :, 5) + Phi(3) * B(:, :, 6);
    B_2 = Phi(2) * B(:, :, 2) + Phi(1) * B(:, :, 3);
    B_3 = B(:, :, 1);

    %Calcolo INTEGRALI

    %Integrale di 1/r
    ris5 = THETA1tilde * I(1) * eye(3,3);
    
    %Integrale di r_i*r_j/r^3
    ris6 = (I(2) * B_1) + (I(3) * B_2) + (THETA1tilde * I(4) * B_3);    
else
    %Caso zeta = 0 (triangoli complanari)
    
    %Calcolo integrali in rho
    I(1) = RS;
    I(2) = RS;
    
    %Calcolo integrali in theta
    Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
    Phi(4) = -1/2 * ((THETA1tilde * cos_gamma) - (sin(THETA1tilde) * cos(gamma-THETA2tilde)));
    Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));

    %STEP 4: CALCOLO della MATRICE B_1 nel caso zeta=0

    %Calcolo MATRICE B_1 (B_2 e B_3 non utilizzate) 
    B_1 = Phi(5) * B(:, :, 4) + Phi(4) * B(:, :, 5) + Phi(3) * B(:, :, 6);

    %Calcolo INTEGRALI

    %INTEGRALE di 1/r
    ris5 = THETA1tilde * I(1) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^3
    ris6 = I(2) * B_1;
end

return
end