function G = BEMenerg_coeffG_TR_P(zeta, THETA1, THETA2, infoTri2D)
% INPUT
% - 
% - 
% OUTPUT
% - ... Milroy


%% CALCOLO COSTANTI ANGOLARI RICORRENTI ed ESTRAZIONE INFORMAZIONI GEOMETRICHE TRIANGOLO

%Estrazione LUNGHEZZA primo lato
a = infoTri2D.a;

%Estrazione SENO ANGOLI
sin_beta = infoTri2D.sin_beta;
sin_alpha = infoTri2D.sin_alpha;

%Estrazione COSENO ANGOLI
cos_alpha = infoTri2D.cos_alpha;
cos_beta  = infoTri2D.cos_beta;
cos_gamma = infoTri2D.cos_gamma;

%Estrazione AMPIEZZA ANGOLI
alpha = infoTri2D.alpha;
beta  = infoTri2D.beta;
gamma = infoTri2D.gamma;

%Calcolo delle COSTANTI RICORRENTI negli INTEGRALI in THETA
THETA1tilde = THETA2 - THETA1;
THETA2tilde = THETA2 + THETA1;

%Calcolo delle COSTANTI RICORRENTI negli INTEGRALI in RHO (vedi Milroy)
F = a * sin_beta;

R1 = F / sin(THETA1 + beta);
R2 = F / sin(THETA2 + beta);

FRAC_zeta_R1 = zeta / R1;
FRAC_zeta_R2 = zeta / R2;
ADD_zexp2_Fexp2 = zeta^2 + F^2;

%% CALCOLO COEFFICIENTI INTERMEDI

%Calcolo integrali Phi
Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
Phi(4) = -1/2 * (-(sin(THETA1tilde) * cos(gamma - THETA2tilde)) + (THETA1tilde * cos_gamma));
Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));

%Calcolo dei COEFFICIENTI g_i (vedi Milroy)
g(1) = sqrt(1 + FRAC_zeta_R1^2);
g(2) = sqrt(1 + FRAC_zeta_R2^2);
g(3) = 1/2 * (log((g(1) + cos(THETA1 + beta)) * (g(2) - cos(THETA2 + beta))) ...
             - log((g(1) - cos(THETA1 + beta)) * (g(2) + cos(THETA2 + beta))));
g(4) = log((1 / FRAC_zeta_R1) + sqrt(1 + (1 / FRAC_zeta_R1)^2));
g(5) = log((1 / FRAC_zeta_R2) + sqrt(1 + (1 / FRAC_zeta_R2)^2));
g(7) = pi - acos((-zeta * cos(THETA2 + beta)) / sqrt(ADD_zexp2_Fexp2)) ...
            -acos((zeta * cos(THETA1 + beta)) / sqrt(ADD_zexp2_Fexp2));
g(8) = F^2 / ADD_zexp2_Fexp2;

%% CALCOLO COEFFICIENTI G_i

%Calcolo degli INTEGRALI ANALITICI G_i (vedi Milroy)
G(1) = (zeta * (-THETA1tilde + g(7))) + (F * g(3));

G(4) = 1 / zeta * (THETA1tilde - g(7));

G(5) = cos(gamma - THETA2) * g(5) ...
        -cos(gamma - THETA1) * g(4) ...
        - cos_alpha * g(3);

G(6) = - cos(THETA2) * g(5) ...
        + cos(THETA1) * g(4) ...
        - cos_beta * g(3);

G(7) = zeta * (g(7) - 2*Phi(5)) + F*(- g(1) * cos(alpha - gamma + THETA1) ...
                                     + g(2) * cos(alpha - gamma + THETA2) ...
                                     + g(3) * sin_alpha^2);

G(8) = -zeta * ((cos_gamma * g(7)) + 2*Phi(4)) ...
        + F * (- sin_alpha * sin_beta * g(3) ...
               + g(1) * cos(THETA1 + alpha) ...
               - g(2) * cos(THETA2 + alpha));

G(9) = zeta * (g(7) - 2*Phi(3)) ...
        + F * (sin_beta^2 * g(3) ...
               + g(1) * cos(THETA1 - beta) ...
               - g(2) * cos(THETA2 - beta));

G(14) = 1 / (3*zeta^2) * (1 / zeta * (THETA1tilde - g(7)) ...
        + (F / ADD_zexp2_Fexp2 * (cos(THETA1 + beta) / g(1) ...
                                   - cos(THETA2 + beta) / g(2))));

G(15) = 1 / (3*zeta^2) * ((sin_alpha * sin(THETA2 + beta) ...
                            - g(8) * cos_alpha * cos(THETA2 + beta)) / g(2) ...
                        + (- sin_alpha * sin(THETA1 + beta) ...
                            + g(8) * cos_alpha * cos(THETA1 + beta)) / g(1));

G(16) = 1 / (3*zeta^2) * ((sin_beta * sin(THETA1+beta) ...
                            + g(8) * cos_beta * cos(THETA1 + beta)) / g(1) ...
                        + (- sin_beta * sin(THETA2 + beta) ...
                            - g(8) * cos_beta * cos(THETA2 + beta)) / g(2));

G(17) = (1 / (3*F) * ((- (sin(THETA2 + beta))^2 * cos(alpha - gamma + THETA2) ...
                        + g(8) * (cos_alpha)^2 * cos(THETA2 + beta)) / g(2) ...
                     -(-(sin(THETA1 + beta))^2 * cos(alpha - gamma + THETA1) ...
                        + g(8) * (cos_alpha)^2 * cos(THETA1 + beta)) / g(1))) ...
        - (1 / (3*zeta) * (g(7) - 2*Phi(5)));

G(18) = (1 / (3*F) * (((sin(THETA2 + beta))^2 * cos(alpha + THETA2) ...
                        + g(8) * (cos_alpha * cos_beta * cos(THETA2 + beta))) / g(2) ...
                     -((sin(THETA1 + beta))^2 * cos(alpha + THETA1) ...
                        + g(8) * (cos_alpha * cos_beta * cos(THETA1 + beta))) / g(1))) ...
        + (1 / (3*zeta) * ((g(7) * cos_gamma) + 2*Phi(4)));
     
G(19) = (1 / (3*F) * (((sin(THETA2 + beta))^2 * cos(beta - THETA2) ...
                        + g(8) * (cos_beta^2 * cos(THETA2 + beta))) / g(2) ...
                    -((sin(THETA1 + beta))^2 * cos(beta - THETA1) ...
                        + g(8) * (cos_beta^2 * cos(THETA1 + beta))) / g(1))) ...
         - (1 / (3*zeta) * (g(7) - 2*Phi(3)));  
     
end