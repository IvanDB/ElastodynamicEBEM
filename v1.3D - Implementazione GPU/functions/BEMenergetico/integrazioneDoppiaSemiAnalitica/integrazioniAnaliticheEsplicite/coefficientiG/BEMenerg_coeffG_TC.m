function G = BEMenerg_coeffG_TC(zeta, RS, rhoS, THETA1, THETA2, infoTri2D)
% INPUT
% - 
% - 
% OUTPUT
% - 

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

%Calcolo integrali Phi
Phi(1) = 2 * sin(THETA2tilde / 2) * sin(THETA1tilde / 2);
Phi(2) = 2 * sin(gamma - (THETA2tilde/2)) * sin(THETA1tilde / 2);
Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
Phi(4) = -1/2 * (-(sin(THETA1tilde) * cos(gamma - THETA2tilde)) + (THETA1tilde * cos_gamma));
Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));

%% CALCOLO COEFFICIENTI

%Controllo in funzione di zeta

%Caso z = 0 (triangoli complanari)

%Calcolo coefficienti g_i (vedi Milroy)
g_3 = 1/2 * (log((1 + cos(THETA1 + beta)) * (1 - cos(THETA2 + beta))) ...
            - log((1 - cos(THETA1 + beta)) * (1 + cos(THETA2 + beta))));

% Calcolo COEFFICIENTI G_i (vedi Milroy)
G(1) = (F * g_3) - (RS * THETA1tilde);

G(4) = 1 / F * (cos(THETA2 + beta) ...
                - cos(THETA1 + beta)) ...
       + 1 / RS * THETA1tilde;

G(7) = -RS * Phi(5) + F * (cos(alpha - gamma + THETA2) ...
                           - cos(alpha - gamma + THETA1) ...
                           + sin_alpha^2 * g_3);

G(8) = -RS * Phi(4) + F * (cos(THETA1 + alpha) ...
                           - cos(THETA2 + alpha) ...
                           - sin_alpha * sin_beta * g_3);

G(9) = -RS * Phi(3) + F * (cos(beta - THETA1) ...
                           - cos(beta - THETA2) ...
                           + sin_beta^2 * g_3);

G(17) = 1 / (3*F) * (((1 + cos_alpha^2) * cos(THETA2 + beta) ...
                       - (sin(THETA2 + beta)^2) * cos(alpha - gamma + THETA2)) ...
                    -((1 + cos_alpha^2) * cos(THETA1 + beta) ...
                       - (sin(THETA1 + beta)^2) * cos(alpha - gamma + THETA1))) ...
        + 1 / RS * Phi(5);
  
G(18) = 1 / (3*F) * ((cos(THETA2 + beta) * ((cos_alpha * cos_beta) + cos(alpha + beta)) ...
                        + (sin(THETA2 + beta)^2) * cos(alpha + THETA2)) ...
                    -(cos(THETA1 + beta) * ((cos_alpha * cos_beta) + cos(alpha + beta)) ...
                        +(sin(THETA1 + beta)^2) * cos(alpha + THETA1))) ...
       + 1 / RS * Phi(4);
     
G(19) = 1 / (3*F) * (((1 + cos_beta^2) * cos(THETA2 + beta) ...
                        + (sin(THETA2 + beta)^2) * cos(THETA2 - beta)) ...
                    -((1 + cos_beta^2) * cos(THETA1 + beta) ...
                        + (sin(THETA1 + beta)^2) * cos(beta - THETA1)))...
       + 1 / RS * Phi(3);

return
end