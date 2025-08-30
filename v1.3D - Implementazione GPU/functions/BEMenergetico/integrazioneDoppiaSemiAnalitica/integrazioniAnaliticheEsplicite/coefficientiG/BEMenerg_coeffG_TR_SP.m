function G = BEMenerg_coeffG_TR_SP(zeta, THETA1, THETA2, infoTri2D)
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

%% CALCOLO COEFFICIENTI

%Controllo in funzione di zeta

%Caso z = 0 (triangoli complanari)

%Calcolo coefficienti g_i (vedi Milroy)
g_3 = 1/2 * (log((1 + cos(THETA1 + beta)) * (1 - cos(THETA2 + beta))) ...
            - log((1 - cos(THETA1 + beta)) * (1 + cos(THETA2 + beta))));

% Calcolo COEFFICIENTI G_i (vedi Milroy)
G(1) = F * g_3;

G(7) = F * (- cos(alpha - gamma + THETA1) + cos(alpha - gamma + THETA2) ...
            + sin_alpha^2 * g_3);

G(8) = F * (- (sin_alpha * sin_beta * g_3) ...
            + cos(THETA1 + alpha) ...
            - cos(THETA2 + alpha));

G(9) = F * (sin_beta^2 * g_3 ...
            + cos(THETA1 - beta) ...
            - cos(THETA2 - beta));

return
end