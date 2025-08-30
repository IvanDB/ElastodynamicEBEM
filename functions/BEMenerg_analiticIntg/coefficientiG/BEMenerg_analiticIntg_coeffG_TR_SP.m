function G = BEMenerg_analiticIntg_coeffG_TR_SP(zeta, THETA1, THETA2, infoTri2D)
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
if(zeta > 1.0e-6)
    %Caso z > 0 (triangoli non complanari)
    
    %Calcolo parametri e coefficienti intemedi
    R1 = F / sin(THETA1 + beta);
    R2 = F / sin(THETA2 + beta);
    
    FRAC_zeta_R1 = zeta / R1;
    FRAC_zeta_R2 = zeta / R2;
    ADD_zexp2_Fexp2 = zeta^2 + F^2;

    %Calcolo integrali Phi
    Phi(3) = 1/2 * (THETA1tilde - (cos(THETA2tilde) * sin(THETA1tilde)));
    Phi(4) = -1/2 * (-(sin(THETA1tilde) * cos(gamma - THETA2tilde)) + (THETA1tilde * cos_gamma));
    Phi(5) = 1/2 * (THETA1tilde - (cos(2*gamma - THETA2tilde) * sin(THETA1tilde)));
    
    %Calcolo coefficienti g_i (vedi Milroy)
    g(1) = sqrt(1 + FRAC_zeta_R1^2);
    g(2) = sqrt(1 + FRAC_zeta_R2^2);
    g(3) = 1/2 * (log((g(1) + cos(THETA1 + beta)) * (g(2) - cos(THETA2 + beta))) ...
                 - log((g(1) - cos(THETA1 + beta)) * (g(2) + cos(THETA2 + beta))));
    g(4) = log((1 / FRAC_zeta_R1) + sqrt(1 + (1 / FRAC_zeta_R1)^2));
    g(5) = log((1 / FRAC_zeta_R2) + sqrt(1 + (1 / FRAC_zeta_R2)^2));
    g(7) = pi - acos((-zeta * cos(THETA2 + beta)) / sqrt(ADD_zexp2_Fexp2)) ...
                -acos((zeta * cos(THETA1 + beta)) / sqrt(ADD_zexp2_Fexp2));

    % Calcolo COEFFICIENTI G_i (vedi Milroy)
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
else
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
end
return
end