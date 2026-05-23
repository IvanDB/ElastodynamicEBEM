function [ris5, ris6] = BEMenerg_analiticIntg_TR_SP(zeta, THETA1, THETA2, infoTri2D, B)
% INPUT
% - 
% - 
% OUTPUT
% - 

%% CALCOLO VALORI SINGOLI NUCLEI INTEGRALI 

%Controllo in funzione di zeta
if (zeta > 1.0e-6)
    %Caso zeta > 0 (triangoli non complanari)

    %Calcolo coeffincienti G_i
    G = BEMenerg_analiticIntg_coeffG_TR_SP(zeta, THETA1, THETA2, infoTri2D);

    %Calcolo INTEGRALI

    %INTEGRALE di 1/r
    ris5 = G(1) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^3
    ris6 = G(4) * B(:, :, 1) + G(5) * B(:, :, 2) + G(6) * B(:, :, 3) ...
            +G(7) * B(:, :, 4) + G(8) * B(:, :, 5) + G(9) * B(:, :, 6);
else
    %Caso zeta = 0 (triangoli complanari)

    %Calcolo coeffincienti G_i
    G = BEMenerg_analiticIntg_coeffG_TR_SP(0, THETA1, THETA2, infoTri2D);
    
    %Calcolo INTEGRALI

    %INTEGRALE di 1/r
    ris5 = G(1) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^3 
    ris6 = G(7) * B(:, :, 4) + G(8) * B(:, :, 5) + G(9) * B(:, :, 6); 
end

return
end