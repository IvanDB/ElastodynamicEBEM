function [ris1, ris2, ris3, ris4] = BEMenerg_esplicitIntg_TC(zeta, RS, rhoS, THETA1, THETA2, infoTri2D, B)
% INPUT
% - 
% - 
% OUTPUT
% - 

%Controllo in funzione di zeta
if (zeta > 1.0e-6)
    %Caso zeta > 0 (triangoli non complanari)

    %Calcolo coefficienti G_i necessari per integrazione analitica
    G = BEMenerg_coeffG_TC(zeta, RS, rhoS, THETA1, THETA2, infoTri2D);
    
    %Calcolo INTEGRALI

    %INTEGRALE di 1/r
    ris1 = G(1) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^3
    ris2 = G(4) * B(:, :, 1) + G(5) * B(:, :, 2) + G(6) * B(:, :, 3) ...
            + G(7) * B(:, :, 4) + G(8) * B(:, :, 5) + G(9)*B(:, :, 6); 

    %INTEGRALE di 1/r^3 
    ris3 = G(4) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^5
    ris4 = G(14) * B(:, :, 1) + G(15)*B(:, :, 2) + G(16)*B(:, :, 3) ...
            + G(17) * B(:, :,4) + G(18) * B(:, :, 5) + G(19) * B(:, :, 6); 

else
    %Caso zeta = 0 (triangoli complanari)
    G = BEMenerg_coeffG_TC(0, RS, rhoS, THETA1, THETA2, infoTri2D);
    
    %Calcolo INTEGRALI

    %INTEGRALE di 1/r
    ris1 = G(1) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^3
    ris2 = G(7) * B(:, :, 4) + G(8) * B(:, :, 5) + G(9) * B(:, :, 6);

    %INTEGRALE di 1/r^3
    ris3 = G(4) * eye(3, 3);

    %INTEGRALE di r_i*r_j/r^5
    ris4 = G(17) * B(:, :, 4) + G(18) * B(:, :, 5) + G(19) * B(:, :, 6);

end
return
end