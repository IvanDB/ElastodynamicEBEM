function B = BEMenerg_coreGPUdiag_coeffB(infoTri2D)
% INPUT
% - 
% - 
% OUTPUT
% - B: Array di 6 matrici 3x3 contenenti i coefficienti da applicare agli
%     integrali delle funzioni T(theta) e U(theta) nel passaggi di integrazione analitica 

%% ESTRAZIONE INFORMAZIONE TRIANGOLO

%Estrazione SENO e COSENO ANGOLO gamma 
cos_gamma = infoTri2D.cos_gamma;
sin_gamma = infoTri2D.sin_gamma;

%% CALCOLO della MATRICE 3x3 di rotazione

%Inizializzazione matrice 3x3
A = zeros(3, 3);

%Calcolo prima riga
A(1, 1) = 1;

%Calcolo seconda riga
A(2, 1 : 2) = [cos_gamma sin_gamma];

%Calcolo terza riga
A(3, 3) = 0;

%% CALCOLO MATRICI B_i

%Matrice B_1
B(:, :, 1) = zeros(3);
%Matrice B_2
B(:, :, 2) = zeros(3);
%Matrice B_3
B(:, :, 3) = zeros(3);
%Matrice B_4
B(:, :, 4) = (A(1, :)' * A(1, :)) / sin_gamma^2;
%Matrice B_5
B(:, :, 5) = ((A(1, :)' * A(2, :)) + (A(2, :)' * A(1, :))) / sin_gamma^2;
%Matrice B_6
B(:, :, 6) = (A(2, :)' * A(2, :)) / sin_gamma^2;

return
end
