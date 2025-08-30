function children = preparazioneTriangoloDiagonale(sP, T3, SdR) 
%INPUT
% - 
% -
%
%OUTPUT
% - 
        
%% ESTRAZIONE INFORMAZIONI NUOVO SISTEMA DI RIFERIMENTO

%Estrazione LUNGHEZZA VETTORE P2 - P1
L1 = SdR.L1;

%Estrazione VERSORE e1
e1 = SdR.e1;

%Estrazione LUNGHEZZA VETTORE P3 - P1
L3 = SdR.L3;

%Estrazione SENO e COSENO ANGOLO tra VERSORI e1 e e3
cos_delta = SdR.cos_e1_e3;
sin_delta = SdR.sin_e1_e3;

%Estrazione VERSORE e2
e2 = SdR.e2;

%% CALCOLO INFORMAZIONI TRIANGOLO DI CAMPO NEL NUOVO SISTEMA DI RIFERIMENTO

%Calcolo VETTORE DISTANZA PUNTO P1 e l'ORIGINE sp_plane
diff = T3(1, :) - sP;

%Inizializzazione matrice 3x2 tri2D contenente le coordinate dei 
% vertici del triangolo di campo nel nuovo sistema di riferimento.
tri2D = zeros(3, 2);

%Calcolo COORDINATE VERTICI TRIANGOLO di CAMPO nel NUOVO SdR
tri2D(1, 1) = diff * e1';
tri2D(1, 2) = diff * e2';

tri2D(2, :) = tri2D(1, :) + [L1 0];
tri2D(3, :) = tri2D(1, :) + L3 * [cos_delta sin_delta];

%% CALCOLO COORDINATE TRIANGOLI FIGLI

%Inizializzazione array 3D 3x3x2 contenente le coordinate dei vertici dei
% triangoli figli nel piano del triangolo di campo
children = zeros(3, 3, 2);

%Calcolo COORDINATE dei VERTICI dei TRIANGOLI FIGLI nel nuovo SdR

%Primo triangolo figlio: O - P1 - P2
children(1, 2, :) = tri2D(1, :);
children(1, 3, :) = tri2D(2, :);

%Secondo triangolo figlio: O - P2 - P3
children(2, 2, :) = tri2D(2, :);  
children(2, 3, :) = tri2D(3, :);

%Terzo triangolo figlio: O - P3 - P1
children(3, 2, :) = tri2D(3, :);
children(3, 3, :) = tri2D(1, :);

return