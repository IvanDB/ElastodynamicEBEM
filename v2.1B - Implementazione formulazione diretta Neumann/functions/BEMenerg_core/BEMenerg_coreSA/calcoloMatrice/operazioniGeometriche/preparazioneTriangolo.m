function [zeta, children, d_MIN, c, sign_prod] = preparazioneTriangolo(sp, T3, vnF, Sistema_rif) 
%INPUT
% - 
% -
%
%OUTPUT
% - 

%% CALCOLO PROIEZIONE PUNTO SORGENTE NEL PIANO DEL TRIANGOLO DI CAMPO

%Calcolo PROIEZIONE PUNTO SORGENTE sul PIANO del TRIANGOLO di CAMPO

param_s = vnF * (sp - T3(1, :))';
sp_plane = sp - param_s*vnF;

%Calcolo SEGNO PRODOTTO SCALARE tra VERSORE NORMALE al piano e 
% VETTORE DISTANZA tra il PUNTO SORGENTE sp e la sua PROIEZIONE sp_plane 
sign_prod = (sp - sp_plane) * vnF';
if abs(sign_prod) <= 1.0e-6 %Pulizia numerica
    sign_prod = 0;
end
sign_prod = sign(sign_prod);

%Calcolo DISTANZA PUNTO SORGENTE - PIANO TRIANGOLO di CAMPO
zeta = abs(param_s);
if zeta < 1.0e-6 %Pulizia numerica
    zeta = 0;
end
        
%% ESTRAZIONE INFORMAZIONI NUOVO SISTEMA DI RIFERIMENTO

%Estrazione LUNGHEZZA VETTORE P2 - P1
L1 = Sistema_rif.L1;

%Estrazione VERSORE e1
e1 = Sistema_rif.e1;

%Estrazione LUNGHEZZA VETTORE P3 - P1
L3 = Sistema_rif.L3;

%Estrazione SENO e COSENO ANGOLO tra VERSORI e1 e e3
cos_delta = Sistema_rif.cos_e1_e3;
sin_delta = Sistema_rif.sin_e1_e3;

%Estrazione VERSORE e2
e2 = Sistema_rif.e2;

%% CALCOLO INFORMAZIONI TRIANGOLO DI CAMPO NEL NUOVO SISTEMA DI RIFERIMENTO

%Calcolo VETTORE DISTANZA PUNTO P1 e l'ORIGINE sp_plane
diff = T3(1, :) - sp_plane;

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

%% CALCOLO COEFFICIENTI (non normalizzati) TRIANGOLI FIGLI 

%COEFFICIENTE PRIMO TRIANGOLO FIGLIO (= det(tri2D_1 tri2D_2]))
c(1) = (tri2D(1, 1) * tri2D(2, 2)) - (tri2D(1, 2) * tri2D(2, 1));

%COEFFICIENTE SECONDO TRIANGOLO FIGLIO (= det(tri2D_2 tri2D_3]))
c(2) = (tri2D(2, 1) * tri2D(3, 2)) - (tri2D(2, 2) * tri2D(3, 1));

%COEFFICIENTE TERZO TRIANGOLO FIGLIO (= det(tri2D_3 tri2D_1]))
c(3) = (tri2D(3, 1) * tri2D(1, 2)) - (tri2D(3, 2) * tri2D(1, 1));

%% CALCOLO DISTANZA TRIANGOLO CAMPO da PROIEZIONE NODO (origine nuovo SdR)

%Calcolo valore DOPPIO AREA CON SEGNO TRIANGOLO di CAMPO CORRENTE
A_doppia = (tri2D(2, 1) - tri2D(1, 1)) * (tri2D(3, 2) - tri2D(1, 2));
            %- (tri2D(2, 2) - tri2D(1, 2)) * (tri2D(3, 1) - tri2D(1, 1));
            % termine inutile (?)

%Calcolo SECONDA COORDINATA BARICENTRICA
eta_spPlane(2) = c(3) / A_doppia;

%Calcolo TERZA COORDINATA BARICENTRICA
eta_spPlane(3) = c(1) / A_doppia;

%Calcolo PRIMA COORDINATA BARICENTRICA
eta_spPlane(1) = 1 - eta_spPlane(2) - eta_spPlane(3);

if(all(-1.0e-6 <= eta_spPlane) && all(eta_spPlane <= 1)) 
        d_MIN = 0;
else 
        d_MIN = calcoloDistanzaTriangolo(tri2D);
end
return