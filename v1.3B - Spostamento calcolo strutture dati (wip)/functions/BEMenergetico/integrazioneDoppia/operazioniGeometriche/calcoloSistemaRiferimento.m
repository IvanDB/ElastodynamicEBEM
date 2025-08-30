function sisRif = calcoloSistemaRiferimento(T3, vnF)

% INPUT
%
%
% OUTPUT
% - 

%% CALCOLO SISTEMA RIFERIMENTO

%Calcolo vettore P2 - P1 (lato P1P2 triangolo di campo corrente)
vl1 = T3(2, :) - T3(1, :);  

%Calcolo lunghezza del vettore P2 - P1
sisRif.L1 = sqrt(sum(vl1.^2));

%Calcolo versore e1 (parallelo vettore P2 - P1)
sisRif.e1 = vl1/sisRif.L1;

%----------------------------------------------------------------------
%Calcolo vettore P3 - P1 (lato P1P3 triangolo di campo corrente)
vl3 = T3(3, :) - T3(1, :);  

%Calcolo lunghezza del vettore P3 - P1
sisRif.L3 = sqrt(sum(vl3.^2));

%Calcolo versore e3 (parallelo vettore P2 - P1)
sisRif.e3 = vl3/sisRif.L3;

%----------------------------------------------------------------------
%Calcolo COSENO ANGOLO compreso tra versori e1 e e3
sisRif.cos_e1_e3 = sisRif.e3 * sisRif.e1';
%Calcolo ANGOLO compreso tra versori e1 e e3
sisRif.delta = acos(sisRif.cos_e1_e3);
%Calcolo SENO ANGOLO compreso tra versori e1 e e3
sisRif.sin_e1_e3 = sin(sisRif.delta);

%----------------------------------------------------------------------
%Calcolo versore e2 (perpendicolare ad e1 nel piano nel piano)
% dato da PRODOTTO VETTORIALE tra VERSORE NORMALE vn e VERSORE e_1 
sisRif.e2 = cross(vnF, sisRif.e1);
sisRif.e2 = sisRif.e2 / norm(sisRif.e2, 2);

return

