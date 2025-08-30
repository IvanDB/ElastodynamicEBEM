function uXT = BEMenerg_postProcA_singlePoint(pbParam, domainMesh, density, constValues, x, t)
% INPUT
% - 
% - 
% OUTPUT
% - 

%Allocazione matrice 3x1 contente le componenti del vettore u(x, t) da
% calcolare
uXT = zeros(3, 1);

%Controllo teorico sul valore temporale
if t <= 0
    return
end
    
%Calcolo PASSO di DISCRETIZZAZIONE TEMPORALE
dt = pbParam.Tfin / pbParam.nT;

%Calcolio massimo indice nHat temporale necessario per il calcolo
nHat = ceil(t ./ dt) - 1;

%Estrazione NUMERO TRIANGOLI DISCRETIZZAZIONE SPAZIALE
numTriangles = domainMesh.numberTriangles;

%Inizializzazione prima differenza temporale (tra t e t0 = 0)
diffT = t; %thk0 = t - t0 = t
    
%Allocazione matrice 3 x 3*N_triangles costituita numTriangles
% blocchetti 3x3
matrixTnext = zeros(3, 3*numTriangles);

%Calcolo elementi della matrice matrix0
matrixTcurr = BEMenerg_postProcA_calcMatrix(pbParam, domainMesh, constValues, x, diffT);

%Ciclo sull'indice che individua i sottointervalli temporali
for indTemp = 1 : nHat
        %Calcolo differenza temporale tra t e indTemp*dt
        diffT = t - (indTemp * dt);

        %Estrazione valori della densità nel sottointervallo temporale  
        % di indice indTemp
        alpha = density(:, indTemp);

        %Calcolo elementi matrice relativa al successivo
        % sottointervallo temporale
        matrixTnext = BEMenerg_postProcA_calcMatrix(pbParam, domainMesh, constValues, x, diffT);
         
        %Aggiornamento vettore soluzione 
        uXT = uXT + ((matrixTcurr - matrixTnext) * alpha);

        %Salvataggio matrice intervallo precedente
        matrixTcurr = matrixTnext;
end
    
%Estrazione valori della densità nel sottointervallo temporale  
% di indice n_hat
alpha = density(:, nHat + 1);

%Aggiornamento finale vettore soluzione
uXT = uXT + (matrixTcurr * alpha);

return