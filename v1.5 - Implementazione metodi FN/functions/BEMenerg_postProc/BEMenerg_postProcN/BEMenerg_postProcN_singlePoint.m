function uXT = BEMenerg_postProcN_singlePoint(pbParam, domainMesh, density, methodInfo, constValues, x, t)
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
%%
%Calcolo PASSO di DISCRETIZZAZIONE TEMPORALE
deltaT = pbParam.Tfin / pbParam.nT;

%Calcolio massimo indice nHat temporale necessario per il calcolo
nHat = ceil(t ./ deltaT) - 1;

%Estrazione NUMERO TRIANGOLI DISCRETIZZAZIONE SPAZIALE
numTriang = domainMesh.numberTriangles;

%Inizializzazione prima differenza temporale (tra t e t0 = 0)
diffT = t; %thk0 = t - t0 = t
    
%Calcolo elementi della matrice matrix0
matrixTcurr = BEMenerg_postProcN_calcMatrix(pbParam, domainMesh, methodInfo, constValues, x, diffT);

%Ciclo sull'indice che individua i sottointervalli temporali
for indTemp = 1 : nHat
    %Calcolo differenza temporale tra t e indTemp*dt
    diffT = t - (indTemp * deltaT);

    %Estrazione valori della densità nel sottointervallo temporale  
    % di indice indTemp
    alpha = density(:, indTemp);

    %Calcolo elementi matrice relativa al successivo
    % sottointervallo temporale
    matrixTnext = BEMenerg_postProcN_calcMatrix(pbParam, domainMesh, methodInfo, constValues, x, diffT);
     
    %Aggiornamento vettore soluzione 
    uXT = uXT + ((matrixTcurr - matrixTnext) * alpha);

    %Salvataggio matrice intervallo precedente
    matrixTcurr = matrixTnext;
end
    
%Estrazione valori della densità nel sottointervallo temporale  
% di indice nHat
alpha = density(:, nHat + 1);

%Aggiornamento finale vettore soluzione
uXT = uXT + (matrixTcurr * alpha);
return