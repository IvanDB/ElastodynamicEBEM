function [matrix, isnotzero] = BEMenerg_calcMatrixBlock_FN(matrix, hk, dt, pb_param, domainMesh, gha, ghw, dgn, dgw, ghcn, ghcw, typeQuad)
%INPUT
%   - matrix: matrice nulla su cui salvare il risultato
%   - hk: intero contente l'indice della matrice da calcolare
%   - dt: double contente l'ampiezza della discetizzazione temporare
%   - pb_param: struct contenente i parametri del problema
%   - domainMesh: struct contente le informazioni sulla mesh
%   - gha: matrice contenente i nodi di quadratura di Gauss-Hammer
%   - ghw: matrice contenente i pesi di quadratura di Gauss-Hammer
%
% OUTPUT
%   - matrix: matrice contenente il risultato
%   - isnotzero: flag segnalante se il risultato è la matrice nulla

%% SETUP VARIABILI come ARRAY di CELLE
%NUMERO di TRIANGOLI della MESH
nTriangles = domainMesh.number_triangles;

%Trasformazione matrice matrix in array di N_triangles x N_triangles celle
% ciascuna di dimensione 3x3
matrix = mat2cell(matrix, 3*ones(1, nTriangles), 3*ones(1, nTriangles));

%Calcolo dimensioni array di celle matrix
sz = size(matrix);

%Allocazione cell-array contenente i valori di isok
cell_isok = cell(nTriangles^2, 1);

%Inizializzazione cell-array contenente i valori di isok
cell_isok(1:nTriangles^2, 1) = {0};

%% PARALLELIZZAZIONE DOPPIO CICLO SUI TRIANGOLI SORGENTE E DI CAMPO

%Procedimento element-by-element
parfor ind = 1 : nTriangles^2

    %Calcolo indici di riga e colonna relativi all'indice lineare ind
    [indS, indF] = ind2sub(sz, ind);
         
    %Estrapolazione COORDINATE BARICENTRO triangolo campo corrente
    cS = domainMesh.center(indS, :);
    
    %Estrapolazione MASSIMA LUNGHEZZA LATI triangolo sorgente corrente
    maxS = domainMesh.maxL(indS);

    %Estrapolazione COORDINATE BARICENTRO triangolo campo corrente
    cF = domainMesh.center(indF, :);
    
    %Estrapolazione MASSIMA LUNGHEZZA LATI triangolo campo corrente
    maxF = domainMesh.maxL(indF);
    
    %Calcolo VETTORE DISTANZA BARICENTRI triangoli di campo e sorgente
    vD = cF - cS;
    
    %Calcolo DISTANZA MINIMA BARICENTRI triangoli di campo e sorgente (?)
    distmin = sqrt(sum(vD.^2)) - maxF - maxS;
    
    %Calcolo DISTANZA MASSIMA BARICENTRI triangoli di campo e sorgente (?)
    distmax = sqrt(sum(vD.^2)) + maxF + maxS;
    
    %Controllo condizione teorica sull'elemento matriciale 3x3 corrente
    if(((hk-1)*pb_param.velS*dt < distmax) && ((hk+1)*pb_param.velP*dt >= distmin))
        %Aggiornamento del parametro iok locale
        cell_isok(ind, 1) = {1};
        
        %Controllo caso specifico diagonale
        if(indS == indF)
            matrix{ind} = BEMenerg_calcMatrxSubBlockDiag_SA(indS, pb_param, domainMesh, hk, dt, gha, ghw);
        else
            switch typeQuad
                case "DGH"
                    matrix{ind} = BEMenerg_calcMatrixSubBlock_DGH(indS, indF, pb_param, domainMesh, hk, dt, gha, ghw);
                case "DQM"
                    matrix{ind} = BEMenerg_calcMatrixSubBlock_DQM(indS, indF, pb_param, domainMesh, hk, dt, gha, ghw, dgn, dgw);
                case "DQC"
                    matrix{ind} = BEMenerg_calcMatrixSubBlock_DQC(indS, indF, pb_param, domainMesh, hk, dt, gha, ghw, ghcn, ghcw);
                case "DQC3"
                    matrix{ind} = BEMenerg_calcMatrixSubBlock_DQC(indS, indF, pb_param, domainMesh, hk, dt, gha, ghw, ghcn, ghcw);
            end
        end
    end
end

%% FORMATTAZIONE RISULTATI

%Trasformazione da array di celle a matrice standard
matrix = cell2mat(matrix);

%Pulizia numerica
matrix(abs(matrix) < 1.0e-14) = 0;

%Calcolo PARAMETRO iok
isnotzero = max(cell2mat(cell_isok));
return