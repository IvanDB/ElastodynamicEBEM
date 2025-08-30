function [matrix, isnotzero] = BEMenerg_coreMXG2D_calcMatrixBlock(matrix, methodInfo, indTemp, deltaT, pbParam, domainMesh, constValues)
%INPUT
%   - matrix: matrice nulla su cui salvare il risultato
%   - methodInfo: struct contenente le informazioni sul metodo
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


%% SETUP VARIABILI in ARRAY di CELLE

%Estrazione NUMERO di TRIANGOLI della MESH
numT = domainMesh.numberTriangles;

%Trasformazione matrice matrix in array di N_triangles x N_triangles celle
% ciascuna di dimensione 3x3
matrix = mat2cell(matrix, 3*ones(1, numT), 3*ones(1, numT));

%Calcolo dimensioni array di celle matrix
sz = size(matrix);

%Allocazione cell-array contenente i valori di isok
cell_isnotzero = cell(numT^2, 1);

%Inizializzazione cell-array contenente i valori di isok
cell_isnotzero(1 : numT^2, 1) = {0};

%% PARALLELIZZAZIONE DOPPIO CICLO SUI TRIANGOLI SORGENTE E DI CAMPO

%Procedimento element-by-element
parfor ind = 1 : numT^2

    %Calcolo indici di riga e colonna relativi all'indice lineare ind
    [indS, indF] = ind2sub(sz, ind);

    %Estrapolazione COORDINATE BARICENTRO triangoli sorgente e campo
    %   correnti
    cS = domainMesh.center(indS, :);
    cF = domainMesh.center(indF, :);

    %Estrapolazione MASSIMA LUNGHEZZA LATI triangoli sorgente e campo
    %   correnti
    maxS = domainMesh.maxL(indS);
    maxF = domainMesh.maxL(indF);
    
    %Calcolo VETTORE DISTANZA BARICENTRI triangoli campo e sorgente
    vettDist = cF - cS;
    
    %Stima DISTANZE MINIMA e MASSIMA triangoli campo e sorgente
    distMin = sqrt(sum(vettDist.^2)) - maxF - maxS;
    distMax = sqrt(sum(vettDist.^2)) + maxF + maxS;
            
    %Controllo condizione teorica sull'elemento matriciale 3x3 corrente
    if(((indTemp-1) * pbParam.velS * deltaT < distMax) && ((indTemp + 1) * pbParam.velP * deltaT >= distMin))
        %Aggiornamento del parametro iok locale
        
        cell_isnotzero(ind, 1) = {1};
        %Estrapolazione COORDINATE VERTICI triangolo campo corrente
        nodeF = domainMesh.triangles(indF, 1:3);
        TF = domainMesh.coordinates(nodeF, :);
        
        %Estrapolazione VERSORE NORMALE all'elemento corrente (triangolo di campo)
        vnF = domainMesh.normal(indF, :);

        if(indS == indF)
            matrixSubBlock = BEMenerg_coreMXdiag_calcSubBlock(pbParam, methodInfo, TF, vnF, indTemp, deltaT, constValues{indS})
        else
            matrixSubBlock = BEMenerg_coreMXG2D_calcMatrixSubBlock(pbParam, methodInfo, indTemp, deltaT, ...
                                                                constValues{indS}.GHnodes, constValues{indS}.GHweights, ...
                                                                constValues{indF}.G2Dnodes, constValues{indF}.G2Dweights);
        end
                
        %Pulizia "SPORCIZIA NUMERICA" dal blocchetto matriciale corrente
        matrixSubBlock(abs(matrixSubBlock) < 1.0e-14) = 0;
        
        %Posizionamento blocco matriciale corrente in matrix
        matrix{ind} = matrixSubBlock;
     end 
end

%% FORMATTAZIONE RISULTATI

%Trasformazione da array di celle a matrice standard
matrix = cell2mat(matrix);

%Calcolo PARAMETRO iok
isnotzero = max(cell2mat(cell_isnotzero));
return