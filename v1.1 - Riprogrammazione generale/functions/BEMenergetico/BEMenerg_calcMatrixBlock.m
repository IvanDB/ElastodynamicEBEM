function [matrix, isnotzero] = BEMenerg_calcMatrixBlock(matrix, hk, dt, pb_param, domainMesh, gha, ghw)
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


%% SETUP VARIABILI in ARRAY di CELLE

%Estrazione NUMERO di TRIANGOLI della MESH
N_triangles = domainMesh.number_triangles;

%Trasformazione matrice matrix in array di N_triangles x N_triangles celle
% ciascuna di dimensione 3x3
matrix = mat2cell(matrix, 3*ones(1, N_triangles), 3*ones(1, N_triangles));

%Calcolo dimensioni array di celle matrix
sz = size(matrix);

%Allocazione cell-array contenente i valori di isok
cell_isok = cell(N_triangles^2, 1);

%Inizializzazione cell-array contenente i valori di isok
cell_isok(1:N_triangles^2, 1) = {0};

%% PARALLELIZZAZIONE DOPPIO CICLO SUI TRIANGOLI SORGENTE E DI CAMPO

%Procedimento element-by-element
parfor ind = 1 : N_triangles^2

    %Calcolo indici di riga e colonna relativi all'indice lineare ind
    [indS, indF] = ind2sub(sz, ind);

    %Estrapolazione INDICI VERTICI triangolo sorgente corrente
    nodeS = domainMesh.triangles(indS, 1:3);
    
    %Estrapolazione COORDINATE VERTICI triangolo sorgente corrente
    TS = domainMesh.coordinates(nodeS, :);
    
    %Estrapolazione COORDINATE BARICENTRO triangolo sorgente corrente
    cS = domainMesh.center(indS, :);
    
    %Estrapolazione AREA triangolo sorgente corrente
    areaS = domainMesh.area(indS);
    
    %Estrapolazione MASSIMA LUNGHEZZA LATI triangolo sorgente corrente
    maxS = domainMesh.maxL(indS);

    %Estrapolazione COORDINATE BARICENTRO triangolo campo corrente
    cF = domainMesh.center(indF, :);
    
    %Estrapolazione MASSIMA LUNGHEZZA LATI triangolo campo corrente
    maxF = domainMesh.maxL(indF);
    
    %Calcolo VETTORE DISTANZA BARICENTRI triangoli di campo e sorgente
    cF = cF - cS;
    
    %Calcolo DISTANZA MINIMA BARICENTRI triangoli di campo e sorgente (?)
    distmin = sqrt(sum(cF.^2)) - maxF - maxS;
    
    %Calcolo DISTANZA MASSIMA BARICENTRI triangoli di campo e sorgente (?)
    distmax = sqrt(sum(cF.^2)) + maxF + maxS;
            
    %Controllo condizione teorica sull'elemento matriciale 3x3 corrente
    if(((hk-1)*pb_param.velS*dt < distmax) && ((hk+1)*pb_param.velP*dt >= distmin))
        
        %Aggiornamento del parametro iok locale
        cell_isok(ind, 1) = {1};
        
        %Indici dei vertici dell'elemento corrente (triangolo di campo)
        nodeF = domainMesh.triangles(indF, 1:3);

        %Coordinate dei vertici dell'elemento corrente (triangolo di campo)
        TF = domainMesh.coordinates(nodeF, :);
        
        %VERSORE NORMALE all'elemento corrente (triangolo di campo)
        vnF = domainMesh.normal(indF, :);
        
        %Calcolo INTEGRALE DOPPIO sul TRIANGOLO SORGENTE (int. esterna)
        % e sul TRIANGOLO di CAMPO (int. interna)
        ris = BEMenerg_calcMatrixSubBlock(pb_param, TS, areaS, TF, vnF, hk, dt, gha, ghw);
                
        %Pulizia "SPORCIZIA NUMERICA" dal blocchetto matriciale corrente
        ris(abs(ris) < 1.0e-14) = 0;
        
        %Posizionamento blocco matriciale corrente in matrix
        matrix{ind} = matrix{ind} + ris;
     end 
end

%% FORMATTAZIONE RISULTATI

%Trasformazione da array di celle a matrice standard
matrix = cell2mat(matrix);

%Calcolo PARAMETRO iok
isnotzero = max(cell2mat(cell_isok));
return