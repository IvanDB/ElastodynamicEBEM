function rhs = BEMenerg_core_calcTnBlock(rhs, methodInfo, indTemp, deltaT, pbParam, domainMesh, constValues)
%INPUT
%   - rhs: vettore nullo su cui memorizzare le componenti
%   - methodInfo: struct contenente le informazioni sul metodo
%   - indTemp: intero contente l'indice della matrice da calcolare
%   - deltaT: double contente l'ampiezza della discretizzazione temporale
%   - pbParam: struct contenente i parametri del problema
%   - domainMesh: struct contente le informazioni sulla mesh
%   - gha: matrice contenente i nodi di quadratura di Gauss-Hammer
%   - ghw: matrice contenente i pesi di quadratura di Gauss-Hammer
% OUTPUT
%   - rhs: vettore con le componenti del termine noto


%% SETUP VARIABILI in ARRAY di CELLE

%Estrazione NUMERO di TRIANGOLI della MESH
numTriangles = domainMesh.numberTriangles;

%Trasformazione vettore rhs in array di N_triangles x 1 celle
% ciascuna di dimensione 3x1
rhs = mat2cell(rhs, 3*ones(1, numTriangles), 1);

%Inizializzazione numero nodi Gauss-Hammer
nGH = methodInfo.numNodiExt;

%Calcolo istanti temporale iniziale e finale
tInz = deltaT * indTemp;
tFin = deltaT * (indTemp + 1);

%% CICLO SUI TRIANGOLI SORGENTE E DI CAMPO

%Procedimento element-by-element
for indT = 1 : numTriangles    
    %Estrapolazione INDICE TIPO DATO al BORDO triangolo corrente
    indT_RHS = domainMesh.triangles(indT, 4);

    %Calcolo VETTORE TERMINE NOTO LOCALE
    ris = zeros(3, 1);
    
    %Ciclo sui nodi di GH del triangolo corrente
    for indGH = 1 : nGH
        %Estrazione coordinate nodo corrente
        ghT = constValues{indT}.GHnodes{indGH};
    
        %Calcolo termine corrispondente al nodo corrente
        risTemp = BEMenerg_core_calcNucleoTn(pbParam, ghT, tInz, tFin, indT_RHS);
        
        %Aggiornamento risultato complessivo
        ris = ris + risTemp * constValues{indT}.GHweights{indGH};
    end

    %Pulizia "SPORCIZIA NUMERICA" dal vettore locale
    ris(abs(ris) < 1.0e-14) = 0;

    %Posizionamento vettore locale in rhs
    rhs{indT} = ris;
end

%% FORMATTAZIONE RISULTATI
%Trasformazione da array di celle a vettore standard
rhs = cell2mat(rhs);

return