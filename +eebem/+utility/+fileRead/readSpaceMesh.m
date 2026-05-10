function domainMesh = readSpaceMesh(basePath, domainType, lev)
% INPUT 
%   - domain_type: stringa contenente il TIPO di DOMINIO
%   - lev: intero contente il LIVELLO di RAFFINAMENTO
%
% OUTPUT:
%   -  domainMesh: struct contente i dettagli della mesh
 
%% APERTURA del FILE MESH
%Apertura del file mesh
meshFilePath = fullfile(basePath, "mesh", domainType, domainType + "_" + lev + ".mesh");
meshFile = fopen(meshFilePath, 'r');
if meshFile == -1
    error("Impossibile aprire file di mesh")
end

%Lettura dell'intestazione del file mesh
fgets(meshFile);
fgets(meshFile);
fgets(meshFile);
fgets(meshFile);

%% LETTURA dei NODI della MESH
%Lettura linea commentata introduzione sezione dei nodi
fgets(meshFile);

%Lettura NUMERO di NODI
domainMesh.numVertices = sscanf(fgets(meshFile), "%d");

%Allocazione MATRICE contenente le COORDINATE dei NODI
domainMesh.coordinates = zeros(domainMesh.numVertices, 4);

%Ogni riga della matrice è relativa ad un determinato nodo: i primi 
% tre elementi della riga contengono le coordinate del nodo, mentre
% l"ultimo elemento indica a quale faccia appartiene il nodo
 
%FORMATO in cui salvare le COORDINATE dei NODI
formatSpec = "%f";

%DIMENSIONE struttura in cui salvare i dati
sizeSpec = [4 4*domainMesh.numVertices];

%Lettura dei nodi della mesh
domainMesh.coordinates = fscanf(meshFile, formatSpec, sizeSpec)';

%Eliminazione quarta colonna (inutile)
domainMesh.coordinates(:, 4) = [];

%% LETTURA dei TRIANGOLI della MESH
%Lettura linea commentata introduzione sezione degli elementi di bordo
fgets(meshFile);

%Lettura NUMERO di TRIANGOLI
domainMesh.numTriangles = sscanf(fgets(meshFile), "%d");

%Allocazione MATRICE contenente le INCIDENZE dei TRIANGOLI
domainMesh.triangles = zeros(domainMesh.numTriangles, 4);

%Ogni riga è relativa ad un triangolo: i primi tre elementi contengono 
% l"indice dei nodi che costituiscono i tre vertici del triangolo mentre
% l"ultimo elemento identifica il tipo di condizione al bordo assegnata 
% sulla porzione di bordo in cui si trova quel triangolo

%Formato dei dati
formatSpec = "%f";

%DIMENSIONE struttura in cui salvare i dati
sizeSpec = [4 4*domainMesh.numTriangles];

%Lettura incidenze dei triangoli e relativo tipo di dato al bordo 
domainMesh.triangles = fscanf(meshFile, formatSpec, sizeSpec)';

%L"i-esima riga di questa matrice contiene nei primi 3 elementi le incidenze 
% dei tre vertici dell"i-esimo triangolo e nel quarto elemento l"indice
% relativo al tipo di dato al bordo da assegnare al triangolo

% Assegnazione forzata dell"indice di ciascun triangolo
ind = 4;
domainMesh.triangles(:, 4) = ind * ones(domainMesh.numTriangles, 1);

if(domainType == "DesCop-sphere")   %Move to a intern/extern problem dedicated flag
    domainMesh.triangles(:, [2 3]) = domainMesh.triangles(:, [3 2]);
end

%% CHIUSURA DEL FILE MESH
%Lettura dell"ultima riga
fgets(meshFile);

%CHIUSURA del FILE mesh
fclose(meshFile);

%% ESTRAPOLAZIONE INFORMAZIONI UTILI sui TRIANGOLI della MESH
%Creazione matrice contenente le coordinate in sequenza di tutti i primi
% vertici, tutti i secondi e tutti i terzi
coordVertTriangles = domainMesh.coordinates(domainMesh.triangles(:, 1:3), :);

%Calcolo dei vettori P1 - P3 = P3P1 al variare dei triangoli
vg1 = coordVertTriangles(1 : domainMesh.numTriangles, :) - ...
    coordVertTriangles(2*domainMesh.numTriangles+1 : end, :);

%Calcolo dei vettori P2 - P3 = P3P2 al variare dei triangoli
vg2 = coordVertTriangles(domainMesh.numTriangles+1 : 2*domainMesh.numTriangles, :) - ...
    coordVertTriangles(2*domainMesh.numTriangles+1 : end, :);

%Calcolo dei vettori P1 - P2 = P2P1 al variare dei triangoli
vg3 = coordVertTriangles(1 : domainMesh.numTriangles, :) - ...
    coordVertTriangles(domainMesh.numTriangles+1 : 2*domainMesh.numTriangles, :);

%Allocazione matrice contenente le componenti dei vettore normali ai triangoli della mesh
domainMesh.normal = zeros(domainMesh.numTriangles, 3);

%Allocazione vettore contenente le aree dei triangoli della mesh
domainMesh.area = zeros(domainMesh.numTriangles, 1);

%Calcolo COMPONENTI dei VETTORI NORMALI (n = vg1 x vg2)
domainMesh.normal(:, 1) =  vg1(:, 2) .* vg2(:, 3) - vg2(:, 2) .* vg1(:, 3);
domainMesh.normal(:, 2) = -vg1(:, 1) .* vg2(:, 3) + vg2(:, 1) .* vg1(:, 3);
domainMesh.normal(:, 3) =  vg1(:, 1) .* vg2(:, 2) - vg2(:, 1) .* vg1(:, 2);

%Calcolo AREE (A = |n|_2 / 2)
domainMesh.area = sqrt(sum(domainMesh.normal.^2, 2)) / 2;

%Normalizzazione VETTORI NORMALI ai triangoli della mesh
domainMesh.normal = domainMesh.normal ./ (2*domainMesh.area);

%Allocazione matrice contenente le componenti dei baricentri dei triangoli della mesh
domainMesh.center = zeros(domainMesh.numTriangles, 3);

%Calcolo componenti dei baricentri (b = P1 + P2 + P3 / 3)
domainMesh.center = (coordVertTriangles(1 : domainMesh.numTriangles, :) + ...
    coordVertTriangles(domainMesh.numTriangles + 1 : 2*domainMesh.numTriangles, :) + ...
    coordVertTriangles(2*domainMesh.numTriangles + 1 : 3*domainMesh.numTriangles, :)) / 3;

%Calcolo lunghezza dei lati dei triangoli della mesh (l_i = |vg_i|_2)
len = [sqrt(sum(vg3.^2, 2)), sqrt(sum(vg2.^2, 2)), sqrt(sum(vg1.^2, 2))];

%Allocazione vettore contenente le massime lunghezze dei lati dei triangoli della mesh
domainMesh.maxL = zeros(domainMesh.numTriangles, 1);

%Calcolo massime lunghezze dei lati dei triangoli della mesh
domainMesh.maxL = max(len, [], 2);

%Allocazione matrice contenente le componenti dei curl dei
%triangoli della mesh
domainMesh.curl = zeros(3, 3, domainMesh.numTriangles);

%Calcolo delle componenti dei curl
domainMesh.curl(:, 1, :) =  (vg2 ./ (2*domainMesh.area))';
domainMesh.curl(:, 2, :) = -(vg1 ./ (2*domainMesh.area))';
domainMesh.curl(:, 3, :) =  (vg3 ./ (2*domainMesh.area))';

%Vertex matrix
domainMesh.indSMmatrix = zeros(domainMesh.numTriangles, domainMesh.numVertices, 'int32');
for indS = 1 : domainMesh.numVertices
    for indM = 1 : domainMesh.numTriangles
        [~, domainMesh.indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
    end
end

%Calcolo informazioni utili sulla mesh
domainMesh = calcParamMesh(domainMesh);
return
end

function domainMesh = calcParamMesh(domainMesh)
    %Inizializzazione valori globali
    lMin = Inf;
    lMax = -Inf;
    
    % Ciclo sui triangoli della mesh
    for indT = 1 : domainMesh.numTriangles
        %Estrazione incidenze vertici triangolo corrente
        incidenze = domainMesh.triangles(indT, 1:3);
        
        %Estrazione coordinate vertici triangolo corrente
        verts(1, :) = domainMesh.coordinates(incidenze(1), :);
        verts(2, :) = domainMesh.coordinates(incidenze(2), :);
        verts(3, :) = domainMesh.coordinates(incidenze(3), :);      
    
        % Calcolo LUNGHEZZA lati triangolo corrente
        l1 = norm(verts(1, :) - verts(2, :));
        l2 = norm(verts(1, :) - verts(3, :));
        l3 = norm(verts(2, :) - verts(3, :));
    
        % Calcolo LUNGHEZZA MASSIMA e LUNGHEZZA MINIMA
        lMaxCurr = max(max(l1, l2), l3);
        lMinCurr = min(min(l1, l2), l3);
       
        %Aggiornamento valori globali
        if lMaxCurr > lMax 
  	        lMax = lMaxCurr;
        end
        if lMinCurr < lMin 
            lMin = lMinCurr;
        end
    end
    
    %Salvataggio valori
    domainMesh.lMin = lMin;
    domainMesh.lMax = lMax;
    return
end