% Script per eseguire i vari test sui blocchetti per il W

import eebem.*
format longG
warning off


basePath = fileparts(mfilename("fullpath"));

% %Set values
% if ~exist('glbIndexFigures', 'var')
%     glbIndexFigures = 0;
% end

problemFileName = "input_barH1-symm_lev1.txt";

pbParam = utility.fileRead.readInputFile(basePath, problemFileName);

domainMesh = utility.fileRead.readSpaceMesh(basePath, pbParam.domainType, pbParam.lev);
% glbIndexFigures = utility.plots.plotMesh(domainMesh, glbIndexFigures);


numNodes = domainMesh.numVertices;
numTriangles = domainMesh.numTriangles;

% Qui Creo Le matrici che contengono le informazioni sulla mesh che mi
% servono per fare il test 

% Faccio il ciclo su ogni triangolo e ne guardo i vertici e aggiungo 1 al
% numero di triangoli che toccano quel vertice
counts = zeros(numNodes, 1);
for t = 1:numTriangles
    for i = 1:3
        nodo = domainMesh.triangles(t, i);
        counts(nodo) = counts(nodo) + 1;
    end
end

% Trovo il numero massimo di triangoli connessi a un singolo nodo
maxNumTriangles = max(counts);

TriangPerNodes = zeros(numNodes, maxNumTriangles);
indSMatrix = zeros(numNodes, maxNumTriangles);

% Vettore per tenere traccia di "quale colonna" stiamo riempiendo per ogni nodo
currentCol = zeros(numNodes, 1);

% 
for t = 1:numTriangles
    for i = 1:3 % I tre vertici del triangolo

        nodo = domainMesh.triangles(t, i); % ID Globale del nodo

        % sposto la colonna dove scrivere per questo nodo (ne devo tenere traccia per il singolo nodo: tutte le voltr
        % che lo trovo devo spostarmi di una colonna)
        currentCol(nodo) = currentCol(nodo) + 1;
        position = currentCol(nodo);

        % Creo la matrice che mi dice per ogni nodo (riga) quali sono gli indici dei triangoli di cui è vertice 
        TriangPerNodes(nodo, position) = t;


        % Con questa matrice so anche che vertice è all'interno di quale
        % triangolo e c'è corrispondenza tra TriangPerNodes e indSMatrix,
        % nel senso che fissato il nodo (riga) la prima colonna di
        % TriangPerNodes ci dice l'indice del primo triangolo che ha tale
        % nodo come vertice, e la prima colonna di indSMatrix ci dice la
        % posizione come vertice di quel nodo in quello tesso triangolo!
        indSMatrix(nodo, position) = i; 
    end
end

% scelta dei nodi e dell'istante temporale (che blocchetto della matriciona voglio calcolare)
outerNode = 1;
innerNode = 1;
timeInstant = 1;
quadData = core.setupCore(10);
constData = core.calcConstValues(domainMesh, quadData);
tic;
testBlock = TestingForWKernel(outerNode, innerNode, timeInstant, indSMatrix,...
    TriangPerNodes, pbParam, domainMesh, quadData, constData, maxNumTriangles);
toc
tic;
testBlock2 = Copy_of_TestingForWKernel(outerNode, innerNode, timeInstant, indSMatrix,...
    TriangPerNodes, pbParam, domainMesh, quadData, constData, maxNumTriangles);
toc
disp(testBlock);
disp(testBlock2);