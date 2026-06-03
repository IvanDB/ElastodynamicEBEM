% Script che fa la matrice intera del W usando le funzioncine di test che
% generano i singoli blocchetti 3x3

delete(gcp("nocreate"));
parInfo = parpool("Threads");

import eebem.*
problemFileName = "input_barH1-symm_lev0.txt";

pbParam = utility.fileRead.readInputFile(".", problemFileName);
domainMesh = utility.fileRead.readSpaceMesh(".", pbParam.domainType, pbParam.lev);
quadData = core.setupCore(10);
constData = core.calcConstValues(domainMesh, quadData);

numNodes = domainMesh.numVertices;
numTriangles = domainMesh.numTriangles;

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
for t = 1:numTriangles
    for i = 1:3 
        nodo = domainMesh.triangles(t, i); 
        currentCol(nodo) = currentCol(nodo) + 1;
        position = currentCol(nodo);
        TriangPerNodes(nodo, position) = t;
        indSMatrix(nodo, position) = i; 
    end
end
[~, ~, numBlocksW] = core.calcNumMatrixBlocks(pbParam, domainMesh);
    disp(numBlocksW)
for indT = 1 : numBlocksW
    disp(indT)
    for indNodeExt = 1 : numNodes
    disp(indNodeExt)
        parfor indNodeInt = 1 : numNodes
            OutputMatrix{indNodeExt,indNodeInt,indT} = Copy_of_TestingForWKernel(indNodeExt, indNodeInt, indT, indSMatrix,...
                    TriangPerNodes, pbParam, domainMesh, quadData, constData, maxNumTriangles);
        end
    end
end
save("W_matrix.mat", "OutputMatrix");
%OutputMatrix{:,:,numBlocksW+1:pbParam.nT} = zeros(3,3);
matrixIGamma = 1/deltaT*kron((domainMesh.indSMmatrix > 0) .* domainMesh.area ./ 6, eye(3));
gV = core.calcBoundDataNeumann(pbParam, domainMesh);

displacement = zeros(3*domainMesh.numVertices, pbParam.nT);

matrixSist = cell2mat(OutputMatrix(:, :, 1));

for currInd = 1 : pbParam.nT
    rhs = matrixIGamma*gV{currInd};

    
    endInd = min(currInd, numBlocksW);

    for indMat = 2 : endInd
        rhs = rhs - cell2mat(OutputMatrix(:, :, indMat)) * displacement(:, currInd - indMat + 1);
    end

    displacement(:, currInd) = matrixSist \ rhs;
end
glbIndexFigures = utility.plots.plotLinear(".", pbParam, domainMesh, displacement, 0);