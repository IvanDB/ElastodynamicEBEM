function matrixK = calcMatrixK(matrixSpecs, nGPU, basePath, pbParam, domainMesh, quadData, constValues)
%CALCMATRIXV Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    matrixSpecs struct
    nGPU        (1, 1) double {mustBeInteger, mustBePositive}
    basePath    (1, 1) string
    pbParam     struct
    domainMesh  struct
    quadData    struct
    constValues cell
end

arguments (Output)
    matrixK cell
end

import eebem.core.*
import eebem.utility.*

%Allocazione array contente i blocchi matriciali
matrixK = cell(matrixSpecs.maxNumBlocksPerIter, matrixSpecs.numIter);

%Ciclo sulle iterazioni necessarie (limite memoria GPU)
for indIter = 1 : matrixSpecs.numIter
    %Calcolo numero blocchi iterazione corrente
    numBlocksThisIter = matrixSpecs.offsets_full(indIter + 1) - matrixSpecs.offsets_full(indIter);

    %Avvio computazione GPU
    spmd(nGPU)
        indGPU = spmdIndex;
        globIdx = nGPU * (indIter - 1) + indGPU;
        numBlockThisLaunch = matrixSpecs.offsets_sing(globIdx + 1) - matrixSpecs.offsets_sing(globIdx);

        matrixOutMulti = [];
        
        if(numBlockThisLaunch > 0)
            gpuID = gpuDevice(spmdIndex);

            srcPath = fullfile(basePath, "+eebem", "+core", "kernelsCUDA", "kernelK.cu");
            ptxPath = fullfile(basePath, "buildDir", "kernelK.ptx");

            kernelK = parallel.gpu.CUDAKernel(ptxPath, srcPath);
            kernelK.GridSize = [matrixSpecs.blockSizes2D numBlockThisLaunch];
            kernelK.ThreadBlockSize = [quadData.methodSpecs.numSRint quadData.methodSpecs.numGHint 1];
            kernelK.SharedMemorySize = quadData.methodSpecs.numINT * 9 * 8;

            gpuInputArrays = copyArrayK(domainMesh, quadData, constValues);

            matrixOutMulti = launchCUDAKernel(gpuID, kernelK, pbParam.deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*pbParam.deltaT, ...
                                                gpuInputArrays.stdEXTw, gpuInputArrays.stdEXTnx, gpuInputArrays.stdEXTny, gpuInputArrays.stdEXTnz, quadData.methodSpecs.numEXT, ...
                                                gpuInputArrays.stdINTw, gpuInputArrays.stdINTnx, gpuInputArrays.stdINTny, gpuInputArrays.stdINTnz, ...
                                                gpuInputArrays.vertsT, gpuInputArrays.areeT, gpuInputArrays.normT, gpuInputArrays.indSMmatrix, gpuInputArrays.matCoeff, gpuInputArrays.vetCoeff, ...
                                                matrixSpecs.offsets_sing(globIdx), matrixSpecs.numBlocks, gpuInputArrays.nodesMesh, domainMesh.lMax);
        end
    end

    %Allocazione array contente i componenti singolari dei blocchi di questa iterazione
    matrixSubBlocksSING = zeros(3, 3, 3, matrixSpecs.blockSizes2D(1), numBlocksThisIter);

    %Avvio computazione CPU
    for indTemp = 1 : numBlocksThisIter
        %Check non sforamento numero di blocchi necessario?
        if matrixSpecs.offsets_full(indIter) + indTemp > matrixSpecs.numBlocks
            continue
        end

        %Set istante temporale
        istTemp = matrixSpecs.offsets_full(indIter) + indTemp - 1;
        
        %Ciclo sull'indice di riga
        parfor indBlock = 1 : matrixSpecs.blockSizes2D(1)
            for indV = 1 : 3
                matrixSubBlocksSING(:, :, indV, indBlock, indTemp) = calcSingSubBlockK(pbParam, domainMesh, quadData.methodSpecs, constValues{indBlock}, istTemp, indBlock, indV);
            end
        end
    end

    matrixOut = cell(nGPU, 1);
    for indGPU = 1 : nGPU
        matrixOut{indGPU} = gather(matrixOutMulti{indGPU});
    end
    matrixOut = vertcat(matrixOut{:});

    %Reshape matrice output GPU
    matrixOut = reshape(matrixOut, [matrixSpecs.blockNumRows matrixSpecs.blockNumCols numBlocksThisIter]);

    %Aggiunta componenti blocchi singolari calcolate su CPU e salvataggio in memoria
    parfor indTemp = 1 : numBlocksThisIter
        %Check non sforamento numero di blocchi necessario?
        if matrixSpecs.offsets_full(indIter) + indTemp > matrixSpecs.numBlocks
            continue
        end

        %Selezione singolo blocco matriciale
        matrixK{indTemp, indIter} = squeeze(matrixOut(:, :, indTemp));

        %Aggiunta sottoblocchi singolari
        for indBlock = 1 : matrixSpecs.blockSizes2D(1)
            [~, indCols] = ismember([1 2 3], domainMesh.indSMmatrix(indBlock, :));
            for indV = 1 : 3
                indR = 3 * (indBlock - 1);
                indC = 3 * (indCols(indV) - 1);
                matrixK{indTemp, indIter}(indR + (1:3), indC + (1:3)) = matrixK{indTemp, indIter}(indR + (1:3), indC + (1:3)) + matrixSubBlocksSING(:, :, indV, indBlock, indTemp);
            end
        end

        %Trasformazione in sparse matrix
        matrixK{indTemp, indIter} = sparse(matrixK{indTemp, indIter});
    end
end

%Reshape matrici salvate
matrixK = reshape(matrixK, [numel(matrixK), 1]);
matrixK = matrixK(1 : matrixSpecs.numBlocks);

% AGGIUNTA BLOCCHI NULLI
for ind = (matrixSpecs.numBlocks + 1) : pbParam.nT
    matrixK{ind} = sparse(zeros(matrixSpecs.blockNumRows, matrixSpecs.blockNumCols));
end

end


function gpuInputArrays = copyArrayK(domainMesh, quadData, constValues)
numT = domainMesh.numTriangles;
numV = domainMesh.numVertices;

%Nodi e pesi GH integrazione esterna
gpuInputArrays.stdEXTw = gpuArray(kron(ones(quadData.methodSpecs.numSRext, 1), quadData.EXTw));
gpuInputArrays.stdEXTnx = gpuArray(quadData.EXTn(:, 1));
gpuInputArrays.stdEXTny = gpuArray(quadData.EXTn(:, 2));
gpuInputArrays.stdEXTnz = gpuArray(quadData.EXTn(:, 3));

%Nodi e pesi GHC integrazione interna
gpuInputArrays.stdINTw = gpuArray(quadData.INTw);
gpuInputArrays.stdINTnx = gpuArray(quadData.INTn(:, 1));
gpuInputArrays.stdINTny = gpuArray(quadData.INTn(:, 2));
gpuInputArrays.stdINTnz = gpuArray(quadData.INTn(:, 3));

%Aree, vertici e vettori normali dei triangoli della mesh spaziale
gpuInputArrays.areeT = gpuArray(domainMesh.area);
gpuInputArrays.vertsT = zeros(9 * domainMesh.numTriangles, 1, 'gpuArray');
gpuInputArrays.normT = zeros(3*numT, 1, 'gpuArray');
for indTemp = 1 : numT
    gpuInputArrays.vertsT(9*(indTemp-1) + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp, 1:3), :), [9 1]);
    gpuInputArrays.normT(3*(indTemp-1) + (1:3), 1) = domainMesh.normal(indTemp, :)';
end

%Matrix and vector base function coefficients
gpuInputArrays.matCoeff = zeros(9*numT, 1, 'gpuArray');
gpuInputArrays.vetCoeff = zeros(3*numT, 1, 'gpuArray');

for indTemp = 1 : numT
    gpuInputArrays.matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    gpuInputArrays.vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end

%Vertex index matrix
gpuInputArrays.indSMmatrix = gpuArray(reshape(domainMesh.indSMmatrix, [numV * numT, 1]));

%Vertex coordinates
gpuInputArrays.nodesMesh = gpuArray(reshape(domainMesh.coordinates', [3*numV 1]));
end