function matrixSavedV = calcMatrixV(matrixSpecs, nGPU, basePath, pbParam, domainMesh, quadData, constValues)
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
    matrixSavedV cell
end

import eebem.core.*

matrixSavedV = cell(matrixSpecs.maxNumBlocksPerIter, matrixSpecs.numIter);

% %Ciclo sulle iterazioni necessarie (limite memoria GPU)
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
            gpuID = gpuDevice(indGPU);

            srcPath = fullfile(basePath, "+eebem", "+core", "kernelsCUDA", "kernelV.cu");
            ptxPath = fullfile(basePath, "buildDir", "kernelV.ptx");

            kernelV = parallel.gpu.CUDAKernel(ptxPath, srcPath);
            kernelV.GridSize = [matrixSpecs.blockSizes2D numBlockThisLaunch];
            kernelV.ThreadBlockSize = [quadData.methodSpecs.numSRint quadData.methodSpecs.numGHint 1];
            kernelV.SharedMemorySize = quadData.methodSpecs.numINT * 9 * 8;

            gpuInputArrays = copyArrayV(domainMesh, quadData);

            matrixOutMulti = launchCUDAKernel(gpuID, kernelV, pbParam.deltaT, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho, ...
                                                gpuInputArrays.stdEXTw, gpuInputArrays.stdEXTnx, gpuInputArrays.stdEXTny, gpuInputArrays.stdEXTnz, quadData.methodSpecs.numEXT, ...
                                                gpuInputArrays.stdINTw, gpuInputArrays.stdINTnx, gpuInputArrays.stdINTny, gpuInputArrays.stdINTnz, ...
                                                gpuInputArrays.vertsT, gpuInputArrays.areeT, matrixSpecs.offsets_sing(globIdx), matrixSpecs.numBlocks);
        end
    end

    %Allocazione array contente i blocchi diagonali di questa iterazione
    matrixSingSubBlocks = zeros(3, 3, matrixSpecs.blockSizes2D(1), numBlocksThisIter);

    %Avvio computazione CPU
    for indTemp = 1 : numBlocksThisIter
        %Check non sforamento numero di blocchi necessario
        if matrixSpecs.offsets_full(indIter) + indTemp > matrixSpecs.numBlocks
            continue
        end

        %Set istante temporale
        istTemp = matrixSpecs.offsets_full(indIter) + indTemp - 1;

        %Ciclo sull'indice dei sottoblocchi
        parfor indBlock = 1 : matrixSpecs.blockSizes2D(1)
            matrixSingSubBlocks(:, :, indBlock, indTemp) = calcSingSubBlockV(pbParam, quadData.methodSpecs, constValues{indBlock}, istTemp);
        end
    end

    %Attesa completamento operazioni GPU
    matrixOut = cell(nGPU, 1);
    for indGPU = 1 : nGPU
        matrixOut{indGPU} = gather(matrixOutMulti{indGPU});
    end
    matrixOut = vertcat(matrixOut{:});

    %Reshape matrice output GPU
    matrixOut = reshape(matrixOut, [matrixSpecs.blockNumRows matrixSpecs.blockNumCols numBlocksThisIter]);

    %Inserimento blocchi diagonali calcolati su CPU
    parfor indTemp = 1 : numBlocksThisIter
        %Check non sforamento numero di blocchi necessario?
        if matrixSpecs.offsets_full(indIter) + indTemp > matrixSpecs.numBlocks
            continue
        end

        %Selezione singolo blocco matriciale
        matrixSavedV{indTemp, indIter} = squeeze(matrixOut(:, :, indTemp));

        %Inserimento sottoblocchi diagonali
        for indBlock = 1 : matrixSpecs.blockSizes2D(1)
            indRC = 3 * (indBlock - 1);
            matrixSavedV{indTemp, indIter}(indRC + (1:3), indRC + (1:3)) = matrixSingSubBlocks(:, :, indBlock, indTemp);
        end

        %Trasformazione in sparse matrix
        matrixSavedV{indTemp, indIter} = sparse(matrixSavedV{indTemp, indIter});
    end
end

%Reshape matrici salvate
matrixSavedV = reshape(matrixSavedV, [numel(matrixSavedV), 1]);
matrixSavedV = matrixSavedV(1 : matrixSpecs.numBlocks);

% AGGIUNTA BLOCCHI NULLI
for ind = (matrixSpecs.numBlocks + 1) : pbParam.nT
    matrixSavedV{ind} = sparse(zeros(matrixSpecs.blockNumRows, matrixSpecs.blockNumCols));
end

end


function gpuInputArrays = copyArrayV(domainMesh, quadData)
numT = domainMesh.numberTriangles;

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
gpuInputArrays.vertsT = zeros(9 * domainMesh.numberTriangles, 1, 'gpuArray');
gpuInputArrays.normT = zeros(3*numT, 1, 'gpuArray');
for indTemp = 1 : numT
    gpuInputArrays.vertsT(9*(indTemp-1) + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp, 1:3), :), [9 1]);
    gpuInputArrays.normT(3*(indTemp-1) + (1:3), 1) = domainMesh.normal(indTemp, :)';
end
end