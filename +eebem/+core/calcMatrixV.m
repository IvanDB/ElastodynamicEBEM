function matrixV = calcMatrixV(matrixSpecs, nGPU, basePath, pbParam, domainMesh, quadData, constValues)
%CALCMATRIXV  Assemble the block-Toeplitz single-layer (V) BEM matrix on the GPU.
%   MATRIXV = CALCMATRIXV(MATRIXSPECS, NGPU, BASEPATH, PBPARAM, DOMAINMESH, QUADDATA, CONSTVALUES)
%   computes the sequence of sparse matrix blocks {V_0, V_1, ..., V_{numBlocks-1}} of
%   the discrete single-layer operator, tested and discretized with piecewise-constant
%   basis functions on the mesh triangles.
%   The regular part of every block is evaluated in parallel on NGPU GPU devices
%   by the CUDA kernel "kernelV.cu"; the singular (self-triangle) contribution
%   is corrected afterwards on the CPU via CALCSINGSUBBLOCKV.
%   Work is split across MATRIXSPECS.numIter iterations to respect the available GPU memory. 
%   Blocks beyond numBlocks (up to PBPARAM.nT) are returned as
%   zero sparse matrices, since the kernel support vanishes there.
%
%   Input arguments:
%       MATRIXSPECS - (struct) block sizes/offsets/iteration plan, see CALCMATRIXSPECS.
%       NGPU        - (positive integer) number of GPU devices to use.
%       BASEPATH    - (string) project root, used to locate
%                     "+core/kernelsCUDA/kernelV.cu" and its compiled PTX.
%       PBPARAM     - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH  - (struct) triangulated boundary mesh, see READSPACEMESH.
%       QUADDATA    - (struct) quadrature nodes/weights/METHODSPECS, see GENERATEQUADDATA.
%       CONSTVALUES - (cell) per-triangle data from CALCCONSTVALUES.
%
%   Output arguments:
%       MATRIXV - (cell, PBPARAM.nT x 1) sparse (3*numTriangles x 3*numTriangles) matrices.
%
%   Notes:
%       Requires one or more available CUDA-capable GPUs
%       and a compiled "kernelV.ptx" (see AUTOBUILD).
%
%   See also CALCMATRIXSPECS, CALCSINGSUBBLOCKV, CALCMATRIXK, CALCMATRIXW, TIMEMARCHINGID, TIMEMARCHINGDD

arguments (Input)
    matrixSpecs (1, 1) struct
    nGPU        (1, 1) double {mustBeInteger, mustBePositive}
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
    constValues (:, 1) cell
end

arguments (Output)
    matrixV (:, 1) cell
end

import eebem.core.*
import eebem.utility.*

matrixV = cell(matrixSpecs.maxNumBlocksPerIter, matrixSpecs.numIter);

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
            matrixSingSubBlocks(:, :, indBlock, indTemp) = calcSingSubBlockV(pbParam, quadData.methodSpecs, constValues{indBlock}, quadData.G1Dn, quadData.G1Dw, istTemp);
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
        matrixV{indTemp, indIter} = squeeze(matrixOut(:, :, indTemp));

        %Inserimento sottoblocchi diagonali
        for indBlock = 1 : matrixSpecs.blockSizes2D(1)
            indRC = 3 * (indBlock - 1);
            matrixV{indTemp, indIter}(indRC + (1:3), indRC + (1:3)) = matrixSingSubBlocks(:, :, indBlock, indTemp);
        end

        %Trasformazione in sparse matrix
        matrixV{indTemp, indIter} = sparse(matrixV{indTemp, indIter});
    end
end

%Reshape matrici salvate
matrixV = reshape(matrixV, [numel(matrixV), 1]);
matrixV = matrixV(1 : matrixSpecs.numBlocks);

% AGGIUNTA BLOCCHI NULLI
for ind = (matrixSpecs.numBlocks + 1) : pbParam.nT
    matrixV{ind} = sparse(zeros(matrixSpecs.blockNumRows, matrixSpecs.blockNumCols));
end

end


function gpuInputArrays = copyArrayV(domainMesh, quadData)
%Copy the quadrature nodes/weights and mesh geometry (vertices,
%areas) needed by "kernelV.cu" onto the GPU as gpuArray inputs.
numT = domainMesh.numTriangles;

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
end