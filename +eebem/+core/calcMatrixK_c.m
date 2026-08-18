function matrixK_c = calcMatrixK_c(matrixSpecs, nGPU, basePath, pbParam, domainMesh, quadData, constValues)
%CALCMATRIXK_C  Assemble the block-Toeplitz double-layer (K) matrix in triangle-to-triangle collocation form.
%   MATRIXK_C = CALCMATRIXK_C(MATRIXSPECS, NGPU, BASEPATH, PBPARAM, DOMAINMESH, QUADDATA, CONSTVALUES)
%   computes the same double-layer operator as CALCMATRIXK, but discretized with both test
%   and trial space collocated on triangles (rather than trial on vertices), for use by the
%   "DNc" variant of the direct Neumann formulation. Each block is the sum of an "internal"
%   domain contribution and a "boundary" (near-edge) correction, evaluated on the GPU by the
%   CUDA kernels "kernelKinternal.cu" and "kernelKboundary.cu" respectively; the singular
%   self-triangle contribution is corrected on the CPU via CALCSINGSUBBLOCKK_C.
%   Work is split across MATRIXSPECS.numIter iterations to respect the available GPU memory.
%
%   Input arguments:
%       MATRIXSPECS - (struct) block sizes/offsets/iteration plan, see CALCMATRIXSPECS.
%       NGPU        - (positive integer) number of GPU devices to use.
%       BASEPATH    - (string) project root, used to locate
%                     "+core/kernelsCUDA/kernelKinternal.cu" and
%                     "kernelKboundary.cu" and their compiled PTX.
%       PBPARAM     - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH  - (struct) triangulated boundary mesh, see READSPACEMESH.
%       QUADDATA    - (struct) quadrature nodes/weights/METHODSPECS, see GENERATEQUADDATA.
%       CONSTVALUES - (cell) per-triangle data from CALCCONSTVALUES.
%
%   Output arguments:
%       MATRIXK_C - (cell, PBPARAM.nT x 1) sparse
%                   (3*numTriangles x 3*numTriangles) matrices.
%
%   Notes:
%       Requires one or more available CUDA-capable GPUs and compiled
%       "kernelKinternal.ptx"/"kernelKboundary.ptx" (see AUTOBUILD).
%
%   See also CALCMATRIXK, CALCSINGSUBBLOCKK_C, TIMEMARCHINGDN_C

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
    matrixK_c (:, 1) cell
end

import eebem.core.*
import eebem.utility.*

%Matrix cell array
matrixK_c = cell(matrixSpecs.maxNumBlocksPerIter, matrixSpecs.numIter);

%Loop over number of iterations (GPU memory limit)
for iterIdx = 1 : matrixSpecs.numIter
    %Number of blocks for the current iteration
    numBlocksThisIter = matrixSpecs.offsets_full(iterIdx + 1) - matrixSpecs.offsets_full(iterIdx);

    %Start first GPU work
    spmd(nGPU)
        gpuIdx = spmdIndex;
        globIdx = nGPU * (iterIdx - 1) + gpuIdx;
        numBlockThisLaunch = matrixSpecs.offsets_sing(globIdx + 1) - matrixSpecs.offsets_sing(globIdx);

        matrixOutMultiInt = [];
        
        if(numBlockThisLaunch > 0)
            gpuID = gpuDevice(spmdIndex);

            srcPath = fullfile(basePath, "+eebem", "+core", "kernelsCUDA", "kernelKinternal.cu");
            ptxPath = fullfile(basePath, "buildDir", "kernelKinternal.ptx");

            kernelK = parallel.gpu.CUDAKernel(ptxPath, srcPath);
            kernelK.GridSize = [matrixSpecs.blockSizes2D numBlockThisLaunch];
            kernelK.ThreadBlockSize = [quadData.methodSpecs.numSRint quadData.methodSpecs.numGHint 1];
            kernelK.SharedMemorySize = quadData.methodSpecs.numINT * 9 * 8;

            gpuInputArrays = copyArrayK_c(domainMesh, quadData, constValues);

            matrixOutMultiInt = launchCUDAKernel(gpuID, kernelK, pbParam.deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*pbParam.deltaT, ...
                                                gpuInputArrays.stdEXTw, gpuInputArrays.stdEXTnx, gpuInputArrays.stdEXTny, gpuInputArrays.stdEXTnz, quadData.methodSpecs.numEXT, ...
                                                gpuInputArrays.stdINTw, gpuInputArrays.stdINTnx, gpuInputArrays.stdINTny, gpuInputArrays.stdINTnz, ...
                                                gpuInputArrays.vertsT, gpuInputArrays.areeT, gpuInputArrays.normT, ...
                                                matrixSpecs.offsets_sing(globIdx), matrixSpecs.numBlocks, domainMesh.lMax);
        end
    end

    %Singular component array
    matrixSubBlocksSING = zeros(3, 3, matrixSpecs.blockSizes2D(1), numBlocksThisIter);

    %Start CPU work
    for blockIdx = 1 : numBlocksThisIter
        %Check (redundant?)
        if matrixSpecs.offsets_full(iterIdx) + blockIdx > matrixSpecs.numBlocks
            continue
        end

        %Set time instant
        istTemp = matrixSpecs.offsets_full(iterIdx) + blockIdx - 1;
        
        %Loop over block-row index
        parfor rowIdx = 1 : matrixSpecs.blockSizes2D(1)
            matrixSubBlocksSING(:, :, rowIdx, blockIdx) = calcSingSubBlockK_c(pbParam, domainMesh, quadData.methodSpecs, constValues{rowIdx}, quadData.G1Dn, quadData.G1Dw, istTemp, rowIdx);
        end
    end

    %Gather and reshape first GPU output
    matrixOutInt = cell(nGPU, 1);
    for gpuIdx = 1 : nGPU
        matrixOutInt{gpuIdx} = gather(matrixOutMultiInt{gpuIdx});
    end
    matrixOutInt = vertcat(matrixOutInt{:});

    %Start second GPU work
    spmd(nGPU)
        gpuIdx = spmdIndex;
        globIdx = nGPU * (iterIdx - 1) + gpuIdx;
        numBlockThisLaunch = matrixSpecs.offsets_sing(globIdx + 1) - matrixSpecs.offsets_sing(globIdx);

        matrixOutMultiBound = [];
        
        if(numBlockThisLaunch > 0)
            gpuID = gpuDevice(spmdIndex);

            srcPath = fullfile(basePath, "+eebem", "+core", "kernelsCUDA", "kernelKboundary.cu");
            ptxPath = fullfile(basePath, "buildDir", "kernelKboundary.ptx");

            kernelK = parallel.gpu.CUDAKernel(ptxPath, srcPath);
            kernelK.GridSize = [matrixSpecs.blockSizes2D numBlockThisLaunch];
            kernelK.ThreadBlockSize = [quadData.methodSpecs.numBOUND 1 1];
            kernelK.SharedMemorySize = quadData.methodSpecs.numBOUND * 9 * 8;

            gpuInputArrays = copyArrayK_c(domainMesh, quadData, constValues);

            matrixOutMultiBound = launchCUDAKernel(gpuID, kernelK, pbParam.deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*pbParam.deltaT, ...
                                                gpuInputArrays.stdEXTw, gpuInputArrays.stdEXTnx, gpuInputArrays.stdEXTny, gpuInputArrays.stdEXTnz, quadData.methodSpecs.numEXT, ...
                                                gpuInputArrays.stdINTw, gpuInputArrays.stdINTnx, gpuInputArrays.stdINTny, gpuInputArrays.stdINTnz, ...
                                                gpuInputArrays.vertsT, gpuInputArrays.areeT, gpuInputArrays.normT, ...
                                                matrixSpecs.offsets_sing(globIdx), matrixSpecs.numBlocks, domainMesh.lMax);
        end
    end

    %Gather and reshape first GPU output
    matrixOutBound = cell(nGPU, 1);
    for gpuIdx = 1 : nGPU
        matrixOutBound{gpuIdx} = gather(matrixOutMultiBound{gpuIdx});
    end
    matrixOutBound = vertcat(matrixOutBound{:});

    %Combine and reshape GPU outputs
    matrixOut = matrixOutInt + matrixOutBound;
    matrixOut = reshape(matrixOut, [matrixSpecs.blockNumRows matrixSpecs.blockNumCols numBlocksThisIter]);

    %CPU-GPU matrix assembly
    parfor blockIdx = 1 : numBlocksThisIter
        %Check (redundant?)
        if matrixSpecs.offsets_full(iterIdx) + blockIdx > matrixSpecs.numBlocks
            continue
        end

        %Matrix block extraction
        matrixK_c{blockIdx, iterIdx} = squeeze(matrixOut(:, :, blockIdx));

        %Add singular blocks
        for rowIdx = 1 : matrixSpecs.blockSizes2D(1)
            indRC = 3 * (rowIdx - 1);
            matrixK_c{blockIdx, iterIdx}(indRC + (1:3), indRC + (1:3)) = matrixK_c{blockIdx, iterIdx}(indRC + (1:3), indRC + (1:3)) + matrixSubBlocksSING(:, :, rowIdx, blockIdx);
        end

        %Sparse matrix store
        matrixK_c{blockIdx, iterIdx} = sparse(matrixK_c{blockIdx, iterIdx});
    end
end

%Array reshape and truncation
matrixK_c = reshape(matrixK_c, [numel(matrixK_c), 1]);
matrixK_c = matrixK_c(1 : matrixSpecs.numBlocks);

%Fill null blocks
for ind = (matrixSpecs.numBlocks + 1) : pbParam.nT
    matrixK_c{ind} = sparse(zeros(matrixSpecs.blockNumRows, matrixSpecs.blockNumCols));
end

end


function gpuInputArrays = copyArrayK_c(domainMesh, quadData, constValues)
%Copy the quadrature nodes/weights, mesh geometry and basis-function coefficients
%needed by "kernelKinternal.cu"/"kernelKboundary.cu" onto the GPU as gpuArray inputs.
numT = domainMesh.numTriangles;
numV = domainMesh.numVertices;

%External quadrature nodes/weights
gpuInputArrays.stdEXTw = gpuArray(kron(ones(quadData.methodSpecs.numSRext, 1), quadData.EXTw));
gpuInputArrays.stdEXTnx = gpuArray(quadData.EXTn(:, 1));
gpuInputArrays.stdEXTny = gpuArray(quadData.EXTn(:, 2));
gpuInputArrays.stdEXTnz = gpuArray(quadData.EXTn(:, 3));

%Internal quadrature nodes/weights
gpuInputArrays.stdINTw = gpuArray(quadData.INTw);
gpuInputArrays.stdINTnx = gpuArray(quadData.INTn(:, 1));
gpuInputArrays.stdINTny = gpuArray(quadData.INTn(:, 2));
gpuInputArrays.stdINTnz = gpuArray(quadData.INTn(:, 3));

%Mesh triangles data
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

%Mesh vertex index matrix
gpuInputArrays.indSMmatrix = gpuArray(reshape(domainMesh.indSMmatrix, [numV * numT, 1]));

%Mesh vertex coordinates
gpuInputArrays.nodesMesh = gpuArray(reshape(domainMesh.coordinates', [3*numV 1]));
end