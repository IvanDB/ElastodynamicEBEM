function matrixW = calcMatrixW(matrixSpecs, nGPU, basePath, pbParam, domainMesh, quadData, constValues)
%CALCMATRIXW  Assemble the block-Toeplitz hypersingular (W) BEM matrix on the GPU.
%   MATRIXW = CALCMATRIXW(MATRIXSPECS, NGPU, BASEPATH, PBPARAM, DOMAINMESH, QUADDATA, CONSTVALUES)
%   computes the sequence of sparse matrix blocks {W_0, W_1, ..., W_{numBlocks-1}} of
%   the discrete hypersingular operator, tested and discretized with piecewise-linear
%   (vertex-based, continuous) basis functions. 
%   The regular part of every block is evaluated on NGPU GPU devices
%   by the CUDA kernel "kernelW.cu"; vertex pairs that share a triangle
%   get an additional singular correction computed on the CPU via CALCSINGSUBBLOCKW. 
%   Work is split across MATRIXSPECS.numIter iterations to respect the available GPU memory. 
%   Blocks beyond numBlocks (up to PBPARAM.nT) are returned as
%   zero sparse matrices, since the kernel support vanishes there.
%
%   Input arguments:
%       MATRIXSPECS - (struct) block sizes/offsets/iteration plan, see CALCMATRIXSPECS.
%       NGPU        - (positive integer) number of GPU devices to use.
%       BASEPATH    - (string) project root, used to locate
%                     "+core/kernelsCUDA/kernelW.cu" and its compiled PTX.
%       PBPARAM     - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH  - (struct) triangulated boundary mesh, see READSPACEMESH.
%       QUADDATA    - (struct) quadrature nodes/weights/METHODSPECS, see GENERATEQUADDATA.
%       CONSTVALUES - (cell) per-triangle data from CALCCONSTVALUES.
%
%   Output arguments:
%       MATRIXW - (cell, PBPARAM.nT x 1) sparse (3*numVertices x 3*numVertices) matrices.
%
%   Notes:
%       Requires one or more available CUDA-capable GPUs
%       and a compiled "kernelW.ptx" (see AUTOBUILD).
%
%   See also CALCMATRIXSPECS, CALCSINGSUBBLOCKW, CALCMATRIXV, CALCMATRIXK, TIMEMARCHINGIN
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
    matrixW cell
end

import eebem.core.*
import eebem.utility.*

%Allocazione array contente i blocchi matriciali
matrixW = cell(matrixSpecs.maxNumBlocksPerIter, matrixSpecs.numIter);
% Create mesh node adiacency matrix for singular integrations
isSingularPair = nodeAdj(domainMesh.triangles, matrixSpecs.blockSizes2D(1));
for indIter = 1 : matrixSpecs.numIter
    %Calcolo numero blocchi iterazione corrente
    numBlocksThisIter = matrixSpecs.offsets_full(indIter + 1) - matrixSpecs.offsets_full(indIter);

    %Avvio computazione GPU
    % spmd(nGPU)
        indGPU = 1; %spmdIndex;

        globIdx = nGPU * (indIter - 1) + indGPU;
        numBlockThisLaunch = matrixSpecs.offsets_sing(globIdx + 1) - matrixSpecs.offsets_sing(globIdx);

        matrixOutMulti = [];
        
        if(numBlockThisLaunch > 0)
            gpuID = gpuDevice(indGPU); % gpuDevice(spmdIndex);

            srcPath = fullfile(basePath, "+eebem", "+core", "kernelsCUDA", "kernelW.cu");
            ptxPath = fullfile(basePath, "buildDir", "kernelW.ptx");

            kernelW = parallel.gpu.CUDAKernel(ptxPath, srcPath);
            kernelW.GridSize = [matrixSpecs.blockSizes2D numBlockThisLaunch];
            kernelW.ThreadBlockSize = [quadData.methodSpecs.numSRint quadData.methodSpecs.numGHint 1];
            kernelW.SharedMemorySize = quadData.methodSpecs.numINT * 9 * 8;

            [gpuInputArrays, indSMmatrix, TriangPerNodes, maxTrianglesPerNode] = copyArrayW(domainMesh, quadData, constValues);
	    
            matrixOutMulti = launchCUDAKernel(gpuID, kernelW, pbParam.deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, pi, ...
                                                gpuInputArrays.stdEXTw, gpuInputArrays.stdEXTnx, gpuInputArrays.stdEXTny, gpuInputArrays.stdEXTnz, quadData.methodSpecs.numEXT, ...
                                                gpuInputArrays.stdINTw, gpuInputArrays.stdINTnx, gpuInputArrays.stdINTny, gpuInputArrays.stdINTnz, ...
                                                gpuInputArrays.vertsT, gpuInputArrays.areeT, gpuInputArrays.normT, gpuInputArrays.indSMmatrix, gpuInputArrays.matCoeff, gpuInputArrays.vetCoeff, ...
                                                matrixSpecs.offsets_sing(globIdx), matrixSpecs.numBlocks, gpuInputArrays.nodesMesh, domainMesh.lMax, maxTrianglesPerNode, gpuInputArrays.trianglesPerNode);
        end
    %end

    %Allocating matrix containig singular sub blocks
    matrixSubBlocksSING = zeros(matrixSpecs.blockNumRows, matrixSpecs.blockNumCols, numBlocksThisIter);
    
    %Avvio computazione CPU
    for indTemp = 1 : numBlocksThisIter
        %Check non sforamento numero di blocchi necessario?
        if matrixSpecs.offsets_full(indIter) + indTemp > matrixSpecs.numBlocks
            continue
        end

        %Set istante temporale
        istTemp = matrixSpecs.offsets_full(indIter) + indTemp - 1;
        
        % 2D block containing all singular contributions fot this time
        % istant block
        %currentMatrixSing = zeros(matrixSpecs.blockNumRows, matrixSpecs.blockNumCols);

        % temporary Cell array to save Rows of 2D Block
        tempMatrixRows = cell(matrixSpecs.blockSizes2D(1), 1);
        parfor indNodeExt = 1 : matrixSpecs.blockSizes2D(1)
            % row slice of current singular block
            tempMatrixRowSing = zeros(3, matrixSpecs.blockNumCols);

            for indNodeInt = 1 : matrixSpecs.blockSizes2D(2)

                if isSingularPair(indNodeExt, indNodeInt)
                    subBlockSingW = calcSingSubBlockW(pbParam, domainMesh, quadData, constValues, indSMmatrix, TriangPerNodes, maxTrianglesPerNode, istTemp, indNodeExt, indNodeInt);
                    % Filling row of blocks
                    tempMatrixRowSing(:, (3 * (indNodeInt-1) + 1) : (3 * indNodeInt)) = subBlockSingW;
                end
            end
            % Putting Row in the corresponding cell
            tempMatrixRows{indNodeExt} = tempMatrixRowSing;
        end
        % stacking all rows to form 2D block
        currentMatrixSing = vertcat(tempMatrixRows{:});
        % placing block at current instant
        matrixSubBlocksSING(:, :, indTemp) = currentMatrixSing;
    end
    

    % matrixOut = cell(nGPU, 1);
    % for indGPU = 1 : nGPU
    %     matrixOut{indGPU} = gather(matrixOutMulti{indGPU});
    % end
    % matrixOut = vertcat(matrixOut{:});
    
    matrixOut = gather(matrixOutMulti);
    
    %Reshape matrice output GPU
    matrixOut = reshape(matrixOut, [matrixSpecs.blockNumRows matrixSpecs.blockNumCols numBlocksThisIter]);
    % Adding singular sub blocks
    matrixOut = matrixOut + matrixSubBlocksSING;

    % Loops to get the same matrix stucture of other kernels
    parfor indTemp = 1 : numBlocksThisIter
        if matrixSpecs.offsets_full(indIter) + indTemp > matrixSpecs.numBlocks
            continue
        end
        
        % Squeeze out the block (indTemp)
        matrixW_temp = squeeze(matrixOut(:, :, indTemp));
        
        % Turn into sparse matrix
        matrixW{indTemp, indIter} = sparse(matrixW_temp);
    end
    disp("Iteration " + num2str(indIter)+ " out of "+num2str(matrixSpecs.numIter)+ " completed");
end 

% Reshape matrixW (numBlocksThisIter x numIter) into 1D cell array
matrixW = reshape(matrixW, [numel(matrixW), 1]);

% Taglio eventuali blocchi vuoti extra allocati
matrixW = matrixW(1 : matrixSpecs.numBlocks);

% add null blocks
for ind = (matrixSpecs.numBlocks + 1) : pbParam.nT
    matrixW{ind} = sparse(zeros(matrixSpecs.blockNumRows, matrixSpecs.blockNumCols));
end
end

function [gpuInputArrays, indSMmatrix, TriangPerNodes, maxTrianglesPerNode] = copyArrayW(domainMesh, quadData, constValues)
%Build the vertex-to-triangle incidence tables (TRIANGPERNODES, INDSMMATRIX,
%0-padded to MAXTRIANGLESPERNODE columns) and copy the quadrature nodes/weights
%and mesh geometry needed by "kernelW.cu" onto the GPU as gpuArray inputs.

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
maxTrianglesPerNode = max(counts);
TriangPerNodes = zeros(numNodes, maxTrianglesPerNode, 'int32');
indSMmatrix = zeros(numNodes, maxTrianglesPerNode, 'int32');
% Vettore per tenere traccia di "quale colonna" stiamo riempiendo per ogni nodo
currentCol = zeros(numNodes, 1);
for t = 1:numTriangles
    for i = 1:3 
        nodo = domainMesh.triangles(t, i); 
        currentCol(nodo) = currentCol(nodo) + 1;
        position = currentCol(nodo);
        TriangPerNodes(nodo, position) = t;
        indSMmatrix(nodo, position) = i; 
    end
end

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
gpuInputArrays.normT = zeros(3*numTriangles, 1, 'gpuArray');
for indTemp = 1 : numTriangles
    gpuInputArrays.vertsT(9*(indTemp-1) + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp, 1:3), :), [9 1]);
    gpuInputArrays.normT(3*(indTemp-1) + (1:3), 1) = domainMesh.normal(indTemp, :)';
end

%Matrix and vector base function coefficients
gpuInputArrays.matCoeff = zeros(9*numTriangles, 1, 'gpuArray');
gpuInputArrays.vetCoeff = zeros(3*numTriangles, 1, 'gpuArray');

for indTemp = 1 : numTriangles
    gpuInputArrays.matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    gpuInputArrays.vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end

%Vertex index matrix
gpuInputArrays.indSMmatrix = zeros(numNodes * maxTrianglesPerNode, 1, 'int32');
gpuInputArrays.indSMmatrix = gpuArray(reshape(indSMmatrix, [numNodes * maxTrianglesPerNode, 1]));

%Triangles per node matrix
gpuInputArrays.trianglesPerNode = zeros(numNodes * maxTrianglesPerNode, 1, 'int32');
gpuInputArrays.trianglesPerNode = gpuArray(reshape(TriangPerNodes, [numNodes * maxTrianglesPerNode, 1]));

%Vertex coordinates
gpuInputArrays.nodesMesh = gpuArray(reshape(domainMesh.coordinates', [3*numNodes 1]));
end

function adjMatrix = nodeAdj(triangleAdj, numNodes)
    %Build the vertex-to-vertex adjacency matrix (including self- adjacency)
    %from the mesh triangle connectivity: ADJMATRIX(i,j) is true iff vertices i
    %and j share at least one triangle (or i == j). Used by CALCMATRIXW to flag
    %which vertex pairs need the singular correction from CALCSINGSUBBLOCKW.
    adj = sparse(triangleAdj(:,1), triangleAdj(:,2), 1, numNodes, numNodes) + ...
    sparse(triangleAdj(:,2), triangleAdj(:,3), 1, numNodes, numNodes) + ...
    sparse(triangleAdj(:,3), triangleAdj(:,1), 1, numNodes, numNodes);
          
    adj = adj + adj';
    adjMatrix = (adj > 0) | speye(numNodes, numNodes);
end