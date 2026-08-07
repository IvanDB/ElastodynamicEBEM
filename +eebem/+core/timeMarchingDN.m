function displacement = timeMarchingDN(basePath, pbParam, domainMesh, quadData, fullFileNames)
%TIMEMARCHINGDN  Solve the direct-Neumann (DN) energetic BEM formulation by block time-marching.
%   DISPLACEMENT = TIMEMARCHINGDN(BASEPATH, PBPARAM, DOMAINMESH, QUADDATA, FULLFILENAMES)
%   computes the unknown surface displacement U (defined at the mesh vertices) that
%   solves, for every discrete time step n = 1:PBPARAM.nT, the block-triangular system
%
%       (K_0 + IGAMMA) * DISPLACEMENT(:,n) = BETAV{n} + IGAMMA * DISPLACEMENT(:,n-1)
%                                             - sum_{k=2}^{n} K_{k-1} * DISPLACEMENT(:,n-k+1),
%
%   where K is the double-layer block-Toeplitz matrix (CALCMATRIXK), IGAMMA is a
%   diagonal "jump" mass-like term built from the mesh incidence/area data, and BETAV is
%   the outer-datum load history from the exact Neumann datum (CALCBETAV).
%   The system matrix (K_0 + IGAMMA) is factorized once (via backslash)
%   and reused, in factorized form, at every time step. 
%   The solution is saved to FULLFILENAMES.outFullFilename and, optionally,
%   intermediate matrices/vectors to FULLFILENAMES.tmpFullFilename.
%
%   Input arguments:
%       BASEPATH      - (string) project root.
%       PBPARAM       - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH    - (struct) triangulated boundary mesh, see READSPACEMESH.
%       QUADDATA      - (struct) quadrature nodes/weights/METHODSPECS, see GENERATEQUADDATA.
%       FULLFILENAMES - (struct) output file paths, see GENERATEFILENAMES.
%
%   Output arguments:
%       DISPLACEMENT - (3*numVertices x nT double) the unknown
%                      displacement at every time step.
%
%   Notes:
%       Requires one or more available CUDA-capable
%       GPUs (used by CALCMATRIXV and CALCMATRIXK).
%
%   See also CALCMATRIXK, CALCMATRIXV, CALCBETAV, TIMEMARCHINGDN_C

arguments (Input)
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
    fullFileNames (1, 1) struct
end

import eebem.core.*

%GPUs initialization
nGPU = gpuDeviceCount("available");
gpuIDs = gpuDevice(1 : nGPU);
reset(gpuIDs);
avMem = min([gpuIDs.AvailableMemory]);

% Matrix calculations
[numBlocksV, numBlocksK, ~] = calcNumMatrixBlocks(pbParam, domainMesh);

constValues = calcConstValues(domainMesh, quadData);

blockSizesV = [domainMesh.numTriangles, domainMesh.numTriangles];
matrixSpecsV = calcMatrixSpecs(nGPU, avMem, blockSizesV, numBlocksV);
matrixV = calcMatrixV(matrixSpecsV, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

blockSizesK = [domainMesh.numTriangles, domainMesh.numVertices];
matrixSpecsK = calcMatrixSpecs(nGPU, avMem, blockSizesK, numBlocksK);
matrixK = calcMatrixK(matrixSpecsK, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

matrixIGamma = kron((domainMesh.indSMmatrix > 0) .* domainMesh.area ./ 6, eye(3));

% Datum vectors calculations
betaV = calcBetaV(pbParam, domainMesh, matrixV, basePath);

% Time-marching process
displacement = zeros(3*domainMesh.numVertices, pbParam.nT);

matrixSist = matrixK{1} + matrixIGamma;

for currInd = 1 : pbParam.nT
    rhs = betaV{currInd};

    if(currInd > 1)
        rhs = rhs + matrixIGamma * displacement(:, currInd - 1);
    end
    endInd = min(currInd, numBlocksK);

    for indMat = 2 : endInd
        rhs = rhs - matrixK{indMat} * displacement(:, currInd - indMat + 1);
    end

    displacement(:, currInd) = matrixSist \ rhs;
end

%Save on disk
tmpFlag = true;
if(tmpFlag)
    save(fullFileNames.tmpFullFilename, 'matrixV', 'matrixK', 'matrixIGamma', 'betaV', '-v7.3');
end
save(fullFileNames.outFullFilename, 'displacement', '-v7.3');
return
end