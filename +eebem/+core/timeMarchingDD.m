function traction = timeMarchingDD(basePath, pbParam, domainMesh, quadData, fullFileNames)
%TIMEMARCHINGDD  Solve the direct-Dirichlet (DD) energetic BEM formulation by block time-marching.
%   TRACTION = TIMEMARCHINGDD(BASEPATH, PBPARAM, DOMAINMESH, QUADDATA,
%   FULLFILENAMES) computes the unknown surface traction P that solves, for
%   every discrete time step n = 1:PBPARAM.nT, the block-triangular system
%
%       V_0 * TRACTION(:,n) = (BETAI{n}/2 + BETAK{n}) - sum_{k=2}^{n} V_{k-1} * TRACTION(:,n-k+1),
%
%   where V and K are the single- and double-layer block-Toeplitz matrices (CALCMATRIXV,
%   CALCMATRIXK) and BETAI/BETAK are the load histories from the exact Dirichlet datum
%   (CALCBETAI, CALCBETAK). The diagonal block V_0 is LU-factorized once and reused at
%   every time step. Intermediate matrices/vectors are optionally saved to
%   FULLFILENAMES.tmpFullFilename, and the solution to FULLFILENAMES.outFullFilename.
%
%   Input arguments:
%       BASEPATH      - (string) project root.
%       PBPARAM       - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH    - (struct) triangulated boundary mesh, see READSPACEMESH.
%       QUADDATA      - (struct) quadrature nodes/weights/METHODSPECS, see GENERATEQUADDATA.
%       FULLFILENAMES - (struct) output file paths, see GENERATEFILENAMES.
%
%   Output arguments:
%       TRACTION - (3*numTriangles x nT double) the unknown traction at every time step.
%
%   Notes:
%       Requires one or more available CUDA-capable
%       GPUs (used by CALCMATRIXV and CALCMATRIXK).
%
%   See also CALCMATRIXV, CALCMATRIXK, CALCBETAI, CALCBETAK, TIMEMARCHINGID

arguments
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

% Datum vectors calculations
betaI = calcBetaI(pbParam, domainMesh, constValues, quadData.methodSpecs, basePath);
betaK = calcBetaK(pbParam, domainMesh, matrixK, basePath);

% Time-marching process
traction = zeros(3*domainMesh.numTriangles, pbParam.nT);

[L, U, P] = lu(matrixV{1});

for currInd = 1 : pbParam.nT
    rhs = betaI{currInd} ./ 2 + betaK{currInd}; 

    endInd = min(currInd, numBlocksV);

    for indMat = 2 : endInd
        rhs = rhs - matrixV{indMat} * traction(:, currInd - indMat + 1);  
    end
  
    traction(:, currInd) = U\(L\(P*rhs));
end

%Save on disk
tmpFlag = true;
if(tmpFlag)
    save(fullFileNames.tmpFullFilename, 'matrixV', 'matrixK', 'betaI', 'betaK', '-v7.3');
end
save(fullFileNames.outFullFilename, 'traction', '-v7.3');
return
end