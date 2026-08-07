function density = timeMarchingID(basePath, pbParam, domainMesh, quadData, fullFileNames)
%TIMEMARCHINGID  Solve the indirect-Dirichlet (ID) energetic BEM formulation by block time-marching.
%   DENSITY = TIMEMARCHINGID(BASEPATH, PBPARAM, DOMAINMESH, QUADDATA,
%   FULLFILENAMES) computes the unknown surface density PHI that solves, for
%   every discrete time step n = 1:PBPARAM.nT, the block-triangular system
%
%       V_0 * DENSITY(:,n) = BETAI{n} - sum_{k=2}^{n} V_{k-1} * DENSITY(:,n-k+1),
%
%   where V is the single-layer block-Toeplitz matrix (CALCMATRIXV) and BETAI is the
%   outer-datum load history (CALCBETAI). The diagonal block V_0 is LU-factorized once
%   and reused at every time step. Intermediate matrices/vectors are optionally saved to
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
%       DENSITY - (3*numTriangles x nT double) the unknown density at every time step.
%
%   Notes:
%       Requires one or more available CUDA-capable GPUs (used by CALCMATRIXV)
%       and, if PBPARAM.lambda + PBPARAM.mu == 0, is the only formulation
%       currently supported by MAIN for that (incompressible-limit) case.
%
%   See also CALCMATRIXV, CALCBETAI, CALCMATRIXSPECS, TIMEMARCHINGDD

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
[numBlocksV, ~, ~] = calcNumMatrixBlocks(pbParam, domainMesh);

constValues = calcConstValues(domainMesh, quadData);

blockSizesV = [domainMesh.numTriangles, domainMesh.numTriangles];
matrixSpecsV = calcMatrixSpecs(nGPU, avMem, blockSizesV, numBlocksV);
matrixV = calcMatrixV(matrixSpecsV, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

% Datum vectors calculations
betaI = calcBetaI(pbParam, domainMesh, constValues, quadData.methodSpecs, basePath);

% Time-marching process
density = zeros(3*domainMesh.numTriangles, pbParam.nT);

[L, U, P] = lu(matrixV{1});

for currInd = 1 : pbParam.nT
    rhs = betaI{currInd}; 

    endInd = min(currInd, numBlocksV);
    for indMat = 2 : endInd
        rhs = rhs - matrixV{indMat} * density(:, currInd - indMat + 1);  
    end
  
    density(:, currInd) = U\(L\(P*rhs));
end

%Save on disk
tmpFlag = true;
if(tmpFlag)
    save(fullFileNames.tmpFullFilename, 'matrixV', 'betaI', '-v7.3');
end
save(fullFileNames.outFullFilename, 'density', '-v7.3');
return
end