function density = timeMarchingIN(basePath, pbParam, domainMesh, quadData, fullFileNames)
%TIMEMARCHINGIN  Solve the indirect-Neumann (IN) energetic BEM formulation by block time-marching.
%   DENSITY = TIMEMARCHINGIN(BASEPATH, PBPARAM, DOMAINMESH, QUADDATA, FULLFILENAMES)
%   computes the unknown surface density PHI (defined at the mesh vertices) that
%   solves, for every discrete time step n = 1:PBPARAM.nT, the block-triangular system
%
%       W_0 * DENSITY(:,n) = IGAMMA * GV{n} - sum_{k=2}^{n} W_{k-1} * DENSITY(:,n-k+1),
%
%   where W is the hypersingular block-Toeplitz matrix (CALCMATRIXW), IGAMMA is a
%   diagonal mass-like term built from the mesh incidence/area data, and GV is the
%   exact Neumann datum sampled at the mesh vertices (CALCBOUNDDATANEUMANN). The
%   diagonal block W_0 is LU-factorized once and reused at every time step.
%   Intermediate matrices/vectors are optionally saved to
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
%       DENSITY - (3*numVertices x nT double) the unknown density at every time step.
%
%   Notes:
%       Requires one or more available CUDA-capable GPUs (used by CALCMATRIXW). Unlike
%       the other four formulations, this datum is sampled directly via
%       CALCBOUNDDATANEUMANN rather than convolved through a matrix (compare CALCBETAV).
%
%   See also CALCMATRIXW, CALCBOUNDDATANEUMANN, TIMEMARCHINGID

arguments
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
    fullFileNames (1, 1) struct
end

import eebem.core.*

%GPUs inizialization
nGPU = gpuDeviceCount("available");
gpuIDs = gpuDevice(1 : nGPU);
reset(gpuIDs);
avMem = min([gpuIDs.AvailableMemory]);

% Matrix calculations
[~, ~, numBlocksW] = calcNumMatrixBlocks(pbParam, domainMesh);

constValues = calcConstValues(domainMesh, quadData);

blockSizesW = [domainMesh.numVertices, domainMesh.numVertices];
matrixSpecsW = calcMatrixSpecs(nGPU, avMem, blockSizesW, numBlocksW);
matrixW = calcMatrixW(matrixSpecsW, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

% Datum vectors calculations
matrixIGamma = kron((domainMesh.indSMmatrix > 0) .* domainMesh.area ./ 3, eye(3))';
gV = calcBoundDataNeumann(pbParam, domainMesh);

% Time-marching process
density = zeros(3*domainMesh.numVertices, pbParam.nT);

[L, U, P] = lu(matrixW{1});

for currInd = 1 : pbParam.nT
    rhs = matrixIGamma*gV{currInd}; 

    endInd = min(currInd, numBlocksW);
    for indMat = 2 : endInd
        rhs = rhs - matrixW{indMat} * density(:, currInd - indMat + 1);  
    end
  
    density(:, currInd) = U\(L\(P*rhs));
end

%Save on disk
tmpFlag = true;
if(tmpFlag)
    save(fullFileNames.tmpFullFilename, 'matrixW', 'gV', '-v7.3');
end
save(fullFileNames.outFullFilename, 'density', '-v7.3');
return
end