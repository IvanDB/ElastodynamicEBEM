function displacement = timeMarchingDN_c(basePath, pbParam, domainMesh, quadData, fullFileNames)
%TIMEMARCHINGDN_C  Solve the triangle-collocated direct-Neumann (DNc) energetic BEM formulation by block time-marching.
%   DISPLACEMENT = TIMEMARCHINGDN_C(BASEPATH, PBPARAM, DOMAINMESH, QUADDATA, FULLFILENAMES)
%   is the triangle-collocated counterpart of TIMEMARCHINGDN: it computes the unknown
%   surface displacement U (here defined per-triangle rather than per-vertex) using the
%   collocated double-layer matrix CALCMATRIXK_C in place of CALCMATRIXK, and a
%   correspondingly redefined diagonal "jump" term IGAMMA. The time-marching recursion,
%   factorization strategy and output files are otherwise identical to TIMEMARCHINGDN.
%
%   Input arguments:
%       BASEPATH      - (string) project root.
%       PBPARAM       - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH    - (struct) triangulated boundary mesh, see READSPACEMESH.
%       QUADDATA      - (struct) quadrature nodes/weights/METHODSPECS, see GENERATEQUADDATA.
%       FULLFILENAMES - (struct) output file paths, see GENERATEFILENAMES.
%
%   Output arguments:
%       DISPLACEMENT - (3*numTriangles x nT double) the unknown
%                      displacement at every time step.
%
%   Notes:
%       Requires one or more available CUDA-capable
%       GPUs (used by CALCMATRIXV and CALCMATRIXK_C).
%
%   See also CALCMATRIXK_C, CALCMATRIXV, CALCBETAV, TIMEMARCHINGDN

arguments (Input)
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
    fullFileNames (1, 1) struct
end

arguments (Output)
    displacement (:, :) double
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

blockSizesK = [domainMesh.numTriangles, domainMesh.numTriangles];
matrixSpecsK = calcMatrixSpecs(nGPU, avMem, blockSizesK, numBlocksK);
matrixK = calcMatrixK_c(matrixSpecsK, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

matrixIGamma = kron(eye(domainMesh.numTriangles) .* domainMesh.area ./ 2, eye(3));

% Datum vectors calculations
betaV = calcBetaV(pbParam, domainMesh, matrixV, basePath);

% Time-marching process
displacement = zeros(3*domainMesh.numTriangles, pbParam.nT);

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