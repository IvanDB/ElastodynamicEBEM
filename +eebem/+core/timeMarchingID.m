function density = timeMarchingID(basePath, pbParam, domainMesh, quadData)
arguments
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
end

import eebem.core.*

%Save paths
tmpPath = fullfile(basePath, "tempData", "ID_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);
outPath = fullfile(basePath, "outputData", "ID_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);

%GPUs inizialization
nGPU = gpuDeviceCount("available");
gpuIDs = gpuDevice(1 : nGPU);
reset(gpuIDs);
avMem = min([gpuIDs.AvailableMemory]);

% Matrix calculations
[numBlocksV, ~] = calcNumMatrixBlocks(pbParam, domainMesh);

constValues = calcConstValues(domainMesh, quadData);

blockSizesV = [domainMesh.numTriangles, domainMesh.numTriangles];
matrixSpecsV = calcMatrixSpecs(nGPU, avMem, blockSizesV, numBlocksV);
matrixV = calcMatrixV(matrixSpecsV, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

% Datum vectors calculations
betaI = calcBetaI(pbParam, domainMesh, constValues, quadData.methodSpecs);

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
    save(tmpPath + "_matrix", 'matrixV', 'betaI');
end
save(outPath + "_density", 'density');
return