function solution = timeMarchingDD(basePath, pbParam, domainMesh, quadData)
arguments
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
end

import eebem.core.*

%Save paths
tmpPath = fullfile(basePath, "tempData", "DD_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);
outPath = fullfile(basePath, "outputData", "DD_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);

%GPUs inizialization
nGPU = gpuDeviceCount("available");
gpuIDs = gpuDevice(1 : nGPU);
reset(gpuIDs);
avMem = min([gpuIDs.AvailableMemory]);

% Matrix calculations
[numBlocksV, numBlocksK] = calcNumMatrixBlocks(pbParam, domainMesh);

constValues = calcConstValues(domainMesh, quadData);

blockSizesV = [domainMesh.numberTriangles, domainMesh.numberTriangles];
matrixSpecsV = calcMatrixSpecs(nGPU, avMem, blockSizesV, numBlocksV);
matrixV = calcMatrixV(matrixSpecsV, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

blockSizesK = [domainMesh.numberTriangles, domainMesh.number_nodes];
matrixSpecsK = calcMatrixSpecs(nGPU, avMem, blockSizesK, numBlocksK);
matrixK = calcMatrixK(matrixSpecsK, nGPU, basePath, pbParam, domainMesh, quadData, constValues);

% Datum vectors calculations
betaI = calcBetaI(pbParam, domainMesh, constValues, quadData.methodSpecs);
betaK = calcBetaK(pbParam, domainMesh, matrixK);

% Time-marching process
solution = zeros(3*domainMesh.numberTriangles, pbParam.nT);

[L, U, P] = lu(matrixV{1});

for currInd = 1 : pbParam.nT
    rhs = betaI{currInd} ./ 2 + betaK{currInd}; 

    endInd = min(currInd, numBlocksV);

    for indMat = 2 : endInd
        rhs = rhs - matrixV{indMat} * solution(:, currInd - indMat + 1);  
    end
  
    solution(:, currInd) = U\(L\(P*rhs));
end

%Save on disk
tmpFlag = true;
if(tmpFlag)
    save(tmpPath + "_matrix", 'matrixV', 'matrixK', 'betaI', 'betaK');
end
save(outPath + "_traction", 'solution');
return