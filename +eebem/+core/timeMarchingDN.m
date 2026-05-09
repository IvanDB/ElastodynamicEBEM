function solution = timeMarchingDN(basePath, pbParam, domainMesh, quadData)
arguments (Input)
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
end

import eebem.core.*

%Save paths
tmpPath = fullfile(basePath, "tempData", "DN_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);
outPath = fullfile(basePath, "outputData", "DN_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);

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

matrixIGamma = kron((domainMesh.indSMmatrix > 0) .* domainMesh.area ./ 6, eye(3));

% Datum vectors calculations
betaV = calcBetaV(pbParam, domainMesh, matrixV);

% Time-marching process
solution = zeros(3*domainMesh.number_nodes, pbParam.nT);

matrixSist = matrixK{1} + matrixIGamma;

for currInd = 1 : pbParam.nT
    rhs = betaV{currInd};

    if(currInd > 1)
        rhs = rhs + matrixIGamma * solution(:, currInd - 1);
    end
    endInd = min(currInd, numBlocksK);

    for indMat = 2 : endInd
        rhs = rhs - matrixK{indMat} * solution(:, currInd - indMat + 1);
    end

    solution(:, currInd) = matrixSist \ rhs;
end

%Save on disk
tmpFlag = true;
if(tmpFlag)
    save(tmpPath + "_matrix", 'matrixV', 'matrixK', 'matrixIGamma', 'betaV');
end
save(outPath + "_displacement", 'solution');

return