function density = timeMarchingIN(basePath, pbParam, domainMesh, quadData)
arguments
    basePath    (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadData    (1, 1) struct
end

import eebem.core.*

%Save paths
tmpPath = fullfile(basePath, "tempData", "IN_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);
outPath = fullfile(basePath, "outputData", "IN_" + pbParam.domainType + pbParam.lev + quadData.methodSpecs.stringID);

%GPUs inizialization
nGPU = gpuDeviceCount("available");
gpuIDs = gpuDevice(1 : nGPU);
reset(gpuIDs);
avMem = min([gpuIDs.AvailableMemory]);

% Matrix calculations
[~, ~, numBlocksW] = calcNumMatrixBlocks(pbParam, domainMesh); %% Need to modify this function so that it also
% computes the number of blocks for matrix W (DONE!)

constValues = calcConstValues(domainMesh, quadData);

blockSizesW = [domainMesh.numVertices, domainMesh.numVertices];
matrixSpecsW = calcMatrixSpecs(nGPU, avMem, blockSizesW, numBlocksW); %Need to modify so it also computes the matrix specs for W
matrixW = calcMatrixW(matrixSpecsW, nGPU, basePath, pbParam, domainMesh, quadData, constValues); %% Need to create this function

% Datum vectors calculations
betaI = calcBetaI(pbParam, domainMesh, constValues, quadData.methodSpecs); %Need to modify this function so it computes the right datum


% Time-marching process
density = zeros(3*domainMesh.numVertices, pbParam.nT);

[L, U, P] = lu(matrixW{1});

for currInd = 1 : pbParam.nT
    rhs = betaI{currInd}; 

    endInd = min(currInd, numBlocksW);
    for indMat = 2 : endInd
        rhs = rhs - matrixW{indMat} * density(:, currInd - indMat + 1);  
    end
  
    density(:, currInd) = U\(L\(P*rhs));
end

%Save on disk
tmpFlag = true;
if(tmpFlag)
    save(tmpPath + "_matrix", 'matrixW', 'betaI');
end
save(outPath + "_density", 'density');
return

end