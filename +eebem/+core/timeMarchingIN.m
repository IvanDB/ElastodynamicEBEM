function density = timeMarchingIN(basePath, pbParam, domainMesh, quadData, fullFileNames)
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