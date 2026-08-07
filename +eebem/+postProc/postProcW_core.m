function solution = postProcW_core(basePath, pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw)
import eebem.postProc.*
% Compute dimentions of solution matrix
xSize = size(xVal, 1);
tSize = length(tVal);
solRaw = cell(xSize * tSize, 1); % "raw" solution initialization


%% Building Kernel and mesh info
currentFolder = fileparts(mfilename('fullpath'));
ptxPath = fullfile(currentFolder, "kernelP.ptx");
cuPath = fullfile(currentFolder, "kernelP.cu");
% kernelP = parallel.gpu.CUDAKernel(ptxPath, cuPath);
kernelPConst = parallel.pool.Constant(@() parallel.gpu.CUDAKernel(ptxPath, cuPath));
constValues = postProc_constData(domainMesh);
indSMmatrix = zeros(domainMesh.numTriangles, domainMesh.numVertices, 'int32');
for indS = 1 : domainMesh.numVertices
    for indM = 1 : domainMesh.numTriangles
        [~, indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
    end
end
indSMmatrixGPU = gpuArray(reshape(indSMmatrix, [domainMesh.numVertices*domainMesh.numTriangles, 1]));
indSMBytes = 4 * domainMesh.numTriangles * domainMesh.numVertices; 

vertsT = zeros(9*domainMesh.numTriangles, 1, "gpuArray");
for indTemp = 0 : (domainMesh.numTriangles - 1)
    vertsT(9*indTemp + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp+1, 1:3), :), [9 1]);
end
vertsTBytes = 8 * 9 * domainMesh.numTriangles;

areasT = gpuArray(domainMesh.area);
areasTBytes = 8 * domainMesh.numTriangles;

normT = zeros(3*domainMesh.numTriangles, 1, 'gpuArray');
for indTemp = 0 : (domainMesh.numTriangles - 1)
    normT(3*indTemp + (1:3), 1) = domainMesh.normal(indTemp+1, :)';
end

normTBytes = 8 * 3 * domainMesh.numTriangles;


matCoeff = zeros(9*domainMesh.numTriangles, 1, 'gpuArray');
vetCoeff = zeros(3*domainMesh.numTriangles, 1, 'gpuArray');
for indTemp = 1 : domainMesh.numTriangles
    matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end

matCoeffBytes = 8 * 9 * domainMesh.numTriangles;   
vetCoeffBytes = 8 * 3 * domainMesh.numTriangles;

stdPPw = gpuArray(PPw);
stdPPnx = gpuArray(PPn(:, 1));
stdPPny = gpuArray(PPn(:, 2));
stdPPnz = gpuArray(PPn(:, 3));

numNodesPP = size(PPn, 1);   % = postProcInfo.numSubRegionPP * postProcInfo.numNodeSingPP
stdPPwBytes  = 8 * numNodesPP;   % stdPPw
stdPPnxBytes = 8 * numNodesPP;   % stdPPnx
stdPPnyBytes = 8 * numNodesPP;   % stdPPny
stdPPnzBytes = 8 * numNodesPP;

%% GPU memory check
gpuInfo = gpuDevice;
fprintf("GPU: %s | TotalMemory=%.2f GB | AvailableMemory(baseline)=%.2f GB | numWorkers=%d\n", ...
    gpuInfo.Name, gpuInfo.TotalMemory/1e9, gpuInfo.AvailableMemory/1e9, gcp().NumWorkers);

%% Memory calcs
constMemPerWorker = indSMBytes + vertsTBytes + areasTBytes + normTBytes ...
                   + matCoeffBytes + vetCoeffBytes ...
                   + stdPPwBytes + stdPPnxBytes + stdPPnyBytes + stdPPnzBytes;
numWorkers = gcp().NumWorkers;
margin = 0.35; % Extra margin for memory
memPerIter = gpuInfo.TotalMemory*(1-margin) - numWorkers*constMemPerWorker;
totalPointBytes = 8 * (9 * pbParam.nT * domainMesh.numVertices + (pbParam.nT + 2) + 3 * pbParam.nT + 3); %number of bytes occupied by P(x,t), u(x,t) (x,t)
maxConcurrent = floor(memPerIter / totalPointBytes);
batchSize = min(maxConcurrent, numWorkers); % Number of points per parallel iteration
numBatches = ceil(numPoints / batchSize);
fprintf("Total number of GPU iterations: %.1f ; Number of points per iteration: %.1f \n", numBatches, batchSize);
%% Computing solution
for b = 1 : numBatches
    idxRange = (b-1)*batchSize + 1 : min(b*batchSize, numPoints);
    parfor ind = idxRange
        [indX, indT] = ind2sub([xSize, tSize], ind);
        solRaw{ind} = postProcW_calc(pbParam, domainMesh, density, methodInfo, xVal(indX, :), tVal(indT), stdPPnx, stdPPny,stdPPnz, stdPPw, kernelPConst.Value, matCoeff, vetCoeff, indSMmatrixGPU, vertsT, areasT, normT);
    end
    fprintf("Iteration %.1f completed \n", b);
end

solution = reshape(solRaw, [xSize, tSize]);


%Salvataggio campo vettoriale calcolato
outPath = fullfile(basePath, "solutionData", "IN_" + pbParam.domainType + pbParam.lev + "FN_16-1_64-3_256_256");
save(outPath + "_solution", 'solution');


return