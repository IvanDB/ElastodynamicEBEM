addpath(genpath("./functions"));

%Avvio pool parallela
delete(gcp("nocreate"));
parInfo = parpool("Threads");
clc

pbParam = readInputFile("input_barH1_small.txt");
domainMesh = readSpaceMesh(pbParam.domainType, pbParam.lev);

[methodInfo, GHn, GHw, GHCn, GHCw, DIAGn, DIAGw] = BEMenerg_setupCore(convertStringsToChars("FN 19 16 3 64"));


gpuID = gpuDevice;
memAv = gpuID.AvailableMemory;
numTriang = domainMesh.numberTriangles;
numNodes = domainMesh.number_nodes;

constValues = BEMenerg_dir_calcCostantData(methodInfo, domainMesh, GHn, GHw, GHCn, GHCw, DIAGn, DIAGw);

deltaT = pbParam.Tfin / pbParam.nT;

%Nodi e pesi GH integrazione esterna
stdGHw = gpuArray(GHw);
stdGHnx = gpuArray(GHn(:, 1));
stdGHny = gpuArray(GHn(:, 2));
stdGHnz = gpuArray(GHn(:, 3));

%Nodi e pesi GHC integrazione interna
stdGHCw = gpuArray(GHCw);
stdGHCnx = gpuArray(GHCn(:, 1));
stdGHCny = gpuArray(GHCn(:, 2));
stdGHCnz = gpuArray(GHCn(:, 3));

%Vertici e aree dei triangoli della mesh spaziale
vertsT = zeros(9*numTriang, 1, 'gpuArray');
for indTemp = 0 : (numTriang - 1)
    vertsT(9*indTemp + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp+1, 1:3), :), [9 1]);
end
areeT = gpuArray(domainMesh.area);

normT = zeros(3*numTriang, 1, 'gpuArray');
for indTemp = 0 : (numTriang - 1)
    normT(3*indTemp + (1:3), 1) = domainMesh.normal(indTemp+1, :)';
end

indSMmatrix = zeros(numTriang, numNodes, 'int32');
for indS = 1 : numNodes
    for indM = 1 : numTriang
        [~, indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
    end
end
indSMmatrixGPU = gpuArray(reshape(indSMmatrix, [numNodes*numTriang, 1]));

matCoeff = zeros(9*domainMesh.numberTriangles, 1, 'gpuArray');
vetCoeff = zeros(3*domainMesh.numberTriangles, 1, 'gpuArray');
for indTemp = 1 : numTriang
    matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end

nodesMesh = gpuArray(reshape(domainMesh.coordinates', [3*numNodes 1]));
maxLen = max(domainMesh.maxL);



numBlocksK = min(BEMenerg_dir_calcNumBlocks(pbParam, domainMesh) + 1, pbParam.nT);
blockSizeK = 81 * numTriang * numTriang * 8;

numBlocksPerIterK = numBlocksK;

numIterK = ceil(numBlocksK / numBlocksPerIterK);
offSetsK = numBlocksPerIterK * (0 : (numIterK-1));

matrixSavedK = cell(numBlocksPerIterK, numIterK);


%% INTER

% Setup kernel interno
kernelKint = parallel.gpu.CUDAKernel("kernelKlinslacInternal.ptx", "kernelKlinslacInternal.cu");
kernelKint.GridSize = [numTriang numTriang numBlocksPerIterK];
blockXint = methodInfo.numSubRegion;
blockYint = methodInfo.numNodiSing;
kernelKint.ThreadBlockSize = [blockXint blockYint 1];
kernelKint.SharedMemorySize = blockXint * blockYint * 81 * 8;

%Avvio computazione GPU
matrixOutInt = BEMenerg_dir_launchKernelK_uConstSpace(gpuID, kernelKint, deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*deltaT, ...
                                                stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, ...
                                                vertsT, areeT, normT, offSetsK(1), numBlocksK, numTriang, numBlocksPerIterK, maxLen, matCoeff, vetCoeff);
%Attesa completamento operazioni GPU
matrixOutInt = gather(matrixOutInt);
wait(gpuID);

%% BOUND

% Setup kernel bordo
kernelKbound = parallel.gpu.CUDAKernel("kernelKlinslacBoundary.ptx", "kernelKlinslacBoundary.cu");
kernelKbound.GridSize = [numTriang numTriang numBlocksPerIterK];
methodInfo.blockXbound = 64; %mettere come parametro in methodInfo fuori;
kernelKbound.ThreadBlockSize = [methodInfo.blockXbound 1 1];
kernelKbound.SharedMemorySize = methodInfo.blockXbound * 81 * 8;


matrixOutBound = BEMenerg_dir_launchKernelK_uConstSpace(gpuID, kernelKbound, deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*deltaT, ...
                                                stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, ...
                                                vertsT, areeT, normT, offSetsK(1), numBlocksK, numTriang, numBlocksPerIterK, maxLen, matCoeff, vetCoeff);

%Attesa completamento operazioni GPU
matrixOutBound = gather(matrixOutBound);
wait(gpuID);

 