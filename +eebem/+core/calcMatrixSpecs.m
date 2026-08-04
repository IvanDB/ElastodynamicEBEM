function matrixSpecs = calcMatrixSpecs(nGPU, avMem, blockSizes, numBlocks)
%CALCOFFSETS Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    nGPU        (1, 1) double {mustBeInteger, mustBePositive}
    avMem       (1, 1) double {mustBeInteger, mustBePositive}
    blockSizes  (1, 2) double {mustBeInteger, mustBePositive}
    numBlocks   (1, 1) double {mustBeInteger, mustBePositive}
end

matrixSpecs.blockSizes2D = blockSizes;
matrixSpecs.blockNumRows = 3 * blockSizes(1);
matrixSpecs.blockNumCols = 3 * blockSizes(2);
matrixSpecs.blockMemSize = 9 * prod(blockSizes) * 8;

GPUMemorySafetyFactor = 2.5;
maxNumBlocksInMemory = floor(avMem / (GPUMemorySafetyFactor * matrixSpecs.blockMemSize));
assert(maxNumBlocksInMemory > 0, "Insufficiente GPU memory for the required problem.")

maxNumBlocksPerIter = maxNumBlocksInMemory * nGPU; % maximum number of blocks that can be computed in each iteration


numIter_ful = floor(numBlocks / maxNumBlocksPerIter); % Number of "complete" iterations
numIter_tot = ceil(numBlocks / maxNumBlocksPerIter); % Number of total iterations, including the last one, which might not be complete
lastIterSize = mod(numBlocks, maxNumBlocksPerIter); % Number of blocks for the last iteration


blkSize_F = floor(lastIterSize / nGPU);
blkSize_I = ceil(lastIterSize / nGPU);

nBlocks_F = mod(-lastIterSize, nGPU);
nBlocks_I = (nGPU - nBlocks_F) * (lastIterSize ~= 0);


offsetsSingleDevices = 0 : maxNumBlocksInMemory : (maxNumBlocksPerIter * numIter_ful);

offsets_bf = offsetsSingleDevices(end) + blkSize_I * (1 : nBlocks_I);
offsetsSingleDevices = [offsetsSingleDevices, offsets_bf];

offsets_ef = offsetsSingleDevices(end) + blkSize_F * (1 : nBlocks_F);
offsetsSingleDevices = [offsetsSingleDevices, offsets_ef];

matrixSpecs.numBlocks = numBlocks;
matrixSpecs.maxNumBlocksPerIter = maxNumBlocksPerIter;
matrixSpecs.numIter = numIter_tot;

matrixSpecs.offsets_sing = offsetsSingleDevices;
matrixSpecs.offsets_full = offsetsSingleDevices(1 : nGPU : end);
end