function matrixSpecs = calcMatrixSpecs(nGPU, avMem, blockSizes, numBlocks)
%CALCMATRIXSPECS  Plan how a block-Toeplitz BEM matrix is split across GPUs and iterations.
%   MATRIXSPECS = CALCMATRIXSPECS(NGPU, AVMEM, BLOCKSIZES, NUMBLOCKS) works out,
%   given the amount of free GPU memory AVMEM and the number of available devices NGPU,
%   how many (3*BLOCKSIZES(1)) x (3*BLOCKSIZES(2)) matrix blocks can be held in
%   memory at once, and derives the per-device offsets needed to split the NUMBLOCKS
%   blocks of a CALCMATRIXV/CALCMATRIXK/CALCMATRIXK_C/CALCMATRIXW computation into a
%   sequence of GPU launches ("iterations"), each as large as memory allows.
%
%   Input arguments:
%       NGPU       - (positive integer) number of available GPU devices.
%       AVMEM      - (positive integer) bytes of free memory on the smallest available GPU.
%       BLOCKSIZES - (1x2 positive integer) [numRows, numCols] size of a matrix block before
%                    the 3x multiplicity from the vector-valued (3D elastodynamic) kernel.
%       NUMBLOCKS  - (positive integer) total number of time-blocks to compute.
%
%   Output arguments:
%       MATRIXSPECS - (struct) with fields blockSizes2D, blockNumRows, blockNumCols,
%                     blockMemSize, numBlocks, maxNumBlocksPerIter, numIter,
%                     offsets_sing (per GPU launch) and offsets_full (per iteration).
%
%   Notes:
%       Asserts if AVMEM cannot hold even a single block
%       (using a safety factor of 2.5x the raw block size).
%
%   See also CALCMATRIXV, CALCMATRIXK, CALCMATRIXK_C, CALCMATRIXW, CALCNUMMATRIXBLOCKS
arguments (Input)
    nGPU        (1, 1) double {mustBeInteger, mustBePositive}
    avMem       (1, 1) double {mustBeInteger, mustBePositive}
    blockSizes  (1, 2) double {mustBeInteger, mustBePositive}
    numBlocks   (1, 1) double {mustBeInteger, mustBePositive}
end

arguments (Output)
    matrixSpecs (1, 1) struct
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