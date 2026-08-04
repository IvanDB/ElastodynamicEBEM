function [numBlocksV, numBlocksK, numBlocksW] = calcNumMatrixBlocks(pbParam, domainMesh)
%CALCNUMMATRIXBLOCKS Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    pbParam    (1, 1) struct
    domainMesh (1, 1) struct
end

arguments (Output)
    numBlocksV (1, 1) {mustBeInteger, mustBePositive}
    numBlocksK (1, 1) {mustBeInteger, mustBePositive}
    numBlocksW (1, 1) {mustBeInteger, mustBePositive}
end

numV = domainMesh.numVertices;
deltaT = pbParam.deltaT;

dists = zeros(numV, numV);
parfor indExt = 1 : numV
    for indInt = 1 : numV
        vettDistC = domainMesh.coordinates(indExt, :) - domainMesh.coordinates(indInt, :);
        dists(indExt, indInt) = sqrt(sum(vettDistC.^2));
    end
end

maxDist = max(dists, [], "all");

maxIndTempV = ceil(maxDist / (pbParam.velS * deltaT)) + 1;
maxIndTempK = ceil(maxDist / (pbParam.velS * deltaT)) + 2;

numBlocksV = min(pbParam.nT, maxIndTempV);
numBlocksK = min(pbParam.nT, maxIndTempK);
numBlocksW = numBlocksK;
end