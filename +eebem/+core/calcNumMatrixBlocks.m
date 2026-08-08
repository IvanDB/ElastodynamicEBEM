function [numBlocksV, numBlocksK, numBlocksW] = calcNumMatrixBlocks(pbParam, domainMesh)
%CALCNUMMATRIXBLOCKS  Compute how many time-blocks each BEM operator needs, from the mesh diameter and wave speeds.
%   [NUMBLOCKSV, NUMBLOCKSK, NUMBLOCKSW] = CALCNUMMATRIXBLOCKS(PBPARAM, DOMAINMESH)
%   exploits the finite speed of propagation of the elastodynamic kernels:
%   a block-Toeplitz matrix entry at time-lag k*deltaT is zero once k*velS*deltaT
%   exceeds the mesh diameter (plus a small safety margin), so only the first
%   NUMBLOCKSV/K/W blocks (out of the PBPARAM.nT that a naive implementation would need)
%   are ever non-zero and have to be actually assembled by CALCMATRIXV/K/W.
%
%   Input arguments:
%       PBPARAM    - (struct) physical/time-discretization parameters
%                    (velS, deltaT, nT), see READINPUTFILE.
%       DOMAINMESH - (struct) triangulated boundary mesh (coordinates), see READSPACEMESH.
%
%   Output arguments:
%       NUMBLOCKSV - (positive integer) blocks needed for the single-layer operator V.
%       NUMBLOCKSK - (positive integer) blocks needed for the
%                    double-layer operator K (one more than V).
%       NUMBLOCKSW - (positive integer) blocks needed for the
%                    hypersingular operator W (equal to NUMBLOCKSK).
%
%   Notes:
%       The all-pairs vertex distance matrix is computed with a PARFOR loop and
%       has O(numVertices^2) cost; this can dominate runtime on very fine meshes.
%
%   See also CALCMATRIXSPECS, CALCMATRIXV, CALCMATRIXK, CALCMATRIXW
arguments (Input)
    pbParam    (1, 1) struct
    domainMesh (1, 1) struct
end

arguments (Output)
    numBlocksV (1, 1) double {mustBeInteger, mustBePositive}
    numBlocksK (1, 1) double {mustBeInteger, mustBePositive}
    numBlocksW (1, 1) double {mustBeInteger, mustBePositive}
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