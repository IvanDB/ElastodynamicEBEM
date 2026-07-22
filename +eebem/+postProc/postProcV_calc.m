function uXT = postProcV_calc(pbParam, domainMesh, density, methodInfo, x, t, PPn, PPw)
import eebem.postProc.*
%% GPU SETUP 

gpuID = gpuDevice;

% no negative time
if t <= 0
    uXT = zeros(3, 1);
    return
end

% Delta t
deltaT = pbParam.Tfin / pbParam.nT;

% n hat index: after this all matrix P blocks are null 
nHat = ceil(t ./ deltaT);

% mesh info
numT = domainMesh.numTriangles;


%Constant data (cell array containing vet coeff and mat coeff for eac triangle)
constValues = postProc_constData(domainMesh);

%% SETUP GPU VARIABLES
% input arrays
matrixG = zeros(3 * nHat * 3 * numT, 1, "gpuArray");

% temporal difference vector (\Delta^t_0, \Delta^t_1, ..., \Delta^t_{\hat{n}_t+1}) 
%(i have all the time instants needed to compute matrix P in this vector)
diffTemp = gpuArray(t - ((0 : (nHat+1)) .* deltaT));

% source point coordinates
sourcePoint = gpuArray(x);

% quadrature nodes and weights
stdPPw = gpuArray(PPw);
stdPPnx = gpuArray(PPn(:, 1));
stdPPny = gpuArray(PPn(:, 2));
stdPPnz = gpuArray(PPn(:, 3));


numNodesPerThread = methodInfo.numNodeSingPP;

%mesh triangles vertexes and areas
vertsT = zeros(9*numT, 1, "gpuArray");
for indTemp = 0 : (numT - 1)
    vertsT(9*indTemp + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp+1, 1:3), :), [9 1]);
end
areasT = gpuArray(domainMesh.area);

%Centri e maxLengths dei triangoli della mesh spaziale
centerT = gpuArray(reshape(domainMesh.center', [3*numT 1]));
maxLenT = gpuArray(domainMesh.maxL);

% normal unict vectors
normT = zeros(3*numT, 1, 'gpuArray');
for indTemp = 0 : (numT - 1)
    normT(3*indTemp + (1:3), 1) = domainMesh.normal(indTemp+1, :)';
end


matCoeff = zeros(9*domainMesh.numTriangles, 1, 'gpuArray');
vetCoeff = zeros(3*domainMesh.numTriangles, 1, 'gpuArray');
for indTemp = 1 : numT
    matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end

%% Matrix G kernel Setup
currentFolder = fileparts(mfilename('fullpath'));
ptxPath = fullfile(currentFolder, "kernelG.ptx");
cuPath = fullfile(currentFolder, "kernelG.cu");

kernelG = parallel.gpu.CUDAKernel(ptxPath, cuPath);

kernelG.GridSize = [numT nHat 1];
blockX = methodInfo.numSubRegionPP;
kernelG.ThreadBlockSize = [blockX 1 1];
kernelG.SharedMemorySize = blockX * 9 * 8;

%% COMPUTAZIONE MATRICE G
wait(gpuID);
matrixG = feval(kernelG, matrixG, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho,  ...
                    sourcePoint, diffTemp, ...
                    stdPPw, stdPPnx, stdPPny, stdPPnz, numNodesPerThread, ...
                    vertsT, areasT, centerT, maxLenT);
wait(gpuID);
matrixG = reshape(matrixG, [3*nHat 3*numT]);
%% SOMMA PER uXT
uXTtemp = zeros(3, nHat, "gpuArray");
for indTemp = 1 : nHat
    matrixGSlice = matrixG(3*(indTemp-1) + (1:3), :);

    uXTtemp(:, indTemp) = matrixGSlice * density(:, indTemp);
end
uXT = sum(uXTtemp, 2);
uXT = gather(uXT);
end