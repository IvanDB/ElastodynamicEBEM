function uXT = calcMatrixP(pbParam, domainMesh, density, methodInfo, x, t, PPn, PPw)
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
numS = domainMesh.numVertices;

%Constant data (cell array containing vet coeff and mat coeff for eac triangle)
constValues = calcConstData(domainMesh);

%% SETUP GPU VARIABLES
% input arrays
matrixP = zeros(3 * nHat * 3 * numS, 1, "gpuArray");

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


% normal unict vectors
normT = zeros(3*numT, 1, 'gpuArray');
for indTemp = 0 : (numT - 1)
    normT(3*indTemp + (1:3), 1) = domainMesh.normal(indTemp+1, :)';
end

% mesh topology matrix

indSMmatrix = zeros(numT, numS, 'int32');
for indS = 1 : numS
    for indM = 1 : numT
        [~, indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
    end
end
indSMmatrixGPU = gpuArray(reshape(indSMmatrix, [numS*numT, 1]));

matCoeff = zeros(9*domainMesh.numTriangles, 1, 'gpuArray');
vetCoeff = zeros(3*domainMesh.numTriangles, 1, 'gpuArray');
for indTemp = 1 : numT
    matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end

%% Matrix P kernel Setup
srcPath = fullfile(".", "+eebem", "+core", "kernelsCUDA", "kernelP.cu");
ptxPath = fullfile(".", "buildDir", "kernelP.ptx");

kernelP = parallel.gpu.CUDAKernel(ptxPath, srcPath);

kernelP.GridSize = [numS nHat 1];
blockX = methodInfo.numSubRegionPP;
kernelP.ThreadBlockSize = [blockX 1 1];
kernelP.SharedMemorySize = blockX * 9 * 8;

%% COMPUTAZIONE MATRICE P
wait(gpuID);
matrixP = feval(kernelP, matrixP, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, 4*pi*pbParam.rho*deltaT, numT, ...
                    sourcePoint, diffTemp, ...
                    stdPPw, stdPPnx, stdPPny, stdPPnz, numNodesPerThread, ...
                    vertsT, areasT, normT, indSMmatrixGPU, matCoeff, vetCoeff);

wait(gpuID);
matrixP = reshape(matrixP, [3*nHat 3*numS]);

%% SOMMA PER uXT
uXTtemp = zeros(3, nHat, "gpuArray");
for indTemp = 1 : nHat
    matrixPSlice = matrixP(3*(indTemp-1) + (1:3), :);

    uXTtemp(:, indTemp) = - matrixPSlice * density(:, indTemp);
end
uXT = sum(uXTtemp, 2);
uXT = gather(uXT);
return