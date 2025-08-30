function uXT = BEMenerg_ppN_calc(pbParam, domainMesh, density, gV, methodInfo, x, t, PPn, PPw)
% INPUT
% - 
% - 
% OUTPUT
% - 
%% INIZIALIZZAZIONE PARAMETRI
gpuID = gpuDevice;
%Controllo teorico sul valore temporale
if t <= 0
    uXT = zeros(3, 1);
    return
end

%Calcolo PASSO di DISCRETIZZAZIONE TEMPORALE
deltaT = pbParam.Tfin / pbParam.nT;

%Calcolo massimo indice temporale nHat necessario per il calcolo
nHat = ceil(t ./ deltaT);

%Estrazione NUMERO TRIANGOLI DISCRETIZZAZIONE SPAZIALE
numT = domainMesh.numberTriangles;

%Estrazione NUMERO VERTICI DISCRETIZZAZIONE SPAZIALE
numS = domainMesh.number_nodes;

%Calcolo valori costanti
constValues = BEMenerg_ppN_calcCostantData(domainMesh);

%% SETUP VARIABILI PER COMPUTAZIONE GPU
%Array di input
matrixG = zeros(3 * nHat * 3 * numT, 1, "gpuArray");
matrixP = zeros(3 * nHat * 3 * numS, 1, "gpuArray");

%Vettore delle differenze temporali
diffTemp = gpuArray(t - ((0 : (nHat+1)) .* deltaT));

%Coordinate del punto sorgente
sourcePoint = gpuArray(x);

%Nodi e pesi GHC integrazione
stdPPw = gpuArray(PPw);
stdPPnx = gpuArray(PPn(:, 1));
stdPPny = gpuArray(PPn(:, 2));
stdPPnz = gpuArray(PPn(:, 3));

%Parametro di ripetizione nodi interni
numNodiPerThread = methodInfo.numNodiSingPP;

%Vertici e aree dei triangoli della mesh spaziale
vertsT = zeros(9*numT, 1, "gpuArray");
for indTemp = 0 : (numT - 1)
    vertsT(9*indTemp + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp+1, 1:3), :), [9 1]);
end
areeT = gpuArray(domainMesh.area);

%Centri e maxLengths dei triangoli della mesh spaziale
centerT = gpuArray(reshape(domainMesh.center', [3*numT 1]));
maxLenT = gpuArray(domainMesh.maxL);

normT = zeros(3*numT, 1, 'gpuArray');
for indTemp = 0 : (numT - 1)
    normT(3*indTemp + (1:3), 1) = domainMesh.normal(indTemp+1, :)';
end

indSMmatrix = zeros(numT, numS, 'int32');
for indS = 1 : numS
    for indM = 1 : numT
        [~, indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
    end
end
indSMmatrixGPU = gpuArray(reshape(indSMmatrix, [numS*numT, 1]));

matCoeff = zeros(9*domainMesh.numberTriangles, 1, 'gpuArray');
vetCoeff = zeros(3*domainMesh.numberTriangles, 1, 'gpuArray');
for indTemp = 1 : numT
    matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end
%% SETUP KERNEL MATRICE G
kernelG = parallel.gpu.CUDAKernel("kernelNEWG.ptx", "kernelNEWG.cu");

kernelG.GridSize = [numT nHat 1];

blockX = methodInfo.numSubRegionPP;

kernelG.ThreadBlockSize = [blockX 1 1];
kernelG.SharedMemorySize = blockX * 9 * 8;

%% AVVIO COMPUTAZIONE G
wait(gpuID);
matrixG = feval(kernelG, matrixG, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho,  ...
                    sourcePoint, diffTemp, ...
                    stdPPw, stdPPnx, stdPPny, stdPPnz, numNodiPerThread, ...
                    vertsT, areeT, centerT, maxLenT);
wait(gpuID);
matrixG = reshape(matrixG, [3*nHat 3*numT]);

%% SETUP KERNEL MATRICE P
kernelP = parallel.gpu.CUDAKernel("kernelNEWP.ptx", "kernelNEWP.cu");
% kernelP = parallel.gpu.CUDAKernel("kernelNEWPusingK.ptx", "kernelNEWPusingK.cu");

kernelP.GridSize = [numS nHat 1];

blockX = methodInfo.numSubRegionPP;

kernelP.ThreadBlockSize = [blockX 1 1];
kernelP.SharedMemorySize = blockX * 9 * 8;

%% AVVIO COMPUTAZIONE P
wait(gpuID);
matrixP = feval(kernelP, matrixP, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, 4*pi*pbParam.rho*deltaT, numT, ...
                    sourcePoint, diffTemp, ...
                    stdPPw, stdPPnx, stdPPny, stdPPnz, numNodiPerThread, ...
                    vertsT, areeT, normT, indSMmatrixGPU, matCoeff, vetCoeff);
% matrixP = feval(kernelP, matrixP, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*deltaT, numT, ...
%                     sourcePoint, diffTemp, ...
%                     stdPPw, stdPPnx, stdPPny, stdPPnz, numNodiPerThread, ...
%                     vertsT, areeT, normT, indSMmatrixGPU, matCoeff, vetCoeff); 

wait(gpuID);
matrixP = reshape(matrixP, [3*nHat 3*numS]);

%% SOMMA PER uXT
uXTtemp = zeros(3, nHat, "gpuArray");
for indTemp = 1 : nHat
    matrixGSlice = matrixG(3*(indTemp-1) + (1:3), :);
    matrixPSlice = matrixP(3*(indTemp-1) + (1:3), :);

    uXTtemp(:, indTemp) = matrixGSlice * gV{indTemp} - matrixPSlice * density(:, indTemp);% 
end
uXT = sum(uXTtemp, 2);
uXT = gather(uXT);
return