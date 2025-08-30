function uXT = BEMenerg_postProcGPU_calcSP(pbParam, domainMesh, density, methodInfo, x, t, PPn, PPw)
% INPUT
% - 
% - 
% OUTPUT
% - 
%% INIZIALIZZAZIONE PARAMETRI
%Controllo teorico sul valore temporale
if t <= 0
    uXT = zeros(3, 1);
    return
end

%Calcolo PASSO di DISCRETIZZAZIONE TEMPORALE
deltaT = pbParam.Tfin / pbParam.nT;

%Calcolo massimo indice temporale nHat necessario per il calcolo
nHat = ceil(t ./ deltaT) - 1;

%Numero istanti temporali richiesi
numIst = nHat + 1;

%Estrazione NUMERO TRIANGOLI DISCRETIZZAZIONE SPAZIALE
numTriang = domainMesh.numberTriangles;


%% SETUP VARIABILI PER COMPUTAZIONE GPU
%Array di input
matrixTotal = zeros(3 * numIst * 3 * numTriang, 1, "gpuArray");

%Vettore delle differenze temporali
diffTemp = gpuArray(t - ((0 : nHat) .* deltaT));

%Coordinate del punto sorgente
sourcePoint = x;

%Nodi e pesi GHC integrazione
stdPPw = gpuArray(PPw);
stdPPnx = gpuArray(PPn(:, 1));
stdPPny = gpuArray(PPn(:, 2));
stdPPnz = gpuArray(PPn(:, 3));

%Parametro di ripetizione nodi interni
numNodiPerThread = methodInfo.numNodiSingPP;

%Vertici e aree dei triangoli della mesh spaziale
vertsT = zeros(9*numTriang, 1, "gpuArray");
for indTemp = 0 : (numTriang - 1)
    vertsT(9*indTemp + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp+1, 1:3), :), [9 1]);
end
areeT = gpuArray(domainMesh.area);

%Centri e maxLengths dei triangoli della mesh spaziale
centerT = gpuArray(reshape(domainMesh.center', [3*numTriang 1]));
maxLenT = gpuArray(domainMesh.maxL);

%% SETUP KERNEL GPU
kernel = parallel.gpu.CUDAKernel("kernelSP.ptx", "kernelSP.cu");

kernel.GridSize = [numTriang numIst 1];

blockX = methodInfo.numSubRegionPP;

kernel.ThreadBlockSize = [blockX 1 1];
kernel.SharedMemorySize = blockX * 9 * 8;

%% AVVIO COMPUTAZIONE ASINCRONA GPU
wait(gpuDevice);
matrixTotal = feval(kernel, matrixTotal, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho,  ...
                    sourcePoint, diffTemp, ...
                    stdPPw, stdPPnx, stdPPny, stdPPnz, numNodiPerThread, ...
                    vertsT, areeT, centerT, maxLenT);

wait(gpuDevice);
matrixTotal = reshape(matrixTotal, [3*numIst 3*numTriang]);

%% SOMMA PER uXT
uXTtemp = zeros(3, nHat, "gpuArray");
for indTemp = 1 : nHat
    matrixSlicePos = matrixTotal(3*(indTemp-1) + (1:3), :);
    matrixSliceNeg = matrixTotal(3*indTemp + (1:3), :);
    uXTtemp(:, indTemp) = (matrixSlicePos - matrixSliceNeg) * density(:, indTemp);
end

uXT = sum(uXTtemp, 2) + (matrixTotal(3*nHat + (1:3), :) * density(:, nHat+1));
uXT = gather(uXT);
return