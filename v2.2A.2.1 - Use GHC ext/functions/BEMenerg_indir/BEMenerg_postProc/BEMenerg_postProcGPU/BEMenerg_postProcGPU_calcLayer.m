function uValues = BEMenerg_postProcGPU_calcLayer(pbParam, domainMesh, density, methodInfo, xVal, t, PPn, PPw)
% INPUT
% - 
% - 
% OUTPUT
% - 
%% INIZIALIZZAZIONE PARAMETRI

xSize = size(xVal, 1);

%Controllo teorico sul valore temporale
if t <= 0
    uValues = cell(xSize, 1);
    for i = 1 : xSize
        uValues{i} = zeros(3, 1);
    end
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
uRaw = zeros(3 * xSize, 1, "gpuArray");

%Vettore delle differenze temporali
diffTemp = gpuArray(t - ((0 : nHat) .* deltaT));

%Coordinate dei punti sorgente
sourcePoints = gpuArray(reshape(xVal', [3*xSize 1]));

%Nodi e pesi GHC integrazione
stdPPw = gpuArray(PPw);
stdPPnx = gpuArray(PPn(:, 1));
stdPPny = gpuArray(PPn(:, 2));
stdPPnz = gpuArray(PPn(:, 3));

%Vertici e aree dei triangoli della mesh spaziale
vertsT = zeros(9*numTriang, 1, "gpuArray");
for indTemp = 0 : (numTriang - 1)
    vertsT(9*indTemp + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp+1, 1:3), :), [9 1]);
end
areeT = gpuArray(domainMesh.area);

%Centri e maxLengths dei triangoli della mesh spaziale
centerT = gpuArray(reshape(domainMesh.center', [3*numTriang 1]));
maxLenT = gpuArray(domainMesh.maxL);

densityLin = gpuArray(reshape(density(:, 1:numIst)', [3*numTriang*numIst 1]));

%% SETUP KERNEL GPU
kernel = parallel.gpu.CUDAKernel("kernelLayer.ptx", "kernelLayer.cu");

kernel.GridSize = [xSize 1 1];

kernel.ThreadBlockSize = [numIst 1 1];
kernel.SharedMemorySize = numIst * 3 * 8;

%% AVVIO COMPUTAZIONE ASINCRONA GPU
wait(gpuDevice);
uRaw = feval(kernel, uRaw, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho,  ...
                    sourcePoints, diffTemp, ...
                    numTriang, methodInfo.numSubRegionPP, methodInfo.numNodiSingPP, ...
                    stdPPw, stdPPnx, stdPPny, stdPPnz, ...
                    vertsT, areeT, centerT, maxLenT, densityLin);
wait(gpuDevice);
uRaw = gather(uRaw);

uRaw = reshape(uRaw, [3, xSize]);
uValues = mat2cell(uRaw, 3, ones(xSize, 1));
end