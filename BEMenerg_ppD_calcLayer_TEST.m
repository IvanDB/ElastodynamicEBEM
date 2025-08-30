function [uValues, uValuesG, uValuesP] = BEMenerg_ppD_calcLayer_TEST(pbParam, domainMesh, density, gK, methodInfo, xVal, t, PPn, PPw)
% INPUT
% - 
% - 
% OUTPUT
% - 
%% INIZIALIZZAZIONE PARAMETRI
gpuID = gpuDevice;

xSize = size(xVal, 1);
%Controllo teorico sul valore temporale
if t <= 0
    uValues = cell(xSize, 1);
    parfor i = 1 : xSize
        uValues{i} = zeros(3, 1);
    end
    return
end

%Calcolo PASSO di DISCRETIZZAZIONE TEMPORALE
deltaT = pbParam.Tfin / pbParam.nT;

%Calcolo massimo indice temporale nHat necessario per il calcolo
nHat = ceil(t ./ deltaT);

%Estrazione NUMERO TRIANGOLI DISCRETIZZAZIONE SPAZIALE
numTriang = domainMesh.numberTriangles;

%Estrazione NUMERO VERTICI DISCRETIZZAZIONE SPAZIALE
numNodes = domainMesh.number_nodes;

%Calcolo valori costanti
constValues = BEMenerg_pp_calcCostantData(domainMesh);

%% SETUP VARIABILI PER COMPUTAZIONE GPU
%Array di input
uRawG = zeros(3 * xSize, 1, "gpuArray");
uRawP = zeros(3 * xSize, 1, "gpuArray");

%Vettore delle differenze temporali
diffTemp = gpuArray(t - ((0 : (nHat+1)) .* deltaT));

%Coordinate del punto sorgente
sourcePoints = gpuArray(reshape(xVal', [3*xSize 1]));
gKlinear = gpuArray(reshape(cell2mat(gK(1 : nHat)'), [3*numNodes*nHat 1]));
densityLin = gpuArray(reshape(density(:, 1 : nHat), [3*numTriang*nHat 1]));

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

normT = zeros(3*numTriang, 1, 'gpuArray');
for indTemp = 0 : (numTriang - 1)
    normT(3*indTemp + (1:3), 1) = domainMesh.normal(indTemp+1, :)';
end

indSMmatrix = zeros(numTriang, numNodes, 'int32');
numTonS = zeros(numNodes, 1, 'int32');
indTonS = zeros(numTriang, numNodes, 'int32');
for indS = 1 : numNodes
    for indM = 1 : numTriang
        [~, indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
        if indSMmatrix(indM, indS) ~= 0 
            numTonS(indS) = numTonS(indS) + 1;
            indTonS(numTonS(indS), indS) = indM;
        end
    end
end

indSMmatrixGPU = gpuArray(reshape(indSMmatrix, [numNodes*numTriang, 1]));
numTonS = gpuArray(numTonS);
indTonS = gpuArray(reshape(indTonS, [numNodes*numTriang, 1]));

matCoeff = zeros(9*domainMesh.numberTriangles, 1, 'gpuArray');
vetCoeff = zeros(3*domainMesh.numberTriangles, 1, 'gpuArray');
for indTemp = 1 : numTriang
    matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end
%% SETUP KERNEL MATRICE G
kernelG = parallel.gpu.CUDAKernel("kernelG_layer.ptx", "kernelG_layer.cu");

kernelG.GridSize = [xSize 1];
kernelG.ThreadBlockSize = [nHat 1 1];
kernelG.SharedMemorySize = nHat * 3 * 8;

%% COMPUTAZIONE MATRICE G
wait(gpuID);
time = tic;

uRawG = feval(kernelG, uRawG, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho,  ...
                sourcePoints, diffTemp, densityLin, ...
                methodInfo.numSubRegionPP, methodInfo.numNodiSingPP, stdPPw, stdPPnx, stdPPny, stdPPnz, ...
                numTriang, vertsT, areeT, centerT, maxLenT);
wait(gpuID);
uRawG = reshape(gather(uRawG), [3 xSize]);

disp("G: " + toc(time));
%% SETUP KERNEL MATRICE P
kernelP = parallel.gpu.CUDAKernel("kernelP_layerMOD2.ptx", "kernelP_layerMOD2.cu");

kernelP.GridSize = [xSize 1];
kernelP.ThreadBlockSize = [nHat 1 1];
kernelP.SharedMemorySize = nHat * 3 * 8;

%% COMPUTAZIONE MATRICE P
uRawPtest2 = zeros(3 * xSize, 1, "gpuArray");
wait(gpuID);
time = tic;
uRawPtest2 = feval(kernelP, uRawPtest2, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, 4*pi*pbParam.rho*deltaT, ...
                sourcePoints, diffTemp, gKlinear, ...
                methodInfo.numSubRegionPP, methodInfo.numNodiSingPP, stdPPw, stdPPnx, stdPPny, stdPPnz, ...
                numTriang, vertsT, areeT, normT, centerT, maxLenT, ...
                numNodes, indSMmatrixGPU, matCoeff, vetCoeff, numTonS, indTonS);

wait(gpuID);
uRawPtest2 = reshape(gather(uRawPtest2), [3 xSize]);

disp("Kmod2: " + toc(time));

% %% SETUP KERNEL MATRICE P
% kernelP = parallel.gpu.CUDAKernel("kernelP_layerMOD.ptx", "kernelP_layerMOD.cu");
% 
% kernelP.GridSize = [xSize 1];
% kernelP.ThreadBlockSize = [nHat 1 1];
% kernelP.SharedMemorySize = nHat * 3 * 8;
% 
% %% COMPUTAZIONE MATRICE P
% wait(gpuID);
% time = tic;
% uRawP = feval(kernelP, uRawP, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, 4*pi*pbParam.rho*deltaT, ...
%                 sourcePoints, diffTemp, gKlinear, ...
%                 methodInfo.numSubRegionPP, methodInfo.numNodiSingPP, stdPPw, stdPPnx, stdPPny, stdPPnz, ...
%                 numTriang, vertsT, areeT, normT, centerT, maxLenT, ...
%                 numNodes, indSMmatrixGPU, matCoeff, vetCoeff);
% 
% wait(gpuID);
% uRawP = reshape(gather(uRawP), [3 xSize]);
% 
% disp("Kmod: " + toc(time));
% 
% %% SETUP KERNEL MATRICE P
% kernelP = parallel.gpu.CUDAKernel("kernelP_layer.ptx", "kernelP_layer.cu");
% 
% kernelP.GridSize = [xSize 1];
% kernelP.ThreadBlockSize = [nHat 1 1];
% kernelP.SharedMemorySize = nHat * 3 * 8;
% 
% %% COMPUTAZIONE MATRICE P
% uRawPtest = zeros(3 * xSize, 1, "gpuArray");
% wait(gpuID);
% time = tic;
% uRawPtest = feval(kernelP, uRawPtest, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, 4*pi*pbParam.rho*deltaT, ...
%                 sourcePoints, diffTemp, gKlinear, ...
%                 methodInfo.numSubRegionPP, methodInfo.numNodiSingPP, stdPPw, stdPPnx, stdPPny, stdPPnz, ...
%                 numTriang, vertsT, areeT, normT, centerT, maxLenT, ...
%                 numNodes, indSMmatrixGPU, matCoeff, vetCoeff);
% 
% wait(gpuID);
% uRawPtest = reshape(gather(uRawPtest), [3 xSize]);
% 
% disp("Kstd: " + toc(time));
% 
% disp("error1 " + max(abs(uRawP - uRawPtest), [], "all"));
% disp("error2 " + max(abs(uRawP - uRawPtest2), [], "all"));

%% SOMMA FINALE
uValues = mat2cell(uRawG - uRawPtest2, 3, ones(xSize, 1));
uValuesG = mat2cell(uRawG, 3, ones(xSize, 1));
uValuesP = mat2cell(uRawPtest2, 3, ones(xSize, 1));
disp(" ");
return