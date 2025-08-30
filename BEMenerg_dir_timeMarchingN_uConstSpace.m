function density = BEMenerg_dir_timeMarchingN_uConstSpace(pbParam, domainMesh, methodInfo, GHn, GHw, GHCn, GHCw, DIAGn, DIAGw)
%

%% SETUP INIZIALE
%Inizializzazione GPU
gpuID = gpuDevice;
reset(gpuID);
memAv = gpuID.AvailableMemory;

%Dimensioni problema
numBlocksV = BEMenerg_dir_calcNumBlocks(pbParam, domainMesh);
numTriang = domainMesh.numberTriangles;
numNodes = domainMesh.number_nodes;

blockSizeV = 81 * (numTriang .^ 2) * 8;   %81 * Mx^2 double (8 byte l'uno)

%Calcolo massimo numero di iterazioni parallelizzabili sulla GPU
fattoreCaricoMemoriaGPU = 2.5;                    %Considero un fattore 1.5 per tenere conto dei dati di input        
maxBlockInMemoryV = floor(memAv / (fattoreCaricoMemoriaGPU*blockSizeV));        
if maxBlockInMemoryV == 0
    error("Memoria GPU insufficiente")
end

%Calcolo numero blocchi per singola iterazione
numBlocksPerIterV = min(maxBlockInMemoryV, numBlocksV);
%Calcolo numero iterazioni
numIterV = ceil(numBlocksV / numBlocksPerIterV);
%Calcolo offsets
offSetsV = numBlocksPerIterV * (0 : (numIterV-1));

%Path per il salvataggio delle matrici e della soluzione
% tmpPath = strcat("./tempData/", pbParam.domainType, num2str(pbParam.lev), "_GPU_", ...
%    num2str(methodInfo.numNodiExt), "_", num2str(methodInfo.numSubRegion), "_", num2str(methodInfo.numNodiSing));
outPath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_NEW_", ...
    num2str(methodInfo.numNodiExt), "_", num2str(methodInfo.numSubRegion), "_", num2str(methodInfo.numNodiSing), "_", ...
    num2str(methodInfo.numNodiDiag));


%% SETUP VARIABILI LOCALI PER COMPUTAZIONE GPU

%Calcolo valori costanti
constValues = BEMenerg_dir_calcCostantData(methodInfo, domainMesh, GHn, GHw, GHCn, GHCw, DIAGn, DIAGw);

%Parametro temporale
deltaT = pbParam.Tfin / pbParam.nT;

%Nodi e pesi GH integrazione esterna
stdGHw = gpuArray(GHw);
stdGHnx = gpuArray(GHn(:, 1));
stdGHny = gpuArray(GHn(:, 2));
stdGHnz = gpuArray(GHn(:, 3));

%Nodi e pesi GHC integrazione interna
stdGHCw = gpuArray(GHCw);
stdGHCnx = gpuArray(GHCn(:, 1));
stdGHCny = gpuArray(GHCn(:, 2));
stdGHCnz = gpuArray(GHCn(:, 3));

%Vertici e aree dei triangoli della mesh spaziale
vertsT = zeros(9*numTriang, 1, 'gpuArray');
for indTemp = 0 : (numTriang - 1)
    vertsT(9*indTemp + (1:9), 1) = reshape(domainMesh.coordinates(domainMesh.triangles(indTemp+1, 1:3), :), [9 1]);
end
areeT = gpuArray(domainMesh.area);

normT = zeros(3*numTriang, 1, 'gpuArray');
for indTemp = 0 : (numTriang - 1)
    normT(3*indTemp + (1:3), 1) = domainMesh.normal(indTemp+1, :)';
end

indSMmatrix = zeros(numTriang, numNodes, 'int32');
for indS = 1 : numNodes
    for indM = 1 : numTriang
        [~, indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
    end
end
indSMmatrixGPU = gpuArray(reshape(indSMmatrix, [numNodes*numTriang, 1]));

matCoeff = zeros(9*domainMesh.numberTriangles, 1, 'gpuArray');
vetCoeff = zeros(3*domainMesh.numberTriangles, 1, 'gpuArray');
for indTemp = 1 : numTriang
    matCoeff(9*(indTemp-1) + (1:9), 1) = reshape(constValues{indTemp}.matCoeff, [9 1]);
    vetCoeff(3*(indTemp-1) + (1:3), 1) = constValues{indTemp}.vetCoeff;
end

nodesMesh = gpuArray(reshape(domainMesh.coordinates', [3*numNodes 1]));
maxLen = max(domainMesh.maxL);

%% COMPUTAZIONE BLOCCHI MATRICIALI V

%Allocazione array contente i blocchi matriciali
matrixSavedV = cell(numBlocksPerIterV, numIterV);

% Setup kernel
kernelV = parallel.gpu.CUDAKernel("kernelVlinslac.ptx", "kernelVlinslac.cu");
kernelV.GridSize = [numTriang numTriang numBlocksPerIterV];
blockXint = methodInfo.numSubRegion;
blockYint = methodInfo.numNodiSing;
kernelV.ThreadBlockSize = [blockXint blockYint 1];
kernelV.SharedMemorySize = blockXint * blockYint * 27 * 8;

%Ciclo sulle iterazioni necessarie (limite memoria GPU)
for indIter = 1 : numIterV
    %Avvio computazione GPU
    matrixOutInt = BEMenerg_dir_launchKernelV(gpuID, kernelV, deltaT, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho, ...
                                                stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, ...
                                                vertsT, areeT, offSetsV(indIter), numBlocksV, ...
                                                numTriang, numBlocksPerIterV, maxLen, matCoeff, vetCoeff);

    %Allocazione array contente i blocchi diagonali di questa iterazione
    matrixSubBlocksDIAG = zeros(9, 3, numTriang, numBlocksPerIterV);

    %Avvio computazione CPU
    for indTemp = 1 : numBlocksPerIterV
        %Check non sforamento numero di blocchi necessario
        if offSetsV(indIter) + indTemp > numBlocksV
            continue
        end

        %Set istante temporale
        istTemp = offSetsV(indIter) + indTemp - 1;

        %Ciclo sull'indice dei sottoblocchi
        parfor indBlock = 1 : numTriang
            %Calcolo blocco mediante procedura semianalitica
            matrixSubBlocksDIAG(:, :, indBlock, indTemp) = BEMenerg_dir_calcSubBlockDiagV(pbParam, methodInfo, istTemp, deltaT, constValues{indBlock}, DIAGn, DIAGw);
        end
    end
    %Attesa completamento operazioni GPU
    matrixOutInt = gather(matrixOutInt);
    wait(gpuID);

    %Reshape matrice output GPU
    matrixOutInt = reshape(matrixOutInt, [9*numTriang 3*numTriang numBlocksPerIterV]);

    %Inserimento blocchi diagonali calcolati su CPU e salvataggio in
    %memoria
    parfor indTemp = 1 : numBlocksPerIterV
        %Check non sforamento numero di blocchi necessario
        if offSetsV(indIter) + indTemp > numBlocksV
            continue
        end

        %Selezione singolo blocco matriciale
        matrixSavedV{indTemp, indIter} = squeeze(matrixOutInt(:, :, indTemp));

        %Inserimento sottoblocchi diagonali
        for indBlock = 1 : numTriang
            indR = 9 * (indBlock - 1);
            indC = 3 * (indBlock - 1);
            matrixSavedV{indTemp, indIter}(indR + (1:9), indC + (1:3)) = matrixSubBlocksDIAG(:, :, indBlock, indTemp);
        end

        %Trasformazione in sparse matrix
        matrixSavedV{indTemp, indIter} = sparse(matrixSavedV{indTemp, indIter});
    end
    %Reset array GPU
    wait(gpuID);
end

%Reshape matrici salvate
matrixSavedV = reshape(matrixSavedV, [numIterV*numBlocksPerIterV, 1]);
matrixSavedV = matrixSavedV(1 : numBlocksV);

% AGGIUNTA BLOCCHI NULLI
for ind = (numBlocksV + 1) : pbParam.nT
    matrixSavedV{ind} = sparse(zeros(9*numTriang, 3*numTriang));
end

%Calcolo betaV
betaV = BEMenerg_dir_calcBetaVgpu(deltaT, pbParam, domainMesh, matrixSavedV);

%% Calcolo blocchi matrice K

%Setup valori
numBlocksK = min(numBlocksV + 1, pbParam.nT);
blockSizeK = 81 * numTriang * numTriang * 8;   %81 Mx^2 double (8 byte l'uno)

%Calcolo massimo numero di iterazioni parallelizzabili sulla GPU
fattoreCaricoMemoriaGPU = 2.5;                    %Considero un fattore 1.5 per tenere conto dei dati di input        
maxBlockInMemoryK = floor(memAv / (fattoreCaricoMemoriaGPU*blockSizeK));        
if maxBlockInMemoryK == 0
    error("Memoria GPU insufficiente")
end

%Calcolo numero blocchi per singola iterazione
numBlocksPerIterK = min(maxBlockInMemoryK, numBlocksK);
%Calcolo numero iterazioni
numIterK = ceil(numBlocksK / numBlocksPerIterK);
%Calcolo offsets
offSetsK = numBlocksPerIterK * (0 : (numIterK-1));
%Allocazione array contente i blocchi matriciali
matrixSavedK = cell(numBlocksPerIterK, numIterK);

% Setup kernel interno
kernelKint = parallel.gpu.CUDAKernel("kernelKlinslacInternal.ptx", "kernelKlinslacInternal.cu");
kernelKint.GridSize = [numTriang numTriang numBlocksPerIterK];
blockXint = methodInfo.numSubRegion;
blockYint = 1; % methodInfo.numNodiSing;
kernelKint.ThreadBlockSize = [blockXint blockYint 1];
kernelKint.SharedMemorySize = blockXint * blockYint * 81 * 8;

% Setup kernel bordo
kernelKbound = parallel.gpu.CUDAKernel("kernelKlinslacBoundary.ptx", "kernelKlinslacBoundary.cu");
kernelKbound.GridSize = [numTriang numTriang numBlocksPerIterK];
methodInfo.blockXbound = 64; %mettere come parametro in methodInfo fuori;
kernelKbound.ThreadBlockSize = [methodInfo.blockXbound 1 1];
kernelKbound.SharedMemorySize = methodInfo.blockXbound * 81 * 8;

%Ciclo sulle iterazioni necessarie (limite memoria GPU)
for indIter = 1 : numIterK
    %Calcolo contributi di bordo
    matrixOutBound = BEMenerg_dir_launchKernelK_uConstSpace(gpuID, kernelKbound, deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*deltaT, ...
                                                    stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                    stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, ...
                                                    vertsT, areeT, normT, offSetsK(indIter), numBlocksK, numTriang, numBlocksPerIterK, maxLen, matCoeff, vetCoeff);

    %Attesa completamento operazioni GPU
    matrixOutBound = gather(matrixOutBound);
    wait(gpuID);
    
    %Avvio computazione GPU
    matrixOutInt = BEMenerg_dir_launchKernelKinternal_uConstSpace(gpuID, kernelKint, deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*deltaT, ...
                                                    stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                    stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, methodInfo.numNodiSing, ...
                                                    vertsT, areeT, normT, offSetsK(indIter), numBlocksK, numTriang, numBlocksPerIterK, maxLen, matCoeff, vetCoeff);

    %Allocazione array contente i componenti singolari dei blocchi di questa iterazione
    matrixSubBlocksSING = zeros(9, 9, numTriang, numBlocksPerIterK);

    %Avvio computazione CPU
    for indTemp = 1 : numBlocksPerIterK
        %Check non sforamento numero di blocchi necessario
        if offSetsK(indIter) + indTemp > numBlocksK
            continue
        end 

        %Set istante temporale
        indTempCurr = offSetsK(indIter) + indTemp - 1;

        %Ciclo sull'indice di riga
        parfor indRow = 1 : numTriang
            matrixSubBlocksSING(:, :, indRow, indTemp) = BEMenerg_dir_calcSubBlockSingK_uConstSpace(methodInfo, deltaT, pbParam, domainMesh, DIAGn, DIAGw, indTempCurr, indRow, constValues{indRow}, 0);
        end
    end

    %Attesa completamento operazioni GPU
    matrixOutInt = gather(matrixOutInt);
    wait(gpuID);

    %Somma contributi interni e di bordo
    matrixOut = matrixOutInt + matrixOutBound;

    %Reshape matrice output GPU
    matrixOut = reshape(matrixOut, [9*numTriang 9*numTriang numBlocksPerIterK]);

    %Aggiunta componenti blocchi singolari calcolate su CPU e salvataggio in memoria
    parfor indTemp = 1 : numBlocksPerIterK
        %Check non sforamento numero di blocchi necessario
        if offSetsK(indIter) + indTemp > numBlocksK
            continue
        end

        %Selezione singolo blocco matriciale
        matrixSavedK{indTemp, indIter} = squeeze(matrixOut(:, :, indTemp));

        %Aggiunta sottoblocchi singolari
        for indRow = 1 : numTriang
            matrixSavedK{indTemp, indIter}(9*(indRow-1) + (1:9), 9*(indRow-1) + (1:9)) = matrixSavedK{indTemp, indIter}(9*(indRow-1) + (1:9), 9*(indRow-1) + (1:9)) + matrixSubBlocksSING(:, :, indRow, indTemp);
        end

        %Trasformazione in sparse matrix
        matrixSavedK{indTemp, indIter} = sparse(matrixSavedK{indTemp, indIter});
    end
    %Reset array GPU
    wait(gpuID);
end

%Reshape matrici salvate
matrixSavedK = reshape(matrixSavedK, [numIterK*numBlocksPerIterK, 1]);
matrixSavedK = matrixSavedK(1 : numBlocksK);

% AGGIUNTA BLOCCHI NULLI
for ind = (numBlocksK + 1) : pbParam.nT
    matrixSavedK{ind} = sparse(zeros(9*numTriang, 9*numTriang));
end

%% Calcolo matrice I
baseMat = (ones(3) + diag(diag(ones(3)))) ./ 24;
matrixIGamma = kron(kron(eye(numTriang) .* domainMesh.area, baseMat) , eye(3));

%% RISOLUZIONE SISTEMA LINEARE
%Fattorizzazione matrice E0
matrixSist = matrixSavedK{1} + matrixIGamma;
%[L, U, P] = lu(matrixSist);

%Allocazione array contente la densità incognita
density = zeros(9*numTriang, pbParam.nT);

%Procedimento di time-marching
for currInd = 1 : pbParam.nT
    tnf = betaV{currInd};%; 

    if(currInd > 1)
        tnf = tnf + matrixIGamma * density(:, currInd - 1);
    end
    endInd = min(currInd, numBlocksK);

    for indMat = 2 : endInd
        tnf = tnf - matrixSavedK{indMat} * density(:, currInd - indMat + 1);  
    end
  
    density(:, currInd) = matrixSist \ tnf;
end
reset(gpuID);

%Salvataggio temporaneo della variabile soluzione
filePath = strcat(outPath, '_density');
save(filePath, 'density');


return
