function density = BEMenerg_coreFN_timeMarchingTEST(pbParam, domainMesh, methodInfo, GHn, GHw, GHCn, GHCw, DIAGn, DIAGw)
%

%% SETUP INIZIALE
%Inizializzazione GPU
gpuID = gpuDevice;
reset(gpuID);
memAv = gpuID.AvailableMemory;

%Dimensioni problema
numBlocksV = BEMenerg_coreFN_calcNumBlocks(pbParam, domainMesh);
numTriang = domainMesh.numberTriangles;
numNodes = domainMesh.number_nodes;

blockSizeV = 9 * (numTriang .^ 2) * 8;   %9Nx^2 double (8 byte l'uno)

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
constValues = BEMenerg_coreFN_calcCostantData(methodInfo, domainMesh, GHn, GHw, GHCn, GHCw, DIAGn, DIAGw);

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
kernelV = parallel.gpu.CUDAKernel("kernelGHC.ptx", "kernelGHC.cu");
kernelV.GridSize = [numTriang numTriang numBlocksPerIterV];
blockX = methodInfo.numSubRegion;
blockY = methodInfo.numNodiSing;
kernelV.ThreadBlockSize = [blockX blockY 1];
kernelV.SharedMemorySize = blockX * blockY * 9 * 8;

% %Ciclo sulle iterazioni necessarie (limite memoria GPU)
for indIter = 1 : numIterV
    %Avvio computazione GPU
    matrixOut = BEMenerg_coreFN_launchKernel(gpuID, kernelV, deltaT, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho, ...
                                                stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, ...
                                                vertsT, areeT, offSetsV(indIter), numBlocksV, ...
                                                numTriang, numBlocksPerIterV);

    %Allocazione array contente i blocchi diagonali di questa iterazione
    matrixSubBlocksDIAG = zeros(3, 3, numTriang, numBlocksPerIterV);

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
            %Estrazione dati 
            node = domainMesh.triangles(indBlock, 1:3);
            verts = domainMesh.coordinates(node, :);
            vn = domainMesh.normal(indBlock, :);

            %Calcolo blocco mediante procedura semianalitica
            matrixSubBlocksDIAG(:, :, indBlock, indTemp) = BEMenerg_coreFN_calcSubBlockDiag(pbParam, methodInfo, istTemp, deltaT, constValues{indBlock}, DIAGn, DIAGw);
        end
    end
    %Attesa completamento operazioni GPU
    matrixOut = gather(matrixOut);
    wait(gpuID);

    %Reshape matrice output GPU
    matrixOut = reshape(matrixOut, [3*numTriang 3*numTriang numBlocksPerIterV]);

    %Inserimento blocchi diagonali calcolati su CPU e salvataggio in
    %memoria
    parfor indTemp = 1 : numBlocksPerIterV
        %Check non sforamento numero di blocchi necessario
        if offSetsV(indIter) + indTemp > numBlocksV
            continue
        end

        %Selezione singolo blocco matriciale
        matrixSavedV{indTemp, indIter} = squeeze(matrixOut(:, :, indTemp));

        %Inserimento sottoblocchi diagonali
        for indBlock = 1 : numTriang
            indRC = 3 * (indBlock - 1);
            matrixSavedV{indTemp, indIter}(indRC + (1:3), indRC + (1:3)) = matrixSubBlocksDIAG(:, :, indBlock, indTemp);
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
    matrixSavedV{ind} = sparse(zeros(3*numTriang, 3*numTriang));
end

%Calcolo betaV
%betaV = BEMenerg_core_calcBetaVgpu(deltaT, pbParam, domainMesh, matrixSavedV);


%% Calcolo blocchi matrice K

%Setup valori
numBlocksK = min(numBlocksV + 1, pbParam.nT);
blockSizeK = 9 * numTriang * 3 * numTriang * 8;   %9NxNn double (8 byte l'uno)

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

% Setup kernel
kernelK = parallel.gpu.CUDAKernel("kernelKlinFree.ptx", "kernelKlinFree.cu");
kernelK.GridSize = [numTriang numTriang numBlocksPerIterK];
blockX = methodInfo.numSubRegion;
blockY = methodInfo.numNodiSing;
kernelK.ThreadBlockSize = [blockX blockY 1];
kernelK.SharedMemorySize = blockX * blockY * 27 * 8;

%Ciclo sulle iterazioni necessarie (limite memoria GPU)
for indIter = 1 : numIterK
    %Avvio computazione GPU
    matrixOut = BEMenerg_coreFN_launchKernelK(gpuID, kernelK, deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*deltaT, ...
                                                    stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                    stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, ...
                                                    vertsT, areeT, normT, indSMmatrixGPU, matCoeff, vetCoeff, ...
                                                    offSetsK(indIter), numBlocksK, numTriang, numNodes, numBlocksPerIterK, nodesMesh, maxLen);
    
    %Allocazione array contente i componenti singolari dei blocchi di questa iterazione
    matrixSubBlocksSING = zeros(3, 3, 3, numTriang, numBlocksPerIterK);
    
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
            for indV = 1 : 3
                matrixSubBlocksSING(:, :, indV, indRow, indTemp) = BEMenerg_core_calcKSubBlockSing(methodInfo, deltaT, pbParam, domainMesh, DIAGn, DIAGw, indTempCurr, indRow, constValues{indRow}, indV);
            end
        end
    end

    %Attesa completamento operazioni GPU
    matrixOut = gather(matrixOut);
    wait(gpuID);

    %Reshape matrice output GPU
    matrixOut = reshape(matrixOut, [3*numTriang 9*numTriang numBlocksPerIterK]);
    
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
            for indVertCurr = 1 : 3
                indCol = 3*(indRow-1) + indVertCurr;
                matrixSavedK{indTemp, indIter}(3*(indRow-1) + (1:3), 3*(indCol-1) + (1:3)) = matrixSavedK{indTemp, indIter}(3*(indRow-1) + (1:3), 3*(indCol-1) + (1:3)) + matrixSubBlocksSING(:, :, indVertCurr, indRow, indTemp);
            end
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
    matrixSavedK{ind} = sparse(zeros(3*numTriang, 9*numTriang));
end

%% Calcolo matrice I
matrixIGamma = kron(eye(numTriang), [1, 1, 1]);
matrixIGamma = kron(matrixIGamma .* domainMesh.area ./ 6, eye(3));

%% Calcolo vettori betaI*

h = 1 * strcmp(pbParam.domainType, "barH1") + 3*strcmp(pbParam.domainType, "barH3");
velP = pbParam.velP;
hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
k_val = (0 : hat_k)';
tildeP = @(x, t) sum((-1).^k_val .* (((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                + ((velP.*t - 2.*h.*(k_val+1) + (h - x(3))) > 0)));
gV = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gVcurr = cell(domainMesh.numberTriangles, 1);
    for indT = 1 : domainMesh.numberTriangles
        gVcurr{indT} = [0; 0; 2*pbParam.mu * tildeP(domainMesh.center(indT, :), indTemp*deltaT) * domainMesh.normal(indT, 3)];
    end
    gV{indTemp} = cell2mat(gVcurr);
end

matrixIGammaStar = matrixIGamma';
for indTemp = 1 : pbParam.nT
    betaIstar{indTemp} = matrixIGammaStar * gV{indTemp};
end

tildeU = @(x, t) sum((-1).^k_val .* ( (velP.*t - 2.*h.*k_val - (h-x(3))) .* ((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                        - (velP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((velP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
gK = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gKcurr = cell(3*domainMesh.numberTriangles, 1);
    for indT = 1 : domainMesh.numberTriangles
        gKcurr{3*(indT-1) + 1} = [0; 0; tildeU(domainMesh.coordinates(domainMesh.triangles(indT, 1), :), indTemp*deltaT)];
        gKcurr{3*(indT-1) + 2} = [0; 0; tildeU(domainMesh.coordinates(domainMesh.triangles(indT, 2), :), indTemp*deltaT)];
        gKcurr{3*(indT-1) + 3} = [0; 0; tildeU(domainMesh.coordinates(domainMesh.triangles(indT, 3), :), indTemp*deltaT)];
    end
    gK{indTemp} = cell2mat(gKcurr);
end

%% CHECK matrix K
rhs = zeros(3*numTriang, pbParam.nT);
%Calcolo elementi del RHS
parfor indTemp = 1 : pbParam.nT
    istTemp = indTemp - 1 ;
    rhs(:, indTemp) = BEMenerg_core_calcTnBlock(rhs(:, indTemp), methodInfo, istTemp, deltaT, pbParam, domainMesh, constValues);
end

betaK = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaK{indTemp} = zeros(3*domainMesh.numberTriangles, 1);
    for j = 1 : indTemp
        betaK{indTemp} = betaK{indTemp} + matrixSavedK{indTemp - j + 1} * gK{j};
    end
end

%Fattorizzazione matrice E0
matrixSist = matrixSavedV{1};
[L, U, P] = lu(matrixSist);

%Allocazione array contente la densità incognita
density = zeros(3*numTriang, pbParam.nT);

%Procedimento di time-marching
for currInd = 1 : pbParam.nT
    tnf = rhs(:, currInd) ./ 2 + betaK{currInd};%; 

    endInd = min(currInd, numBlocksV);

    for indMat = 2 : endInd
        tnf = tnf - matrixSavedV{indMat} * density(:, currInd - indMat + 1);  
    end
  
    density(:, currInd) = U\(L\(P*tnf));
end

%plotDensity(pbParam, domainMesh, density, 100);


%% RISOLUZIONE SISTEMA LINEARE
%Fattorizzazione matrice E0
matrixSist = [matrixSavedK{1} + matrixIGamma, -matrixSavedV{1}; ...
                zeros(9*numTriang, 9*numTriang), matrixIGammaStar];
%[L, U, P] = lu(matrixSist);

%Allocazione array contente la densità incognita
density = zeros(3*(numNodes + numTriang), pbParam.nT);

matDiagI = [matrixIGamma, zeros(3.*numTriang, 3.*numTriang);
            zeros(3.*numTriang, 3.*(numTriang + numNodes))];

%Procedimento di time-marching
for currInd = 1 : pbParam.nT
    tnf = [zeros(3.*numTriang, 1); betaIstar{currInd}]; 

    if(currInd > 1)
        tnf = tnf + matDiagI * density(:, currInd - 1);
    end

    endInd = min(currInd, numBlocksK);

    for indMat = 2 : endInd
        matrixPasso = [matrixSavedK{indMat}, -matrixSavedV{indMat};
                        zeros(3.*numTriang, 3.*(numTriang + numNodes))];
        tnf = tnf - matrixPasso * density(:, currInd - indMat + 1);  
    end
  
    density(:, currInd) = matrixSist \ tnf;
    disp(max(abs(tnf - matrixSist*density(:, currInd))));
    %density(:, currInd) = U\(L\(P*tnf));
end
reset(gpuID);

%Salvataggio temporaneo della variabile soluzione
filePath = strcat(outPath, '_density');
save(filePath, 'density');

return