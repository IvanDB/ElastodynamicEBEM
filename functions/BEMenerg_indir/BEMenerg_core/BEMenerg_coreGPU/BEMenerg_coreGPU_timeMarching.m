function density = BEMenerg_coreGPU_timeMarching(pbParam, domainMesh, methodInfo, GHn, GHw, GHCn, GHCw)
%

%% SETUP INIZIALE
%Inizializzazione GPU
gpuID = gpuDevice;
reset(gpuID);
memAv = gpuID.AvailableMemory;

%Dimensioni problema
numBlocks = BEMenerg_coreGPU_calcNumBlocks(pbParam, domainMesh);
numTriang = domainMesh.numberTriangles;

blockSize = 9 * (numTriang .^ 2) * 8;   %9Nx^2 double (8 byte l'uno)

%Calcolo massimo numero di iterazioni parallelizzabili sulla GPU
fattoreCaricoMemoriaGPU = 1.5;                    %Considero un fattore 1.5 per tenere conto dei dati di input        
maxBlockInMemory = floor(memAv / (fattoreCaricoMemoriaGPU*blockSize));        
if maxBlockInMemory == 0
    error("Memoria GPU insufficiente")
end

%Calcolo numero blocchi per singola iterazione
numBlocksPerIter = min(maxBlockInMemory, numBlocks);
%Calcolo numero iterazioni
numIter = ceil(numBlocks / numBlocksPerIter);
%Calcolo offsets
offSets = numBlocksPerIter * (0 : (numIter-1));

%Path per il salvataggio delle matrici e della soluzione
% tmpPath = strcat("./tempData/", pbParam.domainType, num2str(pbParam.lev), "_GPU_", ...
%    num2str(methodInfo.numNodiExt), "_", num2str(methodInfo.numSubRegion), "_", num2str(methodInfo.numNodiSing));
outPath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_GPU_", ...
    num2str(methodInfo.numNodiExt), "_", num2str(methodInfo.numSubRegion), "_", num2str(methodInfo.numNodiSing));


%% SETUP VARIABILI LOCALI PER COMPUTAZIONE GPU

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

%% COMPUTAZIONE BLOCCHI MATRICIALI
%Calcolo valori costanti per computazione CPU
constValues = BEMenerg_coreGPU_calcCostantData(methodInfo, domainMesh, GHn, GHw);

%Allocazione array contente il termine noto
rhs = zeros(3*numTriang, pbParam.nT);
%Allocazione array contente i blocchi matriciali
matrixSaved = cell(numBlocksPerIter, numIter);

% Setup kernel
kernel = parallel.gpu.CUDAKernel("kernelGHC.ptx", "kernelGHC.cu");
kernel.GridSize = [numTriang numTriang numBlocksPerIter];
blockX = methodInfo.numSubRegion;
blockY = methodInfo.numNodiSing;
kernel.ThreadBlockSize = [blockX blockY 1];
kernel.SharedMemorySize = blockX * blockY * 9 * 8;

%Ciclo sulle iterazioni necessarie (limite memoria GPU)
for indIter = 1 : numIter
    %Avvio computazione GPU
    matrixOut = BEMenerg_coreGPU_launchKernel(gpuID, kernel, deltaT, pbParam.velP, pbParam.velS, 4*pi*pbParam.rho, ...
                                                stdGHw, stdGHnx, stdGHny, stdGHnz, methodInfo.numNodiExt, ...
                                                stdGHCw, stdGHCnx, stdGHCny, stdGHCnz, ...
                                                vertsT, areeT, offSets(indIter), numBlocks, ...
                                                numTriang, numBlocksPerIter);

    %Allocazione array contente i blocchi diagonali di questa iterazione
    matrixSubBlocksSA = zeros(3, 3, numTriang, numBlocksPerIter);

    %Avvio computazione CPU
    for indTemp = 1 : numBlocksPerIter
        %Check non sforamento numero di blocchi necessario
        if offSets(indIter) + indTemp > numBlocks
            continue
        end

        %Set istante temporale
        istTemp = offSets(indIter) + indTemp - 1;

        %Ciclo sull'indice dei sottoblocchi
        parfor indBlock = 1 : numTriang
            %Estrazione dati 
            node = domainMesh.triangles(indBlock, 1:3);
            verts = domainMesh.coordinates(node, :);
            vn = domainMesh.normal(indBlock, :);
            
            %Calcolo blocco mediante procedura semianalitica
            matrixSubBlocksSA(:, :, indBlock, indTemp) = BEMenerg_coreGPUdiag_calcSubBlock(pbParam, methodInfo, verts, vn, istTemp, deltaT, constValues{indBlock});
        end
    end
    %Attesa completamento operazioni GPU
    matrixOut = gather(matrixOut);
    wait(gpuID);

    %Reshape matrice output GPU
    matrixOut = reshape(matrixOut, [3*numTriang 3*numTriang numBlocksPerIter]);

    %Inserimento blocchi diagonali calcolati su CPU e salvataggio in
    %memoria
    parfor indTemp = 1 : numBlocksPerIter
        %Check non sforamento numero di blocchi necessario
        if offSets(indIter) + indTemp > numBlocks
            continue
        end

        %Selezione singolo blocco matriciale
        matrixSaved{indTemp, indIter} = squeeze(matrixOut(:, :, indTemp));
        %Inserimento sottoblocchi diagonali
        for indBlock = 1 : numTriang
            indRC = 3 * (indBlock - 1);
            matrixSaved{indTemp, indIter}(indRC + (1:3), indRC + (1:3)) = matrixSubBlocksSA(:, :, indBlock, indTemp);
        end
        matrixSaved{indTemp, indIter} = sparse(matrixSaved{indTemp, indIter});
    end

    %Reset array GPU
    wait(gpuID);
end

%Calcolo elementi del RHS
parfor indTemp = 1 : pbParam.nT
    istTemp = indTemp - 1 ;
    rhs(:, indTemp) = BEMenerg_core_calcTnBlock(rhs(:, indTemp), methodInfo, istTemp, deltaT, pbParam, domainMesh, constValues);
end

%Reshape matrici salvate
matrixSaved = reshape(matrixSaved, [numIter*numBlocksPerIter, 1]);
matrixSaved = matrixSaved(1 : numBlocks);

%% RISOLUZIONE SISTEMA LINEARE
%Fattorizzazione matrice E0
matrixSist = matrixSaved{1};
[L, U, P] = lu(matrixSist);

%Allocazione array contente la densità incognita
density = zeros(3*numTriang, pbParam.nT);

%Procedimento di time-marching
for currInd = 1 : pbParam.nT
    tnf = rhs(:, currInd);

    endInd = min(currInd, numBlocks);

    for indMat = 2 : endInd
        tnf = tnf - matrixSaved{indMat} * density(:, currInd - indMat + 1);  
    end
  
    density(:, currInd) = U\(L\(P*tnf));
end
reset(gpuID);

%Salvataggio temporaneo della variabile soluzione
filePath = strcat(outPath, '_density');
save(filePath, 'density');
return