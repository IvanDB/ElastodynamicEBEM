function density = BEMenerg_coreSA_timeMarching(pbParam, domainMesh, methodInfo, GHn, GHw)
% INPUT 
%   - pbParam: struct contenente i parametri del problema
%   - domainMesh: struct contente le informazioni sulla mesh
%   - methodInfo: struct contente le informazioni sull metodo di quadratura
%   - GHn: matrice contenente i nodi di quadratura di Gauss-Hammer
%   - GHw: matrice contenente i pesi di quadratura di Gauss-Hammer
%
% OUTPUT:
%   - density: matrice contenente la densità incognita per ogni tempo

%% INIZIALIZZAZIONE COMPONENTI del METODO

%Calcolo del PASSO di DISCRETIZZAZIONE TEMPORALE (deltaT = T_fin / nT)
deltaT = pbParam.Tfin / pbParam.nT;

%Allocazione array SOLUZIONE
density = zeros(3*domainMesh.numberTriangles, pbParam.nT);

%Allocazione array BLOCCO MATRICIALE
matrix = zeros(3*domainMesh.numberTriangles, 3*domainMesh.numberTriangles);

%Allocazione array TERMINE NOTO
rhs = zeros(3*domainMesh.numberTriangles, pbParam.nT);

%Inizializzazione a zero parametro che indica l'indice del SOTTOINTERVALLO 
% TEMPORALE a partire dal quale i successivi BLOCCHI MATRICIALI sono NULLI.     
maxstep = -1;

%Inizializzazione path directory temp e output
tmpPath = strcat("./tempData/", pbParam.domainType, num2str(pbParam.lev), "_SA_", num2str(methodInfo.numNodiExt));
outPath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_SA_", num2str(methodInfo.numNodiExt));

% CALCOLO PARAMETRI COSTANTI
constValues = BEMenerg_coreSA_calcCostantData(methodInfo, domainMesh, GHn, GHw);

%% CALCOLO BLOCCO MATRICIALE E0 e TERMINE NOTO
%Calcolo blocco matriciale E0
[matrix, ~] = BEMenerg_coreSA_calcMatrixBlock(matrix, methodInfo, 0, deltaT, pbParam, domainMesh, constValues);

%Scrittura su file del blocco E0
filePath = strcat(tmpPath, "_matrixE", num2str(0));
save(filePath, 'matrix');

%Fattorizzazione PLU del blocco matriciale E0
[L, U, P] = lu(sparse(matrix));
    
%Calcolo blocchi termine noto
parfor indPos = 1 : pbParam.nT
    indTemp = indPos - 1;
    rhs(:, indPos) = BEMenerg_core_calcTnBlock(rhs(:, indPos), methodInfo, indTemp, deltaT, pbParam, domainMesh, constValues);
end

%% RISOLUZIONE PRIMO ISTANTE 
%Estrazione blocco b0 termine noto
rhsCurr = rhs(:, 1);

%Risoluzione primo sistema
density(:, 1) = U\(L\(P*rhsCurr));

%Inizializzazione file temporaneo per salvataggio variabile soluzione
filePath = strcat(tmpPath, "_soluzione");
save(filePath, 'density');

%% RISOLUZIONE ISTANTI SUCCESSIVI
%Inizializzane indice finale sommatoria per rhs
indFinSom = 0;

%Ciclo (non-parallelizzabile) sugli instanti temporali
for indStep = 2 : pbParam.nT
    %Aggiornamento indice temporale indTemp (= indStep - 1)
    indTemp = indStep - 1;

    %Estrazione blocco termine noto passo corrente
    rhsCurr = rhs(:, indStep);

    %Aggiornamento rhs passo corrente
    for indInner = 1 : indFinSom
        filePath = strcat(tmpPath, "_matrixE", num2str(indInner));
        matrix = sparse(load(filePath).matrix);
        rhsCurr = rhsCurr - matrix*density(:, indStep - indInner);
    end
    
    if maxstep < 0
        %Reset matrice
        matrix = zeros(3*domainMesh.numberTriangles, 3*domainMesh.numberTriangles);

        %Calcolo blocco matriciale passo corrente
        [matrix, isnotzero] = BEMenerg_coreSA_calcMatrixBlock(matrix, methodInfo, indTemp, deltaT, pbParam, domainMesh, constValues);
        
        %Controllo flag isnotzero
        if isnotzero
            %Aggiornamento finale rhs
            rhsCurr = rhsCurr - matrix*density(:, 1);

            %Salvataggio su file matrice passo corrente
            filePath = strcat(tmpPath, "_matrixE", num2str(indTemp));
            save(filePath, 'matrix')
            
            %Aggiornamento indice limite blocchi non nulli calcolati
            indFinSom = indFinSom + 1;
        else
            maxstep = indStep;
        end
    end

    %Risoluzione sistema lineare
    density(:, indStep) = U\(L\(P*rhsCurr));

    %Aggiornamento salvataggio temporaneo della variabile densità
    filePath = strcat(tmpPath, "_soluzione");
    save(filePath, 'density');  
end

%Salvataggio finale densità calcolata
filePath = strcat(outPath, "_density");
save(filePath, 'density'); 
return