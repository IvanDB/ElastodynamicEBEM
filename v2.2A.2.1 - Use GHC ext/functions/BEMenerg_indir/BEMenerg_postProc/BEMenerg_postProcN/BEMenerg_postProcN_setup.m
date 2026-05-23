function campoVett = BEMenerg_postProcN_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, nStd, wStd)

%Calcolo valori costanti
constValues = BEMenerg_postProcN_calcCostantData(domainMesh, methodInfo, nStd, wStd);

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);

%Calcolo dimensioni ed inizializzazione matrice soluzione
solRaw = cell(xSize * tSize, 1);

parfor ind = 1 : numPoints
    [indX, indT] = ind2sub([xSize, tSize], ind);
    solRaw{ind} = BEMenerg_postProcN_singlePoint(pbParam, domainMesh, density, methodInfo, constValues, xVal(indX, :), tVal(indT));
end

campoVett = reshape(solRaw, [xSize, tSize]);

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_N_campoVett");
save(filePath, 'campoVett');   
return