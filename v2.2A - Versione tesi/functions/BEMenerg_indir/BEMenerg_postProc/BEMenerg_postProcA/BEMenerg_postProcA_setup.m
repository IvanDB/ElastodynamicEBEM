function campoVett = BEMenerg_postProcA_setup(pbParam, domainMesh, density, numPoints, xVal, tVal)

%Calcolo valori costanti
constValues = BEMenerg_postProcA_calcCostantData(domainMesh);

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);

%Calcolo dimensioni ed inizializzazione matrice soluzione
solRaw = cell(xSize * tSize, 1);

parfor ind = 1 : numPoints
    [indX, indT] = ind2sub([xSize, tSize], ind);
    solRaw{ind} = BEMenerg_postProcA_singlePoint(pbParam, domainMesh, density, constValues, xVal(indX, :), tVal(indT));
end

campoVett = reshape(solRaw, [xSize, tSize]);

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_A_campoVett");
save(filePath, 'campoVett'); 
return