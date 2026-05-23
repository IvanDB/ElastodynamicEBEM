function campoVett = BEMenerg_postProcGPU_setupSP(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw)

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);

%Calcolo dimensioni ed inizializzazione matrice soluzione
solRaw = cell(xSize * tSize, 1);

parfor ind = 1 : numPoints
    [indX, indT] = ind2sub([xSize, tSize], ind);
    solRaw{ind} = BEMenerg_postProcGPU_calcSP(pbParam, domainMesh, density, methodInfo, xVal(indX, :), tVal(indT), PPn, PPw);
end

campoVett = reshape(solRaw, [xSize, tSize]);

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_GPU_campoVett");
save(filePath, 'campoVett');  

return