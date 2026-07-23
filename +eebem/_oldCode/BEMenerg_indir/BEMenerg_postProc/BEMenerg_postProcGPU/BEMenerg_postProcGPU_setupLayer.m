function campoVett = BEMenerg_postProcGPU_setupLayer(pbParam, domainMesh, density, xVal, tVal, methodInfo, PPn, PPw)

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);

%Calcolo dimensioni ed inizializzazione matrice soluzione
campoVett = cell(xSize, tSize);

for ind = 1 : tSize
    campoVett(:, ind) = BEMenerg_postProcGPU_calcLayer(pbParam, domainMesh, density, methodInfo, xVal, tVal(ind), PPn, PPw);
end

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_GPU_campoVett");
save(filePath, 'campoVett');  

return