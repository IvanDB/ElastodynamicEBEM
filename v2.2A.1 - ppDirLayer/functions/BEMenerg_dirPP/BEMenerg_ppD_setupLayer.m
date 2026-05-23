function campoVett = BEMenerg_ppD_setupLayer(pbParam, domainMesh, density, ~, xVal, tVal, methodInfo, PPn, PPw)

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);

deltaT = pbParam.Tfin / pbParam.nT;

%Calcolo dato al bordo (Dirichlet)
gK = BEMenerg_dir_calcBoundDataDirichlet(deltaT, pbParam, domainMesh);

%Calcolo dimensioni ed inizializzazione matrice soluzione
campoVett = cell(xSize, tSize);

for indT = 1 : tSize
    campoVett(:, indT) = BEMenerg_ppD_calcLayer(pbParam, domainMesh, density, gK, methodInfo, xVal, tVal(indT), PPn, PPw);
end

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_NEW_campoVett");
save(filePath, 'campoVett');  

return