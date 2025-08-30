function [campoVett, campoVettG, campoVettP] = BEMenerg_ppD_setupLayer_TEST(pbParam, domainMesh, density, ~, xVal, tVal, methodInfo, PPn, PPw)

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);

deltaT = pbParam.Tfin / pbParam.nT;

%Calcolo dato al bordo (Dirichlet)
gK = BEMenerg_dir_calcBoundDataDirichlet(deltaT, pbParam, domainMesh);

%Calcolo dimensioni ed inizializzazione matrice soluzione
campoVett = cell(xSize, tSize);
campoVettG = cell(xSize, tSize);
campoVettP = cell(xSize, tSize);

for indT = 1 : tSize
    [campoVett(:, indT), campoVettG(:, indT), campoVettP(:, indT)]  = BEMenerg_ppD_calcLayer_TEST(pbParam, domainMesh, density, gK, methodInfo, xVal, tVal(indT), PPn, PPw);
end

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_NEW_campoVett");
save(filePath, 'campoVett');  
%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_NEW_campoVettG");
save(filePath, 'campoVettG');  
%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_NEW_campoVettP");
save(filePath, 'campoVettP');  
return