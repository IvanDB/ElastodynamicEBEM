function campoVett = BEMenerg_ppN_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw)

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);

deltaT = pbParam.Tfin / pbParam.nT;

%Calcolo dato al bordo (Dirichlet)
gV = BEMenerg_dir_calcBoundDataNeumann(deltaT, pbParam, domainMesh);

%Calcolo dimensioni ed inizializzazione matrice soluzione
solRaw = cell(xSize * tSize, 1);

parfor ind = 1 : numPoints
    [indX, indT] = ind2sub([xSize, tSize], ind);
    solRaw{ind} = BEMenerg_ppN_calc(pbParam, domainMesh, density, gV, methodInfo, xVal(indX, :), tVal(indT), PPn, PPw);
end

campoVett = reshape(solRaw, [xSize, tSize]);

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_NEW_campoVett");
save(filePath, 'campoVett');  

return