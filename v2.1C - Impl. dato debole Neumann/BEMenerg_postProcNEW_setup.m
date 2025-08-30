function campoVett = BEMenerg_postProcNEW_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw)

%Calcolo dimensioni
xSize = size(xVal, 1);
tSize = length(tVal);


%Calcolo dato al bordo (Dirichlet)
h = 1 * strcmp(pbParam.domainType, "barH1") + 3*strcmp(pbParam.domainType, "barH3");
velP = pbParam.velP;
hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
k_val = (0 : hat_k)';
tildeP = @(x, t) sum((-1).^k_val .* (((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                + ((velP.*t - 2.*h.*(k_val+1) + (h - x(3))) > 0)));

deltaT = pbParam.Tfin / pbParam.nT;
gV = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gVcurr = cell(domainMesh.numberTriangles, 1);
    for indT = 1 : domainMesh.numberTriangles
        gVcurr{indT} = [0; 0; 2*pbParam.mu * tildeP(domainMesh.center(indT, :), indTemp*deltaT) * domainMesh.normal(indT, 3)];
    end
    gV{indTemp} = cell2mat(gVcurr);
end

%Calcolo dimensioni ed inizializzazione matrice soluzione
solRaw = cell(xSize * tSize, 1);

parfor ind = 1 : numPoints
    [indX, indT] = ind2sub([xSize, tSize], ind);
    solRaw{ind} = BEMenerg_postProcNEW_calc(pbParam, domainMesh, density, gV, methodInfo, xVal(indX, :), tVal(indT), PPn, PPw);
end

campoVett = reshape(solRaw, [xSize, tSize]);

%Salvataggio campo vettoriale calcolato
filePath = strcat("./outputData/", pbParam.domainType, num2str(pbParam.lev), "_NEW_campoVett");
save(filePath, 'campoVett');  

return