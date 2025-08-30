function betaV = BEMenerg_dirN_calcBetaVgpu(deltaT, pbParam, domainMesh, matrixSavedV)

%Calcolo dato al bordo (Neumann)
gV = BEMenerg_dirN_calcBoundDataNeumann(deltaT, pbParam, domainMesh);

betaV = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaV{indTemp} = zeros(3*domainMesh.numberTriangles, 1);
    for j = 1 : indTemp
        betaV{indTemp} = betaV{indTemp} + matrixSavedV{indTemp - j + 1} * gV{j};
    end
end
return