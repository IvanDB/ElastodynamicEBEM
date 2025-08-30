function betaK = BEMenerg_dir_calcBetaKgpu_uConstSpace(deltaT, pbParam, domainMesh, matrixSavedK)
    
%Calcolo dato al bordo (Dirichlet)
gK = BEMenerg_dir_calcBoundDataDirichlet_uConstSpace(deltaT, pbParam, domainMesh);

%Costruzione vettore betaK
betaK = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaK{indTemp} = zeros(3*domainMesh.numberTriangles, 1);
    for j = 1 : indTemp
        betaK{indTemp} = betaK{indTemp} + matrixSavedK{indTemp - j + 1} * gK{j};
    end
end
return