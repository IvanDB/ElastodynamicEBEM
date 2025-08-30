function betaK = BEMenerg_dir_calcBetaKcpu_uConstSpace(methodInfo, deltaT, pbParam, domainMesh, constValues, DIAGn, DIAGw)
    
%Calcolo dato al bordo (Dirichlet)
gK = BEMenerg_dir_calcBoundDataDirichlet_uConstSpace(deltaT, pbParam, domainMesh);

    
for indTemp = 1 : pbParam.nT
    K{indTemp} = BEMenerg_dir_calcKBlock_uConstSpace(methodInfo, indTemp - 1, deltaT, pbParam, domainMesh, constValues, DIAGn, DIAGw);
end

betaK = cell(1, pbParam.nT);
for indTemp = 1 : pbParam.nT
    betaK{indTemp} = zeros(3*domainMesh.numberTriangles, 1);
    for j = 1 : indTemp
        betaK{indTemp} = betaK{indTemp} + K{indTemp - j + 1} * gK{j};
    end
end
return