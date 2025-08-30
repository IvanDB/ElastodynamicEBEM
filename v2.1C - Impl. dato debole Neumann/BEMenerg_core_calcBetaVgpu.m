function betaV = BEMenerg_core_calcBetaVgpu(deltaT, pbParam, domainMesh, matrixSavedV)
    
%%
h = 1 * strcmp(pbParam.domainType, "barH1") + 3*strcmp(pbParam.domainType, "barH3");
velP = pbParam.velP;
hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
k_val = (0 : hat_k)';
tildeP = @(x, t) sum((-1).^k_val .* (((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                + ((velP.*t - 2.*h.*(k_val+1) + (h - x(3))) > 0)));


gV = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gVcurr = cell(domainMesh.numberTriangles, 1);
    for indT = 1 : domainMesh.numberTriangles
        gVcurr{indT} = [0; 0; 2*pbParam.mu * tildeP(domainMesh.center(indT, :), indTemp*deltaT) * domainMesh.normal(indT, 3)];
        %gVcurr{indT} = [0; 0; 0];
        %if(all(domainMesh.coordinates(domainMesh.triangles(indT, 1:3), 3) > -0.0015) || all(domainMesh.coordinates(domainMesh.triangles(indT, 1:3), 3) < -0.35))
        %if(all(domainMesh.coordinates(domainMesh.triangles(indT, 1:3), 3) == 0) || all(domainMesh.coordinates(domainMesh.triangles(indT, 1:3), 3) == h))
        %    gVcurr{indT}(3) = domainMesh.normal(indT, 3);
        %end
    end
    gV{indTemp} = cell2mat(gVcurr);
end

betaV = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaV{indTemp} = zeros(3*domainMesh.numberTriangles, 1);
    for j = 1 : indTemp
        betaV{indTemp} = betaV{indTemp} + matrixSavedV{indTemp - j + 1} * gV{j};
    end
end
end