function betaK = BEMenerg_dir_calcBetaKcpu(methodInfo, deltaT, pbParam, domainMesh, constValues, DIAGn, DIAGw)
    
%%
    h = 3;
    velP = pbParam.velP;
    hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
    k_val = (0 : hat_k)';
    sol_an = @(x, t) sum((-1).^k_val .* ( (velP.*t - 2.*h.*k_val - (h-x(3))) .* ((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                    - (velP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((velP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
    g = @(x, t) [0; 0; sol_an(x, t)];

    
    for indTemp = 1 : pbParam.nT
        K{indTemp} = BEMenerg_dir_calcKBlock(methodInfo, indTemp - 1, deltaT, pbParam, domainMesh, constValues, DIAGn, DIAGw);

        gKcurr = cell(domainMesh.number_nodes, 1);
        for indVert = 1 : domainMesh.number_nodes
            gKcurr{indVert} = g(domainMesh.coordinates(indVert, :), indTemp*deltaT);
        end
        gK{indTemp} = cell2mat(gKcurr);
    end
    
    betaK = cell(1, pbParam.nT);
    for indTemp = 1 : pbParam.nT
        betaK{indTemp} = zeros(3*domainMesh.numberTriangles, 1);
        for j = 1 : indTemp
            betaK{indTemp} = betaK{indTemp} + K{indTemp - j + 1} * gK{j};
        end
    end
end