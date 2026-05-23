function dirBoundData = BEMenerg_dir_calcBoundDataDirichlet(deltaT, pbParam, domainMesh)
%BEMENERG_CORE_CALCBOUNDDIRICHLET Summary of this function goes here
%   Detailed explanation goes here
%%
switch pbParam.domainType
    case {"barH1", "barH3"}
        h = 1 * strcmp(pbParam.domainType, "barH1") + 3*strcmp(pbParam.domainType, "barH3");
        velP = pbParam.velP;
        hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
        k_val = (0 : hat_k)';
        sol_an = @(x, t) sum((-1).^k_val .* ( (velP.*t - 2.*h.*k_val - (h-x(3))) .* ((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                        - (velP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((velP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
        g = @(x, t) [0; 0; sol_an(x, t)];

    case "sphereWave"
        velP = pbParam.velP;
        ondaIncid = @(x, t) exp(-20 .* ((x(1) - 2 + velP.*t - 0.475).^2));
        g = @(x, t) [ondaIncid(x, t); 0; 0];

    case 'elementoIndustriale'
        g = @(x, t) [0, 0, 0];
    otherwise
        error("Problema non codificato")
end
%%
dirBoundData = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gKcurr = cell(domainMesh.number_nodes, 1);
    parfor indVert = 1 : domainMesh.number_nodes
        gKcurr{indVert} = g(domainMesh.coordinates(indVert, :), indTemp*deltaT);
    end
    dirBoundData{indTemp} = cell2mat(gKcurr);
end

return