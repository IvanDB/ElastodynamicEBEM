function neuBoundData = BEMenerg_dir_calcBoundDataNeumann(deltaT, pbParam, domainMesh)
%BEMENERG_CORE_CALCBOUNDNEUMANN Summary of this function goes here
%   Detailed explanation goes here

switch pbParam.domainType
    case {"barH1", "barH3", "barH3sim"}
        h = 1 * strcmp(pbParam.domainType, "barH1") + 3*strcmp(pbParam.domainType, "barH3") + 3*strcmp(pbParam.domainType, "barH3sim");
        velP = pbParam.velP;
        hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
        k_val = (0 : hat_k)';
        tildeP = @(x, t) sum((-1).^k_val .* (((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                        + ((velP.*t - 2.*h.*(k_val+1) + (h - x(3))) > 0)));
        g = @(x, t, n) [0; 0; 2*pbParam.mu * tildeP(x, t) * n(3)];
    
    case {"elementoIndustriale", "elementoIndustrialeSemplificato"}
	    g = @(x, t, n) [0; 0; ((x(3) > -0.02) - (x(3) < -0.35)) * (abs(n(3)) > 0.5)];

    otherwise
        error("Problema non codificato")
end

neuBoundData = cell(1, pbParam.nT);
for indTemp = 1 : pbParam.nT
    gVcurr = cell(domainMesh.numberTriangles, 1);
    parfor indT = 1 : domainMesh.numberTriangles
        gVcurr{indT} = g(domainMesh.center(indT, :), indTemp*deltaT, domainMesh.normal(indT, :));
    end
    neuBoundData{indTemp} = cell2mat(gVcurr);
end

return