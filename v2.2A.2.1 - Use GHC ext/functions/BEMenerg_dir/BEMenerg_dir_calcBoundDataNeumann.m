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
    
    case "DesCop-cube"
        a = 0.1;
        b = 100;
        cP = pbParam.velP;
        cS = pbParam.velS;
        
        qP = @(r, t) sqrt(a) * (a*b - t + (r / cP));
        qS = @(r, t) sqrt(a) * (a*b - t + (r / cS));
        F = @(r, t) (1 / (2 * a * (r^2))) * (exp(-(qP(r, t)^2)) - exp(-(qS(r, t)^2)) + sqrt(a * pi) * (a*b - t) * (erf(qP(r, t)) - erf(qS(r, t))));
        
        fP = @(r, t) exp(-a * ((t - (r / cP) - (a*b))^2));
        fS = @(r, t) exp(-a * ((t - (r / cS) - (a*b))^2));
        
        fS_t = @(r, t) -a * 2 * (t - (r / cS) - (a*b)) * exp(-a * ((t - (r / cS) - (a*b))^2));
        fP_t = @(r, t) -a * 2 * (t - (r / cP) - (a*b)) * exp(-a * ((t - (r / cP) - (a*b))^2));
        
        d = @(i, j) (i == j);
        x0 = [1, 1, 1];
        
        p_ijk = @(r, t, i, j, k) -6 * (cS^2) * ((5*r(i)*r(j)*r(k) / norm(r)^5) - ((d(i, j)*r(k) + d(i, k)*r(j) + d(j, k)*r(i)) / norm(r)^3)) * F(norm(r), t) ...
                                    + 2 * ((6*r(i)*r(j)*r(k) / norm(r)^5) - ((d(i, j)*r(k) + d(i, k)*r(j) + d(j, k)*r(i)) / norm(r)^3)) * (fS(norm(r), t) - (cS^2 / cP^2)*fP(norm(r), t)) ...
                                    + (2 / cS) * (r(i)*r(j)*r(k) / norm(r)^4) * (fS_t(norm(r), t) - (cS^3 / cP^3)*fP_t(norm(r), t)) ...
                                    - (d(i, j) * r(k) / norm(r)^3) * (1 - 2*(cS^2 / cP^2)) * (fP(norm(r), t) + (norm(r) / cP) * fP_t(norm(r), t)) ...
                                    - ((d(i, k) * r(j) + d(j, k)*r(i)) / norm(r)^3) * (fS(norm(r), t) + (norm(r) / cS) * fS_t(norm(r), t));
        
        
        p_ij = @(r, t, i, j) p_ijk(r, t, i, j, 1) + p_ijk(r, t, i, j, 2) + p_ijk(r, t, i, j, 3);
        
        p_i = @(r, t, n, i) n(1) * p_ij(r, t, i, 1) + n(2) * p_ij(r, t, i, 2) + n(3) * p_ij(r, t, i, 3);
        
        p = @(r, t, n) (1/(4*pi)) .* [p_i(r, t, n, 1); p_i(r, t, n, 2); p_i(r, t, n, 3)];

        g = @(x, t, n) p(x - x0, t, n);

    case "DesCop-sphere"  
        p0 = 143.5;
        a = 1279;
        b = 12792;
        
        g = @(x, t, n) p0 * (exp(-a * t) - exp(-b * t)) .* n';
    case "elementoIndustriale"
        g = @(x, t, n) [0; 0; ((x(3) > -0.002) - (x(3) < -0.35)) * (abs(n(3)) > 0.5)];
    otherwise
        error("Problema non codificato")
end

neuBoundData = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gVcurr = cell(domainMesh.numberTriangles, 1);
    parfor indT = 1 : domainMesh.numberTriangles
        gVcurr{indT} = g(domainMesh.center(indT, :), (indTemp - 0.5)*deltaT, domainMesh.normal(indT, :));
    end
    neuBoundData{indTemp} = cell2mat(gVcurr);
end

return