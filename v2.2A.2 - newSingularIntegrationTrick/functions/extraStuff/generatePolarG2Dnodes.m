function [nodi3D, pesi3D] = generatePolarG2Dnodes(v3D, rInt, rExt, numNodesR, numNodesT)
    l13 = v3D(1, :) - v3D(3, :);
    n13 = norm(l13);
    v13 = l13 / n13;
    l23 = v3D(2, :) - v3D(3, :);
    n23 = norm(l23);
    v23 = l23 / n23;
    n = cross(v13, v23);
    n = n / norm(n);
    p = cross(n, l13);
    p = p / norm(p);

    alpha2D = acos(dot(v13, v23));
    rMax = max(n13, n23);

    rInt = max([rInt, 0]);
    rExt = min([rExt, rMax]);

    [nodiStdT, pesiStdT] = Gauss1D(1, numNodesT, 0, 0);
    [nodiStdR, pesiStdR] = Gauss1D(1, numNodesR, 0, 0);

    nodiT = (nodiStdT' + 1) ./ 2 .* alpha2D;
    pesiT = pesiStdT' ./ 2 .* alpha2D;
    
    nodiR = (((nodiStdR' + 1) ./ 2) .* (rExt - rInt)) + rInt;
    pesiR = pesiStdR' ./ 2 .* (rExt - rInt);
    
    nodi2D = kron(nodiR, [cos(nodiT), sin(nodiT)]);
    pesi2D = kron(pesiR .* nodiR, pesiT);

    nodi3D = v3D(3, :) + nodi2D(:, 1) .* v13 + nodi2D(:, 2) .* p;
    pesi3D = pesi2D;

    flag = checkExt(v3D, nodi3D);
    
    %plotNodes(v3D, nodi3D, ".r", 5, 100);
    nodi3D(flag, :) = [];
    pesi3D(flag, :) = [];
    %plotNodes(v3D, nodi3D, ".b", 2, 100);
end

function flag = checkExt(V, n)
    v0 = V(2, :) - V(1, :);
    v1 = V(3, :) - V(1, :);
    v2 = n - V(1, :);

    d00 = v0 * v0';
    d01 = v0 * v1';
    d11 = v1 * v1';
    d20 = v2 * v0';
    d21 = v2 * v1';

    l2 = (d11*d20 - d01*d21) / (d00 * d11 - d01^2); 
    l3 = (d00*d21 - d01*d20) / (d00 * d11 - d01^2);
    l1 = 1 - l2 - l3;

    flag = ((l1 < 0) | (l1 > 1) | (l2 < 0) | (l2 > 1) | (l3 < 0) | (l3 > 1));
    return
end
