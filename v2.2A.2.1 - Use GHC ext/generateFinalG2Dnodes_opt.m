function [nodi3D, pesi3D] = generateFinalG2Dnodes_opt(v3D, rMin, rInt, rExt, numNodesR, numNodesT)
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

    [nodiStdT, pesiStdT] = Gauss1D(1, numNodesT, 0, 0);
    [nodiStdR, pesiStdR] = Gauss1D(1, numNodesR, 0, 0);

    %Map to real alpha
    nodiT = (nodiStdT' + 1) ./ 2 .* alpha2D;
    pesiT = pesiStdT' ./ 2 .* alpha2D;

    lenSeg = (n13 * n23 * sin(alpha2D)) ./ (n23 .* sin(alpha2D - nodiT) + n13 .* sin(nodiT));
    rMax = min(lenSeg, rExt);

    linV = sort([rMin*ones(numNodesR, 1), rInt*ones(numNodesR, 1), rMax], 2);
    rMax = rMax .* ones(1, 3);
    flags = (linV > rMax);
    linV(flags) = rMax(flags);
    
    intV = linV(:, 1 : 2);
    radV = linV(:, 2 : 3) - linV(:, 1 : 2);
    
    nodiR = kron(intV, ones(1, numNodesR)) + kron(radV, (nodiStdR + 1) ./ 2);
    pesiR = kron(radV, pesiStdR ./ 2);

    nodi2D = nodiR .* reshape([cos(nodiT), sin(nodiT)], numNodesT, 1, 2);
    pesi2D = pesiR .* nodiR .* pesiT;

    nodi2D = reshape(nodi2D, [], 2);
    pesi2D = reshape(pesi2D, [], 1);

    flags = pesi2D == 0;
    nodi2D(flags, :) = [];
    pesi2D(flags) = [];
       
    nodi3D = v3D(3, :) + nodi2D(:, 1) .* v13 + nodi2D(:, 2) .* p;
    pesi3D = pesi2D;

    %plotNodes(v3D, nodi3D, ".r", 5, 100);
end
