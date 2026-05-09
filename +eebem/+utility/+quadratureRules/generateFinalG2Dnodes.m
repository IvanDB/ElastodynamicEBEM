function [nodi3D, pesi3D] = generateFinalG2Dnodes(v3D, rMin, rInt, rExt, numNodes)

import eebem.utility.quadratureRules.*

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

numNodes1D = sqrt(numNodes);
[nodiStd, pesiStd] = Gauss1D(1, numNodes1D, 0, 0);

%Map to real alpha
nodiT = (nodiStd' + 1) ./ 2 .* alpha2D;
pesiT = pesiStd' ./ 2 .* alpha2D;

lenSeg = (n13 * n23 * sin(alpha2D)) ./ (n23 .* sin(alpha2D - nodiT) + n13 .* sin(nodiT));
rMax = min(lenSeg, rExt);

linV = sort([rMin*ones(numNodes1D, 1), rInt*ones(numNodes1D, 1), rMax], 2);
rMax = rMax .* ones(1, 3);
flags = (linV > rMax);
linV(flags) = rMax(flags);

intV = linV(:, 1 : 2);
radV = linV(:, 2 : 3) - linV(:, 1 : 2);

nodiR = kron(intV, ones(1, numNodes1D)) + kron(radV, (nodiStd + 1) ./ 2);
pesiR = kron(radV, pesiStd ./ 2);

nodi2D = nodiR .* reshape([cos(nodiT), sin(nodiT)], numNodes1D, 1, 2);
pesi2D = pesiR .* nodiR .* pesiT;

nodi2D = reshape(nodi2D, [], 2);
pesi2D = reshape(pesi2D, [], 1);

flags = (pesi2D == 0);
nodi2D(flags, :) = [];
pesi2D(flags) = [];
   
nodi3D = v3D(3, :) + nodi2D(:, 1) .* v13 + nodi2D(:, 2) .* p;
pesi3D = pesi2D;
end
