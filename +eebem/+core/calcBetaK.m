function betaK = calcBetaK(pbParam, domainMesh, matrixK, basePath)
%CALCBETAV Summary of this function goes here
%   Detailed explanation goes here
arguments
    pbParam     struct
    domainMesh  struct
    matrixK     (:, 1) cell
    basePath    (1, 1) string = "."
end

import eebem.core.*
gK = calcBoundDataDirichlet(pbParam, domainMesh, basePath);

betaK = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaK{indTemp} = zeros(3*domainMesh.numTriangles, 1);
    for j = 1 : indTemp
        betaK{indTemp} = betaK{indTemp} + matrixK{indTemp - j + 1} * gK{j};
    end
end
return