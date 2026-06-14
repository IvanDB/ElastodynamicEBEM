function betaV = calcBetaV(pbParam, domainMesh, matrixV, basePath)
%CALCBETAV Summary of this function goes here
%   Detailed explanation goes here
arguments
    pbParam     struct
    domainMesh  struct
    matrixV     (:, 1) cell
    basePath    (1, 1) string = "."
end

import eebem.core.*

gV = calcBoundDataNeumann(pbParam, domainMesh, basePath);

betaV = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaV{indTemp} = zeros(3*domainMesh.numTriangles, 1);
    for j = 1 : indTemp
        betaV{indTemp} = betaV{indTemp} + matrixV{indTemp - j + 1} * gV{j};
    end
end
return