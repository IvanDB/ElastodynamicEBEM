function betaV = calcBetaV(pbParam, domainMesh, matrixV)
%CALCBETAV Summary of this function goes here
%   Detailed explanation goes here

import eebem.core.*

gV = calcBoundDataNeumann(pbParam, domainMesh);

betaV = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaV{indTemp} = zeros(3*domainMesh.numberTriangles, 1);
    for j = 1 : indTemp
        betaV{indTemp} = betaV{indTemp} + matrixV{indTemp - j + 1} * gV{j};
    end
end
return