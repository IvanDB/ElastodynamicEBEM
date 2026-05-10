function betaK = calcBetaK(pbParam, domainMesh, matrixSavedK)
%CALCBETAV Summary of this function goes here
%   Detailed explanation goes here

import eebem.core.*
gK = calcBoundDataDirichlet(pbParam, domainMesh);

%Costruzione vettore betaK
betaK = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    betaK{indTemp} = zeros(3*domainMesh.numTriangles, 1);
    for j = 1 : indTemp
        betaK{indTemp} = betaK{indTemp} + matrixSavedK{indTemp - j + 1} * gK{j};
    end
end
return