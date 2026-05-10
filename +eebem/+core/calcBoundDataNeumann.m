function boundDataNeumann = calcBoundDataNeumann(pbParam, domainMesh)
%CALCBOUNDNEUMANN Summary of this function goes here
%   Detailed explanation goes here

import eebem.core.*

g = getDatumHandleNeumann(pbParam);

boundDataNeumann = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gVcurr = cell(domainMesh.numTriangles, 1);
    parfor indT = 1 : domainMesh.numTriangles
        gVcurr{indT} = g(domainMesh.center(indT, :), (indTemp - 0.5)*pbParam.deltaT, domainMesh.normal(indT, :));
    end
    boundDataNeumann{indTemp} = cell2mat(gVcurr);
end

return