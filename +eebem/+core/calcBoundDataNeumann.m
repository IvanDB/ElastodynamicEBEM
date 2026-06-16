function boundDataNeumann = calcBoundDataNeumann(pbParam, domainMesh, basePath)
%CALCBOUNDNEUMANN Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    basePath    (1, 1) string = "."
end

import eebem.core.*

g = getDatumHandleNeumann(pbParam, basePath);

boundDataNeumann = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gVcurr = cell(domainMesh.numTriangles, 1);
    parfor indT = 1 : domainMesh.numTriangles
        gVcurr{indT} = g(domainMesh.center(indT, :), (indTemp - 0.5)*pbParam.deltaT, domainMesh.normal(indT, :));
    end
    boundDataNeumann{indTemp} = cell2mat(gVcurr);
end

return