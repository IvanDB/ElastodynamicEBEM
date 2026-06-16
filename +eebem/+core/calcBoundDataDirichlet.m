function boundDataDirichlet = calcBoundDataDirichlet(pbParam, domainMesh, basePath)
%CALCBOUNDDIRICHLET Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    basePath    (1, 1) string = "."
end

import eebem.core.*

g = getDatumHandleDirichlet(pbParam, basePath);

boundDataDirichlet = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gKcurr = cell(domainMesh.numVertices, 1);
    parfor indVert = 1 : domainMesh.numVertices
        gKcurr{indVert} = g(domainMesh.coordinates(indVert, :), indTemp*pbParam.deltaT);
    end
    boundDataDirichlet{indTemp} = cell2mat(gKcurr);
end

return