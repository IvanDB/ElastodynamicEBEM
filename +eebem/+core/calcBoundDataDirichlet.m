function boundDataDirichlet = calcBoundDataDirichlet(pbParam, domainMesh)
%CALCBOUNDDIRICHLET Summary of this function goes here
%   Detailed explanation goes here

import eebem.core.*

g = getDatumHandleDirichlet(pbParam);

boundDataDirichlet = cell(pbParam.nT, 1);
for indTemp = 1 : pbParam.nT
    gKcurr = cell(domainMesh.number_nodes, 1);
    parfor indVert = 1 : domainMesh.number_nodes
        gKcurr{indVert} = g(domainMesh.coordinates(indVert, :), indTemp*pbParam.deltaT);
    end
    boundDataDirichlet{indTemp} = cell2mat(gKcurr);
end

return