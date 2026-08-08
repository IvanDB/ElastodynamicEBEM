function boundDataDirichlet = calcBoundDataDirichlet(pbParam, domainMesh, basePath)
%CALCBOUNDDATADIRICHLET  Sample the exact Dirichlet datum at every mesh vertex and time step.
%   BOUNDDATADIRICHLET = CALCBOUNDDATADIRICHLET(PBPARAM, DOMAINMESH, BASEPATH)
%   evaluates the datum handle returned by GETDATUMHANDLEDIRICHLET at every vertex of
%   DOMAINMESH and at every discrete time instant t_n = n*deltaT, n = 1:PBPARAM.nT.
%
%   Input arguments:
%       PBPARAM    - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH - (struct) triangulated boundary mesh, see READSPACEMESH.
%       BASEPATH   - (string, optional, default ".") project
%                    root, forwarded to GETDATUMHANDLEDIRICHLET.
%
%   Output arguments:
%       BOUNDDATADIRICHLET - (cell, nT x 1) each entry a
%                            (3*numVertices x 1) double column vector.
%
%   See also GETDATUMHANDLEDIRICHLET, CALCBETAK

arguments (Input)
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    basePath    (1, 1) string = "."
end

arguments (Output)
    boundDataDirichlet (:, 1) cell
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