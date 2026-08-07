function boundDataNeumann = calcBoundDataNeumann(pbParam, domainMesh, basePath)
%CALCBOUNDDATANEUMANN  Sample the exact Neumann datum at every triangle centroid and time step.
%   BOUNDDATANEUMANN = CALCBOUNDDATANEUMANN(PBPARAM, DOMAINMESH, BASEPATH)
%   evaluates the datum handle returned by GETDATUMHANDLENEUMANN at the
%   centroid and outward normal of every triangle of DOMAINMESH, at the
%   mid-point time instants t_n = (n-0.5)*deltaT, n = 1:PBPARAM.nT.
%
%   Input arguments:
%       PBPARAM    - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH - (struct) triangulated boundary mesh, see READSPACEMESH.
%       BASEPATH   - (string, optional, default ".") project
%                    root, forwarded to GETDATUMHANDLENEUMANN.
%
%   Output arguments:
%       BOUNDDATANEUMANN - (cell, nT x 1) each entry a
%                          (3*numTriangles x 1) double column vector.
%
%   See also GETDATUMHANDLENEUMANN, CALCBETAV, TIMEMARCHINGIN
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