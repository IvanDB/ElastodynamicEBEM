function betaV = calcBetaV(pbParam, domainMesh, matrixV, basePath)
%CALCBETAV  Convolve the single-layer block-Toeplitz matrix V with the Neumann datum history.
%   BETAV = CALCBETAV(PBPARAM, DOMAINMESH, MATRIXV, BASEPATH) samples the exact 
%   Neumann (traction) datum at every triangle centroid and half-integer time instant
%   (via CALCBOUNDDATANEUMANN) and forms, for each time step n, the discrete convolution
%
%       BETAV{n} = sum_{j=1}^{n} MATRIXV{n-j+1} * gV{j},
%
%   the right-hand side used by TIMEMARCHINGDN and TIMEMARCHINGDN_C.
%
%   Input arguments:
%       PBPARAM    - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH - (struct) triangulated boundary mesh, see READSPACEMESH.
%       MATRIXV    - (cell, nT x 1) sparse block-Toeplitz
%                    single-layer matrix produced by CALCMATRIXV.
%       BASEPATH   - (string) project root, forwarded to CALCBOUNDDATANEUMANN.
%
%   Output arguments:
%       BETAV - (cell, nT x 1) each entry a (3*numTriangles x 1) double column vector.
%
%   See also CALCMATRIXV, CALCBOUNDDATANEUMANN, TIMEMARCHINGDN, TIMEMARCHINGDN_C

arguments (Input)
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    matrixV     (:, 1) cell
    basePath    (1, 1) string = "."
end

arguments (Output)
    betaV (:, 1) cell
end

import eebem.core.*

gV = calcBoundDataNeumann(pbParam, domainMesh, basePath);

betaV = cell(pbParam.nT, 1);
parfor indTemp = 1 : pbParam.nT
    betaV{indTemp} = zeros(3*domainMesh.numTriangles, 1);
    for j = 1 : indTemp
        betaV{indTemp} = betaV{indTemp} + matrixV{indTemp - j + 1} * gV{j};
    end
end
return