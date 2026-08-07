function betaK = calcBetaK(pbParam, domainMesh, matrixK, basePath)
%CALCBETAK  Convolve the double-layer block-Toeplitz matrix K with the Dirichlet datum history.
%   BETAK = CALCBETAK(PBPARAM, DOMAINMESH, MATRIXK, BASEPATH) samples the exact
%   Dirichlet datum at every mesh vertex and time instant (viaCALCBOUNDDATADIRICHLET) 
%   and forms, for each time step n, the discrete convolution
%
%       BETAK{n} = sum_{j=1}^{n} MATRIXK{n-j+1} * gK{j}.
%
%   The result is combined with BETAI to build the right-hand side used by TIMEMARCHINGDD.
%
%   Input arguments:
%       PBPARAM    - (struct) physical/time-discretization parameters, see READINPUTFILE.
%       DOMAINMESH - (struct) triangulated boundary mesh, see READSPACEMESH.
%       MATRIXK    - (cell, nT x 1) sparse block-Toeplitz
%                    double-layer matrix produced by CALCMATRIXK.
%       BASEPATH   - (string) project root, forwarded to CALCBOUNDDATADIRICHLET.
%
%   Output arguments:
%       BETAK - (cell, nT x 1) each entry a (3*numTriangles x 1) double column vector.
%
%   See also CALCMATRIXK, CALCBOUNDDATADIRICHLET, CALCBETAI, TIMEMARCHINGDD
arguments
    pbParam     struct
    domainMesh  struct
    matrixK     (:, 1) cell
    basePath    (1, 1) string = "."
end

import eebem.core.*
gK = calcBoundDataDirichlet(pbParam, domainMesh, basePath);

betaK = cell(pbParam.nT, 1);
parfor indTemp = 1 : pbParam.nT
    betaK{indTemp} = zeros(3*domainMesh.numTriangles, 1);
    for j = 1 : indTemp
        betaK{indTemp} = betaK{indTemp} + matrixK{indTemp - j + 1} * gK{j};
    end
end
return