function betaI = calcBetaI(pbParam, domainMesh, constValues, methodSpecs, basePath)
%CALCBETAI  Assemble the RHS load history for the indirect Dirichlet formulation.
%   BETAI = CALCBETAI(PBPARAM, DOMAINMESH, CONSTVALUES, METHODSPECS, BASEPATH) builds,
%   for every discrete time step n = 1:PBPARAM.nT, the load vector used as right-hand side
%   by TIMEMARCHINGID (and, together with CALCBETAK, by TIMEMARCHINGDD).
%   For each boundary triangle it integrates, over the interval [(n-1)*deltaT, n*deltaT],
%   the exact Dirichlet datum returned by GETDATUMHANDLEDIRICHLET, using
%   the outer quadrature nodes/weights stored in CONSTVALUES.
%   Entries smaller than 1e-14 in magnitude are truncated to zero.
%
%   Input arguments:
%       PBPARAM     - (struct) physical/time-discretization parameters
%                     (deltaT, nT, domainName, ...), see READINPUTFILE.
%       DOMAINMESH  - (struct) triangulated boundary mesh, see READSPACEMESH.
%       CONSTVALUES - (cell, numTriangles x 1) per-triangle geometric
%                     and quadrature data from CALCCONSTVALUES.
%       METHODSPECS - (struct) quadrature scheme sizes, see GENERATEQUADDATA.
%       BASEPATH    - (string) project root, forwarded to GETDATUMHANDLEDIRICHLET
%                     to load custom "<domainName>_D.m" datum files.
%
%   Output arguments:
%       BETAI - (cell, nT x 1) each entry a (3*numTriangles x 1)
%               double column vector, the load at time step n.
%
%   See also GETDATUMHANDLEDIRICHLET, CALCCONSTVALUES, TIMEMARCHINGID, TIMEMARCHINGDD
arguments
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    constValues (:, 1) cell
    methodSpecs (1, 1) struct
    basePath    (1, 1) string = "."
end

import eebem.core.*

numT = domainMesh.numTriangles;

gI = getDatumHandleDirichlet(pbParam, basePath);

betaI = cell(pbParam.nT, 1);
parfor indTemp = 1 : pbParam.nT
    betaI{indTemp} = cell(numT, 1);

    tInz = pbParam.deltaT * (indTemp - 1);
    tFin = pbParam.deltaT * indTemp;

    for indT = 1 : numT
        locVect = zeros(3, 1);
    
        for indSRext = 1 : methodSpecs.numSRext
            for indGHext = 1 : methodSpecs.numGHext
                indEXT = (indSRext - 1) * methodSpecs.numGHext + indGHext;

                locVect = locVect + (gI(constValues{indT}.EXTn{indEXT}, tFin) - gI(constValues{indT}.EXTn{indEXT}, tInz)) * constValues{indT}.EXTw{indGHext};
            end
        end
    
        locVect(abs(locVect) < 1.0e-14) = 0;
    
        betaI{indTemp}{indT} = locVect;
    end

    betaI{indTemp} = cell2mat(betaI{indTemp});
end

end