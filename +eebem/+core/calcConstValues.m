function constValues = calcConstValues(domainMesh, quadData)
%CALCCONSTVALUES  Pre-compute, once, the per-triangle data reused throughout the time-marching loop.
%   CONSTVALUES = CALCCONSTVALUES(DOMAINMESH, QUADDATA) loops over every triangle
%   of DOMAINMESH (in parallel) and pre-computes quantities that do not depend
%   on the time step and can therefore be cached and reused for the whole simulation: 
%   the affine map coefficients (matCoeff, vetCoeff) of the piecewise-linear basis functions,
%   the outer quadrature nodes/weights mapped onto the triangle (EXTn, EXTw),
%   the inner quadrature nodes/weights mapped onto the triangle (INTn, INTw), and, 
%   for every outer node, the three "child" sub-triangles (childVerts, childArea)
%   obtained by connecting the node to the triangle's edges,
%   used by the singular-kernel integration routines (CALCSINGSUBBLOCKV/K/K_C).
%
%   Input arguments:
%       DOMAINMESH - (struct) triangulated boundary mesh, see READSPACEMESH.
%       QUADDATA   - (struct) quadrature nodes/weights and
%                    METHODSPECS, see GENERATEQUADDATA.
%
%   Output arguments:
%       CONSTVALUES - (cell, numTriangles x 1) one struct per triangle with fields
%                     matCoeff, vetCoeff, EXTn, EXTw, INTn, INTw, childVerts, childArea.
%
%   Notes:
%       Uses a PARFOR loop; a parallel pool speeds this step up but is not required.
%
%   See also eebem.utility.generateQuadData, CALCSINGSUBBLOCKV, CALCSINGSUBBLOCKK
arguments (Input)
    domainMesh (1, 1) struct
    quadData   (1, 1) struct
end

arguments (Output)
    constValues (:, 1) cell
end
numExtN = quadData.methodSpecs.numEXT;
numExtW = quadData.methodSpecs.numGHext;

numIntN = quadData.methodSpecs.numINT;
numIntW = quadData.methodSpecs.numGHint;

constValues = cell(domainMesh.numTriangles, 1);

parfor indT = 1 : domainMesh.numTriangles
    vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
    areaT = domainMesh.area(indT);

    %Valore per operatore Rj su test
    l21 = (vertsT(2, :) - vertsT(1, :))';
    l31 = (vertsT(3, :) - vertsT(1, :))';
    n = cross(l21, l31);
    constValues{indT}.matCoeff = [-1, -1, 0; 1, 0, 0; 0, 1, 0] / [l21, l31, n];
    constValues{indT}.vetCoeff = [1; 0; 0] - constValues{indT}.matCoeff * vertsT(1, :)';

    constValues{indT}.EXTn = cell(numExtN, 1);
    constValues{indT}.EXTw = cell(numExtW, 1);
    for indEXT = 1 : numExtN        
        constValues{indT}.EXTn{indEXT} = quadData.EXTn(indEXT, :) * vertsT;    
    end
    for indEXT = 1 : numExtW
        constValues{indT}.EXTw{indEXT} = areaT * quadData.EXTw(indEXT);
    end

    constValues{indT}.GHCnodes = cell(numIntN, 1);
    for indINT = 1 : numIntN
        constValues{indT}.INTn{indEXT} = quadData.INTn(indEXT, :) * vertsT; 
    end
    for indINT = 1 : numIntW
        constValues{indT}.INTw{indINT} = areaT * quadData.INTw(indINT);
    end

    for indEXT = 1 : numExtN
        singularPoint = constValues{indT}.EXTn{indEXT};

        for indChild = 1 : 3
            vertsChild = [vertsT(indChild, :); vertsT(mod(indChild + 1, 3) + 1, :); singularPoint];
            areaChild = norm(cross(vertsT(indChild, :) - singularPoint, vertsT(mod(indChild + 1, 3) + 1, :) - singularPoint)) ./ 2;
            
            constValues{indT}.childVerts{indEXT, indChild} = vertsChild;
            constValues{indT}.childArea{indEXT, indChild} = areaChild;
        end
    end
end
end