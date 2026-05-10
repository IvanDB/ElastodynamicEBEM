function constValues = calcConstValues(domainMesh, quadData)
%CALCCONSTANTVALUES Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    domainMesh struct
    quadData struct
end

arguments (Output)
    constValues cell
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