function constValues = BEMenerg_dirN_calcCostantData(methodInfo, domainMesh, GHn, GHw, GHCn, GHCw, DIAGn, DIAGw)
    
numTriangle = domainMesh.numberTriangles;

numNodesInt = methodInfo.numNodiInt;
numNodesExt = methodInfo.numNodiExt;
numNodesDiag = methodInfo.numNodiDiag;
numNodesSing = methodInfo.numNodiSing;

constValues = cell(numTriangle, 1);
for indT = 1 : numTriangle
    vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
    areaT = domainMesh.area(indT);


    %Valore per operatore Rj su test
    l21 = (vertsT(2, :) - vertsT(1, :))';
    l31 = (vertsT(3, :) - vertsT(1, :))';
    n = cross(l21, l31);
    constValues{indT}.matCoeff = [-1, -1, 0; 1, 0, 0; 0, 1, 0] / [l21, l31, n];
    constValues{indT}.vetCoeff = [1; 0; 0] - constValues{indT}.matCoeff * vertsT(1, :)';

    constValues{indT}.GHnodes = cell(numNodesExt, 1);
    constValues{indT}.GHweights = cell(numNodesExt, 1);
    for indGH = 1 : numNodesExt        
        constValues{indT}.GHnodes{indGH} = GHn(indGH, :) * vertsT;
        constValues{indT}.GHweights{indGH} = areaT * GHw(indGH);
    end

    constValues{indT}.GHCnodes = cell(numNodesInt, 1);
    for indGHC = 1 : numNodesInt
        nodoMXstd = zeros(1, 3);
        nodoMXstd(1) = GHCn(indGHC, 1);
        nodoMXstd(2) = GHCn(indGHC, 2);
        nodoMXstd(3) = GHCn(indGHC, 3);
        
        constValues{indT}.GHCnodes{indGHC} = nodoMXstd * vertsT;
    end

    constValues{indT}.GHCweights = cell(numNodesSing, 1);
    for indGHC = 1 : numNodesSing
        constValues{indT}.GHCweights{indGHC} = areaT * GHCw(indGHC);
    end 

    for indGH = 1 : numNodesExt
        singularPoint = constValues{indT}.GHnodes{indGH};

        for indChild = 1 : 3
            vertsChild = [vertsT(indChild, :); vertsT(mod(indChild + 1, 3) + 1, :); singularPoint];
            areaChild = norm(cross(vertsT(indChild, :) - singularPoint, vertsT(mod(indChild + 1, 3) + 1, :) - singularPoint)) ./ 2;
            
            %% Minore consumo di memoria - maggior numero di conti
            constValues{indT}.DIAGverts{indGH, indChild} = vertsChild;
            constValues{indT}.DIAGarea{indGH, indChild} = areaChild;
            
            %% Maggior consumo di memoria - minor numero di conti
            % for indDiag = 1 : numNodesDiag
            %     indCurrNode = numNodesDiag*(indChild - 1) + indDiag;
            %     constValues{indT}.DIAGnodes{indGH, indCurrNode} = DIAGn(indDiag, :) * vertsChild;
            %     constValues{indT}.DIAGweights{indGH, indCurrNode} = DIAGw(indDiag) * areaChild;
            % end
        end
    end
end
    
return