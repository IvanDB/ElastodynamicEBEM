function constValues = BEMenerg_coreFN_calcCostantData(methodInfo, domainMesh, GHn, GHw, DIAGn, DIAGw)
    
numTriangle = domainMesh.numberTriangles;

numNodesExt = methodInfo.numNodiExt;
numNodesDiag = methodInfo.numNodiDiag;

constValues = cell(numTriangle, 1);
parfor indT = 1 : numTriangle
    vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
    areaT = domainMesh.area(indT);

    constValues{indT}.GHnodes = cell(numNodesExt, 1);
    constValues{indT}.GHweights = cell(numNodesExt, 1);
    for indGH = 1 : numNodesExt        
        constValues{indT}.GHnodes{indGH} = GHn(indGH, :) * vertsT;
        constValues{indT}.GHweights{indGH} = areaT * GHw(indGH);
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