function constValues = BEMenerg_coreMXGHC_calcCostantData(methodInfo, domainMesh, GHn, GHw, GHCn, GHCw)
    
numTriangle = domainMesh.numberTriangles;
constValues = cell(numTriangle, 1);
parfor indT = 1 : numTriangle
    vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
    vnT = domainMesh.normal(indT, :);
    constValues{indT}.SdR = calcoloSistemaRiferimento(vertsT, vnT);
    areaT = domainMesh.area(indT);

    numNodesExt = methodInfo.numNodiExt;

    constValues{indT}.GHnodes = cell(numNodesExt, 1);
    constValues{indT}.GHweights = cell(numNodesExt, 1)
    for indGHC = 1 : numNodesExt
        nodoGHstd = zeros(1, 3);
        nodoGHstd(1) = GHn(indGHC, 1);
        nodoGHstd(2) = GHn(indGHC, 2);
        nodoGHstd(3) = GHn(indGHC, 3);
        
        constValues{indT}.GHnodes{indGHC} = nodoGHstd * vertsT;

        constValues{indT}.GHweights{indGHC} = areaT * GHw(indGHC);
    end


    numNodesInt = methodInfo.numNodiInt;
    constValues{indT}.GHCnodes = cell(numNodesInt, 1);
    for indGHC = 1 : numNodesInt
        nodoMXstd = zeros(1, 3);
        nodoMXstd(1) = GHCn(indGHC, 1);
        nodoMXstd(2) = GHCn(indGHC, 2);
        nodoMXstd(3) = GHCn(indGHC, 3);
        
        constValues{indT}.GHCnodes{indGHC} = nodoMXstd * vertsT;
    end

    numNodesSing = methodInfo.numNodiSing;
    constValues{indT}.GHCweights = cell(numNodesSing, 1)
    for indGHC = 1 : numNodesSing
        constValues{indT}.GHCweights{indGHC} = areaT * GHCw(indGHC);
    end 
end
    
return