function constValues = BEMenerg_coreMXG2D_calcCostantData(methodInfo, domainMesh, GHn, GHw, G2Dn, G2Dw)
    
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
    constValues{indT}.G2Dnodes = cell(numNodesInt, 1);
    constValues{indT}.G2Dweights = cell(numNodesInt, 1);
    for indGHC = 1 : numNodesInt
        nodoMXstd = zeros(1, 3);
        nodoMXstd(1) = G2Dn(indGHC, 1);
        nodoMXstd(2) = G2Dn(indGHC, 2);
        nodoMXstd(3) = G2Dn(indGHC, 3);
        
        constValues{indT}.G2Dnodes{indGHC} = nodoMXstd * vertsT;
        constValues{indT}.G2Dweights{indGHC} = areaT * G2Dw(indGHC);
    end   
end
    
return