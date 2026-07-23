function constValues = BEMenerg_coreSA_calcCostantData(methodInfo, domainMesh, GHn, GHw)
    
numTriangle = domainMesh.numberTriangles;
constValues = cell(numTriangle, 1);
parfor indT = 1 : numTriangle
    vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
    vnT = domainMesh.normal(indT, :);
    constValues{indT}.SdR = calcoloSistemaRiferimento(vertsT, vnT);
    areaT = domainMesh.area(indT);

    numNodes = methodInfo.numNodiExt;

    constValues{indT}.GHnodes = cell(numNodes, 1);
    constValues{indT}.GHweights = cell(numNodes, 1)
    for indGH = 1 : numNodes
        nodoGHstd = zeros(1, 3);
        nodoGHstd(1) = GHn(indGH, 1);
        nodoGHstd(2) = GHn(indGH, 2);
        nodoGHstd(3) = GHn(indGH, 3);
        
        constValues{indT}.GHnodes{indGH} = nodoGHstd * vertsT;

        constValues{indT}.GHweights{indGH} = areaT * GHw(indGH);
    end
end

return