function constValues = BEMenerg_postProcN_calcCostantData(domainMesh, methodInfo, nStd, wStd)
    
numTriangle = domainMesh.numberTriangles;
constValues = cell(numTriangle, 1);
parfor indT = 1 : numTriangle
    vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
    
    areaT = domainMesh.area(indT);

    numNodes = methodInfo.numNodiPP;

    constValues{indT}.PPnodes = cell(numNodes, 1);
    for indNodo = 1 : numNodes
        nodoPPstd = zeros(1, 3);
        nodoPPstd(1) = nStd(indNodo, 1);
        nodoPPstd(2) = nStd(indNodo, 2);
        nodoPPstd(3) = nStd(indNodo, 3);
        
        constValues{indT}.PPnodes{indNodo} = nodoPPstd * vertsT;
    end

    numNodesSing = methodInfo.numNodiSingPP;
    constValues{indT}.PPweights = cell(numNodesSing, 1)
    for indNodo = 1 : numNodesSing
        constValues{indT}.PPweights{indNodo} = areaT * wStd(indNodo);
    end 

end

return