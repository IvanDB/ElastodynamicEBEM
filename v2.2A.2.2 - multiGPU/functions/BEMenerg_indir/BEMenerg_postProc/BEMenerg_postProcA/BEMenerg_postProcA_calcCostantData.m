function constValues = BEMenerg_postProcA_calcCostantData(domainMesh)
    
numTriangle = domainMesh.numberTriangles;
constValues = cell(numTriangle, 1);
parfor indT = 1 : numTriangle
    constValues{indT}.vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
    constValues{indT}.vnT = domainMesh.normal(indT, :);
    constValues{indT}.SdR = calcoloSistemaRiferimento(constValues{indT}.vertsT, constValues{indT}.vnT);
end

return