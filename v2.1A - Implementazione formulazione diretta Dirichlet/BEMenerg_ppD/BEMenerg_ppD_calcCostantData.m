function constValues = BEMenerg_ppD_calcCostantData(domainMesh)
    
numTriangle = domainMesh.numberTriangles;

constValues = cell(numTriangle, 1);
for indT = 1 : numTriangle
    vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);

    %Valore per operatore Rj su test
    l21 = (vertsT(2, :) - vertsT(1, :))';
    l31 = (vertsT(3, :) - vertsT(1, :))';
    PPn = domainMesh.normal(indT, :)';
    constValues{indT}.matCoeff = [-1, -1, 0; 1, 0, 0; 0, 1, 0] / [l21, l31, PPn];
    constValues{indT}.vetCoeff = [1; 0; 0] - constValues{indT}.matCoeff * vertsT(1, :)';
end
    
return