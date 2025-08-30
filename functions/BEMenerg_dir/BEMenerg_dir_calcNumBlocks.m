function numBlocks = BEMenerg_dir_calcNumBlocks(pbParam, domainMesh)

numV = domainMesh.number_nodes;
deltaT = pbParam.Tfin / pbParam.nT;

dists = zeros(numV, numV);
parfor indExt = 1 : numV
    for indInt = 1 : numV
        vettDistC = domainMesh.coordinates(indExt, :) - domainMesh.coordinates(indInt, :);
        dists(indExt, indInt) = sqrt(sum(vettDistC.^2));
    end
end

maxDist = max(dists, [], "all");

maxIndTemp = ceil(maxDist / (pbParam.velS * deltaT)) + 1;
numBlocks = min(pbParam.nT, maxIndTemp);

end

