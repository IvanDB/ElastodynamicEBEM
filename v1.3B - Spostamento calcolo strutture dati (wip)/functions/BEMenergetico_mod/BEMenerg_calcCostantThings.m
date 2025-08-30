function constValues = BEMenerg_calcCostantThings(pb_param, domainMesh, gha)

    constValues = cell(domainMesh.number_triangles, 1);
    parfor indT = 1 : domainMesh.number_triangles
        vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
        vnT = domainMesh.normal(indT, :);
        constValues{indT}.SdR = calcoloSistemaRiferimento(vertsT, vnT);

        constValues{indT}.GHnodes = cell(pb_param.ngh, 1);
        for indGH = 1 : pb_param.ngh
            nodoGHstd = zeros(1, 3);
            nodoGHstd(1) = gha(pb_param.ngh, indGH, 1);
            nodoGHstd(2) = gha(pb_param.ngh, indGH, 2);
            nodoGHstd(3) = gha(pb_param.ngh, indGH, 3);

            constValues{indT}.GHnodes{indGH} = nodoGHstd * vertsT;
        end
    end

return
end