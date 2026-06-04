function fullFileNames = generateFilenames(basePath, formID, pbParam, domainMesh, quadID)
%GENERATEOUTFILENAMES undefined
%   undefined
arguments
    basePath    (1, 1) string
    formID      (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadID      (1, 1) string
end


baseFileName = pbParam.domainName + "_" + formID + "_mesh=" + domainMesh.name + domainMesh.lev ...
                        + "_bMult=" + pbParam.beta + "_TMult=" + pbParam.tMlt ...
                        + "_quad=" + quadID;

tmpvarName = "matrix";
outvarName = "solution";
fullFileNames.tmpFullFilename = fullfile(basePath, "tempData", baseFileName + "_" + tmpvarName + ".mat");
fullFileNames.outFullFilename = fullfile(basePath, "outputData", baseFileName + "_" + outvarName + ".mat");

end