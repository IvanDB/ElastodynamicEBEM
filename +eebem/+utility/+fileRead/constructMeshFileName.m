function meshFileName = constructMeshFileName(pbParam, meshSpecs)
    arguments
        pbParam             (1, 1) struct
        meshSpecs.meshType  (1, 1) string {mustBeText} = "" 
        meshSpecs.meshLevel (1, 1) double {mustBeInteger, mustBePositive} = 1
    end

    %If no mesh type is provided we use the default one for the problem
    if(meshSpecs.meshType == "")
        meshSpecs.meshType = pbParam.defMeshType;
    end

    meshFileName = pbParam.domainName;
    if(meshSpecs.meshType ~= "")
        meshFileName = meshFileName + "-" + meshSpecs.meshType;
    end
    meshFileName = meshFileName + "_" + meshSpecs.meshLevel + ".mesh";
end