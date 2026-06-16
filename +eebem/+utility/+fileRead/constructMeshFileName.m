function meshFileName = constructMeshFileName(pbParam, meshSpecs)
    arguments
        pbParam             (1, 1) struct
        meshSpecs.meshType  (1, 1) string {mustBeText} = "" 
        meshSpecs.meshLevel (1, 1) double {mustBeInteger, mustBeNonnegative} = 1
    end

    %If no mesh type is provided we use the default one for the problem
    if(meshSpecs.meshType == "")
        meshSpecs.meshType = pbParam.defaultValues.meshType;
    end

    meshFileName = pbParam.meshName;
    if(meshSpecs.meshType ~= "")
        meshFileName = meshFileName + "_" + meshSpecs.meshType;
    end
    meshFileName = meshFileName + "_" + meshSpecs.meshLevel + ".mesh";
end