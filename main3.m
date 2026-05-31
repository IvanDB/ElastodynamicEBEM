%% INIT PHASE
clc
clearvars -except inputStruct glbIndexFigures basePath

import eebem.*
format longG
warning off

%Initialize the workspace
utility.setupWorkspace(inputStruct);

%Build extern functions
utility.autobuild(basePath);

%Start parallel pool
% delete(gcp("nocreate"));
% parInfo = parpool("Processes");

%% SETUP INPUT DATA
problemFileName = utility.fileRead.constructProblemFileName(pbIndex, pbSpecs{:});
pbParam = utility.fileRead.readInputFile(basePath, problemFileName);

meshFileName = utility.fileRead.constructMeshFileName(pbParam, meshSpecs{:});
domainMesh = utility.fileRead.readSpaceMesh(basePath, meshFileName);
glbIndexFigures = utility.plots.plotMesh(domainMesh, glbIndexFigures);

[pbParam, domainMesh] = utility.finalizeParameters(pbParam, domainMesh, timeSpecs{:});

%Check invalid configuration problems -> (Barilli working on it?) 
assert((pbParam.lambda + pbParam.mu ~= 0) || (formSelected == "ID"), "Input error", "Problems with lambda + mu = 0 are not solvable with current implementation of the direct formulations");

%% CORE EXECUTION
% Setup quadrature data
coreQuadData = core.setupCore(quadID, quadSpecs{:});

return
% Main call
switch formSelected
    case "ID"
        assert(coreQuadData.methodSpecs.quadType == "FN", sprintf("Quadrature type (%s) for the ID  formulation is WIP", ...
                                                        coreQuadData.methodSpecs.quadType))

        density = core.timeMarchingID(basePath, pbParam, domainMesh, coreQuadData);

        glbIndexFigures = utility.plots.plotConstant(basePath, pbParam, domainMesh, density, glbIndexFigures);

    case "DD"
        assert(coreQuadData.methodSpecs.quadType == "FN", sprintf("Quadrature type (%s) non available for the selected formulation (%s)", ...
                                                        coreQuadData.methodSpecs.quadType, formSelected))

        traction = core.timeMarchingDD(basePath, pbParam, domainMesh, coreQuadData);

        glbIndexFigures = utility.plots.plotConstant(basePath, pbParam, domainMesh, traction, glbIndexFigures);

    case "DN"
        assert(coreQuadData.methodSpecs.quadType == "FN", sprintf("Quadrature type (%s) non available for the selected formulation (%s)", ...
                                                        coreQuadData.methodSpecs.quadType, formSelected))

        displacement = core.timeMarchingDN(basePath, pbParam, domainMesh, coreQuadData);
        
        glbIndexFigures = utility.plots.plotLinear(basePath, pbParam, domainMesh, displacement, glbIndexFigures);

    case "DNc"
        assert(coreQuadData.methodSpecs.quadType == "FN", sprintf("Quadrature type (%s) non available for the selected formulation (%s)", ...
                                                        coreQuadData.methodSpecs.quadType, formSelected))

        displacement = core.timeMarchingDN_c(basePath, pbParam, domainMesh, coreQuadData);
        
        glbIndexFigures = utility.plots.plotConstant(basePath, pbParam, domainMesh, displacement, glbIndexFigures);

    case "IN"
        assert(false, "Coming soon...")

    otherwise
        error("Unrecognized formulation")
end

%% POST PROCESSING EXECUTION
assert(true, "Post processing in WIP")
return