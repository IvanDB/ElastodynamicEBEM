%% INIT PHASE
clearvars -except inputStruct glbIndexFigures basePath cmd cmdBuffer

import eebem.*
format longG
warning off

%Initialize the workspace
assert(exist('inputStruct', 'var'), "Input error", "An input structure must be provided. Please use the dedicated scripts.")
utility.setupWorkspace(inputStruct);

%Build extern functions
utility.autobuild(basePath, glbFlags.autoBuild);

%Start parallel pool
if(glbFlags.usePool)
    delete(gcp("nocreate"));
    parInfo = parpool("Processes");
end

%% SETUP INPUT DATA
problemFileName = utility.fileRead.constructProblemFileName(pbIndex, pbSpecs{:});
pbParam = utility.fileRead.readInputFile(basePath, problemFileName);

meshFileName = utility.fileRead.constructMeshFileName(pbParam, meshSpecs{:});
domainMesh = utility.fileRead.readSpaceMesh(basePath, meshFileName);

glbIndexFigures = utility.plots.plotMesh(domainMesh, glbIndexFigures, glbFlags);


[pbParam, domainMesh] = utility.finalizeParameters(pbParam, domainMesh, timeSpecs{:});

%Check invalid configuration problems -> (Barilli working on it?) 
assert((pbParam.lambda + pbParam.mu ~= 0) || (formSelected == "ID"), "Input error", "Problems with lambda + mu = 0 are not solvable with current implementation of the direct formulations");

%% CORE EXECUTION
% Setup quadrature data
coreQuadData = utility.generateQuadData(quadID, quadSpecs{:});
assert(coreQuadData.methodSpecs.quadType == "FN", sprintf("Quadrature type (%s) non available for the selected formulation (%s)", ...
                                                                    coreQuadData.methodSpecs.quadType, formSelected))

%Construct full file names
fullFileNames = utility.generateFilenames(basePath, formSelected, pbParam, domainMesh, coreQuadData.methodSpecs.stringID);

% Main call
switch formSelected
    case "ID"
        solution = core.timeMarchingID(basePath, pbParam, domainMesh, coreQuadData, fullFileNames);

    case "DD"
        solution = core.timeMarchingDD(basePath, pbParam, domainMesh, coreQuadData, fullFileNames);

    case "DN"
        solution = core.timeMarchingDN(basePath, pbParam, domainMesh, coreQuadData, fullFileNames);

    case "DNc"
        solution = core.timeMarchingDN_c(basePath, pbParam, domainMesh, coreQuadData, fullFileNames);

    case "IN"
        assert(false, "Coming soon...")

    otherwise
        error("Unrecognized formulation")
end

%Plot
glbIndexFigures = utility.plots.plotSolutions(formSelected, pbParam, domainMesh, solution, glbIndexFigures, glbFlags, basePath);

%% POST PROCESSING EXECUTION
assert(true, "Post processing in WIP")

if(exist("solution", "var"))
    density = solution;
elseif(exist(nomeFile, "file"))
    density = load(nomeFile).density;
else
    error("Error message")
end

[postProcInfo, PPn, PPw] = postProc.setupPost();
[xVal, tVal, iVal, numPoints, ~] = postProc.loadPoints(pbParam);

if numPoints == 0
    return
end

switch formSelected
    case "ID"
        vectorField = postProc.postProcV_core(basePath, pbParam, domainMesh, density, numPoints, xVal, tVal, postProcInfo, PPn, PPw);

    case "IN"
        vectorField = postProc.postProcW_core(basePath, pbParam, domainMesh, density, numPoints, xVal, tVal, postProcInfo, PPn, PPw);

    otherwise
        warning("Post processing is WIP")
end

glbIndexFigures = utility.plots.plotSolution(basePath, pbParam, vectorField, xVal, tVal, iVal, glbIndexFigures, formSelected);
return