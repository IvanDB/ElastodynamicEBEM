%% INIT PHASE
clc
clearvars -except glbIndexFigures indProblem 

import eebem.*
format longG
warning off

%Set values
if ~exist('glbIndexFigures', 'var')
    glbIndexFigures = 0;
end

%Obtain the base path 
basePath = fileparts(mfilename("fullpath"));

%Build functions
utility.autobuild(basePath);

%Start parallel pool
delete(gcp("nocreate"));
parInfo = parpool("Processes");

%% SETUP INPUT DATA
problemFileName = "input_barH1-symm_lev1.txt"; %constructFileName(pbIndex)
pbParam = utility.fileRead.readInputFile(problemFileName);
utility.fileRead.checkImplementation(pbParam);

domainMesh = utility.fileRead.readSpaceMesh(pbParam.domainType, pbParam.lev);
glbIndexFigures = utility.plots.plotMesh(domainMesh, glbIndexFigures);

formSelected = "DD";
%Check invalid configuration problems -> (Barilli working on it?) 
assert((pbParam.lambda + pbParam.mu ~= 0) || (formSelected == "ID"), "Input error", "Problems with lambda + mu = 0 are not solvable with current implementation of the direct formulations");

%% CORE EXECUTION
% Setup quadrature data
coreQuadData = core.setupCore(10);
disp(coreQuadData.methodSpecs.stringID)

% Main call
switch formSelected
    case "ID"
        assert(false, "WIP...")

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

        assert(false, "WIP...")

    case "IN"
        assert(false, "Coming soon...")

    otherwise
        error("Unrecognized formulation")
end

%% POST PROCESSING EXECUTION
assert(true, "Post processing in WIP")
return