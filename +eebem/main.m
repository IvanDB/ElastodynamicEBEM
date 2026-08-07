%MAIN  Top-level script: run one elastodynamic energetic-BEM simulation end to end.
%   Entry point of the ElastodynamicEBEM pipeline. Requires the base
%   workspace variable INPUTSTRUCT (a struct, see SETUPWORKSPACE) to already
%   be assigned before this script is run -- typically by a launcher such as
%   WINSCRIPT, or by assigning it by hand for a single ad hoc run.
%
%   Pipeline: SETUPWORKSPACE parses INPUTSTRUCT into base-workspace variables; 
%   AUTOBUILD (re)compiles the CUDA/MEX kernels if requested; a parallel pool 
%   is started if GLBFLAGS.usePool; the problem and mesh files are read 
%   (READINPUTFILE, READSPACEMESH) and the mesh is optionally plotted (PLOTMESH); 
%   FINALIZEPARAMETERS derives the final time-discretization; GENERATEQUADDATA builds 
%   the quadrature rules and GENERATEFILENAMES the output paths; the formulation 
%   selected by FORMSELECTED ("ID"|"DD"|"DN"|"DNc"|"IN") is solved by the
%   matching TIMEMARCHING* function; PLOTSOLUTIONS renders/exports the result.
%
%   Notes:
%       The "POST PROCESSING" section near the end of the script is still a placeholder
%       (asserts with a "WIP" message) and does not perform any actual post-processing yet.
%
%   See also SETUPWORKSPACE, AUTOBUILD, READINPUTFILE, READSPACEMESH,
%   FINALIZEPARAMETERS, GENERATEQUADDATA, TIMEMARCHINGID, TIMEMARCHINGDD,
%   TIMEMARCHINGDN, TIMEMARCHINGDN_C, TIMEMARCHINGIN, PLOTSOLUTIONS, WINSCRIPT

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
        solution = core.timeMarchingIN(basePath, pbParam, domainMesh, coreQuadData, fullFileNames);

    otherwise
        error("Unrecognized formulation")
end

%Plot
glbIndexFigures = utility.plots.plotSolutions(formSelected, pbParam, domainMesh, solution, glbIndexFigures, glbFlags, basePath);

%% POST PROCESSING EXECUTION
assert(true, "Post processing in WIP")