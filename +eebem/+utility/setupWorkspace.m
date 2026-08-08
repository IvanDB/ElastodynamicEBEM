function setupWorkspace(inputStruct)
%SETUPWORKSPACE  Translate the user input struct into the base-workspace variables MAIN expects.
%   SETUPWORKSPACE(INPUTSTRUCT) reads the fields of INPUTSTRUCT (as assembled by a 
%   launcher script, e.g. WINSCRIPT, or typed directly before calling eebem.main) 
%   and ASSIGNIN's, in the base workspace, the variables MAIN relies on:
%    - pbIndex/pbSpecs (problem selection),
%    - formSelected (formulation ID),
%    - timeSpecs/meshSpecs (name-value cell arrays for FINALIZEPARAMETERS and the mesh file name),
%    - quadID/quadSpecs (quadrature scheme selection), and
%    - glbFlags (boolean configuration switches: usePool, autoBuild, plotFigs, saveFigs, saveTemp).
%   Also initializes glbIndexFigures and basePath in the base workspace
%   if not already present, inferring basePath from the call stack.
%
%   Input arguments:
%       INPUTSTRUCT - (struct) user configuration. Must contain "pbName" and "formID";
%                     may contain any of betaMult, timeMult, meshType, meshLevel, quadID,
%                     quadType, numSRext, numGHext, numSRint, numGHint, numSNGLR,
%                     numBOUND, usePool, autoBuild, plotFigs, saveFigs, saveTemp.
%
%   Notes:
%       Must be called from eebem.main (or another script one call frame above it): it
%       asserts on the call stack depth when it needs to infer basePath. Side effect only:
%       writes directly into the base workspace via ASSIGNIN and never returns a value.
%
%   See also eebem.main, AUTOBUILD, GENERATEQUADDATA

arguments (Input)
    inputStruct (1, 1) struct
end

%Technical variables 
if ~evalin('base', "exist('glbIndexFigures', 'var')")
    assignin('base', 'glbIndexFigures', 0);
end

if ~evalin('base', "exist('basePath', 'var')")
    callStack = dbstack('-completenames');
    assert(length(callStack) > 1, "Internal error", "setupWorkspace shouldn't be called directly but only from main.") 
    [basePath, ~, ~] = fileparts(callStack(2).file);
    assignin('base', 'basePath', basePath);
end

%Field lists
timeFields = ["betaMult", "timeMult"];
meshFields = ["meshType", "meshLevel"];
quadFields = ["quadType", "numSRext", "numGHext", "numSRint", "numGHint", "numSNGLR", "numBOUND"];
confFields = ["usePool", "autoBuild", "plotFigs", "saveFigs", "saveTemp"];

%Problem data
assert(isfield(inputStruct, "pbName"), "Input error", "A problem name (pbName) must be specified in the input structure.")
assignin('base', "pbIndex", 0);

localStruct = struct();
localStruct.pbName = inputStruct.pbName;
assignin('base', "pbSpecs", namedargs2cell(localStruct));

%Formulation data
assert(isfield(inputStruct, "formID"), "Input error", "Formulation ID (formID) must be specified in the input structure.")
assignin('base', "formSelected", inputStruct.formID);

%Time parameters
localStruct = struct();
for fieldName = timeFields
    if isfield(inputStruct, fieldName)
        localStruct.(fieldName) = inputStruct.(fieldName);
    end
end
assignin('base', "timeSpecs", namedargs2cell(localStruct));


%Space mesh data
localStruct = struct();
for fieldName = meshFields
    if isfield(inputStruct, fieldName)
        localStruct.(fieldName) = inputStruct.(fieldName);
    end
end
assignin('base', "meshSpecs", namedargs2cell(localStruct));

%Quadrature data
assignin('base', "quadID", 0);
if isfield(inputStruct, "quadID")
    assignin('base', "quadID", inputStruct.quadID);
end

localStruct = struct();
for fieldName = quadFields
    if isfield(inputStruct, fieldName)
        localStruct.(fieldName) = inputStruct.(fieldName);
    end
end
assignin('base', "quadSpecs", namedargs2cell(localStruct));

%Configuration flags
localStruct = struct();
for fieldName = confFields
    localStruct.(fieldName) = isfield(inputStruct, fieldName);
end
assignin('base', "glbFlags", localStruct);
end