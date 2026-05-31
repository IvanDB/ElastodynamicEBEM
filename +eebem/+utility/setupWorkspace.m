function setupWorkspace(inputStruct)
arguments
    inputStruct (1, 1) struct = struct()
end

if isunix()
    assert(exist('inputStruct', 'var'), "Input error", "On Linux (HPC), an input structure must be provided. Please use the dedicated bash script.")
    setupHPC(inputStruct);
    return
end

if ispc()
    warning("Windows support is WIP");
    return
end

assert(false, "OS different from Unix and Windows not supported yet.")
end

function setupHPC(inputStruct)
%Tecnical variables 
if ~exist("glbIndexFigures", "var")
    assignin('base', 'glbIndexFigures', 0);
end

if ~exist("basePath", "var")
    callStack = dbstack('-completenames');
    assert(length(callStack) > 2, "Internal error", "setupWorkspace shouldn't be called directly but only from main3.") 
    [basePath, ~, ~] = fileparts(callStack(3).file);
    assignin('base', 'basePath', basePath);
end

%Field lists
timeFields = ["betaVal", "TimeMult"];
meshFields = ["meshType", "meshLevel"];
quadFields = ["quadType", "numSRext", "numGHext", "numSRint", "numGHint", "numSNGLR", "numBOUND"];
confFields = ["plotFigs", "saveFigs", "saveTemp"];

%Problem data
assert(isfield(inputStruct, "pbName"), "Input error", "A problem name (pbName) must be specified in the input structure.")
assignin('base', "pbIndex", 0);

localStruct = struct();
localStruct.pbName = inputStruct.pbName;
assignin('base', "pbSpecs", localStruct);

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
assignin('base', "timeSpecs", localStruct);


%Space mesh data
localStruct = struct();
for fieldName = meshFields
    if isfield(inputStruct, fieldName)
        localStruct.(fieldName) = inputStruct.(fieldName);
    end
end
assignin('base', "meshSpecs", localStruct);

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
assignin('base', "quadSpecs", localStruct);

%Configuration flags
localStruct = struct();
for fieldName = confFields
    if isfield(inputStruct, fieldName)
        localStruct.(fieldName) = inputStruct.(fieldName);
    end
end
assignin('base', "glbFlags", localStruct);
end