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
assert(isfield(inputStruct, "pbName"), "Input error", "Problem name (pbName) must be specified in the input structure.")
assignin('base', "pbSpecs.pbName", inputStruct.pbName);

%Formulation data
assert(isfield(inputStruct, "formID"), "Input error", "Formulation ID (formID) must be specified in the input structure.")
assignin('base', "formSelected", inputStruct.formID);

%Time parameters
for fieldName = timeFields
    if isfield(inputStruct, fieldName)
        assignin('base', "timeSpecs." + fieldName, inputStruct.(fieldName));
    end
end

%Space mesh data
for fieldName = meshFields
    if isfield(inputStruct, fieldName)
        assignin('base', "meshSpecs." + fieldName, inputStruct.(fieldName));
    end
end

%Quadrature data
if isfield(inputStruct, "quadID")
    assignin('base', "quadID", inputStruct.quadID);
    assignin('base', "quadSpecs", struct());
else
    assignin('base', "quadID", 0); %Default value
    for fieldName = quadFields
        if isfield(inputStruct, fieldName)
            assignin('base', "quadSpecs." + fieldName, inputStruct.(fieldName));
        end
    end
end

%Configuration flags
for fieldName = confFields
    if isfield(inputStruct, fieldName)
        assignin('base', "glbFlags." + fieldName, inputStruct.(fieldName));
    end
end
end