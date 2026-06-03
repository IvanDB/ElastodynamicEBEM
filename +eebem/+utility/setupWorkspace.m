function setupWorkspace(inputStruct)
arguments
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
timeFields = ["betaVal", "timeMult"];
meshFields = ["meshType", "meshLevel"];
quadFields = ["quadType", "numSRext", "numGHext", "numSRint", "numGHint", "numSNGLR", "numBOUND"];
confFields = ["aBldFlag", "plotFigs", "saveFigs", "saveTemp"];

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