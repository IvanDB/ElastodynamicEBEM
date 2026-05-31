function setupWorkspace(inputStruct)
arguments
    inputStruct (1, 1) struct = struct()
end

if isunix()
    assert(exist('inputStruct', 'var') && isstruct(inputStruct) && inputStruct ~= struct(), ...
        "Input error", "On Linux (HPC), an input structure must be provided. Please use the dedicated bash script.")
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
%FieldList
fieldNames = ["meshType", "meshLevel", "beta", "quadID", "plotFlag", "saveMatx"];

%Set values
for fieldName = fieldNames
    if ~isfield(inputStruct, fieldName)
        assignin('base', fieldName, inputStruct.(fieldName));
    end
end

assignin('base', 'glbIndexFigures', 0);

if ~exist("basePath", "var")
    callStack = dbstack('-completenames');
    assert(length(callStack) > 2, "Internal error", "setupWorkspace shouldn't be called directly but only from main3.") 
    [basePath, ~, ~] = fileparts(callStack(3).file);
    assignin('base', 'basePath', basePath);
end
end