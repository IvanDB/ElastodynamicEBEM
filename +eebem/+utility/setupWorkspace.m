if isunix()
    assert(exist('inputStruct', 'var') && isstruct(inputStruct), "Input error", "On HPC, an input structure must be provided. Please use the dedicated bash script.")
    setupHPC(inputStruct);
end

if ispc()
    assert(false, "WIP");
end

assert(false, "OS different from Unix and Windows not supported yet.")
return

function setupHPC(inputStruct)
%FieldList
fieldNames = ["meshType", "meshLevel", "beta", "quadID", "plotFlag", "saveMatx"];

%Default values
defaultData.meshLevel = 1;
defaultData.meshType = "symm";
defaultData.beta = 1;
defaultData.quadID = 10;

defaultData.plotFlag = true;
defaultData.saveMatx = false;

%Set values
for fieldName = fieldNames
    if isfield(inputStruct, fieldName)
        assignin('base', fieldName, inputStruct.(fieldName));
        continue
    end
    assignin('base', fieldName, defaultData.(fieldName));  
end

if ~exist('glbIndexFigures', 'var')
    glbIndexFigures = 0;
end

if ~exist("basePath", "var")
    callStack = dbstack('-completenames');
    assert(length(callStack) > 2, "Internal error", "setupWorkspace shouldn't be called directly but only from main3.") 
    [basePath, ~, ~] = fileparts(callStack(3).file);
    assignin('base', 'basePath', basePath);
end
end