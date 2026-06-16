clc; clear;

basePath = fileparts(mfilename("fullpath"));
cmdBuffer = generateCommandBuffer(parseConfigFile(basePath));

for cmd = cmdBuffer'
   disp("Running " + cmd)
   eval(cmd + " eebem.main;")
end


function cmdBuffer = generateCommandBuffer(configData)
arguments
    configData (1, 1) struct
end

inputFieldNames = ["formIDs", "formID";         "pbNames", "pbName";        ...
                   "betaMults", "betaMult";     "timeMults", "timeMult";    ...
                   "meshTypes", "meshType";     "meshLvls", "meshLevel";    ...
                   "quadIDs", "quadID";         "quadTypes", "quadType";    ...
                   "numsSRext", "numSRext";     "numsGHext", "numGHext";    ...
                   "numsSRint", "numSRint";     "numsGHint", "numGHint";    ...
                   "numsSNGLR", "numSNGLR";     "numsBOUND", "numBOUND";    ...
                   "poolFlag", "usePool";       "aBldFlag", "autoBuild";    ...
                   "plotFigs", "plotFigs";      "saveFigs", "saveFigs";     ...
                   "saveTemp", "saveTemp"];

quadSpecsFields = ["numsSRext", "numsGHext", "numsSRint", "numsGHint", "numsSNGLR", "numsBOUND"];
isSetQuadID = isfield(configData, "quadIDs");

cmdBuffer = "";

fieldNames = string(fieldnames(configData))';

for fieldName = fieldNames
    if ~ismember(fieldName, inputFieldNames)
        continue
    end

    if isSetQuadID && ismember(fieldName, quadSpecsFields)
        continue
    end

    fieldVals = configData.(fieldName);

    if islogical(fieldVals) && ~fieldVals
        continue
    end

    fieldCmds = "inputStruct." + inputFieldNames(inputFieldNames(:, 1) == fieldName, 2) + "=" + fieldVals + ";";

    [cmdGridF, cmdGridB] = ndgrid(fieldCmds, cmdBuffer);
    cmdBuffer = cmdGridB(:) + " " + cmdGridF(:);
end
end

function configData = parseConfigFile(basePath)
arguments (Input)
    basePath (1, 1) string
end
arguments (Output)
    configData (1, 1) struct
end

configData = struct();

configFileName = "parameterLists.config";
configFilePath = fullfile(basePath, configFileName); 
configFile = fopen(configFilePath);
assert(configFile ~= -1, "Error opening config file")

while configFile
    lineStr = string(fgetl(configFile));
    if lineStr == "-1"
        break
    end
    if lineStr == "" || startsWith(lineStr, "#") || contains(lineStr, "()")
        continue
    end
    if contains(lineStr, "#")
        lineStr = extractBefore(lineStr, "#");
    end

    lineStr = replace(lineStr, "(", "[");
    lineStr = replace(lineStr, ")", "]");

    eval(strrep("configData." + strip(lineStr) + ";", '"', '"""'));
end
fclose(configFile);
end

