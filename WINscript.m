%WINSCRIPT  Batch-launch eebem.main over the full Cartesian product of parameters listed in a config file.
%   Top-level script (not part of the +eebem package). Reads "parameterLists.config"
%   from this file's own directory via the local helper PARSECONFIGFILE, expands
%   every listed parameter combination into one "inputStruct.<field> = <value>; ..."
%   command string per run via the local helper GENERATECOMMANDBUFFER, then EVALs
%   each command followed by "eebem.main;" in turn -- i.e. it runs one full
%   simulation per combination, sequentially, in the base workspace.
%
%   See also eebem.main, eebem.utility.setupWorkspace

clc
clearvars -except glbIndexFigures

basePath = fileparts(mfilename("fullpath"));
cmdBuffer = generateCommandBuffer(parseConfigFile(basePath));

for cmd = cmdBuffer'
    clearvars inputStruct
    disp("Running " + cmd)
    eval(cmd + " eebem.main;")
end


function cmdBuffer = generateCommandBuffer(configData)
%GENERATECOMMANDBUFFER  Expand a parsed config struct into one inputStruct-assignment command string per run.
%   CMDBUFFER = GENERATECOMMANDBUFFER(CONFIGDATA) maps every recognized CONFIGDATA field
%   (via the internal INPUTFIELDNAMES table, e.g. "pbNames" -> "pbName", "betaMults" ->
%   "betaMult", ...) to a list of "inputStruct.<name> = <value>;" assignment snippets,
%   then builds the full Cartesian product of all such lists (via repeated NDGRID) so
%   that CMDBUFFER contains one concatenated command string per parameter combination,
%   ready to be EVAL'd right before "eebem.main;". Boolean fields set to false are
%   skipped entirely (that inputStruct field is left at its SETUPWORKSPACE default).
%   Quadrature-size fields (numsSRext, numsGHext, ...) are skipped whenever CONFIGDATA
%   also defines "quadIDs", since a quadrature preset index takes precedence.
%
%   Input arguments:
%       CONFIGDATA - (struct) parsed config fields, see PARSECONFIGFILE.
%
%   Output arguments:
%       CMDBUFFER - (string array) one command string per parameter combination.
%
%   See also PARSECONFIGFILE, eebem.utility.setupWorkspace
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

onceOnlyFields = ["poolFlag", "aBldFlag"];
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

    if ismember(fieldName, onceOnlyFields)
        continue
    end

    fieldCmd = "inputStruct." + inputFieldNames(inputFieldNames(:, 1) == fieldName, 2) + "=" + fieldVals + ";";

    [cmdGridF, cmdGridB] = ndgrid(fieldCmd, cmdBuffer);
    cmdBuffer = cmdGridB(:) + " " + cmdGridF(:);
end

for onceOnlyField = onceOnlyFields
    assert(ismember(onceOnlyField, inputFieldNames), sprintf("Inconsistent once-only flag: %s", onceOnlyField))

    if ~ismember(onceOnlyField, fieldNames) 
        continue
    end
    fieldVal = configData.(onceOnlyField);
    
    if islogical(fieldVal) && ~fieldVal
        continue
    end

    fieldCmd = "inputStruct." + inputFieldNames(inputFieldNames(:, 1) == onceOnlyField, 2) + "=" + fieldVal + ";";
    cmdBuffer(1) = cmdBuffer(1) + " " + fieldCmd;
end
end

function configData = parseConfigFile(basePath)
%PARSECONFIGFILE  Parse "parameterLists.config" into a struct via literal MATLAB assignment evaluation.
%   CONFIGDATA = PARSECONFIGFILE(BASEPATH) reads BASEPATH/parameterLists.config
%   line by line, skipping blank lines, "#"-prefixed comment lines and lines
%   containing "()", stripping trailing "# ..." comments, rewriting "(" / ")" to
%   "[" / "]" (so MATLAB array literals can be written with parentheses in the
%   config file), and EVALs each resulting line as "configData.<lineContent>;" --
%   i.e. each config line is expected to already be valid MATLAB assignment syntax
%   once its brackets are normalized, e.g. "pbNames = ["barH1" "barH3"];".
%
%   Input arguments:
%       BASEPATH - (string) directory containing "parameterLists.config".
%
%   Output arguments:
%       CONFIGDATA - (struct) one field per assignment found in the config file.
%
%   Notes:
%       Asserts if the config file cannot be opened. Uses EVAL on file content
%       with only light sanitization (comment stripping, bracket substitution);
%       the config file is trusted input, not arbitrary/untrusted text.
%
%   See also GENERATECOMMANDBUFFER
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

