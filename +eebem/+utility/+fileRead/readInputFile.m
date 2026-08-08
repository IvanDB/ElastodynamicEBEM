function pbParam = readInputFile(basePath, problemFileName)
%READINPUTFILE  Parse a problem definition file into the PBPARAM struct.
%   PBPARAM = READINPUTFILE(BASEPATH, PROBLEMFILENAME) reads the fixed-format text
%   file BASEPATH/inputFiles/PROBLEMFILENAME (one value per line, each preceded by a
%   label line that is skipped) and returns a struct with the physical parameters
%   (rho, mu, nu, lambda), the derived P/S wave speeds (velP, velS), the default
%   space/time discretization (defaultValues.timeLimit, defaultValues.numIntvls,
%   STcoupling, defaultValues.meshType), the mesh/domain name (meshName, domainName,
%   the latter from PROBLEMFILENAME itself) and the boundary-integral-equation
%   settings (BIE, BOU), validated via CHECKIMPLEMENTATION before returning.
%
%   Input arguments:
%       BASEPATH        - (string) project root.
%       PROBLEMFILENAME - (string) file name under BASEPATH/inputFiles,
%                         see CONSTRUCTPROBLEMFILENAME.
%
%   Output arguments:
%       PBPARAM - (struct) parsed problem parameters, see the fields listed above.
%
%   Notes:
%       Asserts if the problem file cannot be opened. If PBPARAM. meshName
%       is left blank in the file, it defaults to PBPARAM.domainName.
%
%   See also CONSTRUCTPROBLEMFILENAME, CHECKIMPLEMENTATION, READSPACEMESH

arguments (Input)
    basePath        (1, 1) string
    problemFileName (1, 1) string
end

arguments (Output)
    pbParam (1, 1) struct
end

import eebem.utility.fileRead.*

%Open problem file
problemFilePath = fullfile(basePath, "inputFiles", problemFileName);
problemFile = fopen(problemFilePath, 'r');
assert(problemFile ~= -1, "Error opening problem file")

pbParam = struct();
%Set problem name
pbParam.domainName = extractBefore(problemFileName, ".txt");

%Read physical parameters.
fgets(problemFile);
pbParam.rho = sscanf(fgets(problemFile), '%f');

fgets(problemFile);
pbParam.mu = sscanf(fgets(problemFile), '%f');

fgets(problemFile);
pbParam.nu = sscanf(fgets(problemFile), '%f');

fgets(problemFile);
pbParam.lambda = sscanf(fgets(problemFile), '%f');

%Calc speeds of P and S waves
pbParam.velP = sqrt((pbParam.lambda + 2*pbParam.mu) / pbParam.rho); % c_P = sqrt((lamba + 2*mu) / rho)
pbParam.velS = sqrt(pbParam.mu/pbParam.rho);                        % c_S = sqrt(mu / rho) 
     
%Read default space/time discretization parameters           
fgets(problemFile);       
pbParam.defaultValues.timeLimit = sscanf(fgets(problemFile), '%f');

fgets(problemFile);       
pbParam.defaultValues.numIntvls = sscanf(fgets(problemFile), '%d');

fgets(problemFile);
pbParam.STcoupling = strcmpi(strtrim(fgets(problemFile)), "true");

fgets(problemFile);
pbParam.meshName = sscanf(fgets(problemFile), '%s');
%If the no meshName is provided use domainName
if(pbParam.meshName == "")
    pbParam.meshName = pbParam.domainName;
end

fgets(problemFile);       
pbParam.defaultValues.meshType = sscanf(fgets(problemFile), '%s');

%Other parameters (WIP)
fgets(problemFile);       
pbParam.BIE = sscanf(fgets(problemFile), '%s'); 

fgets(problemFile); 
pbParam.BOU = sscanf(fgets(problemFile), '%s'); 

%Close problem file
fclose(problemFile);

%Check for unsupported configurations - TODO: need refactoring
checkImplementation(pbParam);
return