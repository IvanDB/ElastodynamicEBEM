function pbParam = readInputFile(basePath, problemFileName)
arguments
    basePath        (1, 1) string
    problemFileName (1, 1) string
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
pbParam.defTimeLimit = sscanf(fgets(problemFile), '%f');

fgets(problemFile);       
pbParam.defNumIntvls = sscanf(fgets(problemFile), '%d');

fgets(problemFile);
pbParam.STcoupling = strcmpi(strtrim(fgets(problemFile)), "true");

fgets(problemFile);       
pbParam.defMeshType = sscanf(fgets(problemFile), '%s');

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