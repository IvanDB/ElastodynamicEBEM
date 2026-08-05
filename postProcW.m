%% INIT PHASE
clc
clearvars -except glbIndexFigures problemFileName

import eebem.*
format longG
warning off

%Set values
if ~exist('glbIndexFigures', 'var')
    glbIndexFigures = 0;
end

%Obtain the base path 
basePath = fileparts(mfilename("fullpath"));

%% Start parallel pool
delete(gcp("nocreate"));
parInfo = parpool("Processes");

%% SETUP PROBLEM
% problem data
pbParam = utility.fileRead.readInputFile(basePath, problemFileName);
disp("input file read");
domainMesh = utility.fileRead.readSpaceMesh(basePath, pbParam.domainType, pbParam.lev);
disp("mesh file read")
% quad data
coreQuadData = core.setupCore(10);
%load density
form = "IN";
density = load(fullfile(basePath, "outputData", pbParam.domainType, "IN_" + pbParam.domainType + pbParam.lev + "FN_16-1_64-3_256_256" +"_density.mat")).density;
disp("density loaded");

%% START POST PROCESSING

% set up quadrature points
[postProcInfo, PPn, PPw] = postProc.postProc_setup();
disp("Post-Proc setup done");
% load interest points
[xVal, tVal, iVal, numPoints, typePlot] = postProc.postProc_loadPoints(pbParam);
disp("Post-Proc points loaded");

if numPoints == 0
    return
end
% compute field u

vectorField = postProc.postProcW_core(basePath, pbParam, domainMesh, density, numPoints, xVal, tVal, postProcInfo, PPn, PPw);

disp("Post-Proc computations done");
glbIndexFigures = utility.plots.plotSolution(basePath, pbParam, vectorField, xVal, tVal, iVal, glbIndexFigures, form);
disp("Solution plotted");