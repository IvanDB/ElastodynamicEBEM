%% SETUP WORKSPACE
close all
clear 
clc

format longG
warning off

%Definizione percorsi delle cartelle contenti le functions
addpath(genpath("./functions"));

%Avvio pool parallela
delete(gcp("nocreate"));
parInfo = parpool("Threads");

%Definizione 
glbIndexFigures = 0;

%% LETTURA dei PARAMETRI di INPUT
%Selezione file input
listProblems = ["input_screenUniformLimit_small.txt", "input_screenUniformLimit_mid.txt", "input_screenUniformLimit_large.txt", "input_screenUniformLimit_maxed.txt", ...
               "input_sphereUniform_small.txt", "input_sphereUniform_mid.txt", "input_sphereUniform_large.txt", "input_sphereUniform_maxed.txt", ...
               "input_barH1_small.txt", "input_barH1_mid.txt", "input_barH1_large.txt", "input_barH1_maxed.txt", ...
               "input_barH3_small.txt", "input_barH3_mid.txt", "input_barH3_large.txt", "input_barH3_maxed.txt", ...
               "input_waveOnSphere.txt"
              ];
problemFileName = listProblems(9);

%Lettura del file di input
pbParam = readInputFile(problemFileName);

%Controllo implementazione metodo risolutivo richiesto
checkImplementation(pbParam);

%% LETTURA e PLOT della MESH SPAZIALE
%Lettura del file contenente la MESH SPAZIALE
domainMesh = readSpaceMesh(pbParam.domainType, pbParam.lev);

%Plot MESH SPAZIALE
glbIndexFigures = plotMesh(domainMesh, glbIndexFigures);

%% SCELTA METODO DI CALCOLO e COSTRUZIONE QUADRATURE
[methodInfo, EXTn, EXTw, INTn, INTw] = BEmenerg_setupCore('GPU 12 64 3');

%% TIME-MARCHING
%Calcolo funzione densità incognita
density = BEMenerg_coreGPU_timeMarching(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw);

%Plot della densità ottenuta
glbIndexFigures = plotDensity(pbParam, domainMesh, density, glbIndexFigures);

%% SETUP METODO E QUADRATURA POST-PROCESSING
[methodInfo, PPn, PPw] = BEmenerg_setupPost(methodInfo, 'GPU 16 3');

%% ESECUZIONE POST-PROCESSING
%Recupero punti (x, t) di interesse per il post-processing
[xVal, tVal, iVal, numPoints, typePlot] = BEMenerg_postProc_loadPoints(pbParam);

% Check per verificare sia richiesto il postProcessing
if numPoints == 0
    return
end

%Calcolo componenti campo vettoriale incognito nei punti selezionati
campoVett = BEMenerg_postProcGPU_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);

%Plot soluzione nei punti (x, t) fissati
glbIndexFigures = plotSolution(pbParam, campoVett, xVal, tVal, iVal, typePlot, glbIndexFigures);