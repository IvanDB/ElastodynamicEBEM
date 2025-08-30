%% SETUP WORKSPACE

close all
clearvars -except indProblem indMethodCore indMethodPP

format shortG
warning off

%Impostazione parametri
if ~exist('indProblem', 'var')
    indProblem = 0;
end
if ~exist('indMethodCore', 'var')
    indMethodCore = 27;
end
if ~exist('indMethodPP', 'var')
    indMethodPP = 9;
end

%Definizione percorsi delle cartelle contenti le functions
addpath(genpath("./functions"));

%Avvio pool parallela
delete(gcp("nocreate"));
parInfo = parpool("Threads");

%Definizione 
glbIndexFigures = 50*indProblem;

%% LETTURA dei PARAMETRI di INPUT
%Selezione file input
listProblems = ["input_barH1_mid_base.txt", "input_barH1_mid_delta2.txt", "input_barH1_mid_delta4.txt", ...
                "input_barH3_mid_base.txt", "input_barH3_mid_delta2.txt", "input_barH3_mid_delta4.txt", ...
                "input_barH1_mid_2delta.txt", "input_barH3_mid_2delta.txt", ...
               ];
problemFileName = listProblems(indProblem);

%Lettura del file di input
pbParam = readInputFile(problemFileName);

%Controllo implementazione metodo risolutivo richiesto
checkImplementation(pbParam);

%% LETTURA e PLOT della MESH SPAZIALE
%Lettura del file contenente la MESH SPAZIALE
domainMesh = readSpaceMesh(pbParam.domainType, pbParam.lev);

%Plot MESH SPAZIALE
% glbIndexFigures = plotMesh(domainMesh, glbIndexFigures);

%% SCELTA METODO DI CALCOLO e COSTRUZIONE QUADRATURE
%Selezione metoto core
% - "SA     xx":            integrazione esterna mediante GH a xx nodi e integrazione interna analitica  
% - "MX.G2D xx yy":         integrazione esterna mediante GH a xx nodi e integrazione interna mediante G2D a yy nodi 
% - "MX.GHC xx yy z":       integrazione esterna mediante GH a xx nodi e integrazione interna mediante GHC con yy sottoregioni a z nodi ciascuna 
% - "GPU    xx yy z":       esecuzione su GPU - implementato solo metodo MX.GHC xx yy z
% - "FN     xx yy z kk":    metodo interamente numerico  implementata soltanto combinazione GPU xx yy z (non diag) + MX.G2DC xx kk (diag) 
listMethods = ["SA 03", "SA 07", "SA 12", "SA 19", ...
               "MX.GHC 12 01 12", "MX.GHC 12 01 19", "MX.GHC 19 01 12", "MX.GHC 19 01 19", ...
               "MX.G2D 12 64", "MX.G2D 12 256", "MX.G2D 19 64", "MX.G2D 19 256", ...
               "MX.GHC 12 16 3", "MX.GHC 12 64 3", "MX.GHC 19 16 3", "MX.GHC 19 64 3", ...
               "GPU 12 16 3", "GPU 12 64 3", "GPU 19 16 3", "GPU 19 64 3", ...
               "GPU 12 16 19", "GPU 19 16 19", ...
               "FN 12 64 3 64", "FN 12 64 3 256", "FN 12 64 3 1024", ...
               "FN 19 64 3 64", "FN 19 64 3 256", "FN 19 64 3 1024", ...
               "FN 12 16 3 16", "FN 12 16 3 64", "FN 12 16 3 256", ...
               "FN 19 16 3 16", "FN 19 16 3 64", "FN 19 16 3 256", ...
              ];
methodSelected = convertStringsToChars(listMethods(indMethodCore));

[methodInfo, EXTn, EXTw, INTn, INTw, DIAGn, DIAGw] = BEMenerg_setupCore(methodSelected);

%% TIME-MARCHING
%Calcolo funzione densità incognita
switch methodInfo.typeIntg
    case "SA"
        density = BEMenerg_coreSA_timeMarching(pbParam, domainMesh, methodInfo, EXTn, EXTw);
    case "MX.G2D"
        density = BEMenerg_coreMXG2D_timeMarching(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw);
    case "MX.GHC"
        density = BEMenerg_coreMXGHC_timeMarching(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw);
    case "GPU"
        density = BEMenerg_coreGPU_timeMarching(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw);
    case "FN"
        density = BEMenerg_dirN_timeMarching(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw, DIAGn, DIAGw);
end
%Plot della densità ottenuta
%glbIndexFigures = plotSolFormDirNeumann(pbParam, domainMesh, density, glbIndexFigures);

name = convertStringsToChars(problemFileName);
save("./tempData/" + name(1:(end-4)) + "_density", 'density');

fig = figure(glbIndexFigures);
indR = 27 * strcmp(pbParam.domainType, "barH1") + 111 * strcmp(pbParam.domainType, "barH3");
tVal = (1 : pbParam.nT) .* (pbParam.Tfin ./ pbParam.nT);
plot(tVal, density(indR, :))

print(fig, "./tempData/" + name(1:(end-4)), '-dtiff');

%%
% %% SETUP METODO E QUADRATURA POST-PROCESSING
% %Selezione metoto postProcessing
% % - "A":            integrazione analitica 
% % % - "N.GHC yy z":   integrazione mediante GHC con yy sottoregioni a z nodi ciascuna 
% % % - "GPU   yy z":   esecuzione su GPU - implementato solo metodo N.GHC 
% listPost = ["A", ...
%             "N.GHC 01 12", "N.GHC 01 19", ...
%             "N.GHC 16 3", "N.GHC 64 3", ...
%             "GPU 16 3", "GPU 64 3", ...
%             "NEW 16 3", "NEW 64 3", ...
%            ];
% postSelected = convertStringsToChars(listPost(indMethodPP));
% 
% [methodInfo, PPn, PPw] = BEMenerg_setupPost(methodInfo, postSelected);
% 
% %% ESECUZIONE POST-PROCESSING
% %Recupero punti (x, t) di interesse per il post-processing
% [xVal, tVal, iVal, numPoints, typePlot] = BEMenerg_postProc_loadPoints(pbParam);
% 
% % Check per verificare sia richiesto il postProcessing
% if numPoints == 0
%     return
% end
% 
% %Calcolo componenti campo vettoriale incognito nei punti selezionati
% switch methodInfo.typePost 
%     case "A"
%         campoVett = BEMenerg_postProcA_setup(pbParam, domainMesh, density, numPoints, xVal, tVal);
%     case "N.GHC"
%         campoVett = BEMenerg_postProcN_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
%     case "GPU"
%         switch typePlot
%             case "u(:, t)"
%                 campoVett = BEMenerg_postProcGPU_setupLayer(pbParam, domainMesh, density, xVal, tVal, methodInfo, PPn, PPw);
%             case "u(x, :)"
%                 campoVett = BEMenerg_postProcGPU_setupSP(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
%         end
%     case "NEW"
%         campoVett = BEMenerg_postProcNEW_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
% end
% 
% %Plot soluzione nei punti (x, t) fissati
% glbIndexFigures = plotSolution(pbParam, campoVett, xVal, tVal, iVal, typePlot, glbIndexFigures);