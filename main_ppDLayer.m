%% SETUP WORKSPACE
clc
% close all
clearvars -except indForm indProblem indMethodCore indMethodPP glbIndexFigures

format longG
warning off

%Impostazione parametri
if ~exist('indForm', 'var')
    indForm = 2;
end
if ~exist('indProblem', 'var')
    indProblem = 18;
end
if ~exist('indMethodCore', 'var')
    indMethodCore = 27;
end
if ~exist('indMethodPP', 'var')
    indMethodPP = 7;
end

%Definizione percorsi delle cartelle contenti le functions
addpath(genpath("./functions"));

%Avvio pool parallela
delete(gcp("nocreate"));
parInfo = parpool("Threads");

%Definizione 
%glbIndexFigures = 0;

%% LETTURA dei PARAMETRI di INPUT
%Selezione file input
listProblems = ["input_screenUniformLimit_small.txt", "input_screenUniformLimit_mid.txt", "input_screenUniformLimit_large.txt", "input_screenUniformLimit_maxed.txt", ...
                "input_sphereUniform_small.txt", "input_sphereUniform_mid.txt", "input_sphereUniform_large.txt", "input_sphereUniform_maxed.txt", ...
                "input_barH1_small.txt", "input_barH1_mid.txt", "input_barH1_large.txt", "input_barH1_maxed.txt", ...
                "input_barH3_small.txt", "input_barH3_mid.txt", "input_barH3_large.txt", "input_barH3_maxed.txt", ...
                "input_waveOnSphere.txt", "input_waveOnSphere_refined.txt", "input_waveOnElemInd.txt"
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

%% SCELTA FORMULAZIONE SELEZIONATA
listForm = ["ID", "DD", "DN"];
formSelected = listForm(indForm);

%Check temporaneo problemi
if (pbParam.lambda + pbParam.mu == 0) && (formSelected ~= "ID")
    error("Problema non adatto a questa formulazione")
end

%% SCELTA METODO DI CALCOLO e COSTRUZIONE QUADRATURE
%Selezione metoto core
% - "SA     xx":            integrazione esterna mediante GH a xx nodi e integrazione interna analitica  
% - "MX.G2D xx yy":         integrazione esterna mediante GH a xx nodi e integrazione interna mediante G2D a yy nodi 
% - "MX.GHC xx yy z":       integrazione esterna mediante GH a xx nodi e integrazione interna mediante GHC con yy sottoregioni a z nodi ciascuna 
% - "GPU    xx yy z":       esecuzione su GPU - implementato solo metodo MX.GHC xx yy z
% - "FN     xx yy z kk":    metodo interamente numerico implementata soltanto combinazione GPU xx yy z (non diag) + MX.G2DC xx kk (diag) 

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
switch formSelected
    case "ID"
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
                density = BEMenerg_coreFN_timeMarching(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw, DIAGn, DIAGw);
            otherwise
                error("Metodo non disponibile per la formulazione selezionata")
        end
        %Plot della densità ottenuta
        % glbIndexFigures = plotDensityV(pbParam, domainMesh, density, glbIndexFigures);
    case "DD"
        switch methodInfo.typeIntg
            case "FN"
                density = BEMenerg_dir_timeMarchingD(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw, DIAGn, DIAGw);
            otherwise
                error("Metodo non disponibile per la formulazione selezionata")
        end
        %Plot della densità ottenuta
        % glbIndexFigures = plotDensityV(pbParam, domainMesh, density, glbIndexFigures);
    case "DN"
        switch methodInfo.typeIntg
            case "FN"
                density = BEMenerg_dir_timeMarchingN(pbParam, domainMesh, methodInfo, EXTn, EXTw, INTn, INTw, DIAGn, DIAGw);
            otherwise
                error("Metodo non disponibile per la formulazione selezionata")
        end
        %Plot della densità ottenuta
        % glbIndexFigures = plotDensityU(pbParam, domainMesh, density, glbIndexFigures);
    otherwise
        error("Formulazione non implementata")
end

%% SETUP METODO E QUADRATURA POST-PROCESSING
%Selezione metoto postProcessing
% - "A":            integrazione analitica 
% - "N.GHC yy z":   integrazione mediante GHC con yy sottoregioni a z nodi ciascuna 
% - "GPU   yy z":   esecuzione su GPU - implementato solo metodo N.GHC 
listPost = ["A", ...
            "N.GHC 01 12", "N.GHC 01 19", ...
            "N.GHC 16 3", "N.GHC 64 3", ...
            "GPU 16 3", "GPU 64 3", ...
            "GPU 16 19", ...
           ];
postSelected = convertStringsToChars(listPost(indMethodPP));

[methodInfo, PPn, PPw] = BEMenerg_setupPost(methodInfo, postSelected);

%% ESECUZIONE POST-PROCESSING
%Recupero punti (x, t) di interesse per il post-processing
[xVal, tVal, iVal, numPoints, typePlot] = BEMenerg_postProc_loadPoints(pbParam);

% Check per verificare sia richiesto il postProcessing
if numPoints == 0
    return
end

%Calcolo componenti campo vettoriale incognito nei punti selezionati
switch formSelected
    case "ID"
        switch methodInfo.typePost 
            case "A"
                campoVett = BEMenerg_postProcA_setup(pbParam, domainMesh, density, numPoints, xVal, tVal);
            case "N.GHC"
                campoVett = BEMenerg_postProcN_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
            case "GPU"
                switch typePlot
                    case "u(:, t)"
                        campoVett = BEMenerg_postProcGPU_setupLayer(pbParam, domainMesh, density, xVal, tVal, methodInfo, PPn, PPw);
                    case "u(x, :)"
                        campoVett = BEMenerg_postProcGPU_setupSP(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
                end
            otherwise
                error("Metodo non disponibile per la formulazione selezionata")
        end
    case "DD"
        switch methodInfo.typePost 
            case "GPU"
                time = tic;
                switch typePlot
                    case "u(:, t)"
                        campoVett = BEMenerg_ppD_setupLayer(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
                    case "u(x, :)"
                        campoVett = BEMenerg_ppD_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
                end
                disp(toc(time));
            otherwise
                error("Metodo non disponibile per la formulazione selezionata")
        end
    case "DN"
        switch methodInfo.typePost 
            case "GPU"
                campoVett = BEMenerg_ppN_setup(pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw);
            otherwise
                error("Metodo non disponibile per la formulazione selezionata")
        end
    otherwise
        error("Formulazione non implementata")
end

%Plot soluzione nei punti (x, t) fissati
glbIndexFigures = plotSolution(pbParam, campoVett, xVal, tVal, iVal, typePlot, glbIndexFigures);
% glbIndexFigures = plotSolution(pbParam, campoVett, xVal, tVal, iVal, "u(x, :)", glbIndexFigures);
% glbIndexFigures = plotSolution(pbParam, campoVett, xVal, tVal, iVal, "u(:, t)", glbIndexFigures);