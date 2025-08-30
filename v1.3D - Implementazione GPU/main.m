%% SETUP WORKSPACE
close all
clear 
clc

format long e
warning off

%Definizione percorsi delle cartelle contenti le functions
addpath(genpath("./functions"))

addpath ./sorgentiGPU

%Avvio pool parallela
delete(gcp("nocreate"));
parInfo = parpool("Threads");

%Definizione 
glbIndexFigures = 0;

%% LETTURA dei PARAMETRI di INPUT

%Selezione file input
%input_file_name = 'input_screenTest_3GH_Integral2.txt';
%input_file_name = 'input_screenTest_3GH_doppioIntegral1D.txt';
%input_file_name = 'input_screenTest_3GH_doppioIntegral1DSpezzato.txt';
%input_file_name = 'input_screenTest_7GH_Integral2.txt';
%input_file_name = 'input_screenTest_7GH_doppioIntegral1D.txt';
%input_file_name = 'input_screenTest_7GH_doppioIntegral1DSpezzato.txt';

%input_file_name = 'input_screenUniformLimit_large.txt';
%input_file_name = 'input_screenUniformReal.txt';
%input_file_name = 'input_screenGradedLimit.txt';
%input_file_name = 'input_sphereUniform.txt';
%input_file_name = 'input_sphereNotUniform.txt';
%input_file_name = 'input_barH1.txt';
%input_file_name = 'input_barH3.txt';

%input_file_name = 'input_screenUniformLimit_maxed.txt';
%input_file_name = 'input_sphereUniform_maxed.txt';
%input_file_name = 'input_barH1_maxed.txt';
%input_file_name = 'input_barH3_maxed.txt';
%input_file_name = 'input_waveOnSphere.txt';


input_file_name = 'input_sphereUniform_small.txt';
%input_file_name = 'input_barH1_small.txt';


%input_file_name = 'input_screenUniformLimit_large07.txt';
%input_file_name = 'input_sphereUniform_large07.txt';
%input_file_name = 'input_barH1_large07.txt';
%input_file_name = 'input_barH3_large19.txt';


%Lettura del file di input
pb_param = readInputFile(input_file_name);

%Controllo implementazione metodo risolutivo richiesto
[err_flag, message] = checkImplementation(pb_param);
if err_flag
    error(message)
    return
end

%% LETTURA e PLOT della MESH SPAZIALE

%Lettura del file contenente la MESH SPAZIALE
domainMesh = readSpaceMesh(pb_param.domain_type, pb_param.lev);

%Plot MESH SPAZIALE
%glbIndexFigures = plotMesh(domainMesh, glbIndexFigures);

%% SCELTA METODI DI CALCOLO
%Selezione tipo integrazione
%typeInt = "SA"; %su cluster assegnato mediante script -> commentare riga
typeInt = "FN"; %su cluster assegnato mediante script -> commentare riga

%Selezione tipo quadratura
%typeQuad = "DQM"; %su cluster assegnato mediante script -> commentare riga
%typeQuad = "DGH"; %su cluster assegnato mediante script -> commentare riga
%typeQuad = "DQC"; %su cluster assegnato mediante script -> commentare riga
typeQuad = "DQC3"; %su cluster assegnato mediante script -> commentare riga

%Selezione tipo postProcessing
%typePost = "A"; %su cluster assegnato mediante script -> commentare riga
%typePost = "N"; %su cluster assegnato mediante script -> commentare riga

%% Costruzione NODI e PESI per FORMULE di QUADRATURA 

%Assegnazione massimo numero di nodi e pesi di quadratura
mxghp = 28;   

%Calcolo matrici contenenti nodi e pesi per formule di Gauss-Hammer
[gha, ghw] = GaussHammer_base(mxghp);

%Assegnazione numero nodi doppioGauss1D
if(typeInt == "FN")
    if(typeQuad == "DQM")
        ngF = 8; %su cluster assegnato mediante script -> commentare riga
        pb_param.ngF = ngF;
        pb_param.nghF = 0;
        pb_param.nghc = 0;
        %Calcolo nodi e pesi doppio Gauss1D
        [dgn, dgw] = doppioGauss1D(ngF);
        [ghcn, ghcw] = GaussHammerComposite(1, 1);
    end
    
    %Assegnazione numero nodi GaussHammer interno
    if(typeQuad == "DGH")
        nghF = 12; %su cluster assegnato mediante script -> commentare riga
        pb_param.nghF = nghF;
        pb_param.ngF = 0;
        pb_param.nghc = 0;
        %Calcolo nodi e pesi doppio Gauss1D
        [dgn, dgw] = doppioGauss1D(1);
        [ghcn, ghcw] = GaussHammerComposite(1, 1);
    end

    if(typeQuad == "DQC")
        nghc = 16; %su cluster assegnato mediante script -> commentare riga
        pb_param.nghc = nghc;
        pb_param.nghF = 0;
        pb_param.ngF = 0;
        %Calcolo nodi e pesi doppio Gauss1D
        [dgn, dgw] = doppioGauss1D(1);
        [ghcn, ghcw] = GaussHammerComposite(nghc, 1);
    end

    if(typeQuad == "DQC3")
        nghc = 16; %su cluster assegnato mediante script -> commentare riga
        pb_param.nghc = nghc;
        pb_param.nghF = 0;
        pb_param.ngF = 0;
        %Calcolo nodi e pesi doppio Gauss1D
        [dgn, dgw] = doppioGauss1D(1);
        [ghcn, ghcw] = GaussHammerComposite(nghc, 3);
    end
end

%% TIME-MARCHING

%Calcolo funzione densità incognita
switch typeInt
    case "SA"
        density = BEMenerg_timeMarching_SA(pb_param, domainMesh, gha, ghw);
    case "FN"
        density = BEMenerg_timeMarching_FN(pb_param, domainMesh, gha, ghw, dgn, dgw, ghcn, ghcw, typeQuad);
end

%Plot della densità ottenuta
glbIndexFigures = plotDensity(pb_param, domainMesh, density, glbIndexFigures);

%% POST-PROCESSING
% 
% %Recupero punti (x, t) di interesse per il post-processing
% [x_val, t_val, i_val] = postProcessing_loadSavedPoints(pb_param);
% 
% %Calcolo dimensioni ed inizializzazione matrice soluzione
% x_sz = size(x_val, 1);
% t_sz = length(t_val);
% 
% sol_raw = zeros(x_sz * t_sz, 3);
% sol = zeros(x_sz, t_sz, 3);
% 
% if(x_sz && t_sz)
%     n_points = x_sz * t_sz;
%     parfor ind = 1 : n_points
%         [ind_x, ind_t] = ind2sub([x_sz, t_sz], ind);
%         switch typeInt
%             case "A"
%                 sol_raw(ind, :) = postProcessing_A(pb_param, domainMesh, density, x_val(ind_x, :), t_val(ind_t));
%             case "N"
%                 sol_raw(ind, :) = postProcessing_N(pb_param, domainMesh, density, x_val(ind_x, :), t_val(ind_t), dgn, dgw);
%         end
%     end
%     sol(:, :, 1) = reshape(sol_raw(:, 1), [x_sz, t_sz]);
%     sol(:, :, 2) = reshape(sol_raw(:, 2), [x_sz, t_sz]);
%     sol(:, :, 3) = reshape(sol_raw(:, 3), [x_sz, t_sz]);
% end
% 
% %Plot soluzione nei punti (x, t) fissati
% glbIndexFigures = plotSolution(pb_param, sol, x_val, t_val, i_val, glbIndexFigures);