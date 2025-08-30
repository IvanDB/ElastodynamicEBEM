%% SETUP WORKSPACE

close all
clear 
clc

format long e

warning off

%Definizione percorsi delle cartelle contenti le functions
addpath ./functions
addpath ./functions/input
addpath ./functions/quadrature
addpath ./functions/graphics
addpath ./functions/BEMenergetico
addpath ./functions/BEMenergetico/integrazioneSingola
addpath ./functions/BEMenergetico/integrazioneDoppia
addpath ./functions/BEMenergetico/integrazioneDoppia/operazioniGeometriche
addpath ./functions/BEMenergetico/integrazioneDoppia/integrazioniAnaliticheEsplicite/
addpath ./functions/BEMenergetico/integrazioneDoppia/integrazioniAnaliticheEsplicite/coefficientiB/
addpath ./functions/BEMenergetico/integrazioneDoppia/integrazioniAnaliticheEsplicite/coefficientiG/
addpath ./functions/postProcessing

glb_index_figures = 0;

%% LETTURA dei PARAMETRI di INPUT

%Selezione file input
%input_file_name = 'input_screenTest_3GH_Integral2.txt';
%input_file_name = 'input_screenTest_3GH_doppioIntegral1D.txt';
%input_file_name = 'input_screenTest_3GH_doppioIntegral1DSpezzato.txt';
%input_file_name = 'input_screenTest_7GH_Integral2.txt';
%input_file_name = 'input_screenTest_7GH_doppioIntegral1D.txt';
%input_file_name = 'input_screenTest_7GH_doppioIntegral1DSpezzato.txt';

%input_file_name = 'input_screenUniformLimit.txt';
%input_file_name = 'input_screenUniformReal.txt';
%input_file_name = 'input_screenGradedLimit.txt';
%input_file_name = 'input_sphereUniform.txt';
input_file_name = 'input_sphereUniform_fast.txt';
%input_file_name = 'input_sphereNotUniform.txt';
%input_file_name = 'input_barH1.txt';
%input_file_name = 'input_barH3.txt';

%input_file_name = 'input_screenUniformLimit_maxed.txt';
%input_file_name = 'input_sphereUniform_maxed.txt';
%input_file_name = 'input_barH1_maxed.txt';
%input_file_name = 'input_barH3_maxed.txt';

%input_file_name = 'input_waveOnSphere.txt';

%Lettura del file di input
pb_param = readInputFile(input_file_name);

%Controllo implementazione metodo risolutivo richiesto
[err_flag, message] = checkImplementation(pb_param);
if err_flag
    error(message)
    return
end

%% LETTURA della MESH SPAZIALE

%Lettura del file contenente la MESH SPAZIALE
domainMesh = readSpaceMesh(pb_param.domain_type, pb_param.lev);

%Calcolo parametri caratteristici mesh
%domainMesh = calc_paramMesh(domainMesh);

%Plot MESH SPAZIALE
glb_index_figures = plotMesh(domainMesh, glb_index_figures);
%% Costruzione NODI e PESI per FORMULE di QUADRATURA 

%Assegnazione massimo numero di nodi e pesi di quadratura
mxghp = 28;   

%Calcolo matrici contenenti nodi e pesi per formule di Gauss-Hammer
[gha, ghw] = GaussHammer(mxghp);

%% TIME-MARCHING

%Calcolo funzione densità incognita
density = BEMenerg_timeMarching(pb_param, domainMesh, gha, ghw);

%% PLOT DENSITY
glb_index_figures = plotDensity(pb_param, domainMesh, density, glb_index_figures);

%% CALCOLO VALORI DI ANALISI

phi_1 = density(1:3:end, end);
phi_1_min = min(phi_1);
phi_1_max = max(phi_1);
max_diff_1 = max(abs(phi_1 - ones(length(phi_1), 1)));

phi_2 = density(2:3:end, end);
phi_2_min = min(phi_2);
phi_2_max = max(phi_2);
max_diff_2 = max(abs(phi_2 - ones(length(phi_2), 1)));

phi_3 = density(3:3:end, end);
phi_3_min = min(phi_3);
phi_3_max = max(phi_3);
max_diff_3 = max(abs(phi_3 - ones(length(phi_3), 1)));

%Salvataggio workspace
name = strcat('./output/Workspace_', pb_param.domain_type, '_', num2str(pb_param.lev), ...
    '_Nodi12_T', num2str(pb_param.T_fin), '_N', num2str(pb_param.Nt), '_Delta1_Beta1');
save(name)

%% POST-PROCESSING

%Recupero punti (x, t) di interesse per il post-processing
[x_val, t_val, i_val] = postProcessing_savedPoints(pb_param);

%Calcolo dimensioni ed inizializzazione matrice soluzione
x_sz = size(x_val, 1);
t_sz = length(t_val);

sol_raw = zeros(x_sz * t_sz, 3);
sol = zeros(x_sz, t_sz, 3);

if(x_sz && t_sz)
    n_points = x_sz * t_sz;
    parfor ind = 1 : n_points
        [ind_x, ind_t] = ind2sub([x_sz, t_sz], ind);
        sol_raw(ind, :) = postProcessing(pb_param, domainMesh, density, x_val(ind_x, :), t_val(ind_t));
    end
    sol(:, :, 1) = reshape(sol_raw(:, 1), [x_sz, t_sz]);
    sol(:, :, 2) = reshape(sol_raw(:, 2), [x_sz, t_sz]);
    sol(:, :, 3) = reshape(sol_raw(:, 3), [x_sz, t_sz]);
end

%Plot soluzione nei punti (x, t) fissati
glb_index_figures = plotSolution(pb_param, sol, x_val, t_val, i_val, glb_index_figures);