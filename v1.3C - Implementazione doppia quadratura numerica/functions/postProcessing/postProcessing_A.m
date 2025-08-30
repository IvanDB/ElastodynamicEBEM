function u = postProcessing_A(pb_param, domainMesh, density, x, t)
% INPUT
% - 
% - 
% OUTPUT
% - 

%Allocazione matrice 3x1 contente le componenti del vettore u(x, t) da
% calcolare
u = zeros(3, 1);

%Controllo teorico sul valore temporale
if t > 0
    
    %Calcolo PASSO di DISCRETIZZAZIONE TEMPORALE
    dt = pb_param.T_fin / pb_param.Nt;

    %Calcolio massimo indice n_hat temporale necessario per il calcolo
    n_hat = ceil(t ./ dt) - 1;
    
    %Estrazione NUMERO TRIANGOLI DISCRETIZZAZIONE SPAZIALE
    N_triangles = domainMesh.number_triangles;
    
    %Inizializzazione prima differenza temporale (tra t e t0 = 0)
    DeltaT_next = t; %thk0 = t - t0 = t
        
    %Allocazione matrice 3 x 3*N_triangles costituita N_triangles
    % blocchetti 3x3
    matrix_tcurr = zeros(3, 3*N_triangles);
    matrix_tnext = zeros(3, 3*N_triangles);

    %Calcolo elementi della matrice matrix0
    matrix_tcurr = postProcessing_calcMatrix_A(pb_param, domainMesh, matrix_tcurr, x, DeltaT_next);
    
    %Ciclo sull'indice che individua i sottointervalli temporali
    for ind_temp = 0 : n_hat-1

            %Calcolo differenza temporale tra t e (hk+1)*dt
            DeltaT_next = t - ((ind_temp + 1) * dt);

            %Estrazione valori della densità incognia nel sottointervallo temporale  
            % di indice ind_temp
            alpha = density(:, ind_temp + 1);

            %Pulizia matrice
            matrix_tnext = zeros(3, 3*N_triangles);

            %Calcolo elementi matrice relativa al successivo
            % sottointervallo temporale
            matrix_tnext = postProcessing_calcMatrix_A(pb_param, domainMesh, matrix_tnext, x, DeltaT_next);
             
            %Aggiornamento vettore soluzione 
            u = u + ((matrix_tcurr - matrix_tnext) * alpha);

            %Salvataggio matrice intervallo precedente
            matrix_tcurr = matrix_tnext;
    end
        
    %Estrazione valori della densità nel sottointervallo temporale  
    % di indice n_hat
    alpha = density(:, n_hat + 1);

    %Aggiornamento finael vettore soluzione
    u = u + (matrix_tnext * alpha);
end
return