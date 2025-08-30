function sol = BEMenerg_timeMarching(pb_param, domainMesh, gha, ghw)
% INPUT 
%   - pb_param: struct contenente i parametri del problema
%   - domainMesh: struct contente le informazioni sulla mesh
%   - gha: matrice contenente i nodi di quadratura di Gauss-Hammer
%   - ghw: matrice contenente i pesi di quadratura di Gauss-Hammer
%
% OUTPUT:
%   -  sol: matrice contenente la densit� incognita per ogni tempo

%% INIZIALIZZAZIONE COMPONENTI del METODO

%Calcolo del PASSO di DISCRETIZZAZIONE TEMPORALE (deltaT = T_fin / Nt)
dt = pb_param.T_fin / pb_param.Nt;

%Assegnazione INDICE che individua il blocco E0 della MATRICE di TOEPLITZ E
hk = 0;

%Allocazione matrice contenente la SOLUZIONE
sol = zeros(3*domainMesh.number_triangles, pb_param.Nt);

%Il numero di righe � uguale a tre volte il numero di triangoli della 
% discretizzazione spaziale mentre il numero di colonne � uguale al numero 
% di istanti temporali della discretizzazione temporale.

%Allocazione matrice contenente il BLOCCO MATRICIALE E0 
matrix = zeros(3*domainMesh.number_triangles, 3*domainMesh.number_triangles);

%Il numero di righe e di colonne di questa matrice � pari a tre volte il 
% numero di triangoli della discretizzazione spaziale. 
% Possiamo pensare questa matrice come una matrice formata da 
% domainMesh.number_triangles^2 blocchetti 3x3


%Allocazione array contenente il BLOCCO del TERMINE NOTO b0
rhs = zeros(3*domainMesh.number_triangles, 1);

%Inizializzazione a zero parametro che indica l'indice del SOTTOINTERVALLO 
% TEMPORALE a partire dal quale i successivi BLOCCHI MATRICIALI sono NULLI.     
maxstep = 0;

%% CALCOLO BLOCCO MATRICIALE E0 e BLOCCO VETTORIALE b0 

%Calcolo blocco matriciale E0
[matrix, ~] = BEMenerg_calcMatrixBlock(matrix, hk, dt, pb_param, domainMesh, gha, ghw);

%Scrittura su file del blocco E0
file_name = strcat('./temp/', pb_param.domain_type, '_', num2str(pb_param.lev), '_', 'matrix', '_E', num2str(0), '.txt');
writematrix(matrix, file_name, 'Delimiter', ';')

%Fattorizzazione PLU del blocco matriciale E0
[L, U, P] = lu(sparse(matrix));
    
%Calcolo blocco beta0 termine noto
rhs = BEMenerg_calcTnBlock(rhs, hk, dt, pb_param, domainMesh, gha, ghw);

%% RISOLUZIONE PRIMO ISTANTE 
%Risoluzione primo sistema
sol(:, 1) = U\(L\(P*rhs));

%Inizializzazione file temporaneo per salvattaggio varabiale soluzione
file_name = strcat('./temp/', pb_param.domain_type, '_', num2str(pb_param.lev), '_', 'soluzione');
save(file_name, 'sol');

%% RISOLUZIONE ISTANTI SUCCESSIVI

%Inizializzane indice finale sommatoria per rhs
end_ind = 0;

%Ciclo (non-parallelizzabile) sugli instanti temporali
for ind_step = 2 : pb_param.Nt
    
    %Re-inizializzazione a 0 delle variabili matrix e rhs
    matrix = sparse(zeros(3*domainMesh.number_triangles, 3*domainMesh.number_triangles));
    rhs = zeros(3*domainMesh.number_triangles, 1);

    %Aggiornamento indice temporale hk (= ind_step - 1)
    hk = hk + 1;

    %Calcolo termine noto passo corrente
    rhs = BEMenerg_calcTnBlock(rhs, hk, dt, pb_param, domainMesh, gha, ghw);
    
    
    %Aggiornamento rhs
    for ind_inner = 1 : end_ind
        file_name = strcat('./temp/', pb_param.domain_type, '_', num2str(pb_param.lev), '_', 'matrix', '_E', num2str(ind_inner), '.txt');
        matrix = sparse(readmatrix(file_name));
        rhs = rhs - matrix*sol(:, ind_step-ind_inner);
    end
    
    if(maxstep == 0)
        %Reset matrice (serve?)
        matrix = zeros(3*domainMesh.number_triangles, 3*domainMesh.number_triangles);
        %Calcolo blocco matriciale passo corrente
        [matrix, isnotzero] = BEMenerg_calcMatrixBlock(matrix, hk, dt, pb_param, domainMesh, gha, ghw);
        
        %Controllo flag iszero
        if(isnotzero)
            %Aggiornamento finale rhs
            rhs = rhs - matrix*sol(:, 1);

            %Salvataggio su file matrice passo corrente
            file_name = strcat('./temp/', pb_param.domain_type, '_', num2str(pb_param.lev), '_', 'matrix', '_E', num2str(ind_step-1), '.txt');
            writematrix(matrix, file_name, 'Delimiter', ';')
            
            %Aggiornamento indice limite blocchi non nulli calcolati
            end_ind = end_ind + 1;
        else
            maxstep = ind_step;
        end
    end

    %Risoluzione sistema lineare
    sol(:, ind_step) = U\(L\(P*rhs));

    %Aggiornamento salvataggio temporaneo della variabile soluzione
    file_name = strcat('./temp/', pb_param.domain_type, '_', num2str(pb_param.lev), '_', 'soluzione');
    save(file_name, 'sol');   
end

%Salvataggio su file soluzione completa
file_name = strcat('./output/', pb_param.domain_type, '_', num2str(pb_param.lev), '_', 'soluzione', '.txt');
writematrix(sol, file_name, 'Delimiter', ';')
return