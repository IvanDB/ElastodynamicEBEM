function matrixSubBlock = BEMenerg_calcMatrixSubBlock_DQC(indS, indF, pb_param, domainMesh, hk, dt, gha, ghw, ghcn, ghcw)
%INPUT
% - 
% -
%
%OUTPUT
% - 

%% SETUP
%Numero nodi GH
nghS = pb_param.ngh;

%Coefficienti interni
coef = [1, -2, 1];

%Inizializzazione variabile risultato
matrixSubBlock = zeros(3, 3);

%% ESTRAZIONE INFORMAZIONI DOMINI

%Estrapolazione INDICI VERTICI triangolo sorgente corrente
nodeS = domainMesh.triangles(indS, 1:3);

%Estrapolazione COORDINATE VERTICI triangolo sorgente corrente
TS = domainMesh.coordinates(nodeS, :);

%Estrapolazione AREA triangolo sorgente corrente
areaS = domainMesh.area(indS);

%Estrapolazione INDICI VERTICI triangolo campo corrente
nodeF = domainMesh.triangles(indF, 1:3);

%Estrapolazione COORDINATE VERTICI triangolo campo corrente
TF = domainMesh.coordinates(nodeF, :);

%Estrapolazione AREA triangolo campo corrente
areaF = domainMesh.area(indF);

%% CALCOLO INTEGRALI
%Inizializzazione vettore coordinate nodo sorgente
nodoGHSstd = zeros(1, 3);

for indGHS = 1 : nghS
    %Inizializzazione matrice 3x3 tmp_j
    tmp_j = zeros(3, 3);

    %Estrazione coordinate nodo Gauss-Hammer standard
    nodoGHSstd(1) = gha(nghS, indGHS, 1);
    nodoGHSstd(2) = gha(nghS, indGHS, 2);
    nodoGHSstd(3) = gha(nghS, indGHS, 3);

    %Proiezione nodo corrente sul triangolo sorgente
    nodoGHScurr = nodoGHSstd * TS;

    %Calcolo integrale interno
    ist_temp = hk + [1, 0, -1];

    for var = 1 : 3
        %Controllo necessità di calcolo del sottonucleo
        if(ist_temp(var) > 0)
            %Calcolo parametro temporale del sottonucleo
            curr_t = ist_temp(var) * dt;

            tmp_js = BEMenerg_calcIntgInt_GHcomp(pb_param, nodoGHScurr, TF, areaF, curr_t, ghcn, ghcw);

            tmp_j = tmp_j + (coef(var) .* tmp_js);
        end
    end

    %Somma con peso GH esterno
    matrixSubBlock = matrixSubBlock + tmp_j .* ghw(nghS, indGHS);
end    
%Correttivo triangolo sorgente
matrixSubBlock = matrixSubBlock .* 2 .* areaS;

%Applicazione fattore moltiplicativo
matrixSubBlock = matrixSubBlock ./ (4 .* pi .* pb_param.rho);
return