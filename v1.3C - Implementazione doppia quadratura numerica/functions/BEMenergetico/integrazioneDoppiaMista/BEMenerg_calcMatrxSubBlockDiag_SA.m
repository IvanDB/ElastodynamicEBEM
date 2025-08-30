function matrixSubBlock = BEMenerg_calcMatrxSubBlockDiag_SA(ind, pb_param, domainMesh, hk, dt, gha, ghw)
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

%% ESTRAZIONE INFORMAZIONI DOMINIO

%Estrapolazione INDICI VERTICI triangolo corrente
node = domainMesh.triangles(ind, 1:3);

%Estrapolazione COORDINATE VERTICI triangolo corrente
T = domainMesh.coordinates(node, :);

%Estrapolazione AREA triangolo sorgente corrente
area = domainMesh.area(ind);

%Calcolo VERSORE NORMALE all'elemento corrente
vn = domainMesh.normal(ind, :);

%% CALCOLO SOTTOBLOCCO

ist_temp = hk + [1, 0, -1];

for var = 1 : 3
    %Controllo necessità di calcolo del sottonucleo
    if(ist_temp(var) > 0)
        %Calcolo parametro temporale del sottonucleo
        curr_t = ist_temp(var) * dt;

        %Calcolo dell'integrale relativo al valore temporale curr_t
        intgDoppioParz = BEMenerg_calcIntgExt(pb_param, T, area, T, vn, curr_t, gha, ghw);

        %Aggiornamento pesato del sottoblocco
        matrixSubBlock = matrixSubBlock + (coef(var) .* intgDoppioParz);
    end
end

%Applicazione fattore moltiplicativo
matrixSubBlock = matrixSubBlock ./ (4 .* pi .* pb_param.rho);
return
