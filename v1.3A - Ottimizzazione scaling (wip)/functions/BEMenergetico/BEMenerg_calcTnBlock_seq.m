function rhs = BEMenerg_calcTnBlock_seq(rhs, hk, dt, pb_param, domainMesh, gha, ghw)

%INPUT
% - ...
% - ...
% OUTPUT
% - ...


%% SETUP VARIABILI in ARRAY di CELLE

%Estrazione NUMERO di TRIANGOLI della MESH
N_triangles = domainMesh.number_triangles;

%Trasformazione vettore rhs in array di N_triangles x 1 celle
% ciascuna di dimensione 3x1
rhs = mat2cell(rhs, 3*ones(1, N_triangles), 1);

%Correzione dimensione rhs per errori in parfor (check se ancora presente)
rhs = horzcat(rhs, cell(N_triangles, N_triangles - 1));

%% PARALLELIZZAZIONE CICLO SUI TRIANGOLI SORGENTE E DI CAMPO

%Procedimento element-by-element
for indS = 1 : N_triangles
       
    %---------------------------------------------------------------------
    %Estrapoliamo le INFORMAZIONI utili relative al TRIANGOLO SORGENTE 
    %di INDICE indS
    
    %Estrapolazione INDICI VERTICI triangolo sorgente corrente
    nodeS = domainMesh.triangles(indS, 1:3);
    
    %Estrapolazione COORDINATE VERTICI triangolo sorgente corrente
    TS = domainMesh.coordinates(nodeS, :);
    
    %Estrapolazione AREA triangolo sorgente corrente
    areaS = domainMesh.area(indS);
    
    %Estrapolazione VERSORE NORMALE triangolo sorgente corrente
    vnS = domainMesh.normal(indS, :);
    
    %Estrapolazione ROTORE triangolo sorgente corrente
    curlS = domainMesh.curl(:, :, indS);
    
    %Estrapolazione INDICE TIPO DATO al BORDO triangolo sorgente corrente
    indS_RHS = domainMesh.triangles(indS, 4);

    %Calcolo VETTORE TERMINE NOTO LOCALE
    ris = BEMenerg_calcIntgrGH(pb_param, TS, areaS, indS_RHS, hk, dt, gha, ghw);

    %Pulizia "SPORCIZIA NUMERICA" dal vettore locale
    ris(abs(ris) < 1.0e-14) = 0;

    %Posizionamento vettore locale in rhs
    rhs{indS} = ris;
end


%% FORMATTAZIONE RISULTATI

%Trasformazione da array di celle a vettore standard
rhs = cell2mat(rhs);