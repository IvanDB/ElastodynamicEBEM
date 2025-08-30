function matrix = postProcessing_calcMatrix_N(pb_param, domainMesh, matrix, sp, DeltaT, dgn, dgw)

% INPUT
% - 
% - 
% OUTPUT
% - 

%Estrazione NUMERO TRIANGOLI della MESH
N_triangles = domainMesh.number_triangles;

%%Trasformazione matrice matrix in array di 1 x N_triangles celle
% ciascuna di dimensione 3x3
matrix = mat2cell(matrix, 3, 3 * ones(1, N_triangles));

%Ciclo sul numero di triangoli 
parfor ind = 1 : N_triangles
    
    %Indici dei vertici dell'elemento corrente (triangolo di campo)
    nodeF = domainMesh.triangles(ind, 1:3);

    %Coordinate dei vertici dell'elemento corrente (triangolo di campo)
    vertsT = domainMesh.coordinates(nodeF, :);
    
    %VERSORE NORMALE all'elemento corrente (triangolo di campo)
    areaT = domainMesh.area(ind);
        
    %Calcolo dell'INTEGRALE sul TRIANGOLO SORGENTE
    ris = BEMenerg_calcIntgInt_DG1D(pb_param, sp, vertsT, areaT, DeltaT, dgn, dgw);

    %Applicazione coeffiente costante
    ris = ris / (4 * pi * pb_param.rho);

    %Pulizia "SPORCIZIA NUMERICA" dal blocchetto matriciale corrente
    ris(abs(ris) < 1.0e-14) = 0;
    
    %Posizionamento blocco matriciale corrente in matrix
    matrix{ind} = matrix{ind} + ris;
  
end

%Trasformazione da array di celle a matrice standard
matrix = cell2mat(matrix);
return