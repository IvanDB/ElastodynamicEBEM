function matrix = postProcessing_calcMatrix(pb_param, domainMesh, matrix, sp, DeltaT)

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
    TF = domainMesh.coordinates(nodeF, :);
    
    %VERSORE NORMALE all'elemento corrente (triangolo di campo)
    vnF = domainMesh.normal(ind, :);
        
    %Calcolo dell'INTEGRALE sul TRIANGOLO SORGENTE
    ris = postProcessing_calcMatrixBlock(pb_param, TF, vnF, sp, DeltaT)

    %Pulizia "SPORCIZIA NUMERICA" dal blocchetto matriciale corrente
    ris(abs(ris) < 1.0e-14) = 0;
    
    %Posizionamento blocco matriciale corrente in matrix
    matrix{ind} = matrix{ind} + ris;
  
end

%Trasformazione da array di celle a matrice standard
matrix = cell2mat(matrix);
return