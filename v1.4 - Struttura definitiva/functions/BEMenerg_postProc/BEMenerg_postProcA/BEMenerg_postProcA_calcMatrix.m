function matrix = BEMenerg_postProcA_calcMatrix(pbParam, domainMesh, constValues, sourcePoint, deltaT)
% INPUT
% - 
% - 
% OUTPUT
% - 

%Estrazione NUMERO TRIANGOLI DISCRETIZZAZIONE SPAZIALE
numTriangles = domainMesh.numberTriangles;

%%Trasformazione matrice matrix in array di 1 x N_triangles celle
% ciascuna di dimensione 3x3
matrix = zeros(3, 3*numTriangles);
matrix = mat2cell(matrix, 3, 3 * ones(1, numTriangles));

%Ciclo sul numero di triangoli 
for indT = 1 : numTriangles        
    %Calcolo dell'INTEGRALE sul TRIANGOLO SORGENTE
    ris = BEMenerg_postProcA_calcMatrixBlock(pbParam, constValues{indT}, sourcePoint, deltaT);

    %Pulizia "SPORCIZIA NUMERICA" dal blocchetto matriciale corrente
    ris(abs(ris) < 1.0e-14) = 0;
    
    %Posizionamento blocco matriciale corrente in matrix
    matrix{indT} = ris;
  
end

%Trasformazione da array di celle a matrice standard
matrix = cell2mat(matrix);
return