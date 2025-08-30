function domainMesh = readSpaceMesh(domain_type, lev)

% INPUT 
%   - domain_type: stringa contenente il TIPO di DOMINIO
%   - lev: intero contente il LIVELLO di RAFFINAMENTO
%
% OUTPUT:
%   -  domainMesh: struct contente i dettagli della mesh

   
%% APERTURA del FILE MESH
%Apertura del file mesh
file_mesh_path = strcat('./data/mesh/', domain_type, '/', domain_type, '_', num2str(lev), '.mesh');
mesh_file = fopen(file_mesh_path, 'r');
if mesh_file == -1
    error("Impossibile aprire file di mesh")
end

%Lettura dell'intestazione del file mesh
line = fgets(mesh_file);
line = fgets(mesh_file);
line = fgets(mesh_file);
line = fgets(mesh_file);

%% LETTURA dei NODI della MESH
%Lettura linea commentata introduzione sezione dei nodi
line = fgets(mesh_file);

%Lettura NUMERO di NODI
domainMesh.number_nodes = sscanf(fgets(mesh_file), '%d');

%Allocazione MATRICE contenente le COORDINATE dei NODI
domainMesh.coordinates = zeros(domainMesh.number_nodes, 4);

%Ogni riga della matrice è relativa ad un determinato nodo: i primi 
% tre elementi della riga contengono le coordinate del nodo, mentre
% l'ultimo elemento indica a quale faccia appartiene il nodo

 
%FORMATO in cui salvare le COORDINATE dei NODI
formatSpec = '%f';

%DIMENSIONE struttura in cui salvare i dati
sizeSpec = [4 4*domainMesh.number_nodes];

%Lettura dei nodi della mesh
domainMesh.coordinates = fscanf(mesh_file, formatSpec, sizeSpec)';

%Eliminazione quarta colonna (inutile)
domainMesh.coordinates(:, 4) = [];

%% LETTURA dei TRIANGOLI della MESH
%Lettura linea commentata introduzione sezione degli elementi di bordo
line = fgets(mesh_file);

%Lettura NUMERO di TRIANGOLI
domainMesh.number_triangles = sscanf(fgets(mesh_file), '%d');

%Allocazione MATRICE contenente le INCIDENZE dei TRIANGOLI
domainMesh.triangles = zeros(domainMesh.number_triangles, 4);

%Ogni riga è relativa ad un triangolo: i primi tre elementi contengono 
% l'indice dei nodi che costituiscono i tre vertici del triangolo mentre
% l'ultimo elemento identifica il tipo di condizione al bordo assegnata 
% sulla porzione di bordo in cui si trova quel triangolo

%Formato dei dati
formatSpec = '%f';

%DIMENSIONE struttura in cui salvare i dati
sizeSpec = [4 4*domainMesh.number_triangles];

%Lettura incidenze dei triangoli e relativo tipo di dato al bordo 
domainMesh.triangles = fscanf(mesh_file, formatSpec, sizeSpec)';

%L'i-esima riga di questa matrice contiene nei primi 3 elementi le incidenze 
% dei tre vertici dell'i-esimo triangolo e nel quarto elemento l'indice
% relativo al tipo di dato al bordo da assegnare al triangolo

% Assegnazione forzata dell'indice di ciascun triangolo
ind = 4;
domainMesh.triangles(:, 4) = ind * ones(domainMesh.number_triangles, 1);

%% CHIUSURA DEL FILE MESH
%Lettura dell'ultima riga
line = fgets(mesh_file);

%CHIUSURA del FILE mesh
mesh_file = fclose(mesh_file);

%% ESTRAPOLAZIONE INFORMAZIONI UTILI sui TRIANGOLI della MESH
%Creazione matrice contenente le coordinate in sequenza di tutti i primi
% vertici, tutti i secondi e tutti i terzi
coord_vert_triangles = domainMesh.coordinates(domainMesh.triangles(:, 1:3), :);

%Calcolo dei vettori P1 - P3 = P3P1 al variare dei triangoli
vg1 = coord_vert_triangles(1 : domainMesh.number_triangles, :) - ...
    coord_vert_triangles(2*domainMesh.number_triangles+1 : end, :);

%Calcolo dei vettori P2 - P3 = P3P2 al variare dei triangoli
vg2 = coord_vert_triangles(domainMesh.number_triangles+1 : 2*domainMesh.number_triangles, :) - ...
    coord_vert_triangles(2*domainMesh.number_triangles+1 : end, :);

%Calcolo dei vettori P1 - P2 = P2P1 al variare dei triangoli
vg3 = coord_vert_triangles(1 : domainMesh.number_triangles, :) - ...
    coord_vert_triangles(domainMesh.number_triangles+1 : 2*domainMesh.number_triangles, :);

%Allocazione matrice contenente le componenti dei vettore normali ai triangoli della mesh
domainMesh.normal = zeros(domainMesh.number_triangles, 3);

%Allocazione vettore contenente le aree dei triangoli della mesh
domainMesh.area = zeros(domainMesh.number_triangles, 1);

%Calcolo COMPONENTI dei VETTORI NORMALI (n = vg1 x vg2)
domainMesh.normal(:,1) =   vg1(:,2) .* vg2(:,3) - vg2(:,2) .* vg1(:,3);
domainMesh.normal(:,2) = -(vg1(:,1) .* vg2(:,3) - vg2(:,1) .* vg1(:,3));
domainMesh.normal(:,3) =   vg1(:,1) .* vg2(:,2) - vg2(:,1) .* vg1(:,2);

%Calcolo AREE (A = |n|_2 / 2)
domainMesh.area = sqrt(sum(domainMesh.normal.^2, 2)) / 2;

%Normalizzazione VETTORI NORMALI ai triangoli della mesh
domainMesh.normal = domainMesh.normal ./ (2*domainMesh.area);

%Allocazione matrice contenente le componenti dei baricentri dei triangoli della mesh
domainMesh.center = zeros(domainMesh.number_triangles, 3);

%Calcolo componenti dei baricentri (b = P1 + P2 + P3 / 3)
domainMesh.center = (coord_vert_triangles(1 : domainMesh.number_triangles, :) + ...
    coord_vert_triangles(domainMesh.number_triangles+1 : 2*domainMesh.number_triangles, :) + ...
    coord_vert_triangles(2*domainMesh.number_triangles+1 : 3*domainMesh.number_triangles, :)) / 3;

%Calcolo lunghezza dei lati dei triangoli della mesh (l_i = |vg_i|_2)
len = [sqrt(sum(vg3.^2, 2)) sqrt(sum(vg2.^2, 2)) sqrt(sum(vg1.^2, 2))];

%Allocazione vettore contenente le massime lunghezze dei lati dei triangoli della mesh
domainMesh.maxL = zeros(domainMesh.number_triangles, 1);

%Calcolo massime lunghezze dei lati dei triangoli della mesh
domainMesh.maxL = max(len, [], 2);

%Allocazione matrice contenente le componenti dei curl dei
%triangoli della mesh
domainMesh.curl = zeros(3, 3, domainMesh.number_triangles);

%Calcolo delle componenti dei curl
domainMesh.curl(:, 1, :) =  (vg2 ./ (2*domainMesh.area))';
domainMesh.curl(:, 2, :) = -(vg1 ./ (2*domainMesh.area))';
domainMesh.curl(:, 3, :) =  (vg3 ./ (2*domainMesh.area))';


%Calcolo informazioni utili sulla mesh
domainMesh = calcParamMesh(domainMesh);
return