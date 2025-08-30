function domainMesh=time3D_read_space_mesh(domain_type,lev)
%La function time3D_read_space_mesh() prende in INPUT:
%- la variabile domain_type, contenente il TIPO di DOMINIO
%- la variabile lev, contenente il LIVELLO di RAFFINAMENTO
%e restituisce in OUTPUT la struct domainMesh contenente i cui campi
%descrivono la mesh
    
%NOME del FILE MESH
file_mesh = strcat('./old_function/mesh','/',domain_type,'/',domain_type,'_',num2str(lev),'.mesh');
%'./mesh/cantilever_BEM/cantilever_BEM_0.mesh

%Apertura del file mesh
[mesh_file,~] = fopen(file_mesh,'r');
%Apriamo il file dal nome file_mesh in sola lettura attraverso la
%funzione fopen(), che restituisce in OUTPUT nella variabile 
%mesh_file un numero intero maggiore di 3 che è l'identificativo del
%file aperto

%Lettura dell'intestazione del file mesh
line = fgets(mesh_file);
line = fgets(mesh_file);
line = fgets(mesh_file);
line = fgets(mesh_file);

%A questo punto il programma si trova a dover leggere la quinta riga
%del file

%% LETTURA dei NODI della MESH

%Lettura della linea commentata che introduce la sezione dedicata ai nodi
line = fgets(mesh_file);

%Lettura del NUMERO dei NODI
domainMesh.number_nodes = sscanf(fgets(mesh_file),'%d');

%Allocazione della MATRICE contenente le COORDINATE dei NODI
domainMesh.coordinates = zeros(domainMesh.number_nodes,4);

%Questa matrice ha tante righe quanti sono i nodi e ha 4 colonne.
%Ogni riga della matrice è relativa ad un determinato nodo: i primi 
%tre elementi della riga contengono le coordinate del nodo, mentre
%l'ultimo elemento indica a quale faccia appartiene il nodo

%FORMATO in cui salvare le COORDINATE dei NODI
formatSpec = '%f';

%Struttura in cui salvare i dati
sizeSpec = [4 4*domainMesh.number_nodes];

%Lettura dei nodi della mesh
domainMesh.coordinates = fscanf(mesh_file,formatSpec,sizeSpec)';

%Eliminazione della quarta colonna inutile
domainMesh.coordinates(:,4) = [];

%domainMesh.coordinates(abs(domainMesh.coordinates)<=1.0e-10)=0;
%Shiftiamo in avanti le coordinate dei triangoli per cambiare il piano in
%cui giace la screen quadrata
% for i=1:domainMesh.number_nodes
%     domainMesh.coordinates(i,1:3)=circshift(domainMesh.coordinates(i,1:3),2);
% end

%% LETTURA dei TRIANGOLI della MESH

%Lettura della linea commentata che introduce la sezione dedicata agli
%elementi di bordo
line = fgets(mesh_file);

%Lettura del NUMERO di TRIANGOLI
domainMesh.number_triangles  = sscanf(fgets(mesh_file),'%d');

%Allocazione della MATRICE contenente le INCIDENZE dei TRIANGOLI
domainMesh.triangles = zeros(domainMesh.number_triangles,4);

%domainMesh.triangles è una matrice che ha tante righe quanti sono i
%triangoli della mesh e ha 4 colonne. Ogni riga è relativa ad un
%triangolo, i primi tre elementi contengono l'indice dei nodi che
%costituiscono i tre vertici del triangolo mentre l'ultimo elemento 
%è un indice che identifica il tipo di condizione al bordo assegnata 
%sulla porzione di bordo in cui si trova quel triangolo


%Formato dei dati
formatSpec = '%f';

%Struttura in cui salvare i dati
sizeSpec = [4 4*domainMesh.number_triangles];

%Lettura delle incidenze dei vertici del triangolo corrente e del tipo di
%dato al bordo assegnato a ciascun triangolo
domainMesh.triangles = fscanf(mesh_file,formatSpec,sizeSpec)';

%domainMesh.triangles è una matrice che ha tante righe quanti sono i
%triangoli della discretizzazione e ha 4 colonne. Sull'i-esima riga di
%questa matrice, i primi 3 elementi contenengono le incidenze dei tre
%vertici del triangolo i-esimo mentre il quarto elemento contiene l'indice
%relativo al tipo di dato al bordo da assegnare al triangolo

%Forziamo l'indice da assegnare a ciascun triangolo
ind=4;
domainMesh.triangles(:,4)=ind*ones(domainMesh.number_triangles,1);

%Lettura dell'ultima riga del file mesh
line = fgets(mesh_file);

%CHIUSURA del FILE mesh
mesh_file = fclose(mesh_file);

%% ESTRAPOLO INFORMAZIONI UTILI sui TRIANGOLI della MESH

coord_vert_triangles = domainMesh.coordinates(domainMesh.triangles(:,1:3),:);

%Calcolo il vettore P1-P3=P3P1 (Pi vertice i-esimo del triangolo)
vg1 = coord_vert_triangles(1:domainMesh.number_triangles,:)-...
  coord_vert_triangles(2*domainMesh.number_triangles+1:end,:);
%Calcolo il vettore P2-P3=P3P2 (Pi vertice i-esimo del triangolo)
vg2 = coord_vert_triangles(domainMesh.number_triangles+1:2*domainMesh.number_triangles,:)-...
  coord_vert_triangles(2*domainMesh.number_triangles+1:end,:);
%Calcolo il vettore P1-P2 (Pi vertice i-esimo del triangolo)  
vg3 = coord_vert_triangles(1:domainMesh.number_triangles,:)-...
  coord_vert_triangles(domainMesh.number_triangles+1:2*domainMesh.number_triangles,:);

%Allocazione della matrice contenente le componenti del vettore normale ai
%triangoli della mesh
domainMesh.normal = zeros(domainMesh.number_triangles,3);

%Allocazione del vettore contenente le aree dei triangoli della mesh
domainMesh.area = zeros(domainMesh.number_triangles,1);

%Calcolo delle COMPONENTI del VETTORE NORMALE 
%dato da vg1 x vg2 = (P1-P3) x (P2-P3) = P3P1 x P3P2
domainMesh.normal(:,1) =   vg1(:,2).*vg2(:,3)-vg2(:,2).*vg1(:,3);
domainMesh.normal(:,2) = -(vg1(:,1).*vg2(:,3)-vg2(:,1).*vg1(:,3));
domainMesh.normal(:,3) =   vg1(:,1).*vg2(:,2)-vg2(:,1).*vg1(:,2);

%Calcolo delle AREE
domainMesh.area = sqrt(sum(domainMesh.normal.^2,2))/2;

%Normalizzazione dei vettori normale ai triangoli della mesh
domainMesh.normal = domainMesh.normal./(2*domainMesh.area);

%Allocazione della matrice contenente le componenti dei baricentri dei
%triangoli della mesh
domainMesh.center = zeros(domainMesh.number_triangles,3);
%Calcolo delle componenti dei baricentri
domainMesh.center = (coord_vert_triangles(1:domainMesh.number_triangles,:)+...
coord_vert_triangles(domainMesh.number_triangles+1:2*domainMesh.number_triangles,:)+...
coord_vert_triangles(2*domainMesh.number_triangles+1:3*domainMesh.number_triangles,:))/3;

%Calcolo la lunghezza dei lati dei triangoli della mesh
len = [sqrt(sum(vg3.^2,2)) sqrt(sum(vg2.^2,2)) sqrt(sum(vg1.^2,2))];

%Allocazione del vettore contenente le massime lunghezze dei lati
%dei triangoli della mesh
domainMesh.maxL = zeros(domainMesh.number_triangles,1);

%Massima lunghezza dei lati dei triangoli della mesh
domainMesh.maxL = max(len,[],2);

%Allocazione della matrice contenente le componenti dei curl dei
%triangoli della mesh
domainMesh.curl = zeros(3,3,domainMesh.number_triangles);

%Calcolo delle componenti dei curl
domainMesh.curl(:,1,:) =  (vg2./(2*domainMesh.area))';
domainMesh.curl(:,2,:) = -(vg1./(2*domainMesh.area))';
domainMesh.curl(:,3,:) =  (vg3./(2*domainMesh.area))';

return

