function [matrix] = time3D_syst_post_processing(pb_param,domainMesh,matrix,thk,x)

%La function time3D_syst_post_processing() prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- la struct domainMesh contenente le informazioni sulla mesh
%
%- la matrice matrix che ha 3 righe e tante colonne quanto tre volte il
%  numero di triangoli della mesh. Essa è quindi costituita da tanti 
%  blocchetti 3x3 quanti sono i triangoli della mesh
%
%- la variabile thk (=Delta) contenente la differenza tra l'istante t e un
%  istante della discretizzazione temporale
%
%- il vettore 1x3 x che contiene le coordinate del punto che assume il
%  ruolo di punto sorgente
%
%e restituisce in OUTPUT:
%- la matrice matrix che contiene il risultato dell'integrazione 
%  su tutti i triangoli della mesh ottenuta considerando come punto
%  sorgente il punto x
%
%-------------------------------------------------------------------------

%*************************************************************************

%PARALLELIZZAZIONE DEL DOPPIO CICLO FOR SUI TRIANGOLI SORGENTE 

%Fissiamo il NUMERO di TRIANGOLI della MESH
N_triangles = domainMesh.number_triangles;

%Trasformiamo la matrice matrix in una matrice di celle costituita da
%1 x N_triangles celle ciascuna delle quali contiene una 
%matrice di formato 3x3.
%Al termine di questa function, la cella di indice j conterrà  la 
%matrice 3x3 i cui coefficienti sono il risultato dell'integrazione 
%sul j-esimo triangolo di campo
matrix = mat2cell(matrix,3,3*ones(1,N_triangles));

%Ciclo sul numero di triangoli 
parfor ind=1:N_triangles
%for ind = 1:N_triangles

    indF=ind;   
    %---------------------------------------------------------------------
    %Estrapoliamo le INFORMAZIONI utili relative al TRIANGOLO DI CAMPO 
    %di INDICE indF
    
    %INDICI dei VERTICI dell'elemento corrente (triangolo di campo)
    nodeF = domainMesh.triangles(indF,1:3);
    %nodeF è un array 1x3 che contiene gli indici dei tre vertici del
    %triangolo di campo corrente
    
    %COORDINATE dei VERTICI dell'elemento corrente (triangolo sorgente)
    TF = domainMesh.coordinates(nodeF,:);
    %TF è una matrice 3x3 che contiene le coordinate dei tre vertici del
    %triangolo di campo corrente
    
    %VERSORE NORMALE all'elemento corrente (triangolo di campo)
    vnF = domainMesh.normal(indF,:);
        
    %Calcolo dell'INTEGRALE sul TRIANGOLO SORGENTE
    ris = time3D_dswLE_post_processing(pb_param,TF,vnF,thk,x)

    %Buttiamo via la "SPORCIZIA NUMERICA" dal blocchetto 
    %matriciale locale
    ris(abs(ris)<1.0e-14) = 0;

    %Posizioniamo il blocco matriciale locale nella matrice
    %globale matrix
    matrix{ind} = matrix{ind}+ ris;
  
end

%Trasformiamo la matrice matrix costituita da 1 x N_triangles 
%celle nella corrispondente matrice di formato 3 x 3N_triangles 
%utilizzando la function cell2matrix()
matrix=cell2mat(matrix);

return