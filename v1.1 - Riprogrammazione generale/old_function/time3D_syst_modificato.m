%function [matrix,rhs,iok,mat_screen0] = time3D_syst(matrix,rhs,hk,dt,pb_param,domainMesh,gha,ghw,maxstep)
function [matrix,rhs,iok] = time3D_syst_modificato(matrix,rhs,hk,dt,pb_param,domainMesh,gha,ghw,maxstep)

%La function time3D_syst() prende in INPUT:
%- la matrice matrix in cui memorizzare la matrice del TIME-MARCHING
%
%- il vettore rhs in cui memorizzare il termine noto
%
%- la variabile hk che contiene l'indice del blocco matriciale presente 
%  sulla prima colonna della matrice di Toeplitz che dobbiamo calcolare
%
%- la struct pb_param contenente i parametri del problema
%
%- la struct domainMesh contenente le informazioni sulla mesh
%
%- l'array 3D gha contenente i nodi di quadratura di Gauss-Hammer
%
%- la matrice ghw che contiene i pesi della formula di quadratura di
%  Gauss-Hammer
%
%- la variabile maxstep che contiene l'indice del sottointervallo
%  della discretizzazione temporale a partire dal quale le matrici sono
%  costituite da blocchi matriciali tutti nulli. Se maxstep è uguale a
%  zero, allora significa che non abbiamo ancora incontrato un indice  
%  in corrispondenza del quale il relativo blocco matriciale è nullo.
%
%e restituisce in OUTPUT:
%- la matrice matrix contenente il blocco matriciale di indice hk che ha 
%  tante righe e tante colonne quanto tre volte il numero di triangoli
%  della discretizzazione  
%
%- il vettore rhs che contiene il termine noto corrispondente alla matrice
%  matrix e che ha tanti elementi quanto tre volte il numero dei triangoli
%  della discretizzazione  
%
%- il parametro iok che vale 1 se almeno uno dei blocchi matriciali 3x3
%  relativi all'integrazione su un triangolo sorgente e un triangolo di
%  campo è non nullo, mentre vale zero se tutti i blocchetti sono nulli (e
%  quindi la matrice matrix calcolata è una matrice nulla)
%
%-------------------------------------------------------------------------

%*************************************************************************

%PARALLELIZZAZIONE DEL DOPPIO CICLO FOR SUI TRIANGOLI SORGENTE E DI CAMPO

%Fissiamo il NUMERO di TRIANGOLI della MESH
N_triangles = domainMesh.number_triangles;

%Trasformiamo la matrice matrix in una matrice di celle costituita da
%N_triangles x N_triangles celle ciascuna delle quali contiene una 
%matrice di formato 3x3.
%Al termine di questa function, la cella di indice (i,j) conterrÃ  la 
%matrice 3x3 i cui coefficienti sono il risultato dell'integrazione 
%sull'i-esimo triangolo sorgente e sul j-esimo triangolo di campo
matrix = mat2cell(matrix,3*ones(1,N_triangles),3*ones(1,N_triangles));

%Memorizziamo nel vettore 1x2 sz le dimensioni della matrice di celle
%matrix
sz = size(matrix);

%Trasformiamo il vettore rhs in un vettore costituito da N_triangles x 1
%celle ciascuna delle quali contiene un vettore di formato 3 x 1.
%Al termine di questa function, la cella di indice i conterrÃ  il 
%vettore 3 x 1 i cui coefficienti sono relativi all'integrazione 
%sull'i-esimo triangolo sorgente 
rhs = mat2cell(rhs,3*ones(1,N_triangles),1);
beta = mat2cell(zeros(3*N_triangles,1),3*ones(1,N_triangles),1);
%Aggiungiamo al vettore colonna rhs altre colonne costituite da celle 
%vuote in modo tale che rhs risulti una matrice costituita da 
%N_triangles x N_triangles celle
%N.B. Se non si esegue quest'ultima operazione, si hanno dei problemi con
%l'if presente all'interno del parfor!!!
rhs = horzcat(rhs,cell(N_triangles,N_triangles-1));

beta = horzcat(beta,cell(N_triangles,N_triangles-1));
%mat_screen0=cell(1,N_triangles);

%Allocazione del cell-array contenente tutti i valori di iok
cell_iok = cell(N_triangles^2,1);

%Inizializzazione del cell-array contenente tutti i valori di iok
cell_iok(1:N_triangles^2,1) = {0};

%TECNICA ELEMENT-by-ELEMENT: ciclo sul numero di triangoli al quadrato
%in cui, ad ogni iterazione viene fissato un triangolo sorgente e 
%un triangolo di campo


parfor ind=1:N_triangles^2
%for ind = N_triangles^2
     
    [indS, indF] = ind2sub(sz,ind);
    %La function ind2sub() prende in input il vettore 1x2 sz contenente le
    %dimensioni della matrice e l'indice linerare ind e restituisce in
    %output gli indici di riga indS e di colonna indF a cui corrisponde 
    %l'indice lineare ind
    
    %---------------------------------------------------------------------
    %Estrapoliamo le INFORMAZIONI utili relative al TRIANGOLO SORGENTE 
    %di INDICE indS
    
    %INDICI dei VERTICI dell'elemento corrente (triangolo sorgente)
    nodeS = domainMesh.triangles(indS,1:3);
    %nodeS è un array 1x3 che contiene gli indici dei tre vertici del
    %triangolo sorgente corrente
    
    %COORDINATE dei VERTICI dell'elemento corrente (triangolo sorgente)
    TS = domainMesh.coordinates(nodeS,:);
    %TS è una matrice 3x3 che contiene le coordinate dei tre vertici del
    %triangolo sorgente corrente
    
    %COORDINATE del BARICENTRO dell'elemento corrente (triangolo sorgente)
    cS = domainMesh.center(indS,:);
    
    %AREA dell'elemento corrente (triangolo sorgente)
    areaS = domainMesh.area(indS);
    
    %MASSIMA LUNGHEZZA dei LATI dell'elemento corrente (triangolo sorgente)
    maxS = domainMesh.maxL(indS);
    
    %VERSORE NORMALE all'elemento corrente (triangolo sorgente)
    vnS = domainMesh.normal(indS,:);
    
    %ROTORE dell'elemento corrente (triangolo sorgente)
    curlS = domainMesh.curl(:,:,indS); 
    
    %INDICE che identifica il TIPO di DATO al BORDO relativo 
    %all'ELEMENTO CORRENTE (triangolo sorgente)
    indS_RHS = domainMesh.triangles(indS,4);
    %---------------------------------------------------------------------

    if(maxstep==0) 
        %Se maxstep è UGUALE a 0 allora sappiamo che dobbiamo PROCEDERE 
        %a CALCOLARE il BLOCCO MATRICIALE relativo all'indice hk.
        %IN CASO CONTRARIO sappiamo che maxstep contiene l'indice
        %dell'istante di tempo della discretizzazione temporale
        %relativamente al quale il corrispondente blocco matriciale è
        %nullo. Pertanto NON DOBBIAMO CALCOLARE un nuovo blocco 
        %matriciale del time-marching.
        
        %Estrapoliamo le INFORMAZIONI utili relative al TRIANGOLO di CAMPO
        %di INDICE indF
        
        %COORDINATE del BARICENTRO dell'elemento corrente (triangolo di campo)
        cF = domainMesh.center(indF,:);
        
        %MASSIMA LUNGHEZZA dei LATI dell'elemento corrente (triangolo di campo)
        maxF = domainMesh.maxL(indF);
        
        %VETTORE DISTANZA tra i BARICENTRI del triangolo di campo e del
        %triangolo sorgente
        cF = cF-cS;
        
        %DISTANZA MINIMA tra i BARICENTRI del triangolo di campo e del
        %triangolo sorgente
        distmin = sqrt(sum(cF.^2))-maxF-maxS;
        
        %DISTANZA MASSIMA tra i BARICENTRI del triangolo di campo e del
        %triangolo sorgente
        distmax = sqrt(sum(cF.^2))+maxF+maxS;
                
        %Controlliamo se l'elemento matriciale 3x3 ris corrispondente al 
        %contenente il risultato dell'integrazione sul triangolo
        %sorgente corrente e sul triangolo di campo corrente è non nullo.
        %In caso contrario non dobbiamo eseguire nessuna integrazione
        
        %if(((hk+1)*pb_param.velS*dt>distmin) && ((hk-2)*pb_param.velP*dt<distmax))
        
        if(((hk-1)*pb_param.velS*dt<distmax) && ((hk+1)*pb_param.velP*dt>=distmin))
            
            %Aggiornamento del parametro iok locale
            cell_iok(ind,1) = {1};
            %Poniamo la cella del vettore iok di indice ind  uguale a 1
            %poiché il blocco matriciale 3x3 ris che dobbiamo calcolare è
            %non nullo
            
            %Indici dei vertici dell'elemento corrente (triangolo di campo)
            nodeF = domainMesh.triangles(indF,1:3);
            %Coordinate dei vertici dell'elemento corrente (triangolo di campo)
            TF = domainMesh.coordinates(nodeF,:);
            
            %             fill3(TS(1,:),TS(2,:),TS(3,:),'r')
            %             hold on
            %             fill3(TF(1,:),TF(2,:),TF(3,:),'b')
            %             hold off
            %             pause
            
            %%AREA dell'elemento corrente (triangolo di campo)
            %areaF = domainMesh.area(indF);
            
            %VERSORE NORMALE all'elemento corrente (triangolo di campo)
            vnF = domainMesh.normal(indF,:);
            
            %%ROTORE dell'elemento corrente (triangolo di campo)
            %curlF = domainMesh.curl(:,:,indF);
            %---------------------------------------------------------------------
            
            %Calcolo dell'INTEGRALE DOPPIO sul TRIANGOLO SORGENTE
            %(int. esterna) e sul TRIANGOLO di CAMPO (int. interna)
            ris=zeros(3,3);
            %ris = time3D_doubleLE(pb_param,TS,areaS,TF,vnF,hk,dt,gha,ghw);
                    
            %         if (indS==1 && indF==1)
            %             r1=ris;
            %         end
            %         if (indS==12 && indF==47)
            %             r2=ris;
            %         end
            
            %Buttiamo via la "SPORCIZIA NUMERICA" dal blocchetto 
            %matriciale locale
            ris(abs(ris)<1.0e-14) = 0;
            
            %Posizioniamo il blocco matriciale locale nella matrice
            %globale matrix
            matrix{ind} = matrix{ind}+ ris;

         end %Fine if(((hk-1)*pb_param.velS*dt<distmax) && ((hk+1)*pb_param.velP*dt>=distmin))
        
    end %Fine if(maxstep==0)
    
    %Contributo al termine noto dell'integrazione sul triangolo sorgente
    if (ind<=N_triangles)

        %Calcolo del VETTORE TERMINE NOTO LOCALE
%        [ris_rhs, riga_screen0] = time3D_singleLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,hk,dt,gha,ghw,indS);
        [ris_rhs] = time3D_singleLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,hk,dt,gha,ghw,indS);

        %Buttiamo via la "SPORCIZIA NUMERICA" dal vettore temine noto 
        %locale
        ris_rhs(abs(ris_rhs)<1.0e-14) = 0;
%         riga_screen0(abs(riga_screen0)<1.0e-14)=0;

        %Posizionamento nel vettore termine noto globale
        rhs{ind} = rhs{ind}+ris_rhs;
%         mat_screen0=vertcat(mat_screen0,mat2cell(riga_screen0,3,3*ones(1,N_triangles)));
        
        beta{ind} = beta{ind}+ris_rhs;
    end
    
end

%Trasformiamo la matrice matrix costituita da N_triangles x N_triangles 
%celle nella corrispondente matrice di formato 3N_triangles x 3N_triangles 
%utilizzando la function cell2matrix()
matrix=cell2mat(matrix);

%Ricaviamo il vettore rhs costituito da 3N_triangles x 1 elementi 
%utilizzando la function cell2matrix() a partire dalla matrice di celle 
%rhs
rhs=cell2mat(rhs);

%Salvataggio del termine noto "non aggiornato"
beta=cell2mat(beta);
file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','beta','_',num2str(hk+1),'.txt');
writematrix(beta,file_name,'Delimiter',';')

%Calcolo del PARAMETRO iok
iok = max(cell2mat(cell_iok));
%iok risulta uguale a 1 se almeno uno dei blocchetti matriciali 3x3
%relativi all'integrazione su un triangolo sorgente e un triangolo di
%campo è diverso da zero. In caso contrario esso è uguale a 0 (quando
%questo accade significa che la matrice matrix è la matrice nulla)

iok=1;
%**********************************************************************

% %TECNICA ELEMENT-by-ELEMENT: ciclo sugli elementi per fissare il 
% %triangolo sorgente
% for indS=1:domainMesh.number_triangles 
%     
%     %Estrapoliamo le informazioni utili relative al triangolo sorgente 
%     %di indice indS
%     
%     %INDICI dei VERTICI dell'elemento corrente (triangolo sorgente)
%     nodeS = domainMesh.triangles(indS,1:3);
%     %nodeS Ã¨ un array 1x3 contenente gli indici dei nodi dei suoi vertici
%     
%     %COORDINATE dei VERTICI dell'elemento corrente (triangolo sorgente)
%     TS = domainMesh.coordinates(nodeS,:);
%     %TS Ã¨ una matrice 3x3 contenente le coordinate dei tre vertici 
%     
%     %COORDINATE del BARICENTRO dell'elemento corrente (triangolo sorgente)
%     cS = domainMesh.center(indS,:);
%     
%     %AREA dell'elemento corrente (triangolo sorgente)
%     areaS = domainMesh.area(indS);
%     
%     %MASSIMA LUNGHEZZA dei LATI dell'elemento corrente (triangolo sorgente)
%     maxS = domainMesh.maxL(indS);
%     
%     %VERSORE NORMALE all'elemento corrente (triangolo sorgente)
%     vnS = domainMesh.normal(indS,:);
%     
%     %ROTORE dell'elemento corrente (triangolo sorgente)
%     curlS = domainMesh.curl(:,:,indS);
%         
%     %TECNICA ELEMENT-by-ELEMENT: ciclo sugli elementi per fissare il 
%     %triangolo di campo
%     for indF = 1:domainMesh.number_triangles
%         
%         %Estrapoliamo le informazioni utili relative al triangolo di campo 
%         %di indice indF
%         
%         %COORDINATE del BARICENTRO dell'elemento corrente (triangolo di campo)
%         cF = domainMesh.center(indF,:);
%         
%         %MASSIMA LUNGHEZZA dei LATI dell'elemento corrente (triangolo di campo)
%         maxF = domainMesh.maxL(indF);
%         
%         %VETTORE DISTANZA tra i BARICENTRI del triangolo di campo e del 
%         %triangolo sorgente
%         cF = cF-cS;
%         
%         %DISTANZA MINIMA tra i BARICENTRI del triangolo di campo e del 
%         %triangolo sorgente
%         distmin = sqrt(sum(cF.^2))-maxF-maxS;
%         
%         %DISTANZA MASSIMA tra i BARICENTRI del triangolo di campo e del 
%         %triangolo sorgente
%         distmax = sqrt(sum(cF.^2))+maxF+maxS;
%                 
%         %if(((hk+1)*c*dt>distmin).and.((hk-2)*c*dt<distmax))
%             
%             iok = 1;
%             
%             %Indici dei vertici dell'elemento corrente (triangolo di campo)
%             nodeF = domainMesh.triangles(indF,1:3);
%             %Coordinate dei vertici dell'elemento corrente (triangolo di campo)
%             TF = domainMesh.coordinates(nodeF,:);
%              
% %             fill3(TS(1,:),TS(2,:),TS(3,:),'r')
% %             hold on
% %             fill3(TF(1,:),TF(2,:),TF(3,:),'b')
% %             hold off
% %             pause
%             
%             %AREA dell'elemento corrente (triangolo di campo)
%             areaF = domainMesh.area(indF);
%             
%             %VERSORE NORMALE all'elemento corrente (triangolo di campo)
%             vnF = domainMesh.normal(indF,:);
%             
%             %ROTORE dell'elemento corrente (triangolo di campo)
%             curlF = domainMesh.curl(:,:,indF);    
%             
%             %Calcolo dell'integrale doppio sul triangolo sorgente (int. esterna)
%             %e sul triangolo di campo (int. interna)
%             ris = time3D_doubleLE(pb_param,TS,areaS,TF,vnF,hk,dt,gha,ghw);
% %             ris = time3D_doubleLE(pb_param,TS,vnS,areaS,curlS,TF,vnF,...
% %                 areaF,curlF,hk,dt,gha,ghw);
% 
% %             if (indS==1 && indF==1)
% %                 r1=ris;
% %             end
% %             if (indS==12 && indF==47)
% %                 r2=ris;
% %             end
%             
%             %Posizionamento nella matrice globale
%             matrix(3*(indS-1)+1:3*indS,3*(indF-1)+1:3*indF) = matrix(3*(indS-1)+1:3*indS,3*(indF-1)+1:3*indF)+ris;
%                       
%         %end
%       
%     end %fine ciclo sui triangoli di campo (indF)
%     
%     %Contributo al termine noto dell'integrazione sul triangolo sorgente
%     ris_rhs = time3D_singleLE_rhs(pb_param,TS,vnS,areaS,curlS,hk,dt,gha,ghw);
% 
%      %Buttiamo via la "SPORCIZIA NUMERICA" dal vettore temine noto 
%      %locale
%      ris_rhs(abs(ris_rhs)<1.0e-14) = 0;
%     
%     %Posizionamento nel vettore termine noto globale
%     rhs(3*(indS-1)+1:3*indS,1) = rhs(3*(indS-1)+1:3*indS,1)+ris_rhs;
%     
% end %fine ciclo sui triangoli sorgente (indS)

% diff=abs(matrix-matrix');
% diff(diff<1.0e-10)=0;
% max(max(diff))
    
return