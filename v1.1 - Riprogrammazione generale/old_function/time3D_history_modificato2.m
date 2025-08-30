function sol = time3D_history_modificato2(pb_param,domainMesh,gha,ghw)

%La function time3D_history() prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- la struct domainMesh contenente le informazioni sulla mesh
%
%- l'array 3D gha contenente i nodi di quadratura di Gauss-Hammer
%
%- la matrice ghw che contiene i pesi della formula di quadratura di
%  Gauss-Hammer
%
%------------------------------------------------------------------------

%Calcoliamo il PASSO di DISCRETIZZAZIONE TEMPORALE (rapporto tra l'istante 
%di tempo finale e il numero di sottointervalli temporali)
dt=pb_param.T_fin/pb_param.Nt;

%Fissiamo l'INDICE hk=0 che individua l'indice del blocco E_0 della 
%della MATRICE di TOEPLITZ E che dobbiamo calcolare
%(è l'indice relativo al primo intervallo della discretizzazione temporale)
hk=0;

%Allochiamo la matrice sol in cui andremo a memorizzare la SOLUZIONE 
%del PROBLEMA.
sol = zeros(3*domainMesh.number_triangles,pb_param.Nt);
%Il numero di righe di questa matrice è pari a tre volte il numero
%di triangoli della discretizzazione mentre il numero di colonne è
%uguale al numero di istanti temporali della discretizzazione temporale.
%Quindi per ogni istante temporale, la corrispondente soluzione viene
%memorizzata in una diversa colonna della matrice soluzione sol


%Allochiamo la matrice matrix che contiene il BLOCCO MATRICIALE E0 
%(relativo all'INDICE hk=0) che dobbiamo calcolare
matrix = zeros(3*domainMesh.number_triangles,3*domainMesh.number_triangles);
%Il numero di righe e di colonne di questa matrice è pari a tre volte il 
%numero di triangoli della discretizzazione. 
%Possiamo pensare questa matrice come una matrice formata da 
%domainMesh.number_triangles^2 blocchetti 3x3


%Allochiamo l'array rhs che contiene il BLOCCO del TERMINE NOTO relativo
%all'INDICE hk=0.
%Esso ha tanti elementi quanto tre il volte il numero dei triangoli della
%discretizzazione
rhs = zeros(3*domainMesh.number_triangles,1);

%Inizializziamo a zero il parametro maxstep che indica l'indice 
%del SOTTOINTERVALLO della DISCRETIZZAZIONE TEMPORALE a partire dal 
%quale i successivi BLOCCHI MATRICIALI da calcolare sono NULLI.     
maxstep = 0;

%Calcoliamo i COEFFICIENTI della MATRICE E0 del TIME-MARCHING e del 
%corrispondente TERMINE NOTO
[matrix,rhs,~] = time3D_syst_modificato(matrix,rhs,hk,dt,pb_param,domainMesh,gha,ghw,maxstep);

%Nome del file in cui è stata memorizzata la prima matrice E0
file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(1),'.txt');
%Lettura della prima matrice E0
matrix = sparse(readmatrix(file_name));
file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(1),'.mat');
load(file_name);

% %file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(1));
% %save(file_name,'matrix');
% 

%FATTORIZZAZIONE PLU della matrice E0 del TIME-MARCHING 
[L,U,P] = lu(sparse(matrix));
 
%RISOLUZIONE del primo sistema del TIME-MARCHING. 
%Memorizziamo la soluzione nella prima colonna della matrice sol (che è
%quella relativa al primo istante di tempo della discretizzazione 
%temporale)
sol(:,1) = U\(L\(P*rhs));

%Salviamo la soluzione ottenuta al primo passo
soluzione=sol(:,1);
file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','soluzione','_',num2str(1),'.txt');
writematrix(soluzione,file_name,'Delimiter',';')

% load('screengraded_2_soluzione.mat','sol')
% hk=46;

for ind_step = 2:pb_param.Nt
    %Ciclo su hk=ind_step-1=1,...,N_(Delta t)-1 
    
    %Ciclo sull'INDICE hk+1=ind_step che varia tra 2 e il NUMERO di
    %SOTTOINTERVALLI TEMPORALI pb_param.Nt della discretizzazione 
    %temporale.
    %L'indice ind_step=1 (cioè hk=0) lo abbiamo già considerato in 
    %precedenza prima del ciclo for
    
    %Ad ogni iterazione, il blocco matriciale che dobbiamo calcolare
    %è il blocco di indice hk=ind_step-1=1,...,pb_param.Nt-1 che si trova
    %nella prima colonna della matrice di Toeplitz
    
    %RE-INIZIALIZZAZIONE a 0 del VETTORE TERMINE NOTO
    rhs = zeros(3*domainMesh.number_triangles,1);
    
    if(maxstep==0) 
        %Se la variabile maxstep è uguale a zero, allora significa che non
        %abbiamo ancora incontrato un indice hk in corrispondenza
        %del quale il relativo blocco matriciale è nullo. 
        %Quindi in tal caso possiamo procedere a calcolare un nuovo 
        %blocco matriciale. 
        %Viceversa non dovremo calcolare nessun nuovo blocco matriciale
        
        %AGGIORNIAMO il TERMINE NOTO con le matrici già esistenti che
        %abbiamo memorizzato nelle precedenti iterazioni, ciascuna in un 
        %diverso file di testo
        for ind_inner = 1:ind_step-2
            %Osserviamo che se siamo nell'iterazione in cui ind_step è
            %uguale a 2 e quindi dovremo calcolare la matrice di indice 
            %hk=1, non entriamo in questo ciclo for in quanto non dobbiamo
            %aggiornare il vettore dei termini noti
            
            %Nome del file in cui è memorizzata la matrice di indice 
            %ind_inner+1 (ovvero di indice hk=ind_inner)
            file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(ind_inner+1),'.txt');
            
            %Lettura della matrice corrente dal file file_name
            matrix = sparse(readmatrix(file_name));
            
            %Aggiornamento del termine noto
            %rhs = rhs-matrix*sol(:,ind_inner+1);
            rhs = rhs-matrix*sol(:,ind_step-ind_inner);
        end
        
        %Aggiorniamo l'INDICE TEMPORALE hk dell'ULTIMA MATRICE da  
        %calcolare per poi aggiornare il termine noto
        %(hk individua l'indice del sottointevrallo temporale che 
        %ad ogni interazione viene considerato)
        hk = hk+1;
        
        %Inizializzazione a zero della variabile matrix
        matrix = sparse(zeros(3*domainMesh.number_triangles,3*domainMesh.number_triangles));
        
        %Assemblaggio della matrice corrente e del termine noto corrente
        [matrix,rhs,iok] = time3D_syst_modificato(matrix,rhs,hk,dt,pb_param,domainMesh,gha,ghw,maxstep);
        
        if (ind_step<=10)
            
            %Lettura della matrice E_(ind_step)
            file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(ind_step),'.txt');
            matrix = sparse(readmatrix(file_name));
            
            file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(ind_step),'.mat');
            load(file_name);
        else
            iok=0;
        end
        
        if(iok==0) 
            %Se iok è uguale a zero allora questo significa che l'ultimo
            %blocco matriciale matrix di indice hk=ind_step-1 che è
            %stato calcolato è una matrice nulla. Quindi non è necessario 
            %aggiornare il termine noto e non salviamo nemmeno l'ultima
            %matrice calcolata
                        
            maxstep = ind_step;
            %Assegniamo alla variabile maxstep il valore dell'indice 
            %ind_step in corrispondenza del quale il relativo blocco 
            %matriciale di indice hk=ind_step-1 è uguale a zero
            
        else
            %Se iok è diverso da zero allora questo significa che l'ultimo
            %blocco matriciale matrix di indice hk=ind_step-1 che è
            %appena stato calcolato non è nullo. Allora dobbiamo 
            %aggiornare il termine noto
            
            %AGGIORNAMENTO del TERMINE NOTO
            %rhs = rhs-matrix*sol(:,ind_step-1);
            rhs = rhs-matrix*sol(:,1);
            
%             %Nome del file in cui verrà memorizzata la nuova matrice 
%             %di indice hk=ind_step-1 calcolato a questo passo
%             file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(ind_step),'.txt');
%             
%             %Scrittura sul file file_name dell'ultima matrice corrente
%             writematrix(matrix,file_name,'Delimiter',';')
            
            %file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(ind_step));
            %save(file_name,'matrix');
            
        end %Fine if(iok==0) 
       
    else
        %Se la variabile maxstep è diversa da zero, allora questo
        %significa che contiene l'indice hk=maxstep-1 in corrispondenza 
        %del quale la relativa matrice del time-marching risulta uguale a 
        %zero. Pertanto non è necessario calcolare nuove matrici e 
        %dobbiamo solamente aggiornare il termine noto
        
        %Aggiornamento del termine noto con le matrici già esistenti
        for ind_inner = 2:maxstep-1
            
            %Nome del file in cui è memorizzata la matrice corrente
            file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','matrix','_',num2str(ind_inner),'.txt');
            
            %Lettura della matrice corrente (di indice ind_inner)
            matrix = sparse(readmatrix(file_name));
            
            %Aggiornamento del termine noto
            rhs = rhs-matrix*sol(:,ind_step-ind_inner+1);
            
        end %Fine ciclo for su ind_inner = 2:maxstep-1
        
        %Aggiorniamo l'INDICE TEMPORALE del TERMINE NOTO da 
        %calcolare
        hk = hk+1;
        
        %Assemblaggio del termine noto corrente 
        %N.B. nella subroutine time3D_syst() la matrice matrix di indice 
        %hk=ind_step-1 non viene calcolata in quanto siamo nel caso in cui 
        %maxstep è diverso da zero
        [matrix,rhs,~] = time3D_syst_modificato(matrix,rhs,hk,dt,pb_param,domainMesh,gha,ghw,maxstep);
        
    end %Fine if(maxstep==0) 
 
    %Risolviamo il SISTEMA e memorizziamo la SOLUZIONE (relativa 
    %all'istante temporale di indice ind_step della discretizzazione 
    %temporale) nella colonna di indice ind_step della matrice sol
    sol(:,ind_step) = U\(L\(P*rhs));
    %Nome del file in cui verrà memorizzata la soluzione
    file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','soluzione');
    %Scrittura sul file file_name della soluzione
    save(file_name,'sol');
    
    %Salviamo la soluzione di indice ind_step
    soluzione=sol(:,ind_step);
    file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','soluzione','_',num2str(ind_step),'.txt');
    writematrix(soluzione,file_name,'Delimiter',';')

    
end %Fine ciclo for su ind_step = 2:pb_param.Nt

%Nome del file in cui verrà memorizzata la soluzione
file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','soluzione','.txt');
%Scrittura sul file file_name della soluzione
writematrix(sol,file_name,'Delimiter',';')
% file_name = strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_','soluzione');
% save(file_name,'sol');


return