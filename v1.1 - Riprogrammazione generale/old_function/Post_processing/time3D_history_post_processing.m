function U = time3D_history_post_processing(pb_param,domainMesh,density,x,t)

%La function time3D_history_post_processing() prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- la struct domainMesh contenente le informazioni sulla mesh
%
%- la matrice density che ha tante righe quanto tre volte il numero dei
%  triangoli della mesh e tante colonne quanti sono i sottotintervalli della
%  discretizzazione temporale
%
%- il vettore 3x1 x contenente le coordinate del punto in cui vogliamo
%  calcolare il valore della soluzione
%
%- la variabile t contenente l'istante temporale in cui vogliamo
%  calcolare il valore della soluzione
%------------------------------------------------------------------------

%Allochiamo la matrice 3x1 U in cui l'i-esima riga contiene il valore
%dell'i-esima componente u_i della soluzione u valutata nel punto x
%all'istante di tempo t
U=zeros(3,1);

if t>0
    %Se l'istante di tempo t in cui vogliamo calcolare la soluzione è
    %uguale a zero, allora sappiamo già che in quel caso la soluzione è 
    %nulla e quindi non dobbiamo calcolarla
    
    %Calcoliamo il PASSO di DISCRETIZZAZIONE TEMPORALE (rapporto tra 
    %l'istante di tempo finale e il numero di sottointervalli temporali)
    dt=pb_param.T_fin/pb_param.Nt;

    %Calcoliamo il massimo indice n_hat degli intervalli temporali in
    %corrispondenza dei quali dobbiamo calcolare la soluzione
    n_hat=ceil(t/dt)-1;
    
    %Fissiamo il NUMERO di TRIANGOLI della DISCRETIZZAZIONE SPAZIALE
    N_triangles = domainMesh.number_triangles;
    
    %Definiamo la prima differenza temporale tra gli istanti t e t0=0
    thk0=t; %thk0=t-t0=t 
        
    %Allochiamo la matrice matrix0 di formato 3 x 3*N_triangles che è
    %costituita da tanti blocchetti 3x3 quanti sono i triangoli della 
    %discretizzaizone spaziale  
    matrix0=zeros(3,3*N_triangles);

    %Calcoliamo gli elementi della matrice matrix0
    [matrix0] = time3D_syst_post_processing(pb_param,domainMesh,matrix0,thk0,x);

    for hk=0:n_hat-1
        %Ciclo sull'indice n che varia da 0 a n_hat-1 e che individua i
        %vari sottointervalli temporali della discretizzazione temporale        
                  
            %Calcoliamo la differenza temporale tra gli istanti di
            %tempo t e (hk+1)*dt
            thk1=t-(hk+1)*dt;

            %Memorizziamo nel vettore alpha di formato 3*N_triangles x 1 la 
            %densità ottenuta in corrispondenza del sottointervallo temporale  
            %di indice hk
            alpha=density(:,hk+1);

            %Allochiamo la matrice matrix0 di formato 3 x 3*N_triangles che è
            %costituita da tanti blocchetti 3x3 quanti sono i triangoli della 
            %discretizzaizone spaziale  
            matrix1=zeros(3,3*N_triangles);

            %Calcoliamo gli elementi della matrice matrix1
            [matrix1] = time3D_syst_post_processing(pb_param,domainMesh,matrix1,thk1,x);
             
            %Aggiorniamo la soluzione U
            U=U+(matrix0-matrix1)*alpha;     
    %         for j=1:N_triangles
    %             U=U+(matrix0{j}-matrix1{j})*alpha(3*(j-1)+1:3*j,1);
    %         end

            %Assegnamo alla matrice matrix0 la matrice matrix1 appena 
            %calcolata
            matrix0=matrix1;
            
    end %Fine for hk=0:n_hat-1
        
    %Memorizziamo nel vettore alpha di formato 3*N_triangles x 1 la 
    %densità ottenuta in corrispondenza del sottointervallo temporale di 
    %indice n_hat
    alpha=density(:,n_hat+1);

    %Aggiorniamo la soluzione U
    U=U+matrix0*alpha;
        
end %Fine if t>0

return