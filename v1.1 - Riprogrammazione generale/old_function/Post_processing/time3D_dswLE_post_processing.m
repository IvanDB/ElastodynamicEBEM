function ris=time3D_dswLE_post_processing(pb_param,TF,vnF,thk,sp)

%La function time3D_dswLE prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- la matrice 3x3 TF in cui ogni riga contiene le coordinate di uno dei 
%  tre vertici del triangolo di campo  
%
%- il vettore 1x3 vnF contenente le coordinate del versore normale al 
%  triangolo di campo
%
%- la variabile thk (=delta) che contiene la differenza tra gli istanti 
%  temporali relativa al blocco matriciale che si vuole calcolare 
%
%- il vettore 1x3 sp contenente le coordinate del punto sorgente
%
%e restituisce in OUTPUT:
%- la matrice 3x3 ris contenente il risultato dell'integrazione analitica
%  sul triangolo di campo considerato considernado come punto sorgente il
%  punto x
%
%-------------------------------------------------------------------------

%RAGGIO della SUPERFICIE SFERICA individuata dal FRONTE d'ONDA relativo 
%alle ONDE P
R_P = pb_param.velP*thk; %c_P*Delta

%RAGGIO della SUPERFICIE SFERICA individuata dal FRONTE d'ONDA relativo 
%alle ONDE S
R_S = pb_param.velS*thk; %c_S*Delta


%Ricaviamo le informazioni utili a individuare gli assi del nuovo 
%sistema di riferimento
Sistema_rif=reference_system(TF,vnF);

%Inizializzazione delle matrici 3x3 contenenti i RISULTATI PARZIALI 
%dell'INTEGRAZIONE ANALITICA
rist1 = zeros(3,3); %Integrale di 1/r (Onda P)
rist2 = zeros(3,3); %Integrale di r_i*r_j/r^3 (Onda P)
rist3 = zeros(3,3); %Integrale di 1/r^3 (Onda P)
rist4 = zeros(3,3); %Integrale di r_i*r_j/r^5 (Onda P)
rist5 = zeros(3,3); %Integrale di 1/r (Onda P e S)
rist6 = zeros(3,3); %Integrale di r_i*r_j/r^3 (Onda P e S)
    
    
%Spostiamoci sul PIANO individuato dal TRIANGOLO di CAMPO e fissiamo 
%un SISTEMA di RIFERIMENTO che ha come ORIGINE la PROIEZIONE del
%PUNTO SORGENTE e come ASSI CARTESIANI gli assi memorizzati
%all'interno della struct Sistema_rif
[zeta,tri2D,children,eta_children,c,sign_prod] = prep_triangle(sp,TF,vnF,Sistema_rif);
       
if(R_P>zeta)    

    %Se R_P >zeta allora la sfera di raggio R_P interseca il piano
    %contenente il triangolo di campo in più di un punto e
    %possiamo procedere con lo studio delle diverse regioni di
    %integrazione.
    %In caso contrario, ovvero se R_P<=zeta, le regioni di 
    %integrazione sono vuote oppure collassano in un solo punto e 
    %quindi il risultato degli integrali è uguale a zero.

%         %Calcoliamo la DISTANZA di ciascuno dei TRE VERTICI del 
%         %TRIANGOLO di CAMPO corrente dall'ORIGINE (0,0) del nuovo 
%         %SISTEMA di RIFERIMENTO BIDIMENSIONALE 
%         d1=sqrt(sum(tri2D(1,:).^2));
%         d2=sqrt(sum(tri2D(2,:).^2));
%         d3=sqrt(sum(tri2D(3,:).^2));        
%          
%         %Calcoliamo il MINIMO tra queste tre DISTANZE
%         d_MIN=min([d1 d2 d3]);

    %Controlliamo se la PROIEZIONE del PUNTO SORGENTE cade 
    %all'INTERNO o all'ESTERNO del TRIANGOLO di CAMPO sfruttando 
    %le sue COORDINATE BARICENTRICHE
    %appartiene=0;
    if (min(-1.0e-6<=eta_children)==1)&&(min(eta_children<=1)==1)     
        appartiene=1;
        %La proiezione sp_plane del punto sorgente appartiene al
        %triangolo di campo se le sue coordinate baricentriche sono
        %comprese tra 0 e 1
    else
        appartiene=0;
        d_MIN=distanza_triangolo(tri2D);
    end
    %La variabile appartiene è uguale ad 1 se la proieizione del punto
    %sorgente appartiene al triangolo di campo, mentre è uguale a zero
    %nell'altro caso

    if(R_S>zeta && zeta <=1.06e-6)
        %CONTRIBUTO del FRONTE delle ONDE P e S nel caso in cui il
        %TRIANGOLO SORGENTE e di CAMPO sono COMPLANARI (zeta=0)       

        %In questo caso i RAGGI pR_S e pR_P dei CERCHI che 
        %rappresentano il FRONTE delle ONDE S e P nel PIANO del 
        %TRIANGOLO di CAMPO coincidono, rispettivamente, con i RAGGI
        %R_S e R_P delle SUPERFICI SFERICHE        
        pR_P=R_P; 
        pR_S=R_S;


        if (appartiene==1 || (appartiene==0 && pR_S>d_MIN))                
            %Se appartiene==1 il punto sp_plane appartiene al
            %triangolo di campo e quindi bisogna considerare il 
            %contributo delle ONDE S e P. 

            %Invece se appartiene==0 e pR_S>d_MIN il punto sp_plane
            %non appartiene al triangolo di campo ma il FRONTE 
            %dell'ONDA S INTERSECA il TRIANGOLO di CAMPO (quindi anche
            %l'onda P interseca il triangolo)

            %Calcolo dell'INTEGRALE ANALITICO INTERNO
            [rist1,rist2,rist3,rist4,rist5,rist6] = time3D_intLE_S_P_1(0,R_P,pR_P,R_S,pR_S,children,c,sign_prod);

        elseif (appartiene==0 && (pR_S<=d_MIN && pR_P>d_MIN))
            %Se appartiene==0 il punto sp_plane non appartiene al 
            %triangolo di campo.
            %Inoltre se pR_S<=d_MIN e pR_P>d_MIN SOLO FRONTE 
            %dell'ONDA P INTERSECA il TRIANGOLO di CAMPO (quindi non 
            %dobbiamo considerare il contributo dell'onda S)     

            %Calcolo dell'INTEGRALE ANALITICO INTERNO
            [rist1,rist2,rist3,rist4] = time3D_intLE_S_P_2(0,R_P,pR_P,R_S,pR_S,children,c,sign_prod);

        end %Fine if (appartiene==1 || (appartiene==0 && pR_S>d_MIN))

    elseif(R_S>zeta && zeta >1.06e-6)

        %CONTRIBUTO del FRONTE delle ONDE P e S nel caso in cui il
        %TRIANGOLO SORGENTE e di CAMPO NON sono COMPLANARI (zeta >0)

        %RAGGIO del CERCHIO RAPPRESENTANTE il FRONTE dell'ONDA P nel
        %PIANO del TRIANGOLO di CAMPO
         pR_P = sqrt(R_P^2-zeta^2);

        %RAGGIO del CERCHIO RAPPRESENTANTE il FRONTE dell'ONDA S nel
        %PIANO del TRIANGOLO di CAMPO
        pR_S = sqrt(R_S^2-zeta^2);

        if (appartiene==1 || (appartiene==0 && pR_S>d_MIN))

            %Se appartiene==1 il punto sp_plane appartiene al
            %triangolo di campo e quindi bisogna considerare il 
            %contributo delle ONDE S e P. 

            %Invece se appartiene==0 e pR_S>d_MIN il punto sp_plane
            %non appartiene al triangolo di campo ma il FRONTE 
            %dell'ONDA S INTERSECA il TRIANGOLO di CAMPO (quindi anche
            %l'onda P interseca il triangolo)

            %Calcolo dell'INTEGRALE ANALITICO INTERNO
            [rist1,rist2,rist3,rist4,rist5,rist6] = time3D_intLE_S_P_1(zeta,R_P,pR_P,R_S,pR_S,children,c,sign_prod);

        elseif (appartiene==0 && (pR_S<=d_MIN && pR_P>d_MIN))

            %Se appartiene==0 il punto sp_plane non appartiene al 
            %triangolo di campo.
            %Inoltre se pR_S<=d_MIN e pR_P>d_MIN SOLO FRONTE 
            %dell'ONDA P INTERSECA il TRIANGOLO di CAMPO (quindi non 
            %dobbiamo considerare il contributo dell'onda S)     

            %Calcolo dell'INTEGRALE ANALITICO INTERNO
            [rist1,rist2,rist3,rist4] = time3D_intLE_S_P_2(zeta,R_P,pR_P,R_S,pR_S,children,c,sign_prod);

        end %Fine if (appartiene==1 || (appartiene==0 && pR_S>MIN))

    else
%             disp('qui')
        %Siamo nel caso in cui R_S<=zeta e quindi abbiamo il            
        %CONTRIBUTO del FRONTE D'ONDA della SOLA ONDA P.

        %N.B. In questo caso i triangoli sorgente e di campo sono
        %sicuramente non complanari e quindi zeta>0

        %RAGGIO del CERCHIO RAPPRESENTANTE il FRONTE dell'ONDA P nel
        %PIANO del TRIANGOLO di CAMPO
        pR_P = sqrt(R_P^2-zeta^2);

        if (appartiene==1 || (appartiene==0 && pR_P>d_MIN))

           %Calcoliamo gli integrali analitici solo nei casi in cui:
           %- la proiezione sp_plane del punto sorgente appartiene 
           %  al triangolo di campo (appartiene=1)
           %- la proiezione sp_plane del punto sorgente non 
           %  appartiene al triangolo di campo (appartiene=0) ma il
           %  FRONTE dell'ONDA P INTERSECA il TRIANGOLO di CAMPO
           %  (quindi pR_P è maggiore del minimo tra le DISTANZA dei 
           %  tre VERTICI dall'origine)     

           %Calcolo dell'INTEGRALE ANALITICO INTERNO
           [rist1,rist2,rist3,rist4] = time3D_intLE_P(zeta,R_P,pR_P,children,c,sign_prod);

       end %Fine if (appartiene==1 || (appartiene==0 && pR_P>MIN))

    end %Fine if(R_S>zeta && zeta<=1.0e-6)                 

end%Fine if(R_P>zeta)


rist2=rot_inv_angoli_eulero(Sistema_rif,vnF,rist2);
rist4=rot_inv_angoli_eulero(Sistema_rif,vnF,rist4);
rist6=rot_inv_angoli_eulero(Sistema_rif,vnF,rist6);
    
    
%Memorizziamo nella matrice 3x3 ris il RISULTATO dell'integrazione sul 
%triangolo di campo corrente (integrazione analitica)
ris = ((rist1-rist2)/pb_param.velP^2+thk^2*(-rist3+3*rist4)+((pb_param.velP^2+pb_param.velS^2)*rist5...
    +(pb_param.velP^2-pb_param.velS^2)*rist6)/(pb_param.velP*pb_param.velS)^2)/2; 
    
ris=ris/(4*pi*pb_param.rho);
return