function [RIST]=rot_inv_angoli_eulero(Sistema_rif,vnF,rist)

%La function rot_inv_angoli_eulero() prende in input:
%- la struct Sistema_rif che contiene le informazioni sul sistema di
%  riferimento che abbiamo fissato sul triangolo di campo
%
%- il vettore 3x1 vnF che contiene il versore normale al triangolo di
%  campo
%
%- la matrice 3x3 rist che contiene i risultati degli integrali di
%  r_ir_j/r^3 e r_ir_j/r^5 calcolati sul triangolo di campo nel nuovo 
%  sistema di riferimento che vanno moltiplicati per i coefficienti di 
%  "correzione" derivanti dal cambiamento di sistema di riferimento, 
%  in modo tale da ricondurre questi integrali ad essere integrali
%  calcolati nel sitema di riferimento canonico di R3
%
%e restituisce in OUTPUT:
%- la matrice 3x3 RIST contenente i risultati degli integrali "corretti"
%  in seguito al cambiamento di sistema di riferimento


%Definiamo i versori della terna (X,Y,Z) che individuano il sistema di 
%riferimento fissato sul triangolo di campo 
X=Sistema_rif.e1;
%Y=Sistema_rif.e2;
Z=vnF;

%Definiamo i versori della terna (e1,e2,e3) che individuano il sistema di 
%riferimento canonico di R3
e1=[1 0 0];
%e2=[0 1 0];
e3=[0 0 1];

%Determiniamo i tre angoli di Eulero (phi,omega,psi) che permettono di
%ruotare il sistema di riferimento canonico di R3 fino a portarlo in un
%sistema di riferimento che ha i versori paralleli ai versori che
%individuano il sistema di riferimento fissato sul triangolo di campo

%Calcoliamo la LINEA dei NODI
if (norm(cross(e3,Z),2)<1.0e-10)
    %Se il prodotto vettoriale tra il versore e3 normale al piano
    %(e1,e2) e il versore Z normale al piano (X,Y) è uguale a zero
    %significa che i PIANI (e1,e2) e (X,Y) risultano PARALLELI.
    %Quindi per convenzione si fa coincidere la LINEA dei NODI 
    %con l'ASSE X
    nodi=X;    
    
    %--------------------------------------------------------------------
    %Calcoliamo l'ANGOLO di PRECESSIONE phi (orientato in SENSO ANTIORARIO)
    %tra l'ASSE e1 e la LINEA dei NODI, che può variare in [0,2*pi).
    cos_phi=e1*nodi';
    phi=acos(cos_phi);
    if (cross(e1,nodi)*e3'<-1.0e-10)
        %Se il prodotto vettoriale tra l'asse e1 e la linea dei nodi è 
        %uguale al vettore nullo oppure è un vettore che ha lo stesso
        %verso del versore e3, allora l'angolo compreso tra e1 e la  
        %linea dei nodi appartiene all'intervallo [0,pi] e quindi 
        %phi=acos(cos_phi). Questa situazione si traduce con il richiedere
        %che il prodotto scalare <(e1 x nodi),e3> sia maggiore o uguale
        %di zero.
        %In caso contrario abbiamo che i vettori e1 x nodi e e3 sono 
        %antiparalleli e quindi l'angolo phi è uguale a
        %phi=2*pi-acos(cos_phi)
        phi=2*pi-phi;
    end
    %L'angolo di precessione permette di ruotare il sistema di 
    %riferimento attorno all'asse e3 in modo tale da portare l'asse e1 
    %a coincidere con la linea dei nodi
   
    %--------------------------------------------------------------------
    %Calcoliamo l'ANGOLO di NUTAZIONE omega compreso tra gli ASSI e3 e Z
    %Siamo nel caso in cui i PIANI (e1,e2) e (X,Y) risultano
    %PARALLELI e quindi l'angolo compreso tra e3 e Z è uguale a zero se
    %essi hanno lo stesso verso oppure è uguale a pi se hanno verso 
    %opposto
    if (e3*Z'<0)
        omega=pi;
    else
        omega=0;
    end
    %L'angolo di nutazione permette di ruotare il sistema di riferimento 
    %attorno alla linea dei nodi in modo tale da portare l'asse e3 a 
    %coincidere con l'asse Z
    
    %--------------------------------------------------------------------
    %Se siamo in questo caso, abbiamo che automaticamente l'ANGOLO di 
    %ROTAZIONE PROPRIA psi compreso tra la linea
    %dei nodi e l'asse X è uguale a zero
    psi=0;
    %--------------------------------------------------------------------

else
    %Se il prodotto vettoriale tra il versore e3 normale al piano
    %(e1,e2) e il versore Z normale al piano (X,Y) è maggiore di zero
    %significa che i PIANI (e1,e2) e (X,Y) NON sono PARALLELI.
    %Siccome i piani (e1,e2) e (X,Y) NON sono PARALLELI, determiniamo la
    %LINEA dei NODI attraverso il prodotto vettoriale tre e3 e Z
    nodi=cross(e3,Z)/norm(cross(e3,Z),2);   
    
    %--------------------------------------------------------------------
    %Calcoliamo l'ANGOLO di PRECESSIONE phi (orientato in SENSO ANTIORARIO)
    %tra l'ASSE e1 e la LINEA dei NODI, che può variare in [0,2*pi).
    cos_phi=e1*nodi';
    phi=acos(cos_phi);
    if (cross(e1,nodi)*e3'<-1.0e-10)
        %Se il prodotto vettoriale tra l'asse e1 e la linea dei nodi è 
        %uguale al vettore nullo oppure è un vettore che ha lo stesso
        %verso del versore e3, allora l'angolo compreso tra e1 e la  
        %linea dei nodi appartiene all'intervallo [0,pi] e quindi 
        %phi=acos(cos_phi). Questa situazione si traduce con il richiedere
        %che il prodotto scalare <(e1 x nodi),e3> sia maggiore o uguale
        %di zero.
        %In caso contrario abbiamo che i vettori e1 x nodi e e3 sono 
        %antiparalleli e quindi l'angolo phi è uguale a
        %phi=2*pi-acos(cos_phi)
        phi=2*pi-phi;
    end
    %L'angolo di precessione permette di ruotare il sistema di 
    %riferimento attorno all'asse e3 in modo tale da portare l'asse e1 
    %a coincidere con la linea dei nodi
    
    %--------------------------------------------------------------------
    %Calcoliamo l'ANGOLO di NUTAZIONE omega compreso tra gli ASSI e3 e Z,
    %che appartiene sempre all'intervallo [0,pi]
    cos_omega=e3*Z';
    omega=acos(cos_omega);
    %L'angolo di nutazione permette di ruotare il sistema di riferimento 
    %attorno alla linea dei nodi in modo tale da portare l'asse e3 a 
    %coincidere con l'asse Z
    
    %--------------------------------------------------------------------
    %Calcoliamo l'ANGOLO di ROTAZIONE PROPRIA psi compreso tra la linea
    %dei nodi e l'asse X, che appartiene all'intervallo [0,2*pi)
    cos_psi=nodi*X';
    psi=acos(cos_psi);
    psi=real(psi);
    if (cross(nodi,X)*Z'<-1.0e-10)
        %Se il prodotto vettoriale tra la linea dei nodi e l'asse X è 
        %uguale al vettore nullo oppure è se è un vettore che ha lo
        %stesso verso del versore Z, allora l'angolo compreso tra la  
        %linea dei nodi e l'asse X è minore o uguale di pi e quindi 
        %psi=acos(cos_psia).
        %In caso contrario, psi=2*pi-acos(cos_psi)
        psi=2*pi-psi;

    end
    %L'angolo di rotazione propria permette di ruotare il sistema di 
    %riferimento attorno all'asse Z in modo tale da portare la linea dei 
    %nodi a coincidere con l'asse X
    %--------------------------------------------------------------------

end

%     %Calcoliamo la LINEA dei NODI
%     if (norm(cross(e3,Z),2)<1.0e-10)
%         %Se il prodotto vettoriale tra il versore e3 normale al piano
%         %(e1,e2) e il versore Z normale al piano (X,Y) risulta
%         %uguale a zero significa che i PIANI (e1,e2) e (X,Y) risultano
%         %PARALLELI.
%         %Quindi la LINEA dei NODI coincide con l'ASSE X
%         nodi=X;
%  
%     else
%         %Se i piani (e1,e2) e (X,Y) NON sono PARALLELI determiniamo la
%         %LINEA dei NODI attraverso il prodotto vettoriale tre e3 e Z
%         nodi=cross(e3,Z);        
%         
%     end
%     
%     %--------------------------------------------------------------------
%     %Calcoliamo l'ANGOLO di PRECESSIONE phi (orientato in SENSO ANTIORARIO)
%     %tra l'ASSE e1 e la LINEA dei NODI, che può variare in [0,2*pi).
%     cos_phi=e1*nodi';
%     phi=acos(cos_phi);
%         
%     if (cross(e1,nodi)*e3'<0)
%         %Se il prodotto vettoriale tra l'asse e1 e la linea dei nodi ha lo
%         %stesso verso del versore e3, allora l'angolo compreso tra e1 e la  
%         %linea dei nodi è minore o uguale di pi e quindi 
%         %phi=acos(cos_phi).
%         %In caso contrario, phi=2*pi-acos(cos_phi)
%         phi=2*pi-phi;
%         
%     end
%     %L'angolo di precessione permette di ruotare il sistema di 
%     %riferimento attorno all'asse e3 in modo tale da portare l'asse e1 
%     %a coincidere con la linea dei nodi
%     
%     %--------------------------------------------------------------------
%     %Calcoliamo l'ANGOLO di NUTAZIONE omega compreso tra gli ASSI e3 e Z,
%     %che appartiene sempre all'intervallo [0,pi]
%     cos_omega=e3*Z';
%     omega=acos(cos_omega);
%     %L'angolo di nutazione permette di ruotare il sistema di riferimento 
%     %attorno alla linea dei nodi in modo tale da portare l'asse e3 a 
%     %coincidere con l'asse Z
%     
%     %--------------------------------------------------------------------
% 
%     %Calcoliamo l'ANGOLO di ROTAZIONE PROPRIA psi compreso tra la linea
%     %dei nodi e l'asse X, che appartiene all'intervallo [0,2*pi)
%     
% 
%     cos_psi=nodi*X';
%     psi=acos(cos_psi);
%     psi=real(psi);
% 
%     if (cross(nodi,X)*Z'<0)
%         %Se il prodotto vettoriale tra la linea dei nodi e l'asse X ha lo
%         %stesso verso del versore Z, allora l'angolo compreso tra la  
%         %linea dei nodi e l'asse X è minore o uguale di pi e quindi 
%         %psi=acos(cos_psia).
%         %In caso contrario, psi=2*pi-acos(cos_psi)
%         psi=2*pi-psi;
% 
%     end
%     %L'angolo di rotazione propria permette di ruotare il sistema di 
%     %riferimento attorno all'asse Z in modo tale da portare la linea dei 
%     %nodi a coincidere con l'asse X

    
%Calcoliamo le tre matrici relative agli angoli di Eulero trovati che
%permetteranno poi di effettuare il cambiamento del sistema di riferimento
A_psi=[cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
A_omega=[1 0 0;0 cos(omega) sin(omega); 0 -sin(omega) cos(omega)];
A_phi=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1];

%Calcoliamo il prodotto della tre matrici A_phi, A_omega e A_psi
C=A_psi*A_omega*A_phi;

%Buttiamo via la "SPORCIZIA NUMERICA" dalla matrice 3x3 A
C(abs(C)<1.0e-14)=0;

%Riconduciamo i risultati degli integrali calcolati nel sistema di 
%riferimento fissato sul triangolo di campo ad essere integrali calcolati 
%nel sistema di riferimento canonico di R3 attraverso la matrice C che
%permette di effettuare il cambiamento del sistema di riferimento
RIST=C'*(rist*C);
   
return