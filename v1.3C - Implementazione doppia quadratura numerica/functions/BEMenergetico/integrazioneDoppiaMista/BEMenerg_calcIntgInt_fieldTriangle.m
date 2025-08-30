function [ris1, ris2, ris3, ris4, ris5, ris6] = BEMenerg_calcIntgInt_fieldTriangle(zeta, R_P, rho_P, R_S, rho_S, children, c, sign_prod, type_intg)

%Inizializzazione variabili risultato
ris1 = zeros(3, 3); %integrale di 1/r (onda P)
ris2 = zeros(3, 3); %integrale di r_i*r_j/r^3 (onda P)
ris3 = zeros(3, 3); %integrale di 1/r^3 (onda P)
ris4 = zeros(3, 3); %integrale di r_i*r_j/r^5 (onda P) 
ris5 = zeros(3, 3); %integrale di 1/r (onde S e P)
ris6 = zeros(3, 3); %integrale di r_i*r_j/r^3 (onde S e P)


%CICLO sui 3 TRIANGOLI FIGLI
for ind_child = 1 : 3
    
    %% OPERAZIONI PRELIMINARI

    %Estrazione COORDINATE VERTICI TRIANGOLO FIGLIO corrente
    child = squeeze(children(ind_child, :, :)); 
        
    %Estrazione COEFFICIENTE del TRIANGOLO FIGLIO corrente
    C = c(ind_child);

    %Controllo che il TRIANGOLO FIGLIO sia NON DEGENERE
    if(abs(C) < 1.0e-06) 
        continue
    end

    %Controllo segno del coefficiente ed eventuale scambio di vertici
    % per garantirne la disposizione in SENSO ANTIORARIO
    if(C < 0)
        temp = child(2, :);
        child(2, :) = child(3, :);
        child(3, :) = temp;
    end
    
    %% ROTAZIONE TRIANGOLO FIGLIO CORRENTE

    %Calcolo COORDINATA POLARE individuante il SECONDO VERTICE nel piano
    % (angolo rispetto al semiasse positivo delle ascisse) 
    theta0 = atan2(child(2, 2), child(2, 1));
    
    %Calcolo valore di rotazione antioraria necessario per portare il primo
    % vertice sul semiasse x positivo
    if (theta0 < 0)
        theta0 = abs(theta0);
    else
        theta0 = 2*pi - theta0;
    end
    
    %Calcolo vertici POST_ROTAZIONE e INFORMAZIONI GEOMETRICHE del triangolo figlio corrente
    [tri2D, infoTri2D] = calcoloDatiTriangolo(child);
    
    %% CALCOLO INTEGRALE ANALITICO sul TRIANGOLO FIGLIO
    [ris1t, ris2t, ris3t, ris4t,ris5t, ris6t] = BEMenerg_calcIntgInt_childTriangle(tri2D, infoTri2D, zeta, R_P, rho_P, R_S, rho_S, sign_prod, type_intg);

    %% APPLICAZIONE CORRETTIVO ROTAZIONE INVERSA

    %Calcolmo coseno seno angolo theta0
    cos_theta0 = cos(theta0);
    sin_theta0 = sin(theta0);

    %Costruzione matrice di rotazione inversa (rotazione in senso orario
    % attorno all'asse z)
    A = [cos_theta0 sin_theta0  0;
        -sin_theta0 cos_theta0  0;
         0          0           1];

    %Applicazione rotazione
    ris2t = A*(ris2t*A');
    ris4t = A*(ris4t*A');
    ris6t = A*(ris6t*A');
    
    %% AGGIORNAMENTO RISULTATO TOTALE
    
    %Calcolo coefficiente
    coeff = sign(C);
    
    %Somma dei contributi
    ris1 = ris1 + ris1t*coeff; %Integrale di 1/r (Onda P)
    ris2 = ris2 + ris2t*coeff; %Integrale di r_i*r_j/r^3 (Onda P)
    ris3 = ris3 + ris3t*coeff; %Integrale di 1/r^3 (Onda P)
    ris4 = ris4 + ris4t*coeff; %Integrale di r_i*r_j/r^5 (Onda P)
    ris5 = ris5 + ris5t*coeff; %Integrale di 1/r (Onde S e P)
    ris6 = ris6 + ris6t*coeff; %Integrale di r_i*r_j/r^3 (Onde S e P)
             
end
return