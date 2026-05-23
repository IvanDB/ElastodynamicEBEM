function [ris1, ris2, ris3, ris4, ris5, ris6] = BEMenerg_coreMXdiag_calcIntgInt_childTriangle(tri2D, infoTri2D, rhoP, rhoS)
%INPUT
% - 
%
%OUTPUT
% - 
% - 
zeta = 0;
%% ESTRAZIONE INFORMAZIONI GEOMETRICHE TRIANGOLO
%Estrazione LUNGHEZZA PRIMO LATO a del TRIANGOLO FIGLIO 
a = infoTri2D.a;

%Estrazione seno secondo angolo
sin_beta = infoTri2D.sin_beta;

%Estrazione ampiezza secondo e terzo angolo
beta = infoTri2D.beta;
gamma = infoTri2D.gamma;

%% CALCOLO MATRICI B_i

%Calcolo coefficienti di correzione necessari a seguito della rotazione
matB = BEMenerg_coreMXdiag_coeffB(infoTri2D);

%% CALCOLO INTERSEZIONI RETTA SECONDO LATO TRIANGOLO e CIRCONFERENZE P ed S

%Calcolo VERSORE retta LATO V2V3 TRIANGOLO
vt = (tri2D(3, :) - tri2D(2, :)) / infoTri2D.c;

%Calcolo DELTA equazione di secondo grado derivante dal sistema fra la
% retta e la circonferenza di raggio rho_P
B_mezzi = a * vt(1);
C = a^2 - rhoP^2;
DeltaP = B_mezzi^2 - C;

%Calcolo DELTA equazione di secondo grado derivante dal sistema fra la
% retta e la circonferenza di raggio rho_S
B_mezzi = a * vt(1);
C = a^2 - rhoS^2;
DeltaS = B_mezzi^2 - C;

%Calcolo costanti intermedie
RP = a * sin_beta / rhoP;
RS = a * sin_beta / rhoS;

%% INTEGRAZIONE
%Inizializzazione variabili risultato
ris1 = zeros(3, 3); %integrale di 1/r (onda P)
ris2 = zeros(3, 3); %integrale di r_i*r_j/r^3 (onda P)
ris3 = zeros(3, 3); %integrale di 1/r^3 (onda P)
ris4 = zeros(3, 3); %integrale di r_i*r_j/r^5 (onda P) 
ris5 = zeros(3, 3); %integrale di 1/r (onde S e P)
ris6 = zeros(3, 3); %integrale di r_i*r_j/r^3 (onde S e P)

%Distinzione casi differenti tipologie di intersezione
if(DeltaP <= 1.0e-06) 
    %Caso di nessuna intersezione o intersezioni coicidenti per il fronte P.
    % Dominio composto da una singola corona circolare

    %Calcolo valori angolari definenti i domini (in questo caso solo la
    % prima corona circolare è non degenere)
    THETA1 = 0;
    THETA2 = gamma;
    THETA3 = gamma; %regioni E2 E7 sicuramente degeneri
    THETA4 = gamma; %regione E3 sicuramente degenere
    THETA5 = gamma; %regioni E4 E8 sicuramente degeneri
    THETA6 = gamma; %regioni E5 E6 sicuramente degeneri

elseif ((DeltaP > 1.0e-06) && (DeltaS <= 1.0e-06))
    %Caso di due intersezioni non coicidenti per il fronte P e 
    % nessuna intersezione o intersezioni coicidenti per il fronte S.
    % Dominio composto da due corone circolari ed trapezio/doppio
    % triangolo circolare.
    
    %Calcolo valori angolari assoluti
    Omega1 = asin(RP) - beta;
    Omega2 = pi - asin(RP) - beta;
    
    %Ordinamento valori angolari assoluti
    theta1 = min(Omega1, Omega2);
    theta2 = max(Omega1, Omega2);
    
    %Calcolo valori angolari definenti i domini
    THETA1 = 0;
    THETA2 = min(max(0, theta1), gamma);
    THETA3 = max(0, min(theta2, gamma));
    THETA4 = max(0, min(theta2, gamma)); %regione E3 sicuramente degenere
    THETA5 = max(0, min(theta2, gamma)); %regioni E5 E6 sicuramente degeneri
    THETA6 = gamma;        
else
    %Caso di due intersezioni non coicidenti per il fronte P
    % e due intersezioni non coicidenti per il fronte S.
    % Dominio composto da due corone circolari, due triangoli circolari,
    % un triangolo e quattro settori circolari

    %Calcolo valori angolari assoluti
    Omega1 = asin(RP) - beta;
    Omega2 = pi - asin(RP) - beta;
    Omega3 = asin(RS) - beta;
    Omega4 = pi - asin(RS) - beta;
    
    %Ordinamento valori angolari assoluti
    theta1 = min(Omega1, Omega2);
    theta2 = max(Omega1, Omega2);
    theta3 = min(Omega3, Omega4);
    theta4 = max(Omega3, Omega4);
    
    %Calcolo valori angolari definenti i domini
    THETA1 = 0;
    THETA2 = min(max(0, theta1), gamma);
    THETA3 = min(max(0, theta3), gamma);
    THETA4 = max(0, min(theta4, gamma));
    THETA5 = max(0, min(theta2, gamma));
    THETA6 = gamma;
end


%Controllo che PRIMO SETTORE CIRCOLARE (regione E1) e
% PRIMA CORONA CIRCOLARE (regione E6) siano NON DEGENERI.
if (THETA1 < THETA2)
    %Calcolo integrali su SETTORE CIRCOLARE (regione E1)
    [ris5t, ris6t] = BEMenerg_coreMXdiag_esplicitIntg_SC(rhoS, THETA1, THETA2, infoTri2D, matB);
    %Calcolo integrali su CORONA CIRCOLARE (regione E6)
    [ris1t, ris2t, ris3t, ris4t] = BEMenerg_coreMXdiag_esplicitIntg_CC(rhoS, rhoP, THETA1, THETA2, infoTri2D, matB);

    %Aggiornamento valore integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    ris5 = ris5 + ris5t;
    ris6 = ris6 + ris6t;
end

%Controllo che SECONDO SETTORE CIRCOLARE (regione E2) e
% PRIMO TRIANGOLO CIRCOLARE (regione E7) siano NON DEGENERI.
if (THETA2 < THETA3)
    %Calcolo integrali su SETTORE CIRCOLARE (regione E2)
    [ris5t, ris6t] = BEMenerg_coreMXdiag_esplicitIntg_SC(rhoS, THETA2, THETA3, infoTri2D, matB);

    %Calcolo integrali su TRIANGOLO CIRCOLARE (regione E7)
    [ris1t, ris2t, ris3t, ris4t] = BEMenerg_coreMXdiag_esplicitIntg_TC(rhoS, THETA2, THETA3, infoTri2D, matB);

    %Aggiornamento valore integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    ris5 = ris5 + ris5t;
    ris6 = ris6 + ris6t;
end
   
%Controllo che TRIANGOLO (Regione E3) sia NON DEGENERE.
if (THETA3 < THETA4)
    %Calcolo integrali su TRIANGOLO (Regione E3)
    [ris5t, ris6t] = BEMenerg_coreMXdiag_esplicitIntg_TR(THETA3, THETA4, infoTri2D, matB);
    %Aggiornamento valore integrali globali
    ris5 = ris5 + ris5t;
    ris6 = ris6 + ris6t;
end

%Controllo che TERZO SETTORE CIRCOLARE (regione E4) e
% SECONDO TRIANGOLO CIRCOLARE (regione E8) siano NON DEGENERI.
if (THETA4 < THETA5)
    %Calcolo integrali su TERZO SETTORE CIRCOLARE (regione E4)
    [ris5t, ris6t] = BEMenerg_coreMXdiag_esplicitIntg_SC(rhoS, THETA4, THETA5, infoTri2D, matB);

    %Calcolo integrali sul SECONDO TRIANGOLO CIRCOLARE (regione E8)
    [ris1t, ris2t, ris3t, ris4t] = BEMenerg_coreMXdiag_esplicitIntg_TC(rhoS, THETA4, THETA5, infoTri2D, matB);

    %Aggiornamento valore integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    ris5 = ris5 + ris5t;
    ris6 = ris6 + ris6t;
end
 
%Controllo che QUARTO SETTORE CIRCOLARE (regione E5) e
% SECONDA CORONA CIRCOLARE (regione E9) siano NON DEGENERI.
if (THETA5 < THETA6)
    %Calcolo integrali su SETTORE CIRCOLARE (regione E5)
    [ris5t, ris6t] = BEMenerg_coreMXdiag_esplicitIntg_SC(rhoS, THETA5, THETA6, infoTri2D, matB);

    %Calcolo integrali su CORONA CIRCOLARE (regione E9)
    [ris1t, ris2t, ris3t, ris4t] = BEMenerg_coreMXdiag_esplicitIntg_CC(rhoS, rhoP, THETA5, THETA6, infoTri2D, matB);

    %Aggiorniamo il valore degli integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    ris5 = ris5 + ris5t;
    ris6 = ris6 + ris6t;  
end 

return