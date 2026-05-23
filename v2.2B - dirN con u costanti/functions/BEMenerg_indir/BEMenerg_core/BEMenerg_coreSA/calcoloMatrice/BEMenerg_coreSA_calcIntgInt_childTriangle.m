function [ris1, ris2, ris3, ris4, ris5, ris6] = BEMenerg_coreSA_calcIntgInt_childTriangle(tri2D, infoTri2D, zeta, RP, rhoP, RS, rhoS, sign_prod, typeIntg)
%INPUT
% - 
%
%OUTPUT
% - 
% - 

%% ESTRAZIONE INFORMAZIONI GEOMETRICHE TRIANGOLO
%Estrazione LUNGHEZZA PRIMO LATO a del TRIANGOLO FIGLIO 
a = infoTri2D.a;

%Estrazione seno secondo angolo
sinBeta = infoTri2D.sin_beta;

%Estrazione ampiezza secondo e terzo angolo
beta = infoTri2D.beta;
gamma = infoTri2D.gamma;

%% CALCOLO MATRICI B_i

%Calcolo coefficienti di correzione necessari a seguito della rotazione
matB = BEMenerg_coreSA_coeffB(zeta, infoTri2D, sign_prod);

%% CALCOLO INTERSEZIONI RETTA SECONDO LATO TRIANGOLO e CIRCONFERENZE P ed S

%Calcolo VERSORE retta LATO V2V3 TRIANGOLO
vt = (tri2D(3, :) - tri2D(2, :)) / infoTri2D.c;

%Calcolo DELTA equazione di secondo grado derivante dal sistema fra la
% retta e la circonferenza di raggio rho_P
Bmezzi = a * vt(1);
C = a^2 - rhoP^2;
DeltaP = Bmezzi^2 - C;

%Calcolo DELTA equazione di secondo grado derivante dal sistema fra la
% retta e la circonferenza di raggio rho_S
Bmezzi = a * vt(1);
C = a^2 - rhoS^2;
DeltaS = Bmezzi^2 - C;

%Calcolo costanti intermedie
rP = a * sinBeta / rhoP;
rS = a * sinBeta / rhoS;

%Reset valori S in caso di sole onde P
if(typeIntg == "P")
    rS = 0;
    DeltaS = 0;
end

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
    Omega1 = asin(rP) - beta;
    Omega2 = pi - asin(rP) - beta;
    
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
    Omega1 = asin(rP) - beta;
    Omega2 = pi - asin(rP) - beta;
    Omega3 = asin(rS) - beta;
    Omega4 = pi - asin(rS) - beta;
    
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
    if(typeIntg == "SP_SP")
        [ris5t, ris6t] = BEMenerg_analiticIntg_SC_SP(RS, rhoS, zeta, THETA1, THETA2, infoTri2D, matB);
    end
    
    %Calcolo integrali su CORONA CIRCOLARE (regione E6)
    if(typeIntg == "P")
        [ris1t, ris2t, ris3t, ris4t] = BEMenerg_analiticIntg_SC_P(RP, rhoP, zeta, THETA1, THETA2, infoTri2D, matB);
    else
        [ris1t, ris2t, ris3t, ris4t] = BEMenerg_analiticIntg_CC(RS, rhoS, RP, rhoP, zeta, THETA1, THETA2, infoTri2D, matB);
    end

    %Aggiornamento valore integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    if(typeIntg == "SP_SP")
        ris5 = ris5 + ris5t;
        ris6 = ris6 + ris6t;
    end
end

%Controllo che SECONDO SETTORE CIRCOLARE (regione E2) e
% PRIMO TRIANGOLO CIRCOLARE (regione E7) siano NON DEGENERI.
if (THETA2 < THETA3)
    %Calcolo integrali su SETTORE CIRCOLARE (regione E2)
    if(typeIntg == "SP_SP")
        [ris5t, ris6t] = BEMenerg_analiticIntg_SC_SP(RS, rhoS, zeta, THETA2, THETA3, infoTri2D, matB);
    end

    %Calcolo integrali su TRIANGOLO CIRCOLARE (regione E7)
    if(typeIntg == "P")
        [ris1t, ris2t, ris3t, ris4t] = BEMenerg_analiticIntg_TR_P(zeta, THETA2, THETA5, infoTri2D, matB);
    else
        [ris1t, ris2t, ris3t, ris4t] = BEMenerg_analiticIntg_TC(zeta, RS, rhoS, THETA2, THETA3, infoTri2D, matB);
    end
    %Aggiornamento valore integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    if(typeIntg == "SP_SP")
        ris5 = ris5 + ris5t;
        ris6 = ris6 + ris6t;
    end
end
   
%Controllo che TRIANGOLO (Regione E3) sia NON DEGENERE.
if (THETA3 < THETA4)
    %Calcolo integrali su TRIANGOLO (Regione E3)
    if(typeIntg == "SP_SP")
        [ris5t, ris6t] = BEMenerg_analiticIntg_TR_SP(zeta, THETA3, THETA4, infoTri2D, matB);
    end
    %Aggiornamento valore integrali globali
    if(typeIntg == "SP_SP")
        ris5 = ris5 + ris5t;
        ris6 = ris6 + ris6t;
    end
end

%Controllo che TERZO SETTORE CIRCOLARE (regione E4) e
% SECONDO TRIANGOLO CIRCOLARE (regione E8) siano NON DEGENERI.
if (THETA4 < THETA5)
    %Calcolo integrali su TERZO SETTORE CIRCOLARE (regione E4)
    if(typeIntg == "SP_SP")
        [ris5t, ris6t] = BEMenerg_analiticIntg_SC_SP(RS, rhoS, zeta, THETA4, THETA5, infoTri2D, matB);
    end

    %Calcolo integrali sul SECONDO TRIANGOLO CIRCOLARE (regione E8)
    [ris1t, ris2t, ris3t, ris4t] = BEMenerg_analiticIntg_TC(zeta, RS, rhoS, THETA4, THETA5, infoTri2D, matB);

    %Aggiornamento valore integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    if(typeIntg == "SP_SP")
        ris5 = ris5 + ris5t;
        ris6 = ris6 + ris6t;
    end
end
 
%Controllo che QUARTO SETTORE CIRCOLARE (regione E5) e
% SECONDA CORONA CIRCOLARE (regione E9) siano NON DEGENERI.
if (THETA5 < THETA6)
    %Calcolo integrali su SETTORE CIRCOLARE (regione E5)
    if(typeIntg == "SP_SP")
        [ris5t, ris6t] = BEMenerg_analiticIntg_SC_SP(RS, rhoS, zeta, THETA5, THETA6, infoTri2D, matB);
    end

    %Calcolo integrali su CORONA CIRCOLARE (regione E9)
    if(typeIntg == "P")
        [ris1t, ris2t, ris3t, ris4t] = BEMenerg_analiticIntg_SC_P(RP, rhoP, zeta, THETA5, THETA6, infoTri2D, matB);
    else
        [ris1t, ris2t, ris3t, ris4t] = BEMenerg_analiticIntg_CC(RS, rhoS, RP, rhoP, zeta, THETA5, THETA6, infoTri2D, matB);
    end

    %Aggiorniamo il valore degli integrali globali
    ris1 = ris1 + ris1t;
    ris2 = ris2 + ris2t;
    ris3 = ris3 + ris3t;
    ris4 = ris4 + ris4t;
    if(typeIntg == "SP_SP")
        ris5 = ris5 + ris5t;
        ris6 = ris6 + ris6t;
    end   
end 
return