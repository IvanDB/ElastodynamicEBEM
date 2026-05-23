function ris = BEMenerg_coreMXdiag_calcSubBlock(pbParam, methodInfo, vertsT, vn, indTemp, deltaT, constValues)
%INPUT
% - 
% -
%
%OUTPUT
% - 

%% INIZIALIZZAZIONE VALORI

%Inizializzazione matrice 
ris = zeros(3, 3);

%Inizializzazione coefficienti
coef = [1, -2, 1];

%Estrazione dati
nGH = methodInfo.numNodiExt;

velP = pbParam.velP;
velS = pbParam.velS;

SdR = constValues.SdR;
GHnodes = constValues.GHnodes;
GHweights = constValues.GHweights;

%% CALCOLO INTEGRALI

%Calcolo indici temporali
istTemp = indTemp + [1, 0, -1];

%Ciclo sui tre istanti temporali
for var = 1 : 3
    %Controllo necessità di calcolo del sottonucleo
    if(istTemp(var) <= 0)
        continue
    end

    %Inizializzazione matrici contenenti RISULTATI COMPLESSIVI integrazione 
    % di ciascun tipo di nucleo
    ris_1 = zeros(3, 3); %Integrale di 1/r (Onda P)
    ris_2 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onda P)
    ris_3 = zeros(3, 3); %Integrale di 1/r^3 (Onda P)
    ris_4 = zeros(3, 3); %Integrale di r_i*r_j/r^5 (Onda P)
    ris_5 = zeros(3, 3); %Integrale di 1/r (Onde P e S)
    ris_6 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onde P e S)

    %Calcolo parametro temporale del sottonucleo
    currT = istTemp(var) * deltaT;
    
    %Calcolo RAGGI SUPERFICIE SFERICHE relative alle ONDE P ed S
    rP = velP * currT; % = c_P * Delta
    rS = velS * currT; % = c_S * Delta

    %Ciclo sui nodi di Gauss Hammer del triangolo sorgente
    for indGH = 1 : nGH    
        %Mappatura nodo sul triangolo sorgente 
        nodoGHcurr = GHnodes{indGH};

        %Spostamento su PIANO TRIANGOLO di CAMPO con SdR avente ORIGINE nella
        % PROIEZIONE del NODO corrente e ASSI CARTESIANI salvati in SdR
        children = preparazioneTriangoloDiagonale(nodoGHcurr, vertsT, constValues.SdR);
        
        %Calcolo integrazioni interne sul triangolo di campo
        [rist_1, rist_2, rist_3, rist_4, rist_5, rist_6] = BEMenerg_coreMXdiag_calcIntgInt_fieldTriangle(rP, rS, children);
        
        %ROTAZIONE finale all'indietro    
        rist_2 = rotazioneInversaEulero(SdR, vn, rist_2);
        rist_4 = rotazioneInversaEulero(SdR, vn, rist_4);
        rist_6 = rotazioneInversaEulero(SdR, vn, rist_6);
        
        %Aggiunta contributi parziali di ciascun nucleo
        ris_1 = ris_1 + rist_1 * GHweights{indGH}; %INTEGRALE di 1/r (onda P)    
        ris_2 = ris_2 + rist_2 * GHweights{indGH}; %INTEGRALE di r_i*r_j/r^3  (onda P)  
        ris_3 = ris_3 + rist_3 * GHweights{indGH}; %INTEGRALE di 1/r^3 (onda P)
        ris_4 = ris_4 + rist_4 * GHweights{indGH}; %INTEGRALE di r_i*r_j/r^5 (onda P)
        ris_5 = ris_5 + rist_5 * GHweights{indGH}; %INTEGRALE di 1/r (onda P e S)
        ris_6 = ris_6 + rist_6 * GHweights{indGH}; %INTEGRALE di r_i*r_j/r^3 (onda P e S)
    end
    
    %Combinazione risultati dei nuclei con relativi coefficienti
    rist = ((ris_1 - ris_2)/velP^2 ...
                + currT^2 * (-ris_3 + 3*ris_4) ...
                    + ((velP^2 + velS^2) * ris_5 ...
                    + (velP^2 - velS^2) * ris_6) / (velP * velS)^2) / 2; 
    
    %Aggiornamento risultato complessivo
    ris = ris + (coef(var) .* rist); 
end

%Applicazione coefficiente moltiplicativo comune
ris = ris ./ (4*pi*pbParam.rho);

return