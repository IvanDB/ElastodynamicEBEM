function ris = BEMenerg_coreSA_calcMatrixSubBlock(pbParam, methodInfo, TF, vnF, indTemp, deltaT, SdRF, GHnodesS, GHweightsS)
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

%Inizializzazione numero nodi Gauss-Hammer
nGH = methodInfo.numNodiExt;

%% CALCOLO INTEGRALI

%Calcolo indice temporale
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
    rP = pbParam.velP * currT; % = c_P * Delta
    rS = pbParam.velS * currT; % = c_S * Delta

    %Ciclo sui nodi di Gauss Hammer del triangolo sorgente
    for indGH = 1 : nGH    
        %Mappatura nodo sul triangolo sorgente 
        nodoGHcurr = GHnodesS{indGH};

        %Spostamento su PIANO TRIANGOLO di CAMPO con SdR avente ORIGINE nella
        % PROIEZIONE del NODO corrente e ASSI CARTESIANI salvati in SdR
        [zeta, children, dMIN, c, sign_prod] = preparazioneTriangolo(nodoGHcurr, TF, vnF, SdRF);
        
        %Calcolo rho_P e rho_S
        if(zeta == 0)
            rhoP = rP;
            rhoS = rS;
        elseif(zeta < rS)
            rhoP = sqrt(rP^2 - zeta^2);
            rhoS = sqrt(rS^2 - zeta^2);
        elseif(rS < zeta && zeta < rP)
            rhoP = sqrt(rP^2 - zeta^2);
            rhoS = -1;
        else
            rhoP = -1;
            rhoS = -1;
        end
        
        %Suddivisione casi
        if (dMIN < rhoS)
            typeIntg = "SP_SP";
        elseif (0 < rhoS && rhoS <= dMIN && dMIN < rhoP)
            typeIntg = "SP_P";
        elseif (rhoS <= 0 && dMIN < rhoP)
            typeIntg = "P";
        else
            continue
        end

        %Calcolo integrazioni interne sul triangolo di campo
        [rist_1, rist_2, rist_3, rist_4, rist_5, rist_6] = BEMenerg_coreSA_calcIntgInt_fieldTriangle(zeta, rP, rhoP, rS, rhoS, children, c, sign_prod, typeIntg);
        
        %ROTAZIONE finale all'indietro    
        rist_2 = rotazioneInversaEulero(SdRF, vnF, rist_2);
        rist_4 = rotazioneInversaEulero(SdRF, vnF, rist_4);
        rist_6 = rotazioneInversaEulero(SdRF, vnF, rist_6);
        
        %Aggiunta contributi parziali di ciascun nucleo
        ris_1 = ris_1 + rist_1 * GHweightsS{indGH}; %INTEGRALE di 1/r (onda P)    
        ris_2 = ris_2 + rist_2 * GHweightsS{indGH}; %INTEGRALE di r_i*r_j/r^3  (onda P)  
        ris_3 = ris_3 + rist_3 * GHweightsS{indGH}; %INTEGRALE di 1/r^3 (onda P)
        ris_4 = ris_4 + rist_4 * GHweightsS{indGH}; %INTEGRALE di r_i*r_j/r^5 (onda P)
        ris_5 = ris_5 + rist_5 * GHweightsS{indGH}; %INTEGRALE di 1/r (onda P e S)
        ris_6 = ris_6 + rist_6 * GHweightsS{indGH}; %INTEGRALE di r_i*r_j/r^3 (onda P e S)
    end
    
    %Combinazione risultati dei nuclei con relativi coefficienti
    rist = ((ris_1 - ris_2)/pbParam.velP^2 ...
                + currT^2 * (-ris_3 + 3*ris_4) ...
                    + ((pbParam.velP^2 + pbParam.velS^2) * ris_5 ...
                    + (pbParam.velP^2 - pbParam.velS^2) * ris_6) / (pbParam.velP * pbParam.velS)^2) / 2; 
    
    %Aggiornamento risultato complessivo
    ris = ris + (coef(var) .* rist); 
end

%Applicazione coefficiente moltiplicativo comune
ris = ris ./ (4*pi*pbParam.rho);

return