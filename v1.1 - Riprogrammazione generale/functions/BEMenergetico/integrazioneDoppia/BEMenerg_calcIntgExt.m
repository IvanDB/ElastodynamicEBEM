function ris = BEMenerg_calcIntgExt(pb_param, TS, areaS, TF, vnF, thk, gha, ghw)
%INPUT
% - 
% -
%
%OUTPUT
% - 

%% Inizializzazione variabili

%Inizializzazione numero nodi Gauss-Hammer
ngh = pb_param.ngh;

%Inizializzazione matrici contenenti RISULTATI COMPLESSIVI integrazione 
% di ciascun tipo di nucleo
ris_1 = zeros(3, 3); %Integrale di 1/r (Onda P)
ris_2 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onda P)
ris_3 = zeros(3, 3); %Integrale di 1/r^3 (Onda P)
ris_4 = zeros(3, 3); %Integrale di r_i*r_j/r^5 (Onda P)
ris_5 = zeros(3, 3); %Integrale di 1/r (Onde P e S)
ris_6 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onde P e S)

%Calcolo RAGGIO SUPERFICIE SFERICA relativa alle ONDE P
R_P = pb_param.velP * thk; % = c_P * Delta

%Calcolo RAGGIO SUPERFICIE SFERICA relativa alle ONDE S
R_S = pb_param.velS * thk; % = c_S * Delta

%Calcolo nuovo sistema di riferimento
SdR = calcoloSistemaRiferimento(TF, vnF);

%% Calcolo integrazione esterna
for indGH = 1 : ngh
    %Inizializzazione matrici contenenti RISULTATI PARZIALI integrazione 
    % di ciascun nucleo
    rist_1 = zeros(3, 3); %Integrale di 1/r (Onda P)
    rist_2 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onda P)
    rist_3 = zeros(3, 3); %Integrale di 1/r^3 (Onda P)
    rist_4 = zeros(3, 3); %Integrale di r_i*r_j/r^5 (Onda P)
    rist_5 = zeros(3, 3); %Integrale di 1/r (Onda P e S)
    rist_6 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onda P e S)
    
    %Estrazione COORDINATE NODO GAUSS-HAMMER corrente nel 
    % TRIANGOLO DI RIFERIMENTO
    nodoGHstd(1) = gha(ngh, indGH, 1);
    nodoGHstd(2) = gha(ngh, indGH, 2);
    nodoGHstd(3) = gha(ngh, indGH, 3);
    
    %Mappatura nodo sul triangolo sorgente 
    nodoGHcurr = nodoGHstd * TS; 
    
    %Spostamento su PIANO TRIANGOLO di CAMPO con SdR avente ORIGINE nella
    % PROIEZIONE del NODO corrente e ASSI CARTESIANI salvati in SdR
    [zeta, children, d_MIN, c, sign_prod] = preparazioneTriangolo(nodoGHcurr, TF, vnF, SdR);
    
    %Calcolo rho_P e rho_S
    if(zeta == 0)
        rho_P = R_P;
        rho_S = R_S;
    elseif(zeta < R_S)
        rho_P = sqrt(R_P^2 - zeta^2);
        rho_S = sqrt(R_S^2 - zeta^2);
    elseif(R_S < zeta && zeta < R_P)
        rho_P = sqrt(R_P^2 - zeta^2);
        rho_S = -1;
    else
        rho_S = -1;
        rho_P = -1;
    end
    
    %Suddivisione casi
    if (d_MIN < rho_S)
        [rist_1, rist_2, rist_3, rist_4, rist_5, rist_6] = BEMenerg_calcIntgInt_fieldTriangle(zeta, R_P, rho_P, R_S, rho_S, children, c, sign_prod, "SP_SP");
    elseif (0 < rho_S && rho_S <= d_MIN && d_MIN < rho_P)
        [rist_1, rist_2, rist_3, rist_4, ~, ~] = BEMenerg_calcIntgInt_fieldTriangle(zeta, R_P, rho_P, R_S, rho_S, children, c, sign_prod, "SP_P");
    elseif (rho_S <= 0 && d_MIN < rho_P)
        [rist_1, rist_2, rist_3, rist_4, ~, ~] = BEMenerg_calcIntgInt_fieldTriangle(zeta, R_P, rho_P, R_S, rho_S, children, c, sign_prod, "P");
    end

    %ROTAZIONE finale all'indietro    
    rist_2 = rotazioneInversaEulero(SdR, vnF, rist_2);
    rist_4 = rotazioneInversaEulero(SdR, vnF, rist_4);
    rist_6 = rotazioneInversaEulero(SdR, vnF, rist_6);
    
    %Aggiunta contributi parziali di ciascun nucleo
    ris_1 = ris_1 + rist_1*ghw(ngh, indGH); %INTEGRALE di 1/r (onda P)    
    ris_2 = ris_2 + rist_2*ghw(ngh, indGH); %INTEGRALE di r_i*r_j/r^3  (onda P)  
    ris_3 = ris_3 + rist_3*ghw(ngh, indGH); %INTEGRALE di 1/r^3 (onda P)
    ris_4 = ris_4 + rist_4*ghw(ngh, indGH); %INTEGRALE di r_i*r_j/r^5 (onda P)
    ris_5 = ris_5 + rist_5*ghw(ngh, indGH); %INTEGRALE di 1/r (onda P e S)
    ris_6 = ris_6 + rist_6*ghw(ngh, indGH); %INTEGRALE di r_i*r_j/r^3 (onda P e S)
end 

%Combinazione risultati dei nuclei con relativi coefficienti
ris = ((ris_1 - ris_2)/pb_param.velP^2 ...
       + thk^2 * (-ris_3 + 3*ris_4) ...
       + ((pb_param.velP^2 + pb_param.velS^2) * ris_5 ...
       + (pb_param.velP^2 - pb_param.velS^2) * ris_6) / (pb_param.velP * pb_param.velS)^2) / 2; 

%Applicazione coefficiente comune di scala per i pesi di GH
ris = ris * 2 * areaS;

return