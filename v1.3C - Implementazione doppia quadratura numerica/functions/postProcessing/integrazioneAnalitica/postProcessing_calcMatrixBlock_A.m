function ris = postProcessing_calcMatrixBlock_A(pb_param, TF, vnF, sp, DeltaT)
%INPUT
% - 
% -
%
%OUTPUT
% -

%% Inizializzazione variabili

%Inizializzazione matrici contenenti RISULTATI integrazione 
% di ciascun tipo di nucleo
ris_1 = zeros(3, 3); %Integrale di 1/r (Onda P)
ris_2 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onda P)
ris_3 = zeros(3, 3); %Integrale di 1/r^3 (Onda P)
ris_4 = zeros(3, 3); %Integrale di r_i*r_j/r^5 (Onda P)
ris_5 = zeros(3, 3); %Integrale di 1/r (Onde P e S)
ris_6 = zeros(3, 3); %Integrale di r_i*r_j/r^3 (Onde P e S)

%Calcolo RAGGIO SUPERFICIE SFERICA relativa alle ONDE P
R_P = pb_param.velP * DeltaT; % = c_P * Delta

%Calcolo RAGGIO SUPERFICIE SFERICA relativa alle ONDE S
R_S = pb_param.velS * DeltaT; % = c_S * Delta

%Calcolo nuovo sistema di riferimento
SdR = calcoloSistemaRiferimento(TF, vnF);
    
%% Calcolo integrazione

%Spostamento su PIANO TRIANGOLO di CAMPO con SdR avente ORIGINE nella
% PROIEZIONE del NODO corrente e ASSI CARTESIANI salvati in SdR
[zeta, children, d_MIN, c, sign_prod] = preparazioneTriangolo(sp, TF, vnF, SdR);

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
        [ris_1, ris_2, ris_3, ris_4, ris_5, ris_6] = BEMenerg_calcIntgInt_fieldTriangle(zeta, R_P, rho_P, R_S, rho_S, children, c, sign_prod, "SP_SP");
    elseif (0 < rho_S && rho_S <= d_MIN && d_MIN < rho_P)
        [ris_1, ris_2, ris_3, ris_4, ~, ~] = BEMenerg_calcIntgInt_fieldTriangle(zeta, R_P, rho_P, R_S, rho_S, children, c, sign_prod, "SP_P");
    elseif (rho_S <= 0 && d_MIN < rho_P)
        [ris_1, ris_2, ris_3, ris_4, ~, ~] = BEMenerg_calcIntgInt_fieldTriangle(zeta, R_P, rho_P, R_S, rho_S, children, c, sign_prod, "P");
    end

%ROTAZIONE finale all'indietro    
ris_2 = rotazioneInversaEulero(SdR, vnF, ris_2);
ris_4 = rotazioneInversaEulero(SdR, vnF, ris_4);
ris_6 = rotazioneInversaEulero(SdR, vnF, ris_6);
   
%Combinazione risultati dei nuclei con relativi coefficienti
ris = ((ris_1 - ris_2)/pb_param.velP^2 ...
       + DeltaT^2 * (-ris_3 + 3*ris_4) ...
       + ((pb_param.velP^2 + pb_param.velS^2) * ris_5 ...
       + (pb_param.velP^2 - pb_param.velS^2) * ris_6) / (pb_param.velP * pb_param.velS)^2) / 2; 

%Applicazione coeffiente costante
ris = ris / (4 * pi * pb_param.rho);
return