function ris = BEMenerg_postProcA_calcMatrixBlock(pbParam, constValuesT, sourcePoint, deltaT)
%INPUT
% - 
% -
%
%OUTPUT
% -

%% Inizializzazione variabili
%Inizializzazione matrice 
ris = zeros(3, 3);

%Calcolo RAGGI SUPERFICIE SFERICHE relative alle ONDE P ed S
rP = pbParam.velP * deltaT; % = c_P * Delta
rS = pbParam.velS * deltaT; % = c_S * Delta

%Estrazione informazioni triangolo
SdR = constValuesT.SdR;
vnT = constValuesT.vnT;
vertsT = constValuesT.vertsT;

%% Calcolo integrazione

%Spostamento su PIANO TRIANGOLO di CAMPO con SdR avente ORIGINE nella
% PROIEZIONE del NODO corrente e ASSI CARTESIANI salvati in SdR
[zeta, children, dMIN, c, sign_prod] = preparazioneTriangolo(sourcePoint, vertsT, vnT, SdR);

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
    rhoS = -1;
    rhoP = -1;
end

%Suddivisione casi
if (dMIN < rhoS)
    typeIntg = "SP_SP";
elseif (0 < rhoS && rhoS <= dMIN && dMIN < rhoP)
    typeIntg = "SP_P";
elseif (rhoS <= 0 && dMIN < rhoP)
    typeIntg = "P";
else
    return
end

[ris_1, ris_2, ris_3, ris_4, ris_5, ris_6] = BEMenerg_coreSA_calcIntgInt_fieldTriangle(zeta, rP, rhoP, rS, rhoS, children, c, sign_prod, typeIntg);

%ROTAZIONE finale all'indietro    
ris_2 = rotazioneInversaEulero(SdR, vnT, ris_2);
ris_4 = rotazioneInversaEulero(SdR, vnT, ris_4);
ris_6 = rotazioneInversaEulero(SdR, vnT, ris_6);
   
%Combinazione risultati dei nuclei con relativi coefficienti
ris = ((ris_1 - ris_2)/pbParam.velP^2 ...
       + deltaT^2 * (-ris_3 + 3*ris_4) ...
       + ((pbParam.velP^2 + pbParam.velS^2) * ris_5 ...
       + (pbParam.velP^2 - pbParam.velS^2) * ris_6) / (pbParam.velP * pbParam.velS)^2) / 2; 

%Applicazione coeffiente costante
ris = ris ./ (4 * pi * pbParam.rho);
return