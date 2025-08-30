function ris = BEMenerg_calcIntgrGH(pb_param, TS, areaS, indS_RHS, hk, dt, gha, ghw)
%INPUT
% - ...
% - ...
% OUTPUT
% - ...

%Inizializzazione variabile risultato
ris = zeros(3, 1);

%Inizializzazione numero nodi Gauss-Hammer
ngh = pb_param.ngh;

%Calcolo istanti temporale iniziale e finale
t_inz = hk * dt;
t_fin = (hk+1) * dt;

for indGH = 1 : ngh
    %Estrazione COORDINATE NODO GAUSS-HAMMER sul triangolo di riferimento
    ghR(1) = gha(ngh, indGH, 1);
    ghR(2) = gha(ngh, indGH, 2);
    ghR(3) = gha(ngh, indGH, 3);

    %Mappaggio nodi sul triangolo sorgente
    ghS = ghR * TS;

    %Calcolo termine corrispondente al nodo corrente
    ris_t = BEMenerg_calcNucleoTn(pb_param, ghS, t_inz, t_fin, indS_RHS);
    
    %Aggiornamento risultato complessivo
    ris = ris + ris_t * ghw(ngh, indGH);
end

%Applicazione coefficiente pesi di G-H da rifermiento a sorgente
ris = ris * 2 * areaS;
return