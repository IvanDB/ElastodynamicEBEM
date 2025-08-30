function ris = BEMenerg_calcMatrixSubBlock_mod(pb_param, areaS, TF, vnF, hk, dt, ghw, SdRF, GHnodesS)
%INPUT
% - 
% -
%
%OUTPUT
% - 

%% INIZIALIZZAZIONE VALORI

%Inizializzazione variabile risultato
ris = zeros(3, 3);

%Inizializzazione coefficienti
coef = [1, -2, 1];

%% CALCOLO INTEGRALI

%Calcolo indice temporale
ist_temp = hk + [1, 0, -1];

for var = 1 : 3
    %Controllo necessità di calcolo del sottonucleo
    if(ist_temp(var) > 0)
        %Calcolo parametro temporale del sottonucleo
        curr_t = ist_temp(var) * dt;

        %Calcolo valore integrale per Delta = curr_t
        rist = BEMenerg_calcIntgExt_mod(pb_param, areaS, TF, vnF, curr_t, ghw, SdRF, GHnodesS);
        
        %Aggiornamento risultato complessivo
        ris = ris + (coef(var) .* rist); 
    end
end

%Applicazione coefficiente moltiplicativo comune
ris = ris ./ (4*pi*pb_param.rho);

return