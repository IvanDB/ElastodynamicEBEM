function ris = time3D_doubleLE(pb_param,TS,areaS,TF,vnF,hk,dt,gha,ghw)
% function ris = time3D_doubleLE(pb_param,TS,vnS,areaS,curlS,TF,vnF,...
%     areaF,curlF,hk,dt,gha,ghw)

%La function time3D_doubleLE() prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- la matrice 3x3 TS in cui ogni riga contiene le coordinate di uno dei 
%  tre vertici del triangolo sorgente  
%
%- la matrice 3x3 TF in cui ogni riga contiene le coordinate di uno dei 
%  tre vertici del triangolo di campo  
%
%- il vettore 1x3 vnF contenente le coordinate del versore normale al 
%  triangolo di campo
%
%- la variabile hk che contiene l'indice del blocco matriciale presente 
%  sulla prima colonna della matrice di Toeplitz che dobbiamo calcolare
%
%- la variabile dt che contiene il passo di discretizzazione temporale
%
%- l'array 3D gha contenente i nodi di quadratura di Gauss-Hammer
%
%- la matrice ghw che contiene i pesi della formula di quadratura di
%  Gauss-Hammer
%
%e restituisce in OUTPUT:
%- la matrice 3x3 ris contenente il risultato della doppia integrazione 
%  in spazio relativa al triangolo sorgente corrente (integrazione esterna)
%  e al triangolo sorgente di campo (integrazione interna analitica)
%  in corrispondenza della differenza temporale tkh fissata
%
%-------------------------------------------------------------------------

%INIZIALIZZIAMO il valore dell'integrale
ris = zeros(3,3);

%Numero di nodi e pesi di quadratura per le formule di Gauss-Hammer
%ngaush = 12;
ngaush = 3;
%N.B. Ricordarsi di cambiare il numero di nodi di quadratura anche per 
%quanto riguarda il calcolo del termine noto

% if(hk>10) 
%     ngaush = 7; 
% end
% 
% if(hk>20) 
%     ngaush = 3; 
% end

if(hk>-1)
    %Caso hk=0,1,....,N_(Delta t)-1
    
    %Istante di tempo in cui valutare il nucleo 
    thk=(hk+1)*dt; %thk=(hk+1)*Delta t
    
    %Calcolo del valore dell'integrale
    rist = time3D_dswLE(pb_param,TS,areaS,TF,vnF,thk,ngaush,gha,ghw);
    
    %Aggiorno il contributo locale alla matrice del sistema
    ris = ris+rist; 
end

if(hk>0)
    %Caso hk=1,2,....,N_(Delta t)-1
    
    %Istante di tempo in cui valutare il nucleo
    thk = hk*dt; %thk=hk*Delta t
    
    %Calcolo del valore dell'integrale
    rist = time3D_dswLE(pb_param,TS,areaS,TF,vnF,thk,ngaush,gha,ghw);
    
    %Aggiorno il contributo locale alla matrice del sistema
    ris=ris-2*rist;
end

if(hk>1)
    %Caso hk=2,3,....,N_(Delta t)-1
    
    %Istante di tempo in cui valutare il nucleo
    thk = (hk-1)*dt; %thk=(hk-1)*Delta t
    
    %Calcolo del valore dell'integrale
    rist = time3D_dswLE(pb_param,TS,areaS,TF,vnF,thk,ngaush,gha,ghw);
    
    %Aggiorno il contributo locale alla matrice del sistema
    ris = ris+rist; 
end

%Divido per il coefficiente del nucleo
ris = ris/(4*pi*pb_param.rho);

return