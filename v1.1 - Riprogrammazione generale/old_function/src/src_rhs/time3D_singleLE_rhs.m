%function [ris,riga] = time3D_singleLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,hk,dt,gha,ghw,indS)
function [ris] = time3D_singleLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,hk,dt,gha,ghw,indS)

%La function time3D_singleLE_rhs() prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- la matrice 3x3 TS in cui ogni riga contiene le coordinate di uno dei 
%  tre vertici del triangolo sorgente  
%
%- il vettore 1x3 vnS contenente le coordinate del versore normale al 
%  triangolo sorgente
%
%- la variabile hk che contiene l'indice del blocco matriciale presente 
%  sulla prima colonna della matrice che dobbiamo calcolare
%
%- la variabile dt che contiene il passo di discretizzazione temporale
%
%- l'array 3D gha contenente i nodi di quadratura di Gauss-Hammer
%
%- la matrice ghw che contiene i pesi della formula di quadratura di
%  Gauss-Hammer
%
%e restituisce in OUTPUT:
%- il vettore 3x1 ris che contiene i coefficienti del termine noto 
%  relativo al triangolo sorgente corrente in corrispondenza dell'istante
%  di tempo hk fissato
%
%-------------------------------------------------------------------------

%INIZIALIZZIAMO il valore dell'integrale
ris = zeros(3,1);

%riga=zeros(3,24); screen
%riga=zeros(3,96); bar1

%Numero di nodi e pesi di quadratura per le formule di Gauss-Hammer
%ngaush = 12;
ngaush =7;
%N.B. Ricordarsi di cambiare il numero di nodi di quadratura anche per 
%quanto riguarda il calcolo della matrice
 
% if(hk>10) 
%     ngaush = 7; 
% end
% 
% if(hk>20) 
%     ngaush = 3; 
% end


%PRIMO BLOCCO del TERMINE NOTO del TIME-MARCHING
if(hk>-1)
    
    %Istante di tempo corrente
    thk=(hk+1)*dt;
    
    %Calcolo del valore dell'integrale
    [rist] = time3D_dswLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,thk-dt,thk,ngaush,gha,ghw,indS);
%     [rist,riga_screen0] = time3D_dswLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,thk-dt,thk,ngaush,gha,ghw,indS);

    %Moltiplico il blocco per il COEFFICIENTE DAVANTI al NUCLEO
    ris = ris+rist;
%      riga=riga+riga_screen0;
end


% if(hk>0)
%        thk = hk*dt; %thk=hk*Delta t
%     
%    [rist,riga_screen0] = time3D_dswLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,thk-dt,thk,ngaush,gha,ghw,indS);    
%     riga=riga-2*riga_screen0;
% end
% 
% if(hk>1)
%     thk = (hk-1)*dt; %thk=(hk-1)*Delta t
%     [rist,riga_screen0] = time3D_dswLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,thk-dt,thk,ngaush,gha,ghw,indS);    
%     riga = riga+riga_screen0; 
% end








% if(hk>0) then
%  thk=hk*dt
%  call time3D_dswLE(rist,TS,vnS,areaS,curlS,TF,vnF,areaF,curlF,thk,c,ngaush)
%  ris=ris+3.d0*rist
% endif
% 
% if(hk>1) then
%  thk=(hk-1)*dt
%  call time3D_dswLE(rist,TS,vnS,areaS,curlS,TF,vnF,areaF,curlF,thk,c,ngaush)
%  ris=ris-3.d0*rist
% endif
% 
% if(hk>2) then
%  thk=(hk-2)*dt
%  call time3D_dswLE(rist,TS,vnS,areaS,curlS,TF,vnF,areaF,curlF,thk,c,ngaush)
%  ris=ris+rist
% endif

return