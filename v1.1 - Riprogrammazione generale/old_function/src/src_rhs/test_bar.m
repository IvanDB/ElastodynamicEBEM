function [ris]= test_bar(pb_param,sp,t0,t1)

%La function test_bar() prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- il vettore 1x3 sp contenente le coordinate del punto sorgente
%
%- la variabile t0 contenente il primo istante temporale in cui valutare il
%  dato al bordo
%
%- la variabile t1 secondo il primo istante temporale in cui valutare il
%  dato al bordo
%
%e restituisce in OUTPUT:
%- il vettore 3x1 ris che contiene le componenti del termine noto ottenuto
%  considerando il dato al bordo assegnato sulle facce della barretta
%
%-------------------------------------------------------------------------

%Inizializziamo a zero il vettore ris che conterrà le componenti del
%termine noto
ris=zeros(3,1);

%Istante di tempo finale
T=pb_param.T_fin; 

%Velocità c_P
velc_P=pb_param.velP;
% %Velocità c_S
% velc_S=pb_param.velS;

%Densità di massa
rho=pb_param.rho;
%rho=1/velc_P^2;

%Pressione p=p0H[t] uniforme sulla faccia superiore
p0=1;

%Altezza della barretta
h_min=0;
h_max=1;
h=abs(h_max-h_min);

%Quota del nodo di Gauss-Hammer
z=sp(3);

%Distanza del nodo di quadratura dalla faccia superiore
d=h_max-z;
%d=abs(h_max-z);

%Definizione della funzione greca che mi permetterà di costruire il dato 
%al bordo
greca=@(x,t,kk,c,L)(p0/(rho*c^2))*(heaviside(c*t-x-2*kk*L).*(c*t-x-2*kk*L)-heaviside(c*t+x-(kk+1)*2*L).*(c*t+x-(kk+1)*2*L));

if (abs(z-h_min)>1.0e-10)
    %Se abs(z-h_min)<=1.0e-10 allora siamo sulla faccia inferiore della 
    %barretta e quindi il dato al bordo non deve essere implementato in 
    %quanto è nullo
    kbar=ceil(velc_P*T/(2*h))-1;
    u=zeros(kbar+1,1);
    for k=0:kbar
        u(k+1,1)=-(-1)^k*greca(d,t0,k,velc_P,h)+(-1)^k*greca(d,t1,k,velc_P,h);
    end
    ris(3)=sum(u,1);
end


return



