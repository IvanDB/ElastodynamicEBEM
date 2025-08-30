%function [ris,riga] = time3D_dswLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,t0,t1,ngaush,gha,ghw,indS)
function [ris] = time3D_dswLE_rhs(pb_param,TS,vnS,areaS,curlS,indS_RHS,t0,t1,ngaush,gha,ghw,indS)

%La function time3D_dswLE_rhs() prende in INPUT:
%- la struct pb_param contenente i parametri del problema
%
%- la matrice 3x3 TS in cui ogni riga contiene le coordinate di uno dei 
%  tre vertici del triangolo sorgente  
%
%- la variabile areaS che contiene l'area del triangolo sorgente corrente
%
%- la variabile indS_RHS che contiene l'indice relativo al tipo di dato al
%  bordo assegnato al triangolo sorgente corrente
%
%- la variabile t0 che contiene il primo istante temporale
%
%- la variabile t1 che contiene il secondo istante temporale
%
%- la variabile ngaush che contiene il numero di nodi e pesi di quadratura 
%  da utilizzare per le formule di Gauss-Hammer
%
%- l'array 3D gha contenente i nodi di quadratura di Gauss-Hammer
%
%- la matrice ghw che contiene i pesi della formula di quadratura di
%  Gauss-Hammer
%
%e restituisce in OUTPUT:
%- il vettore 3x1 ris che contiene il risultato dell'integrazione in 
%  spazio sul triangolo sorgente corrente relativamente al calcolo del 
%  termine noto tra gli istanti t0 e t1
%
%-------------------------------------------------------------------------

%Inizializzazione dell'array 3x1 ris che conterrà il risultato
%dell'integrazione in spazio sul triangolo sorgente corrente per il 
%calcolo del termine noto tra gli istanti t0 e t1
ris = zeros(3,1);

%riga=zeros(3,24); screen
%riga=zeros(3,96); bar1

% fill3([TS(1,1) TS(2,1) TS(3,1)],[TS(1,2) TS(2,2) TS(3,2)],[TS(1,3) TS(2,3) TS(3,3)],'g')
% view(2)
% hold on

%CICLO sul NUMERO di NODI di GAUSS-HAMMER all'interno del triangolo 
%sorgente per eseguire l'INTEGRAZIONE ESTERNA (da eliminare?)
for indGH =1:ngaush
    
    %COORDINATE del NODO di GAUSS-HAMMER corrente nel TRIANGOLO DI RIFERIMENTO
    a(1) = gha(ngaush,indGH,1);
    a(2) = gha(ngaush,indGH,2);
    a(3) = gha(ngaush,indGH,3);
    
    %MAPPIAMO il NODO di GAUSS-HAMMER corrente SUL TRIANGOLO SORGENTE
    sp = a*TS;
    
%     hold on
%     plot3(sp(1),sp(2),sp(3),'bo','markerfacecolor','b','markersize',8)
%     hold on
%     xlabel('$x$','interpreter','latex','fontsize' ,14)
%     ylabel('$y$','interpreter','latex','fontsize' ,14)
%     zlabel('$z$','interpreter','latex','fontsize' ,14)

    rist = time3D_intLE_rhs(1,pb_param,sp,t0,t1, 8);
    %rist = time3D_intLE_rhs(1,pb_param,sp,t0,t1,indS_RHS);

%    thk=t1;
%    riga_screen0=time3D_intLE_screen0(pb_param,sp(1),sp(2),thk);
%    rist = time3D_intLE_rhs_screen_x_0(1,pb_param,sp(2),sp(3),t0,t1,indS_RHS);
%    rist = time3D_intLE_rhs_screen_y_0(1,pb_param,sp(1),sp(3),t0,t1,indS_RHS);
%    riga_screen0=time3D_intLE_bar1(pb_param,sp(1),sp(2),sp(3),thk);
   
    %FORMULA di QUADRATURA di GAUSS-HAMMER per l'INTEGRALE 
    %ESTERNO
    %Moltiplichiamo l'integrale analitico trovato in corrispondenza 
    %del nodo di Gauss-Hammer gha(ngaush,indGH,1:3) per il peso 
    %ghw(ngaush,indGH) e sommiamo il risultato ai risultati ottenuti in
    %precedenza

    ris = ris+rist*ghw(ngaush,indGH);
%     riga=riga+riga_screen0*ghw(ngaush,indGH);
    
end %fine ciclo su indGH (nodi di Gauss-Hammer)

%Mappiamo i pesi di gauss Hammer dal triangolo di riferimento al triangolo
%sorgente
ris=ris*2*areaS;
% riga=riga*2*areaS;

return