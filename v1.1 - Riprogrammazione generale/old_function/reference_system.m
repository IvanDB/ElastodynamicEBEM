function Sys_rif=reference_system(T3,vnF)

%La function Sistema_rif() prende in INPUT:
%- la matrice 3x3 T3 in cui l'i-esima riga contiene le 3 componenti
%  dell'i-esimo vertice del triangolo di campo corrente
%
%- il vettore vnF contenente le coordinate del versore normale al 
%  triangolo di campo corrente
%
%e restituisce in OUTPUT:
%- la struct Sys_rif che contiene le informazioni utili a individuare 
%  gli assi del nuovo sistema di riferimento bidimensionale nel piano del
%  triangolo di campo corrente
%
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%VETTORE P_2-P_1 (PRIMO LATO P1P2 del TRIANGOLO di CAMPO) 
vL1 = T3(2,:)-T3(1,:);  
%LUNGHEZZA del VETTORE P_2-P_1
Sys_rif.L1 = sqrt(sum(vL1.^2));
%VERSORE e_1 PARALLELO al VETTORE P_2-P_1
Sys_rif.e1 = vL1/Sys_rif.L1;
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%VETTORE P_3-P_1 (TERZO LATO P1P3 del TRIANGOLO di CAMPO)
vL3 = T3(3,:)-T3(1,:); 
%LUNGHEZZA del VETTORE P_3-P_1
Sys_rif.L3 = sqrt(sum(vL3.^2));
%VERSORE e_3 PARALLELO al VETTORE P_3-P_1
Sys_rif.e3 = vL3/Sys_rif.L3;
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%COSENO dell'ANGOLO compreso TRA e_1 e e_3
Sys_rif.cos_e1_e3 = Sys_rif.e3*Sys_rif.e1';
%ANGOLO compreso TRA i VERSORI e_1 e e_3
Sys_rif.delta = acos(Sys_rif.cos_e1_e3);
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%VERSORE e_2 (che è la direzione perpendicolare al versore e_1 nel piano 
%del triangolo di campo) è dato dal PRODOTTO VETTORIALE tra il VERSORE vn 
%NORMALE al TRIANGOLO di CAMPO e il VERSORE e_1 
Sys_rif.e2 = cross(vnF,Sys_rif.e1);
Sys_rif.e2 = Sys_rif.e2/norm(Sys_rif.e2,2);
%----------------------------------------------------------------------

return

