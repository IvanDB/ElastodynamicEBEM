function [h_max,h_min] = max_edge_size(verts)
% *****************************************************************
% La function max_edge_size() prende in INPUT:
% - la matrice 3x3 verts in cui l'i-esima riga contiene le 
%   coordinate dell'i-esimo vertice del triangolo corrente 
%   della mesh
% e restituisce in OUTPUT:
% - la variabile h_max che contiene la lunghezza massima tra le 
%   lunghezze dei tre lati del triangolo
% - la variabile h_min che contiene la lunghezza minima tra le 
%   lunghezze dei tre lati del triangolo
% *****************************************************************
   
  % Calcoliamo la LUNGHEZZA dei tre lati del triangolo corrente
  l1=norm(verts(1,:)-verts(2,:));
  l2=norm(verts(1,:)-verts(3,:));
  l3=norm(verts(2,:)-verts(3,:));
  
  % Calcolo della LUNGHEZZA MASSIMA tra le lunghezze dei lati 
  % del triangolo
  h_max=max(max(l1,l2),l3);
  
  % Calcolo della LUNGHEZZA MINIMA tra le lunghezze dei lati 
  % del triangolo
  h_min=min(min(l1,l2),l3);

return