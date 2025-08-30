function [nodes, weights] = doppioGauss1D(nNodes)
% INPUT 
%   - nNodes: intero contenente il numero di nodi richiesto
%
% OUTPUT 
%   - nodes: matrice contente le coordinate dei nodi sul
%               triangolo di riferimento
%   - weight: matrice contenente i pesi per ciascun nodo

%% CALCOLO VALORI GAUSS1D
%Numero nodi monodimensionali
nNodi1D = sqrt(nNodes);
if floor(nNodi1D) ~= nNodi1D
    error("Numero di nodi non valido")
end
[nodiStd, pesiStd] = Gauss1D(1, nNodi1D, 0, 0);

%Trasposizione su (0, 1)
nodiStd = 0.5 + (nodiStd ./ 2);
pesiStd = pesiStd ./2;

%% CALCOLO NODI 2D 
% -> x_ij = ((1-x_j)(1-x_i), x_j(1-x_i), x_i)
% Calcolo prima coordinata
nodes(:, :, 1) = (1 - nodiStd) .* (1 - nodiStd)';
% Calcolo seconda coordinata
nodes(:, :, 2) = nodiStd .* (1 - nodiStd)';
% Calcolo seconda coordinata
nodes(:, :, 3) = ones(1, nNodi1D) .* (nodiStd)';

%% CALCOLO PESI
% -> wij = 2A_Rw_iw_j(1-x_i) 
% -> normalizzo su area 1 ->> wij = 2w_iw_j(1-x_i) 
weights = 2 .* pesiStd .* (pesiStd' .* (1 - nodiStd'));


%% RESHAPE in linea
nodes = reshape(nodes, [nNodes , 3]);
weights = reshape(weights, [nNodes , 1]);

end

