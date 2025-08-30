function [nodes, weights] = doppioGauss1D(nNodes)


%% CALCOLO VALORI GAUSS1D
[nodiStd, pesiStd] = Gauss1D(1, nNodes, 0, 0);

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
nodes(:, :, 3) = ones(1, nNodes) .* (nodiStd)';

%% CALCOLO PESI
% -> wij = 2A_Rw_iw_j(1-x_i) 
weights = sqrt(3) .* pesiStd .* (pesiStd' .* (1 - nodiStd'));

end

