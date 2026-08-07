function [nodes, weights] = doppioGauss1D(nNodes)
%DOPPIOGAUSS1D  Build a tensor-product ("double Gauss") quadrature rule on the reference triangle via a Duffy-type collapse.
%   [NODES, WEIGHTS] = DOPPIOGAUSS1D(NNODES) takes a SQRT(NNODES)-point 1D Gauss-Legendre
%   rule (GAUSS1D) on (-1,1), rescales it to (0,1), and collapses the resulting 2D
%   tensor-product grid onto the reference triangle with barycentric coordinates
%   [(1-xj)(1-xi), xj(1-xi), xi], with a matching Jacobian weight correction (the "double
%   Gauss" / Duffy transform commonly used to build triangle rules from 1D rules).
%
%   Input arguments:
%       NNODES - (perfect-square positive integer) total number of nodes/weights requested.
%
%   Output arguments:
%       NODES   - (NNODES x 3 double) barycentric quadrature nodes.
%       WEIGHTS - (NNODES x 1 double) corresponding weights.
%
%   Notes:
%       Errors ("Numero di nodi non valido") if NNODES is not a
%       perfect square. Not currently called elsewhere in the codebase
%       (an alternative to GAUSSHAMMERCOMPOSITE/GAUSSHAMMER_BASE).
%
%   See also GAUSS1D, GAUSSHAMMERCOMPOSITE

import eebem.utility.quadratureRules.*

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

