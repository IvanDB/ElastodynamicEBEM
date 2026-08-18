function [nodes, weights] = doppioGauss1D(nNodes)
%DOPPIOGAUSS1D  Build a tensor-product Gauss quadrature rule on the reference triangle via a Duffy-type collapse.
%   [NODES, WEIGHTS] = DOPPIOGAUSS1D(NNODES) takes a SQRT(NNODES)-point 1D Gauss-Legendre
%   rule (GAUSS1D) on (-1, 1), rescales it to (0, 1), and collapses the resulting 2D
%   tensor-product grid onto the reference triangle with barycentric coordinates
%   [(1 - xj)(1 - xi), xj(1 - xi), xi] and a matching Jacobian weight correction
%   (the Duffy transform commonly used to build triangle rules from 1D rules).
%
%   Input arguments:
%       NNODES - (perfect-square positive integer) total number of nodes/weights requested.
%
%   Output arguments:
%       NODES   - (NNODES x 3 double) barycentric quadrature nodes.
%       WEIGHTS - (NNODES x 1 double) corresponding weights.
%
%   Notes:
%       Not currently called elsewhere in the codebase
%       (an alternative to GENERATEFINALG2DNODES).
%
%   See also GAUSS1D, GAUSSHAMMERCOMPOSITE, GENERATEFINALG2DNODES

arguments (Input)
    nNodes (1, 1) double {mustBeInteger, mustBePositive}
end

arguments (Output)
    nodes   (:, 3) double
    weights (:, 1) double
end

import eebem.utility.quadratureRules.*

%Number of 1D nodes
n1Dnodes = sqrt(nNodes);
assert(floor(n1Dnodes) == n1Dnodes, "Total number of nodes/weights must be a perfect square")

%Generate 1D Gauss-Legendre standard nodes/weights
[stdNodes, stdWeights] = Gauss1D(1, n1Dnodes, 0, 0);

%Map to (0, 1)
stdNodes = 0.5 + (stdNodes ./ 2);
stdWeights = stdWeights ./ 2;

%Duffy nodes transformation x_ij = ((1 - x_j)(1 - x_i), x_j(1 - x_i), x_i)
nodes(:, :, 1) = (1 - stdNodes) .* (1 - stdNodes)';
nodes(:, :, 2) = stdNodes .* (1 - stdNodes)';
nodes(:, :, 3) = ones(1, n1Dnodes) .* (stdNodes)';

%Area-normalized weights wij = 2 w_i w_j(1-x_i) 
weights = 2 .* stdWeights .* (stdWeights' .* (1 - stdNodes'));

%Flattening nodes/weights arrays
nodes = reshape(nodes, [nNodes , 3]);
weights = reshape(weights, [nNodes , 1]);
end