function [nodi3D, pesi3D] = generateFinalG2Dnodes(v3D, rMin, rInt, rExt, quad1D)
%GENERATEFINALG2DNODES  Build quadrature nodes on a triangle intersected with a light-cone annulus.
%   [NODI3D, PESI3D] = GENERATEFINALG2DNODES(V3D, RMIN, RINT, REXT, QUAD1D) integrates
%   over the part of the (2D, embedded-in-3D) triangle V3D that lies within radial
%   distance [RMIN, min(RINT/REXT region boundary, ...)] from vertex V3D(3,:), using a
%   polar- coordinate change of variables: an angular 1D Gauss-Legendre rule over the
%   triangle's opening angle, and a further radial 1D Gauss- Legendre rule split at
%   RMIN/RINT (so nodes can be weighted differently across the P-wave/S-wave light-cone
%   boundary implied by RINT and REXT). This realizes, sub-triangle by sub-triangle,
%   the light-cone-intersected singular integration used by CALCSINGSUBBLOCKV/K/K_C/W.
%
%   Input arguments:
%       V3D    - (3x3 double) the three vertices of the (child) triangle, as rows; the
%                third row is the "singular point" the radial coordinate is centered on.
%       RMIN   - (double) inner radius of the integration annulus.
%       RINT   - (double) S-wave radius (velS * currT) for this time instant.
%       REXT   - (double) P-wave radius (velP * currT) for this time instant.
%       QUAD1D - (name-value) numNodes (nonnegative integer, default 0) number of 1D
%                nodes to generate internally via GAUSS1D if G1Dn/G1Dw are not
%                supplied; G1Dn/G1Dw (double, default []) precomputed 1D
%                Gauss-Legendre nodes/ weights to reuse instead of regenerating them.
%
%   Output arguments:
%       NODI3D - (Mx3 double) quadrature nodes in 3D Cartesian
%                coordinates, with zero-weight nodes already removed.
%       PESI3D - (Mx1 double) corresponding quadrature weights.
%
%   Notes:
%       Either QUAD1D.numNodes must be nonzero, or both
%       QUAD1D.G1Dn and QUAD1D.G1Dw must be supplied (asserted).
%
%   See also GAUSS1D, CALCSINGSUBBLOCKV, CALCSINGSUBBLOCKK,
%   CALCSINGSUBBLOCKK_C, CALCSINGSUBBLOCKW

arguments
    v3D     (3, 3) double
    rMin    (1, 1) double
    rInt    (1, 1) double
    rExt    (1, 1) double
    quad1D.numNodes (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
    quad1D.G1Dn     double = []
    quad1D.G1Dw     double = []
end

import eebem.utility.quadratureRules.*

l13 = v3D(1, :) - v3D(3, :);
n13 = norm(l13);
v13 = l13 / n13;
l23 = v3D(2, :) - v3D(3, :);
n23 = norm(l23);
v23 = l23 / n23;
n = cross(v13, v23);
n = n / norm(n);
p = cross(n, l13);
p = p / norm(p);

alpha2D = acos(dot(v13, v23));

if(quad1D.numNodes ~= 0)
    if(~isempty(quad1D.G1Dn) && ~isempty(quad1D.G1Dw))
        assert((quad1D.numNodes == length(quad1D.G1Dw)) && (quad1D.numNodes == length(quad1D.G1Dn)), "Invalid input")
    else
        [quad1D.G1Dn, quad1D.G1Dw] = Gauss1D(1, quad1D.numNodes, 0, 0);
    end
end

assert(~isempty(quad1D.G1Dn) && ~isempty(quad1D.G1Dw) && (length(quad1D.G1Dw) == length(quad1D.G1Dn)), "Invalid input")
quad1D.numNodes = length(quad1D.G1Dw);

%Map to real alpha
nodiT = (quad1D.G1Dn' + 1) ./ 2 .* alpha2D;
pesiT = quad1D.G1Dw' ./ 2 .* alpha2D;

lenSeg = (n13 * n23 * sin(alpha2D)) ./ (n23 .* sin(alpha2D - nodiT) + n13 .* sin(nodiT));
rMax = min(lenSeg, rExt);

linV = sort([rMin*ones(quad1D.numNodes, 1), rInt*ones(quad1D.numNodes, 1), rMax], 2);
rMax = rMax .* ones(1, 3);
flags = (linV > rMax);
linV(flags) = rMax(flags);

intV = linV(:, 1 : 2);
radV = linV(:, 2 : 3) - linV(:, 1 : 2);

nodiR = kron(intV, ones(1, quad1D.numNodes)) + kron(radV, (quad1D.G1Dn + 1) ./ 2);
pesiR = kron(radV, quad1D.G1Dw ./ 2);

nodi2D = nodiR .* reshape([cos(nodiT), sin(nodiT)], quad1D.numNodes, 1, 2);
pesi2D = pesiR .* nodiR .* pesiT;

nodi2D = reshape(nodi2D, [], 2);
pesi2D = reshape(pesi2D, [], 1);

flags = (pesi2D == 0);
nodi2D(flags, :) = [];
pesi2D(flags) = [];
   
nodi3D = v3D(3, :) + nodi2D(:, 1) .* v13 + nodi2D(:, 2) .* p;
pesi3D = pesi2D;
end
