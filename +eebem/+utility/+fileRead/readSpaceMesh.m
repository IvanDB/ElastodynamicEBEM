function domainMesh = readSpaceMesh(basePath, meshFileName)
%READSPACEMESH  Load a triangular surface mesh file and derive all the per-triangle geometric data the BEM code needs.
%   DOMAINMESH = READSPACEMESH(BASEPATH, MESHFILENAME) parses the mesh file
%   BASEPATH/mesh/<domainName>/MESHFILENAME (vertex coordinates then triangle
%   connectivity, in a fixed custom text format), then computes and stores: outward unit
%   normals and areas per triangle (from the cross product of two edges; triangle vertex
%   order is flipped for the "DesCop-sphere" mesh to keep normals outward- pointing on
%   that interior-scattering geometry), centroids, edge lengths and the maximum edge
%   length per triangle (maxL, used by the light-cone cutoffs in CALCSINGSUBBLOCKK/K_C),
%   the local-index-per- vertex incidence matrix (indSMmatrix), and the mesh-wide
%   minimum/ maximum edge length (lMin/lMax, via the local helper CALCPARAMMESH).
%
%   Input arguments:
%       BASEPATH     - (string) project root.
%       MESHFILENAME - (string) file name, see CONSTRUCTMESHFILENAME.
%
%   Output arguments:
%       DOMAINMESH - (struct) with fields name, lev, numVertices, coordinates, numTriangles,
%                    triangles, normal, area, center, maxL, curl, indSMmatrix, lMin, lMax.
%
%   Notes:
%       Asserts if the mesh file cannot be opened. The mesh-wide min/max edge length
%       computation in the local CALCPARAMMESH helper uses a plain (non-vectorized,
%       non-parallel) loop over all triangles and can be a bottleneck on very fine meshes.
%
%   See also CONSTRUCTMESHFILENAME, READINPUTFILE, CALCCONSTVALUES

arguments
    basePath     (1, 1) string
    meshFileName (1, 1) string
end

import eebem.utility.fileRead.*
%Extract mesh name parts
fileNameParts = split(meshFileName, ["_", "."]);
domainName = fileNameParts(1);
meshName = domainName;
if(length(fileNameParts) > 2)
    meshType = fileNameParts(end-2);
    meshName = meshName + "_" + meshType;
end
meshLev = str2double(fileNameParts(end-1));

%Open mesh file
meshFilePath = fullfile(basePath, "mesh", domainName, meshFileName);
meshFile = fopen(meshFilePath, 'r');
assert(meshFile ~= -1, "Error opening mesh file")

domainMesh = struct();
%Set mesh data
domainMesh.name = meshName;      % TODO: decide convention 
domainMesh.lev = meshLev;

%Skip headers
fgets(meshFile);
fgets(meshFile);
fgets(meshFile);
fgets(meshFile);

%Read mesh verteces
fgets(meshFile);

domainMesh.numVertices = sscanf(fgets(meshFile), "%d");
domainMesh.coordinates = zeros(domainMesh.numVertices, 4);

formatSpec = "%f";
sizeSpec = [4 4*domainMesh.numVertices];

domainMesh.coordinates = fscanf(meshFile, formatSpec, sizeSpec)';

domainMesh.coordinates(:, 4) = [];

%Read mesh triangles
fgets(meshFile);

domainMesh.numTriangles = sscanf(fgets(meshFile), "%d");

domainMesh.triangles = zeros(domainMesh.numTriangles, 4);

formatSpec = "%f";
sizeSpec = [4 4*domainMesh.numTriangles];

domainMesh.triangles = fscanf(meshFile, formatSpec, sizeSpec)';

ind = 4;
domainMesh.triangles(:, 4) = ind * ones(domainMesh.numTriangles, 1);

if(domainMesh.name == "DesCop-sphere")   %Move to a intern/extern problem dedicated flag
    domainMesh.triangles(:, [2 3]) = domainMesh.triangles(:, [3 2]);
end

%Close mesh file
fgets(meshFile);
fclose(meshFile);

%% Calculate useful information about the mesh
coordVertTriangles = domainMesh.coordinates(domainMesh.triangles(:, 1:3), :);

%P_iP_j vectors
vg1 = coordVertTriangles(1 : domainMesh.numTriangles, :) - ...
    coordVertTriangles(2*domainMesh.numTriangles+1 : end, :);

vg2 = coordVertTriangles(domainMesh.numTriangles+1 : 2*domainMesh.numTriangles, :) - ...
    coordVertTriangles(2*domainMesh.numTriangles+1 : end, :);

vg3 = coordVertTriangles(1 : domainMesh.numTriangles, :) - ...
    coordVertTriangles(domainMesh.numTriangles+1 : 2*domainMesh.numTriangles, :);

%Normal vectors (n = vg1 x vg2)
domainMesh.normal = zeros(domainMesh.numTriangles, 3);
domainMesh.normal(:, 1) =  vg1(:, 2) .* vg2(:, 3) - vg2(:, 2) .* vg1(:, 3);
domainMesh.normal(:, 2) = -vg1(:, 1) .* vg2(:, 3) + vg2(:, 1) .* vg1(:, 3);
domainMesh.normal(:, 3) =  vg1(:, 1) .* vg2(:, 2) - vg2(:, 1) .* vg1(:, 2);

%Areas (A = |n| / 2)
domainMesh.area = zeros(domainMesh.numTriangles, 1);
domainMesh.area = sqrt(sum(domainMesh.normal.^2, 2)) / 2;

%Normalization of normal vectors (n = n / (2*A))
domainMesh.normal = domainMesh.normal ./ (2 * domainMesh.area);

%Centroids (b = P1 + P2 + P3 / 3)
domainMesh.center = zeros(domainMesh.numTriangles, 3);
domainMesh.center = (coordVertTriangles(1 : domainMesh.numTriangles, :) + ...
    coordVertTriangles(domainMesh.numTriangles + 1 : 2*domainMesh.numTriangles, :) + ...
    coordVertTriangles(2*domainMesh.numTriangles + 1 : 3*domainMesh.numTriangles, :)) / 3;

%Side lengths (l_i = |vg_i|_2)
len = [sqrt(sum(vg3.^2, 2)), sqrt(sum(vg2.^2, 2)), sqrt(sum(vg1.^2, 2))];

%Max side length (maxL = max(l1, l2, l3))
domainMesh.maxL = zeros(domainMesh.numTriangles, 1);
domainMesh.maxL = max(len, [], 2);

%Curl (curl(:, i) = vg_i' / (2*A))
domainMesh.curl = zeros(3, 3, domainMesh.numTriangles);
domainMesh.curl(:, 1, :) =  (vg2 ./ (2*domainMesh.area))';
domainMesh.curl(:, 2, :) = -(vg1 ./ (2*domainMesh.area))';
domainMesh.curl(:, 3, :) =  (vg3 ./ (2*domainMesh.area))';

%Vertex matrix
indSMmatrix = zeros(domainMesh.numTriangles, domainMesh.numVertices, 'int32');
parfor indM = 1 : domainMesh.numTriangles
    [~, indSMmatrix(indM, :)] = ismember(1 : domainMesh.numVertices, domainMesh.triangles(indM, 1 : 3));
end
domainMesh.indSMmatrix = indSMmatrix;

%Other useful information about the mesh
domainMesh = calcParamMesh(domainMesh); % --> TODO: refactor all section
return
end

function domainMesh = calcParamMesh(domainMesh)
    %Compute the mesh-wide minimum/maximum triangle edge length
    %(DOMAINMESH.lMin/lMax) by looping over every triangle.
    %Initialize values
    lMin = Inf;
    lMax = -Inf;

    for indT = 1 : domainMesh.numTriangles

        incidenze = domainMesh.triangles(indT, 1:3);
        verts(1, :) = domainMesh.coordinates(incidenze(1), :);
        verts(2, :) = domainMesh.coordinates(incidenze(2), :);
        verts(3, :) = domainMesh.coordinates(incidenze(3), :);      
    
        l1 = norm(verts(1, :) - verts(2, :));
        l2 = norm(verts(1, :) - verts(3, :));
        l3 = norm(verts(2, :) - verts(3, :));
    
        lMaxCurr = max(max(l1, l2), l3);
        lMinCurr = min(min(l1, l2), l3);
       
        %Aggiornamento valori globali
        lMax = max(lMax, lMaxCurr);
        lMin = min(lMin, lMinCurr);
    end
    
    %Salvataggio valori
    domainMesh.lMin = lMin;
    domainMesh.lMax = lMax;
    return
end