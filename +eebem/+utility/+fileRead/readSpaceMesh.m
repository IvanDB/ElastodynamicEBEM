function domainMesh = readSpaceMesh(basePath, domainType, inputStruct)
arguments
    basePath    (1, 1) string
    domainType  (1, 1) string
    inputStruct (1, 1) struct
end

%Open mesh file
meshFilePath = fullfile(basePath, "mesh", domainType, domainType + "_" + inputStruct.meshLevel + ".mesh");
meshFile = fopen(meshFilePath, 'r');
if meshFile == -1
    error("Impossibile aprire file di mesh")
end

%Set values of the mesh structure
domainMesh = struct();
domainMesh.domainType = domainType;
domainMesh.meshType = inputStruct.meshType;
domainMesh.meshLevel = inputStruct.meshLevel;

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

if(domainType == "DesCop-sphere")   %Move to a intern/extern problem dedicated flag
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
domainMesh.indSMmatrix = zeros(domainMesh.numTriangles, domainMesh.numVertices, 'int32');
for indS = 1 : domainMesh.numVertices
    for indM = 1 : domainMesh.numTriangles
        [~, domainMesh.indSMmatrix(indM, indS)] = ismember(indS, domainMesh.triangles(indM, 1 : 3));
    end
end

%Other useful information about the mesh
domainMesh = calcParamMesh(domainMesh); % --> TODO: refactor all section
return
end

function domainMesh = calcParamMesh(domainMesh)
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