function GHdata = GaussHammer_extended()
%GAUSSHAMMER_EXTENDED  (Currently non-callable / flagged buggy) load an extended set of Gauss-Hammer rules from file.
%   GHDATA = GAUSSHAMMER_EXTENDED() reads "NodiPesiGiusti.txt" (expected in this same folder)
%   and parses 47 Gauss-Hammer quadrature rules on the reference triangle from it,
%   each as a struct with the node count, barycentric/Cartesian coordinates and weights,
%   meant as a richer alternative to the five hard-coded rules in GAUSSHAMMER_BASE.
%
%   Output arguments:
%       GHDATA - (47x1 cell) one struct per rule, each with (at least) a numPoints field.
%
%   Notes:
%       NOT CURRENTLY USABLE: the file is named "GaussHammer_extended(bugged).m" --
%       the parentheses and space make it an invalid MATLAB identifier, and a MATLAB
%       function file name must exactly match the function it defines
%       ("GaussHammer_extended") to be callable, so this cannot be invoked
%       through the normal package syntax at all in its current location/name. 
%       Not referenced by any other function in the codebase.
%
%   See also GAUSSHAMMER_BASE, GAUSSHAMMERCOMPOSITE

%Read .txt file
fileName = fullfile(fileparts(mfilename("fullpath")), "NodiPesiGiusti.txt");
[inputFile, errmsg] = fopen(fileName, 'r');
assert(inputFile ~= -1, "Error opening file: " + errmsg)

%Setup number of available quadrature rules
numRules = 47;

%Setup reference triangle data
V1T2D = [-1, -sqrt(3)/3, 0];
V1T3D = [1, 0, 0];

versH = [-1, 1, 0];
versV = [-0.5, -0.5, 1];

%Data array initialization
GHdata = cell(numRules, 1);

for i = 1 : numRules
    %Read header line
    lineValues = sscanf(fgets(inputFile), '%d %d');
    numPoints = lineValues(2);

    GHdata{i}.numPoints = numPoints;

    %Variable initialization
    GHdata{i}.nodes = zeros(numPoints, 3);
    GHdata{i}.weights = zeros(numPoints, 1);

    %Nodes/weights read
    for j = 1 : numPoints
        %Read j-th node data
        lineValues = sscanf(fgets(inputFile), '%f %f %f');

        %Extract 2D coordinates
        GHdata{i}.nodes(j, 1 : 2) = lineValues(1 : 2);

        %2D translation and rescaling 
        GHdata{i}.nodes(j, :) = GHdata{i}.nodes(j, :) - V1T2D;
        GHdata{i}.nodes(j, 1) = GHdata{i}.nodes(j, 1) ./ 2;
        GHdata{i}.nodes(j, 2) = GHdata{i}.nodes(j, 2) ./ sqrt(3);

        %3D projection
        GHdata{i}.nodes(j, :) = V1T3D + (GHdata{i}.nodes(j, 1) * versH) + (GHdata{i}.nodes(j, 2) * versV);

        %Extract and normalize weight value
        GHdata{i}.weights(j) = lineValues(3);
        GHdata{i}.weights(j) = GHdata{i}.weights(j) ./ sqrt(3);
    end
    fgets(inputFile);
end

end