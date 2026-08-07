function GHdata = GaussHammer_extended()
%GAUSSHAMMER_EXTENDED  (Currently non-callable / flagged buggy) load an extended set of Gauss-Hammer rules from file.
%   GHDATA = GAUSSHAMMER_EXTENDED() reads "NodiPesiGiusti.txt" (expected in this same
%   folder) and parses 47 Gauss-Hammer quadrature rules on the reference triangle from it,
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
%       ("GaussHammer_extended") to be callable, so this cannot be invoked through the
%       normal eebem.utility.quadratureRules.GaussHammer_extended(...) package syntax
%       at all in its current location/name. The "(bugged)" suffix is the original
%       author's own flag; this documentation pass did not attempt to identify or fix
%       the underlying bug(s), only to describe the intended behavior from reading the
%       source. Not referenced by any other function in the codebase.
%
%   See also GAUSSHAMMER_BASE, GAUSSHAMMERCOMPOSITE

%% Lettura dati da file
fileName = fullfile(fileparts(mfilename("fullpath")), "NodiPesiGiusti.txt");
inputFile = fopen(fileName, 'r');
if ~inputFile
    error("Impossibile aprire file")
end

%Setup numero quadrature disponibili
numRules = 47;

%Setup variabili caratteristiche di triangoli std
V1T2D = [-1, -sqrt(3)/3, 0];
V1T3D = [1, 0, 0];

versH = [-1, 1, 0];
versV = [-0.5, -0.5, 1];


%Inizializzazione array dati GH
GHdata = cell(numRules, 1);

for i = 1 : numRules
    %Lettura riga intestazione
    lineValues = sscanf(fgets(inputFile), '%d %d');
    numPoints = lineValues(2);                        %Estrazione numero di punti
    
    GHdata{i}.numPoints = numPoints;
    
    %Inizializzazione variabili
    GHdata{i}.nodes = zeros(numPoints, 2);
    GHdata{i}.weights = zeros(numPoints, 1);
    
    %Lettura punti e pesi
    for j = 1 : numPoints
        %Lettura riga del j-esimo nodo
        lineValues = sscanf(fgets(inputFile), '%f %f %f');
        
        %Estrazione coordinate 2D del nodo
        GHdata{i}.nodes(j, :) = lineValues([1, 2]);
        
        %Traslazione e riscalamento delle coordinate 2D
        GHdata{i}.nodes(j, :) = GHdata{i}.nodes(j, :) - V1T2D;
        GHdata{i}.nodes(j, 1) = GHdata{i}.nodes(j, 1) ./ 2;
        GHdata{i}.nodes(j, 2) = GHdata{i}.nodes(j, 2) ./ sqrt(3);
        
        %Mappaggio del nodo sul triangolo STD 3D
        GHdata{i}.nodes(j, :) = V1T3D + (GHdata{i}.nodes(j, 1) * versH) + (GHdata{i}.nodes(j, 2) * versV);
        
        %Estrazione e riscalamento del j-esimo peso
        GHdata{i}.weights(j) = lineValues(3);
        GHdata{i}.weights(j) = GHdata{i}.weights(j) ./ (sqrt(3)/2); %Fattore dato per triangolo con area 1/2
    end
    fgets(inputFile);
end

end