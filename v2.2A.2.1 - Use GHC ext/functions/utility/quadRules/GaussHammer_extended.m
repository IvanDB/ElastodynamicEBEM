function GHdata = GaussHammer_extended()



%% Lettura dati da file
fileName = "NodiPesiGiusti.txt";
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