function [nodes, weights] = GaussHammerComposite(nSubPart, nNodGH)
%
%

[nGHstd, wGHstd] = GaussHammer_base(28);
nGHstd = squeeze(nGHstd(nNodGH, 1:nNodGH, :));
wGHstd = squeeze(wGHstd(nNodGH, 1:nNodGH))';

if(nNodGH == 1)
    nGHstd = nGHstd';
end

%Check valore valido
if(sqrt(nSubPart) ~= round(sqrt(nSubPart)))
    error("Valore non valido")
end
nDiv = sqrt(nSubPart);

%Calcolo nodi
nodes = zeros(nSubPart .* nNodGH, 3);

n = 0;
fattLin = 1 / nDiv;
for i = 1 : nDiv
    zMin = (i - 1) .* fattLin;
    zMax = i .* fattLin;

    nTrInRow = 2 * (nDiv - i) + 1;
    
    xVal = [1-zMin, 1-zMin, 1-zMin] - [0, 0, fattLin];
    yVal = [0, 0, 0] - [fattLin, 0, 0];
    for j = 1 : nTrInRow
        if(mod(j, 2) == 1)
            xVal = [xVal(2); xVal(3); xVal(3)];
            yVal = [yVal(2); yVal(3); yVal(3) + fattLin];
            zVal = [zMin; zMax; zMin];
        end
        if(mod(j, 2) == 0)
            xVal = [xVal(2); xVal(3); xVal(3) - fattLin];
            yVal = [yVal(2); yVal(3); yVal(3)];
            zVal = [zMax; zMin; zMax];
        end

        nodes((n .* nNodGH) + (1 : nNodGH), :) = nGHstd * [xVal, yVal, zVal];
        n = n + 1;
    end
end


%Calcolo dei pesi
fatAreaR = sqrt(3); 
weights = (fatAreaR / nSubPart) .* wGHstd;

end

