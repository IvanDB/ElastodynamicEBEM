function [xVal, tVal, iVal, numPoints, typePlot] = BEMenerg_postProc_loadPoints(pbParam)
% INPUT
% - 
% - 
% OUTPUT
% - 

typePlot = "";
%Selezione del caso di studio per componente x
switch pbParam.domainType
    case 'screenTest'
        xVal = [];     %Post-processing non presente
        tVal = [];
        iVal = [];
    case 'screenUniform'
        xVal = [];     %Post-processing non presente
        tVal = [];     
        iVal = 3;  
    case 'screenGraded'
        xVal = [];     %Post-processing non presente
        tVal = [];     
        iVal = [];      
    case 'sphereUniform'
        xVal = [];     %Post-processing non presente
        tVal = [];
        iVal = [];    
    case 'sphereNotUniform'
        xVal = [];     %Post-processing non presente
        tVal = [];
        iVal = [];
    case 'barH1'
        xVal = [0, 0, 0.25;
                0, 0, 0.50;
                0, 0, 0.75];
        tVal = linspace(0, pbParam.Tfin, 100 * pbParam.Tfin);
        iVal = 3; 
        typePlot = "u(x, :)";
    case 'barH3'
        xVal = [0.25, -0.25, 2];
        tVal = linspace(0, pbParam.Tfin, 100 * pbParam.Tfin);
        iVal = 3;     
        typePlot = "u(x, :)";
    case 'sphereWave'
        coordLin = linspace(-3, 3, 150);
        [X2D, Y2D] = meshgrid(coordLin, coordLin);
        normSq2D = X2D.^2 + Y2D.^2;
        X1D = X2D(normSq2D > 1);
        Y1D = Y2D(normSq2D > 1);
        xVal = [X1D, Y1D, zeros(length(X1D), 1)];
        tVal = [1.2676, 1.5211, 1.7746, 2.0282, 2.2817, 2.5352, 2.7887, 3.0423];
        iVal = [1, 2];     
        typePlot = "u(:, t)";
    case 'elementoIndustriale'
        coordLinH = linspace(-0.2, 0.2, 250);
        coordLinV = linspace(-0.4, 0.0, 250);
        [Y2D, Z2D] = meshgrid(coordLinH, coordLinV);
        isExt = (abs(Y2D) > 0.17) .* (Z2D < -0.01);
        Y1D = Y2D(~isExt);
        Z1D = Z2D(~isExt);
        xVal = [zeros(length(Y1D), 1), Y1D, Z1D];
        tVal = [0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0045, 0.0050, 0.0055, 0.0065, 0.0070, 0.0075];
        iVal = [1, 2, 3];     
        typePlot = "u(:, t)";
                
    otherwise
        error('Caso ancora non riportato/Errore nei dati')
end

%Calcolo numero di punti complessivo
numPoints = size(xVal, 1) * length(tVal);

return