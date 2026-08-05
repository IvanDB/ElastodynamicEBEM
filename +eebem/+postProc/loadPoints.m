function [xVal, tVal, iVal, numPoints, typePlot] = loadPoints(pbParam)
% Function loading interest point for post processing
% iVal tells us which component of the solution to plot
typePlot = "";
%Selezione del caso di studio per componente x
switch pbParam.domainType
    case {'barH1','barH1-symm', 'barH1-asym'}
        xVal = [0, 0, 0.25;
                0, 0, 0.50;
                0, 0, 0.75];
        tVal = linspace(0, pbParam.Tfin, 100*pbParam.Tfin);
        iVal = 3;   
        typePlot = "u(x, :)";  
    case 'barH3'
        xVal = [0.25, -0.25, 2];
        tVal = linspace(0, pbParam.Tfin, 100*pbParam.Tfin);
        iVal = 3;     
        typePlot = "u(x, :)";
    case 'DesCop-sphere'
        xVal = [0, 0, 1.003];
        tVal = linspace(0, pbParam.Tfin, pbParam.nT);
        iVal = 3; 
        typePlot = "u(x, :)";
    case {'DesCop-cube-symm', 'DesCop-cube-asym'}
	xVal = [0.3, 0.3, 0.3;
		-0.3, -0.3, -0.3];
	tVal = linspace(0, pbParam.Tfin, 100*pbParam.Tfin);
	iVal = [1, 2, 3];
	typePlot = "u(x, :)";
		
end

% total number of points
numPoints = size(xVal, 1) * length(tVal);

return