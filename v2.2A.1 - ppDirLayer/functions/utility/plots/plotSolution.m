function glbIndexFigures = plotSolution(pbParam, u, xVal, tVal, iVal, typePlot, glbIndexFigures)
%%

%%
%Calcolo dimensioni ed inizializzazione matrice soluzione
switch typePlot
    case "u(x, :)"
        numPlotGroups = size(xVal, 1);
        numValPerPlot = length(tVal);
        uSP = u;

        varLib = "t";
        plot3D = false;
        plotCM = false;

        xLabelText = "Asse $t$";
        yLabelText = "$u(x, t)$";

    case "u(:, t)"
        numPlotGroups = length(tVal);
        numValPerPlot = size(xVal, 1);
        uSP = u';

        varLib = "x";
        plot3D = true;
        plotCM = true;

        xLabelText = "$x_1$";
        yLabelText = "$x_2$";
        zLabelText = "$u$";
        
end

extLimColor = max(abs(cell2mat(u)), [], "all");


for indPG = 1 : numPlotGroups
    for indI = iVal
        %Inizializzazione figura
        glbIndexFigures = glbIndexFigures + 1;
        figure(glbIndexFigures);

        funcValue = zeros(numValPerPlot, 1);
        %Lettura valori funzione
        for indVal = 1 : numValPerPlot
            funcValue(indVal) = uSP{indPG, indVal}(indI);
        end

        %Plot funzione
        if varLib == "t" 
            plot(tVal, funcValue);
        end
        if varLib == "x"
            scatter3(xVal(:, 1), xVal(:, 2), xVal(:, 3), 15, funcValue, ".");
        end

        %Set parametri figura
        xlabel(xLabelText, 'interpreter', 'latex', 'fontsize', 14)
        ylabel(yLabelText, 'interpreter', 'latex', 'fontsize', 14)
        if plot3D
            zlim([-0.1 0.1])
            zlabel(zLabelText, 'interpreter', 'latex', 'fontsize', 14)
        end
        
        if plotCM
            titleText = "$u_" + num2str(indI) + "(: , " + num2str(tVal(indPG)) + ")$";
            title(titleText, 'interpreter', 'latex', 'fontsize', 14)
            colormap jet
            c = colorbar;
            c.Location = 'westoutside';
        end

        axis equal
        grid off
        box off

        %salvataggio figura come file immagine (.tiff) e come file .fig
        mkdir(strcat('./outputPlot/', pbParam.domainType));
        path = strcat('./outputPlot/', pbParam.domainType, "/", pbParam.domainType, '_', num2str(pbParam.lev), '_', "u_" + num2str(indI), "punto_", num2str(indPG));
        print(path, '-djpeg')
    end
end

return