function glbIndexFigures = plotSolution(basePath, pbParam, u, xVal, tVal, iVal, glbIndexFigures, form)
%%

%%
%Calcolo dimensioni ed inizializzazione matrice soluzione
numPlotGroups = size(xVal, 1);
numValPerPlot = length(tVal);
uSP = u;


xLabelText = "$t$ ($s$)";


for indPG = 1 : numPlotGroups
    for indI = iVal
        %Inizializzazione figura
        glbIndexFigures = glbIndexFigures + 1;
        fig = figure(glbIndexFigures);
        hold on;
        funcValue = zeros(numValPerPlot, 1);
        %Lettura valori funzione
        for indVal = 1 : numValPerPlot
            funcValue(indVal) = uSP{indPG, indVal}(indI);
        end

        %Plot funzione 
        plot(tVal, funcValue);
	yLabelText = "$u_{"+num2str(indI)+"}(x_0, t)$";
        %Set parametri figura
        xlabel(xLabelText, 'interpreter', 'latex', 'fontsize', 14)
        ylabel(yLabelText, 'interpreter', 'latex', 'fontsize', 14)

        %axis equal
        grid off
        box off

        folderPath = fullfile(basePath, "solutionPlot", pbParam.domainType);
        if(~exist(folderPath, 'dir'))
            mkdir(folderPath);
        end
        figName = fullfile(folderPath, form + pbParam.domainType + "-" + pbParam.lev + "-point_"+ num2str(indPG) + "-u_"+ num2str(indI)+"_solution");
        %exportgraphics(fig, figName + ".svg", ContentType = "vector")
        exportgraphics(fig, figName + ".jpg", ContentType = "image");
        savefig(fig,figName+".fig");
    end
end

end