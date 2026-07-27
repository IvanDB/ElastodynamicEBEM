function glbIndexFigures = plotSolutions(form, pbParam, domainMesh, density, glbIndexFigures, glbFlags, basePath)
%PLOTSOLUTIONS undefined
%   undefined

arguments
    form            (1, 1) string {mustBeMember(form, ["ID", "DD", "DN", "DNc", "IN"])}
    pbParam         (1, 1) struct
    domainMesh      (1, 1) struct
    density         (:, :) double
    glbIndexFigures (1, 1) double {mustBeInteger, mustBeNonnegative}
    glbFlags        (1, 1) struct
    basePath        (1 ,1) string = "."
end

assert((size(density, 1) == 3*domainMesh.numTriangles) || (size(density, 1) == 3*domainMesh.numVertices), "Value matrix as an invalid size.")

%Early return if both flags are false
if(~glbFlags.plotFigs && ~glbFlags.saveFigs)
    return
end

%Get plots specs
[nDim, tVal, jVal, climCustom] = eebem.utility.plots.getPlotProblemParameter(pbParam);

%Map on triangle if the values are related to the vertices
if(size(density, 1) == 3*domainMesh.numVertices)
    density = kron((domainMesh.indSMmatrix > 0), eye(3)) * density / 3;
end

funChar = [repelem("\phi", ismember(form, "ID")),    ...
           repelem("p", ismember(form, "DD")),       ...
           repelem("u", ismember(form, ["DN", "DNc"]))];

funStrings = dictionary([0 1 2 3], ["|" + funChar + "|", funChar + "_" + [1 2 3]]);
varStrings = dictionary([1 2 3], "(" + ["x_1", "x_1, x_2", "x_1, x_2, x_3"] + ", t)");

%Loop in time instants and vector components
for i = tVal
    for j = jVal
        %Figure initialization
        glbIndexFigures = glbIndexFigures + 1;
        fig = figure(glbIndexFigures);
        fig.Visible = glbFlags.plotFigs;
        
        %Select correct matrix slice
        if(j == 0)
            plotValues = hypot(hypot(density(1:3:end, i), density(2:3:end, i)), density(3:3:end, i));
        else
            plotValues = density(j:3:end, i);
        end

        %Data plot
        trisurf(domainMesh.triangles(:, 1 : 3), domainMesh.coordinates(:, 1), domainMesh.coordinates(:, 2), domainMesh.coordinates(:, 3), plotValues, ...
                 LineStyle = "none");

        %Figure view (default if 3D, planar if 2D)
        view(nDim)

        % Axis labels
        xlabel("$x_1$", Interpreter = 'latex', FontSize = 14, Rotation = 0)
        ylabel("$x_2$", Interpreter = 'latex', FontSize = 14, Rotation = 0)
        if(nDim == 3)
            zlabel("$x_3$", Interpreter = 'latex', FontSize = 14, Rotation = 0)
        end
        ax = gca;
        ax.TickLabelInterpreter = "latex";
    
        % Figure title
        
        figTitle = "$" + funStrings(j) + varStrings(nDim) + "$ with $t = " + pbParam.deltaT * i + "$";
        
        %Graphics settings
        title(figTitle, Interpreter = 'latex', FontSize = 14)
        colormap(jet)
        colorbar
        daspect([1 1 1]);
        if climCustom
            warning("Custom color limit are still WIP")
        end     

        %Export image
        if(glbFlags.saveFigs)
            folderPath = fullfile(basePath, "outputPlot", pbParam.domainName);
            if(~exist(folderPath, 'dir'))
                mkdir(folderPath);
            end
            
            figName = form + "_" + domainMesh.name + "-" + domainMesh.lev + "_" + pbParam.deltaT * i;

            figPath = fullfile(folderPath, figName + ".jpg");
            if(~exist(figPath, "file"))
                exportgraphics(fig, figPath, ContentType = "image")
            end

            figPath = fullfile(folderPath, figName + ".fig");
            if(~exist(figPath, "file"))
                savefig(fig, figPath)
            end
            
            figPath = fullfile(folderPath, figName + ".svg");
            if(ispc && ~exist(figPath, "file"))
                exportgraphics(fig, figPath, ContentType = "vector")
            end
        end

        %Close figure if not visible
        if(~glbFlags.plotFigs)
            close(fig)
        end
    end
end
end