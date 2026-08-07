function glbIndexFigures = plotLinear(basePath, pbParam, domainMesh, values, glbIndexFigures)
%PLOTLINEAR  (Deprecated) plot a piecewise-linear (per-vertex) solution snapshot.
%   GLBINDEXFIGURES = PLOTLINEAR(BASEPATH, PBPARAM, DOMAINMESH, VALUES,
%   GLBINDEXFIGURES) manually builds per-triangle patch coordinates from vertex-based
%   VALUES and plots a single snapshot of a piecewise- linear solution field.
%
%   Input arguments:
%       BASEPATH        - (string) project root.
%       PBPARAM         - (struct) physical/problem parameters, see READINPUTFILE.
%       DOMAINMESH      - (struct) triangulated boundary mesh, see READSPACEMESH.
%       VALUES          - (double) per-vertex solution values.
%       GLBINDEXFIGURES - (nonnegative integer) running figure-number counter.
%
%   Output arguments:
%       GLBINDEXFIGURES - (nonnegative integer) updated figure counter.
%
%   Notes:
%       DEPRECATED: superseded by PLOTSOLUTIONS. As with the sibling
%       "_deprecated_plotConstant.m", the file name does not match the internal function
%       name ("plotLinear") and starts with an underscore, which is not a valid leading
%       character for a MATLAB identifier, so this file likely cannot be invoked through
%       the normal eebem.utility.plots.<name>(...) package syntax at all.
%
%   See also PLOTSOLUTIONS

% Linear vectors with the x-y-z coordinates of all the mesh triangles
X = [(domainMesh.coordinates(domainMesh.triangles(:, 1), 1))'; ...
     (domainMesh.coordinates(domainMesh.triangles(:, 2), 1))'; ...
     (domainMesh.coordinates(domainMesh.triangles(:, 3), 1))'];
 
Y = [(domainMesh.coordinates(domainMesh.triangles(:, 1), 2))'; ...
     (domainMesh.coordinates(domainMesh.triangles(:, 2), 2))'; ...
     (domainMesh.coordinates(domainMesh.triangles(:, 3), 2))']; 
 
Z = [(domainMesh.coordinates(domainMesh.triangles(:, 1), 3))'; ...
     (domainMesh.coordinates(domainMesh.triangles(:, 2), 3))'; ...
     (domainMesh.coordinates(domainMesh.triangles(:, 3), 3))']; 

% Get plots specs
[nDim, tVal, jVal, climCustom] = eebem.utility.plots.getPlotProblemParameter(pbParam);

% Map the values onto the mesh triangles for plotting
densityOnTriangles = kron((domainMesh.indSMmatrix > 0), eye(3)) * values / 3;

for i = tVal
    for j = jVal
        % Figure inizialization
        glbIndexFigures = glbIndexFigures + 1;
        fig = figure(glbIndexFigures);

        % Data plot
        if(nDim == 3)
            if(j == 0)
                normDens = sqrt(densityOnTriangles(1:3:end, i).^2 + densityOnTriangles(2:3:end, i).^2 + densityOnTriangles(3:3:end, i).^2);
                fill3(X, Y, Z, normDens, LineStyle="none");
            else
                fill3(X, Y, Z, densityOnTriangles(j:3:end, i), LineStyle="none");
            end
        elseif(nDim == 2)
            fill(X, Y, densityOnTriangles(j:3:end, i));
        end

        % Axis labels
        xlabel('Asse $x_1$', Interpreter = 'latex', FontSize = 14)
        ylabel('Asse $x_2$', Interpreter = 'latex', FontSize = 14)
        if(nDim == 3)
            zlabel('Asse $x_3$', Interpreter = 'latex', FontSize = 14)
        end
    
        % Figure title
        varStrings = dictionary([1 2 3], ["(x_1, T)", "(x_1, x_2, T)", "(x_1, x_2, x_3, T)"]);
        funStrings = dictionary([0 1 2 3], ["|u|", "u_1", "u_2", "u_3"]);
        
        figTitle = "$" + funStrings(j) + varStrings(nDim) + "$ where $T = " + pbParam.deltaT * i + "$";
        
        %Graphics settings
        title(figTitle, Interpreter = 'latex', FontSize = 14)
        colormap(jet)
        colorbar
        daspect([1 1 1]);
        if climCustom
            if(j == 0)
                clim(max(densityOnTriangles, [], "all"));  %TODO: check and fix
            else
                clim([min(densityOnTriangles(j:3:end, :), [], "all"), -min(densityOnTriangles(j:3:end, :), [], "all")]);
            end
        end
    
        % Export image
        folderPath = fullfile(basePath, "outputPlot", pbParam.domainName);
        if(~exist(folderPath, 'dir'))
            mkdir(folderPath);
        end
        figName = fullfile(folderPath, "L " + pbParam.domainName + "-" + domainMesh.lev + "-" + erase(funStrings(j), "\") + "-t_i=" + i);
        exportgraphics(fig, figName + ".svg", ContentType = "vector")
        exportgraphics(fig, figName + ".jpg", ContentType = "image")
        close(fig);
        disp("done " + figName)
    end
end
end