function glbIndexFigures = plotConstant(basePath, pbParam, domainMesh, density, glbIndexFigures)
    
%Calcolo dei vettori conteneti le coordinate dei triangoli della mesh
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

%Plot componenti della funzione densità all'istante temporale di indice hk
for i = tVal
    for j = jVal
        % Figure inizialization
        glbIndexFigures = glbIndexFigures + 1;
        fig = figure(glbIndexFigures);

        % Data plot
        if(nDim == 3)
            if(j == 0)
                normDens = sqrt(density(1:3:end, i).^2 + density(2:3:end, i).^2 + density(3:3:end, i).^2);
                fill3(X, Y, Z, normDens, LineStyle="none");
            else
                fill3(X, Y, Z, density(j:3:end, i), LineStyle="none");
            end
        elseif(nDim == 2)
            fill(X, Y, density(j:3:end, i));
        end

        % Axis labels
        xlabel('Asse $x_1$', Interpreter = 'latex', FontSize = 14)
        ylabel('Asse $x_2$', Interpreter = 'latex', FontSize = 14)
        if(nDim == 3)
            zlabel('Asse $x_3$', Interpreter = 'latex', FontSize = 14)
        end
    
        % Figure title
        varStrings = dictionary([1 2 3], ["(x_1, T)", "(x_1, x_2, T)", "(x_1, x_2, x_3, T)"]);
        funStrings = dictionary([0 1 2 3], ["|p|", "p_1", "p_2", "p_3"]);
        
        figTitle = "$" + funStrings(j) + varStrings(nDim) + "$ where $T = " + pbParam.deltaT * i + "$";
        
        %Graphics settings
        title(figTitle, Interpreter = 'latex', FontSize = 14)
        colormap(jet)
        colorbar
        daspect([1 1 1]);
        if climCustom
            if(j == 0)
                clim(max(density, [], "all"));
            else
                clim([min(density(j:3:end, :), [], "all"), -min(density(j:3:end, :), [], "all")]);
            end
        end
    
        % Export image
        folderPath = fullfile(basePath, "outputPlot", pbParam.domainName);
        if(~exist(folderPath, 'dir'))
            mkdir(folderPath);
        end
        figName = fullfile(folderPath, pbParam.domainName + "-" + domainMesh.lev + "-" + erase(funStrings(j), "\") + "-t_i=" + i);
        exportgraphics(fig, figName + ".svg", ContentType = "vector")
        exportgraphics(fig, figName + ".jpg", ContentType = "image")
    end
end
end