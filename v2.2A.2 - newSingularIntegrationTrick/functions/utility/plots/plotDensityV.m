function glbIndexFigures = plotDensityV(pbParam, domainMesh, density, glbIndexFigures)
    
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
    
    %Assegnazione parametri in base al problema
    switch pbParam.domainType
        case 'screenTest'
            nDim = 2;
            tVal = pbParam.nT;
            jVal = 1 : 3;
            climCustom = false;
        case 'screenUniform'
            nDim = 2;
            tVal = pbParam.nT;
            jVal = 3;
            climCustom = false;
        case 'screenGraded'
            nDim = 2;
            tVal = pbParam.nT;
            jVal = 3;
            climCustom = false;
        case 'sphereUniform'
            nDim = 3;
            tVal = [pbParam.nT .* (1/3), pbParam.nT * (2/3), (pbParam.nT * (2/3)) + 1, pbParam.nT];
            jVal = 3;
            climCustom = false;
        case 'sphereNotUniform'
            nDim = 3;
            tVal = [pbParam.nT .* (1/3), pbParam.nT * (2/3), (pbParam.nT * (2/3)) + 1, pbParam.nT];
            jVal = 3;
            climCustom = false;
        case 'barH1'
            nDim = 3;
            tVal = pbParam.nT .* (1 : pbParam.Tfin) ./ pbParam.Tfin ;
            jVal = 1 : 3;
            climCustom = false;
        case 'barH3'
            nDim = 3;
            tVal = pbParam.nT .* (1 : pbParam.Tfin) ./ pbParam.Tfin ;
            jVal = 1 : 3;
            climCustom = false;
        case 'sphereWave'
            return
        case 'DesCop-cube'
            nDim = 3;
            tVal = floor(pbParam.nT .* [7.50, 8.25, 9.00, 9.75, 10.50, 11.25] ./ pbParam.Tfin);
            jVal = 0;
            climCustom = true;
        case "DesCop-sphere"
            nDim = 3;
            tVal = 1 : 2 : pbParam.nT;
            jVal = [1, 3];
            climCustom = true;
        case 'elementoIndustriale'
            nDim = 3;
            tVal = 1;
            jVal = 3;
            climCustom = false;
    end

    %Plot componenti della funzione densità all'istante temporale di indice hk
    for i = tVal
        for j = jVal
            %Inizializzazione figura
            glbIndexFigures = glbIndexFigures + 1;
            figure(glbIndexFigures)

            %plot dati
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

            %nomi assi
            xlabel('Asse $x_1$', 'interpreter', 'latex', 'fontsize', 14);
            ylabel('Asse $x_2$', 'interpreter', 'latex', 'fontsize', 14);
            if(nDim == 3)
                zlabel('Asse $x_3$', 'interpreter', 'latex', 'fontsize', 14);
            end
        
            %titolo figura
            if(nDim == 3)
                figTitle = strcat('$\phi_', num2str(j), '(x_1, x_2, x_3, T)$ con $ T=', num2str((pbParam.Tfin / pbParam.nT) * i),'$');
            elseif(nDim == 2)
                figTitle = strcat('$\phi_', num2str(j), '(x_1, x_2, T)$ con $ T=', num2str((pbParam.Tfin / pbParam.nT) * i),'$'); 
            end

            %settaggi grafici
            title(figTitle, 'interpreter','latex', 'fontsize', 14);
            colormap(jet)
            colorbar
            daspect([1 1 1]);
            if climCustom
                if(j == 0)
                    clim([0, max(density, [], "all")]);
                else
                    clim([min(density(j:3:end, :), [], "all"), -min(density(j:3:end, :), [], "all")]);
                end
            end
        
            %salvataggio figura come file immagine (.tiff) e come file .fig
            folderPath = strcat('./outputPlot/', pbParam.domainType);
            mkdir(folderPath);
            name = strcat(folderPath, "/", pbParam.domainType, "_", num2str(pbParam.lev), "_phi", num2str(j), "_intervallo_temporale_", num2str(i));
            print(name, '-dtiff');
            %savefig(strcat(name, '.fig'))
            %print(name, '-depsc2', '-r500')
        end
    end
end