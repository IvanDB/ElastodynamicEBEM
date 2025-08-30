function glb_index_figures = plotDensity(pb_param, domainMesh, density, glb_index_figures)

    %Lettura densità da file
    %density = readmatrix(strcat('./output/', pb_param.domain_type, '_', num2str(pb_param.lev), '_', 'soluzione', '.txt'));
    
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
    switch pb_param.domain_type
        case 'screenTest'
            n_dim = 2;
            t_val = pb_param.Nt;
            j_val = 1 : 3;
        case 'screenUniform'
            n_dim = 2;
            t_val = pb_param.Nt;
            j_val = 3;
        case 'screenGraded'
            n_dim = 2;
            t_val = pb_param.Nt;
            j_val = 3;
        case 'sphereUniform'
            n_dim = 3;
            t_val = [pb_param.Nt .* (1/3), pb_param.Nt * (2/3), (pb_param.Nt * (2/3)) + 1, pb_param.Nt];
            j_val = 3;
        case 'sphereNotUniform'
            n_dim = 3;
            t_val = [pb_param.Nt .* (1/3), pb_param.Nt * (2/3), (pb_param.Nt * (2/3)) + 1, pb_param.Nt];
            j_val = 3;
        case 'barH1'
            n_dim = 3;
            t_val = [pb_param.Nt];
            j_val = 1 : 3;
        case 'barH3'
            n_dim = 3;
            t_val = [pb_param.Nt];
            j_val = 1 : 3;
        case 'sphereWave'
            return
    end

    %Plot componenti della funzione densità all'istante temporale di indice hk
    for i = t_val
        for j = j_val
            %Inizializzazione figura
            glb_index_figures = glb_index_figures + 1;
            figure(glb_index_figures)

            %plot dati
            if(n_dim == 3)
                fill3(X, Y, Z, density(j:3:end, i));
            elseif(n_dim == 2)
                fill(X, Y, density(j:3:end, i));
            end

            %nomi assi
            xlabel('Asse $x_1$', 'interpreter', 'latex', 'fontsize', 14);
            ylabel('Asse $x_2$', 'interpreter', 'latex', 'fontsize', 14);
            if(n_dim == 3)
                zlabel('Asse $x_3$', 'interpreter', 'latex', 'fontsize', 14);
            end
        
            %titolo figura
            if(n_dim == 3)
                title_fig = strcat('$\phi_', num2str(j), '(x_1, x_2, x_3, T)$ con $ T=', num2str((pb_param.T_fin / pb_param.Nt) * i),'$');
            elseif(n_dim == 2)
                title_fig = strcat('$\phi_', num2str(j), '(x_1, x_2, T)$ con $ T=', num2str((pb_param.T_fin / pb_param.Nt) * i),'$'); 
            end

            %settaggi grafici
            title(title_fig, 'interpreter','latex', 'fontsize', 14);
            colormap(jet)
            colorbar
            daspect([1 1 1]);
            %clim([0 3]);
        
            %salvataggio figura come file immagine (.tiff) e come file .fig
            mkdir(strcat('./output/', pb_param.domain_type));
            name = strcat('./output/', pb_param.domain_type, "/", pb_param.domain_type, '_', num2str(pb_param.lev), '_phi', num2str(j), '_intervallo_temporale_', num2str(i));
            print(name, '-dtiff')
            savefig(strcat(name, '.fig'))
            print(name, '-depsc2', '-r500')
        end
    end
end