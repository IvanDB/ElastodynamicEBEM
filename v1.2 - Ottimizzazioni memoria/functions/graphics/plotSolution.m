function glb_index_figures = plotSolution(pb_param, sol, x_val, t_val, i_val, glb_index_figures)

%Calcolo dimensioni ed inizializzazione matrice soluzione
num_plot = size(x_val, 1);

for ind_x = 1 : num_plot
    for ind_i = i_val
        %Inizializzazione figura
        glb_index_figures = glb_index_figures + 1;
        figure(glb_index_figures);
        
        %Plot funzione
        plot(t_val, sol(ind_x, :, ind_i));

        %Set parametri figura
        xlabel('Asse $t$', 'interpreter', 'latex', 'fontsize', 14)
        ylabel('$u$', 'interpreter', 'latex', 'fontsize', 14)
        title_fig = "$u_" + num2str(ind_i) + "(" + num2str(x_val(ind_x, 1)) + ", " + num2str(x_val(ind_x, 2)) + ", " + num2str(x_val(ind_x, 3)) + ", t)$";
        title(title_fig, 'interpreter', 'latex');
        %axis equal
        grid on
        box off

        %salvataggio figura come file immagine (.tiff) e come file .fig
        mkdir(strcat('./output/', pb_param.domain_type));
        path = strcat('./output/', pb_param.domain_type, "/", pb_param.domain_type, '_', num2str(pb_param.lev), '_', "u_" + num2str(ind_i), "punto_", num2str(ind_x));
        print(path, '-djpeg')
        savefig(strcat(path, '.fig'))
        print(path, '-depsc2', '-r500')
    end
end
return
end

