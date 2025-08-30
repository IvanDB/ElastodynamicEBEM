function glb_index_figures = plotMesh(domainMesh, glb_index_figures)

    %Inizializzazione figura
    glb_index_figures = glb_index_figures + 1;
    figure(glb_index_figures)

    % Ciclo sui triangoli della mesh
    for i = 1 : domainMesh.number_triangles
        %Estrazione incidenze vertici triangolo corrente
        incidenze = domainMesh.triangles(i, 1 : 3);
        
        %Estrazione coordinate vertici triangolo corrente
        verts(1, :) = domainMesh.coordinates(incidenze(1), :);
        verts(2, :) = domainMesh.coordinates(incidenze(2), :);
        verts(3, :) = domainMesh.coordinates(incidenze(3), :);      

        %Plot superficie triangolo corrente
        fill3(verts(:, 1), verts(:, 2), verts(:, 3), 'w')
        hold on

        %Plot lati del triangolo corrente
        plot3([verts(1, 1) verts(2, 1)], [verts(1, 2) verts(2, 2)], [verts(1, 3) verts(2, 3)], 'b-')
        plot3([verts(2, 1) verts(3, 1)], [verts(2, 2) verts(3, 2)], [verts(2, 3) verts(3, 3)], 'b-')
        plot3([verts(3, 1) verts(1, 1)], [verts(3, 2) verts(1, 2)], [verts(3, 3) verts(1, 3)], 'b-')
    end
    %Set parametri figura
    xlabel('Asse $x$', 'interpreter', 'latex', 'fontsize', 14)
    ylabel('Asse $y$', 'interpreter', 'latex', 'fontsize', 14)
    zlabel('Asse $z$', 'interpreter', 'latex', 'fontsize', 14)
    title('Plot della mesh', 'interpreter', 'latex', 'fontsize', 14)
    axis equal
    grid on
    box off
return