function glbIndexFigures = plotMesh(domainMesh, glbIndexFigures)

%Inizializzazione figura
glbIndexFigures = glbIndexFigures + 1;
figure(glbIndexFigures)

% Ciclo sui triangoli della mesh
for i = 1 : domainMesh.numberTriangles
    %Estrazione incidenze vertici triangolo corrente
    incidenze = domainMesh.triangles(i, 1 : 3);
    
    verts = zeros(3, 3);
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

%folderPath = strcat('./outputPlot/', pbParam.domainType);
%mkdir(folderPath);
%name = strcat(folderPath, "/", pbParam.domainType, "_", num2str(pbParam.lev), "_mesh");
%print(name, '-dtiff');

return