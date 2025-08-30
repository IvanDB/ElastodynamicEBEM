function time3D_plot(domainMesh)
% *****************************************************************
% La function time3D_plot() prende in INPUT:
% - la struct domainMesh contenente le informazioni sulla mesh
% e costruisce il grafico della MESH considerata, oltre a calcolare 
% il minimo e il masismo tra tutte le lunghezze dei lati dei 
% triangoli della mesh
% *****************************************************************

% Inizializziamo a zero la variabile h_max contenente la lunghezza 
% massima tra tutte le lunghezze dei lati dei triangoli dela mesh   
h_max=0;

% Inizializziamo a zero la variabile h_min contenente la lunghezza 
% minima tra tutte le lunghezze dei lati dei triangoli dela mesh   
h_min=0;
figure(2)
for i=1:domainMesh.number_triangles
    % Ciclo sui triangoli della mesh
    
    % Ricaviamo le incidenze dei vertici del triangolo corrente
    incidenze=domainMesh.triangles(i,1:3);
    
    % Ricaviamo le coordinate dei vertici del triangolo corrente
    v1 = domainMesh.coordinates(incidenze(1),:);
    v2 = domainMesh.coordinates(incidenze(2),:);
    v3 = domainMesh.coordinates(incidenze(3),:);      

    verts(1, :) = v1;
    verts(2, :) = v2;
    verts(3, :) = v3;
       
    % Calcoliamo la lunghezza minima e la lunghezza massima tra le 
    % lunghezze dei tre lati del triangolo corrente
    % [h_size_max,h_size_min]=max_edge_size(verts(1:3,:));
    
    % Aggiorniamo il valore delle variabili h_max e h_min    
    % if h_size_max>h_max 
    %     h_max=h_size_max;
    % end
    % if i==1
    %     h_min=h_size_min;
    % else
    %     if h_size_min<h_min 
    %         h_min=h_size_min;
    %     end
    % end
    
    % Plottiamo i vertici e i lati del triangolo corrente
    fill3(verts(:, 1), verts(:, 2), verts(:, 3), 'w')
    hold on
    plot3([v1(1) v2(1)],[v1(2) v2(2)],[v1(3) v2(3)],'b-')
    hold on
    plot3([v1(1) v3(1)],[v1(2) v3(2)],[v1(3) v3(3)],'b-')
    hold on
    plot3([v3(1) v2(1)],[v3(2) v2(2)],[v3(3) v2(3)],'b-')
    hold on
    xlabel('$x$','interpreter','latex','fontsize' ,14)
    ylabel('$y$','interpreter','latex','fontsize' ,14)
    zlabel('$z$','interpreter','latex','fontsize' ,14)
    xlabel('Asse $x$','interpreter','latex','fontsize',13);
    ylabel('Asse $y$','interpreter','latex','fontsize',13);
    zlabel('Asse $z$','interpreter','latex','fontsize',13);
    title('Plot della mesh','interpreter','latex','fontsize',13);
    axis equal
    grid on
    box off
end

% Lunghezza massima dei lati dei triangoli della mesh
disp('Massimo tra le lunghezza dei lati triangoli della mesh')
disp(h_max)

disp('Minimo tra le lunghezze dei lati triangoli della mesh')
% Lunghezza minima dei lati dei triangoli della mesh
disp(h_min)

return