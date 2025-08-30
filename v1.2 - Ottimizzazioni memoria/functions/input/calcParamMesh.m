function domainMesh = calcParamMesh(domainMesh)

%Inizializzazione valori globali
l_min = Inf;
l_max = -Inf;

% Ciclo sui triangoli della mesh
for i = 1 : domainMesh.number_triangles
    %Estrazione incidenze vertici triangolo corrente
    incidenze = domainMesh.triangles(i, 1:3);
    
    %Estrazione coordinate vertici triangolo corrente
    verts(1, :) = domainMesh.coordinates(incidenze(1), :);
    verts(2, :) = domainMesh.coordinates(incidenze(2), :);
    verts(3, :) = domainMesh.coordinates(incidenze(3), :);      

    % Calcolo LUNGHEZZA lati triangolo corrente
	    l1 = norm(verts(1, :) - verts(2, :));
	    l2 = norm(verts(1, :) - verts(3, :));
	    l3 = norm(verts(2, :) - verts(3, :));
  
	    % Calcolo LUNGHEZZA MASSIMA e LUNGHEZZA MINIMA
	    l_max_curr = max(max(l1, l2), l3);
	    l_min_curr = min(min(l1, l2), l3);
   
    %Aggiornamento valori globali
    if l_max_curr > l_max 
  	    l_max = l_max_curr;
    end
    if l_min_curr < l_min 
        l_min = l_min_curr;
    end
end
%Salvataggio valori
domainMesh.l_min = l_min;
domainMesh.l_max = l_max;
return
end