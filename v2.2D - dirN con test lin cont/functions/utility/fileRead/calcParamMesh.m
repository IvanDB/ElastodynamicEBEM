function domainMesh = calcParamMesh(domainMesh)

%Inizializzazione valori globali
lMin = Inf;
lMax = -Inf;

% Ciclo sui triangoli della mesh
for indT = 1 : domainMesh.numberTriangles
    %Estrazione incidenze vertici triangolo corrente
    incidenze = domainMesh.triangles(indT, 1:3);
    
    %Estrazione coordinate vertici triangolo corrente
    verts(1, :) = domainMesh.coordinates(incidenze(1), :);
    verts(2, :) = domainMesh.coordinates(incidenze(2), :);
    verts(3, :) = domainMesh.coordinates(incidenze(3), :);      

    % Calcolo LUNGHEZZA lati triangolo corrente
    l1 = norm(verts(1, :) - verts(2, :));
    l2 = norm(verts(1, :) - verts(3, :));
    l3 = norm(verts(2, :) - verts(3, :));

    % Calcolo LUNGHEZZA MASSIMA e LUNGHEZZA MINIMA
    lMaxCurr = max(max(l1, l2), l3);
    lMinCurr = min(min(l1, l2), l3);
   
    %Aggiornamento valori globali
    if lMaxCurr > lMax 
  	    lMax = lMaxCurr;
    end
    if lMinCurr < lMin 
        lMin = lMinCurr;
    end
end

%Salvataggio valori
domainMesh.lMin = lMin;
domainMesh.lMax = lMax;
return
end