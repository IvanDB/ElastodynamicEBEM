function KSubBlock = BEMenerg_dir_calcSubBlockSingK_uConstSpace(methodInfo, deltaT, pbParam, domainMesh, DIAGn, DIAGw, indTemp, indM, constValuesCurr, ~)


%Calcolo istanti temporali e relativi coefficienti
istTemp = indTemp + [-2, -1, 0, 1];
coeffTemp = [-1, 3, -3, 1];

%Inizializzazione valore blocco 3x3
KSubBlock = zeros(9, 9);

%Check condizione sul blocco
if(indTemp >= ceil(domainMesh.maxL(indM) / (pbParam.velS * deltaT)) + 2)
    return;
end

%Ciclo sugli istanti temporali
for indZeta = 1 : 4
    %Check controllo temporale
    if (istTemp(indZeta) <= 0)
        continue
    end

    %Calcolo paraemetro temporale       
    diffTemp = istTemp(indZeta) * deltaT;

    %Estrazione dati triangolo corrente
    normInt = domainMesh.normal(indM, :);        

    intgTot = BEMenerg_dir_calcKIntgSing_uConstSpace(methodInfo, pbParam, diffTemp, constValuesCurr, constValuesCurr, normInt, DIAGn, DIAGw, domainMesh.coordinates(domainMesh.triangles(indM, 1:3), :));

    KSubBlock = KSubBlock + coeffTemp(indZeta) .* intgTot;
end

KSubBlock = KSubBlock ./ (4 * pi * deltaT);
end

