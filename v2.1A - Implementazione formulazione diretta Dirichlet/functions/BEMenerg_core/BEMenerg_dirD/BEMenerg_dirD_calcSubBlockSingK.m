function KSubBlock = BEMenerg_dirD_calcSubBlockSingK(methodInfo, deltaT, pbParam, domainMesh, DIAGn, DIAGw, indTemp, indM, constValuesCurr, indVsmCurr)


%Calcolo istanti temporali e relativi coefficienti
istTemp = indTemp + [-2, -1, 0, 1];
coeffTemp = [-1, 3, -3, 1];

%Inizializzazione valore blocco 3x3
KSubBlock = zeros(3, 3);

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
    
    matCoeffCurr = constValuesCurr.matCoeff;
        
    %Calcolo del vettore Vms corrente
    vettVMS = cross(normInt, matCoeffCurr(indVsmCurr, :));

    intgTot = BEMenerg_dirD_calcKIntgSing(methodInfo, pbParam, diffTemp, vettVMS, indVsmCurr, constValuesCurr, constValuesCurr, normInt, DIAGn, DIAGw);

    KSubBlock = KSubBlock + coeffTemp(indZeta) .* intgTot;
end

KSubBlock = KSubBlock ./ (4 * pi * deltaT);
end

