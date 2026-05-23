function KBlock = BEMenerg_dir_calcKBlock_uConstSpace(methodInfo, indTemp, deltaT, pbParam, domainMesh, constValues, DIAGn, DIAGw)

%Estrazioni variabili dei metodi
numT = domainMesh.numberTriangles;

%Inizializzazione matrice di celle del risultato
KBlock = cell(numT, numT);
sizeMatrix = size(KBlock);


%Calcolo istanti temporali e relativi coefficienti
istTemp = indTemp + [-2, -1, 0, 1];
coeffTemp = [-1, 3, -3, 1];

%Doppio ciclo unificato sugli elementi 
for ind = 1 : numT * numT
    %Inizializzazione valore del singolo blocco 3x3
    KBlock{ind} = zeros(3);

    %Calcolo indici di riga e colonna relativi all'indice lineare ind
    [indMtilde, indM] = ind2sub(sizeMatrix, ind);

    %Ciclo sugli istanti temporali
    for indZeta = 1 : 4
        %Check controllo temporale
        if (istTemp(indZeta) <= 0)
            continue
        end
                            
        %Calcolo paraemetro temporale       
        diffTemp = istTemp(indZeta) * deltaT;

        %Inizializzazione valore del singolo blocco 3x3
        intgTotCompl = zeros(3);
        
        %Estrazione valori triangolo sorgente
        constValuesS = constValues{indMtilde};
         
        %Estrazione dati associati del triangolo di campo corrente
        normInt = domainMesh.normal(indM, :);
        constValuesF = constValues{indM};
        
        %Inizializzazione variabile integrale corrente
        intgTotSing = zeros(3);

        %Divisione caso singolare dal caso non singolare
        if indM ~= indMtilde
            intgTotSing = BEMenerg_dir_calcKIntgNotSing_uConstSpace(methodInfo, pbParam, diffTemp, constValuesS, constValuesF, normInt, domainMesh.coordinates(domainMesh.triangles(indM, 1:3), :));
        end
        if indM == indMtilde
            intgTotSing = BEMenerg_dir_calcKIntgSing_uConstSpace(methodInfo, pbParam, diffTemp, constValuesS, constValuesF, normInt, DIAGn, DIAGw, domainMesh.coordinates(domainMesh.triangles(indM, 1:3), :));
        end
        intgTotCompl = intgTotCompl + intgTotSing;

        KBlock{ind} = KBlock{ind} + coeffTemp(indZeta) .* intgTotCompl;
    end
    KBlock{ind} = KBlock{ind} ./ (4 * pi * deltaT) .* (abs(KBlock{ind}) > 10^(-14)); %
end

KBlock = cell2mat(KBlock);
end