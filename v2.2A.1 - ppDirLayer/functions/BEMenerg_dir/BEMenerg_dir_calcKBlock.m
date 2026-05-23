function KBlock = BEMenerg_dir_calcKBlock(methodInfo, indTemp, deltaT, pbParam, domainMesh, constValues, DIAGn, DIAGw)

%Inizializzazione matrice di celle del risultato
KBlock = cell(domainMesh.numberTriangles, domainMesh.number_nodes);
sizeMatrix = size(KBlock);

%Estrazioni variabili dei metodi
numT = domainMesh.numberTriangles;
numS = domainMesh.number_nodes;

%Calcolo istanti temporali e relativi coefficienti
istTemp = indTemp + [-2, -1, 0, 1];
coeffTemp = [-1, 3, -3, 1];

%Doppio ciclo unificato sugli elementi 
for ind = 1 : numT * numS
    %Inizializzazione valore del singolo blocco 3x3
    KBlock{ind} = zeros(3);

    %Calcolo indici di riga e colonna relativi all'indice lineare ind
    [indMtilde, indS] = ind2sub(sizeMatrix, ind);

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
        
        %Ciclo sui possibili triangoli di campo
        for indM = 1 : numT
            %Estrazione indici dei vertici del triangolo 
            indVertsT = domainMesh.triangles(indM, 1 : 3);

            %Estrazione informazioni coppia vertice-triangolo
            [flag, indVsmCurr] = ismember(indS, indVertsT);
            if ~flag
                continue
            end
            
            %Estrazione dati associati del triangolo di campo corrente
            normInt = domainMesh.normal(indM, :);
            constValuesF = constValues{indM};
        
            matCoeffCurr = constValuesF.matCoeff;
        
            %Calcolo del vettore Vms corrente
            vettVMS = cross(normInt, matCoeffCurr(indVsmCurr, :)); %Vms(normInt, matCoeffCurr, indVsmCurr);
        
            % ---> Inserire controlli necessità di calcolo
            
            %Inizializzazione variabile integrale corrente
            intgTotSing = zeros(3);

            %Divisione caso singolare dal caso non singolare
            if indM ~= indMtilde
                intgTotSing = BEMenerg_dir_calcKIntgNotSing(methodInfo, pbParam, diffTemp, vettVMS, indVsmCurr, constValuesS, constValuesF, normInt);
            end
            if indM == indMtilde
                intgTotSing = BEMenerg_dir_calcKIntgSing(methodInfo, pbParam, diffTemp, vettVMS, indVsmCurr, constValuesS, constValuesF, normInt, DIAGn, DIAGw);
            end
            intgTotCompl = intgTotCompl + intgTotSing;
        end
        KBlock{ind} = KBlock{ind} + coeffTemp(indZeta) .* intgTotCompl;
    end
    KBlock{ind} = KBlock{ind} ./ (4 * pi * deltaT); %.* (abs(KBlock{ind}) > 10^(-14))
end

KBlock = cell2mat(KBlock);
end