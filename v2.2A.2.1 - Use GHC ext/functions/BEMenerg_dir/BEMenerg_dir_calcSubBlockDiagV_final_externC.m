function ris = BEMenerg_dir_calcSubBlockDiagV_final_externC(pbParam, methodInfo, indTemp, deltaT, constValuesBlock, DIAGn, DIAGw)
%% INIZIALIZZAZIONE VALORI

%Velocità caratteristiche
cP = pbParam.velP;
cS = pbParam.velS;

%Inizializzazione matrice 
ris = zeros(3, 3);

%Inizializzazione coefficienti
coef = [1, -2, 1];

%Inizializzazione numero nodi 
nGH = methodInfo.numNodiExt;

%% CALCOLO INTEGRALI
%Calcolo indice temporale
istTemp = indTemp + [1, 0, -1];
minDt = istTemp(3) * deltaT; % 0; %

%Ciclo sui tre istanti temporali
for var = 1 : 3
    %Controllo necessità di calcolo del sottonucleo
    if(istTemp(var) <= 0)
        continue
    end
    intgGHG2DC = zeros(3, 3);

    %Calcolo parametro temporale del sottonucleo
    currT = istTemp(var) * deltaT;

    for indGH = 1 : nGH
        nodoExt = constValuesBlock.GHnodes{indGH};
        intgG2DC = zeros(3, 3);

        for indChild = 1 : 3            
            minRad = max(cS * minDt, 0);
            [G2DCnodes, G2DCweights] = generateFinalG2Dnodes_opt(constValuesBlock.DIAGverts{indGH, indChild}, minRad, cS * currT, cP * currT, sqrt(methodInfo.numNodiDiag), sqrt(methodInfo.numNodiDiag));
            intgG2DC = intgG2DC + kernelV_dfor(length(G2DCweights), G2DCnodes, G2DCweights, nodoExt, currT, cP, cS);
        end

        intgGHG2DC = intgGHG2DC + constValuesBlock.GHweights{indGH} .* intgG2DC;

    end
    
    ris = ris + coef(var) .* intgGHG2DC;
end

%Applicazione coefficiente moltiplicativo comune
ris = ris ./ (4*pi*pbParam.rho);

return