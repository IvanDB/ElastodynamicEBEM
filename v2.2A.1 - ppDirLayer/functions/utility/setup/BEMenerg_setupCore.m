function [methodInfo, EXTn, EXTw, INTn, INTw, DIAGn, DIAGw] = BEMenerg_setupCore(methodSelected)

[methodInfo.typeIntg, ~, ~, nextIdx] = sscanf(methodSelected, "%s", 1);
methodSelected = methodSelected(nextIdx : end);


%Selezione numero nodi interni
switch methodInfo.typeIntg
    case "SA"
        %Nodi esterni
        [methodInfo.numNodiExt, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1); 
        methodSelected = methodSelected(nextIdx : end);
        if ~ismember(methodInfo.numNodiExt, [3, 7, 12, 19])
            error("Valori nodi esterni non validi")
        end
        [EXTn, EXTw] = GaussHammerComposite(1, methodInfo.numNodiExt);

        INTn = [];
        INTw = [];
        DIAGn = [];
        DIAGw = [];

    case "MX.G2D"
        %Nodi esterni
        [methodInfo.numNodiExt, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1);      %Numero nodi esterni
        methodSelected = methodSelected(nextIdx : end);        
        if ~ismember(methodInfo.numNodiExt, [3, 7, 12, 19])
            error("Valori nodi esterni non validi")
        end
        [EXTn, EXTw] = GaussHammerComposite(1, methodInfo.numNodiExt);

        %Nodi interni
        [methodInfo.numNodiInt, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1);      %Numero nodi interni
        methodSelected = methodSelected(nextIdx : end);
        if ~ismember(methodInfo.numNodiInt, [4, 16, 64, 256])
            error("Valori nodi interni non validi")
        end
        [INTn, INTw] = doppioGauss1D(methodInfo.numNodiInt);

        DIAGn = [];
        DIAGw = [];
        
    case {"MX.GHC", "GPU"}
        %Nodi esterni
        [methodInfo.numNodiExt, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1); 
        methodSelected = methodSelected(nextIdx : end);
        if ~ismember(methodInfo.numNodiExt, [3, 7, 12, 19])
            error("Valori nodi esterni non validi")
        end
        [EXTn, EXTw] = GaussHammerComposite(1, methodInfo.numNodiExt);

        %Nodi interni
        [methodInfo.numSubRegion, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1);
        methodSelected = methodSelected(nextIdx : end);
        [methodInfo.numNodiSing, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1);
        methodSelected = methodSelected(nextIdx : end);
        if ~ismember(methodInfo.numSubRegion, [1, 4, 16, 64]) || ~ismember(methodInfo.numNodiSing, [1, 3, 7, 12, 19])
            error("Valori nodi interni non validi")
        end
        methodInfo.numNodiInt = methodInfo.numSubRegion * methodInfo.numNodiSing; 
        [INTn, INTw] = GaussHammerComposite(methodInfo.numSubRegion, methodInfo.numNodiSing);

        DIAGn = [];
        DIAGw = [];

    case "FN"
        %Nodi esterni
        [methodInfo.numNodiExt, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1); 
        methodSelected = methodSelected(nextIdx : end);
        if ~ismember(methodInfo.numNodiExt, [3, 7, 12, 19])
            error("Valori nodi esterni non validi")
        end
        [EXTn, EXTw] = GaussHammerComposite(1, methodInfo.numNodiExt);

        %Nodi interni
        [methodInfo.numSubRegion, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1);
        methodSelected = methodSelected(nextIdx : end);
        [methodInfo.numNodiSing, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1);
        methodSelected = methodSelected(nextIdx : end);
        if ~ismember(methodInfo.numSubRegion, [1, 4, 16, 64]) || ~ismember(methodInfo.numNodiSing, [1, 3, 7, 12, 19])
            error("Valori nodi interni non validi")
        end
        methodInfo.numNodiInt = methodInfo.numSubRegion * methodInfo.numNodiSing; 
        [INTn, INTw] = GaussHammerComposite(methodInfo.numSubRegion, methodInfo.numNodiSing);

        %Nodi interni diagonali
        [methodInfo.numNodiDiag, ~, ~, nextIdx] = sscanf(methodSelected, "%d", 1);      %Numero nodi interni
        methodSelected = methodSelected(nextIdx : end);
        if ~ismember(methodInfo.numNodiDiag, [16, 64, 256, 1024])
            error("Valori nodi interni non validi")
        end
        [DIAGn, DIAGw] = doppioGauss1D(methodInfo.numNodiDiag);

    otherwise
        error("Metodo non esistente")
end

end

