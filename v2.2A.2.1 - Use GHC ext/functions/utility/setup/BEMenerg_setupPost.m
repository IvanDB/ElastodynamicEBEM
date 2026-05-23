function [methodInfo, PPn, PPw] = BEMenerg_setupPost(methodInfo, postSelected)

[methodInfo.typePost, ~, ~, nextIdx] = sscanf(postSelected, "%s", 1);
postSelected = postSelected(nextIdx : end);


%Selezione numero nodi interni
switch methodInfo.typePost
    case "A"
        %Nodi esterni
        PPn = [];
        PPw = [];
           
    case {"N.GHC", "GPU"}
        %Nodi interni
        [methodInfo.numSubRegionPP, ~, ~, nextIdx] = sscanf(postSelected, "%d", 1);
        postSelected = postSelected(nextIdx : end);
        [methodInfo.numNodiSingPP, ~, ~, nextIdx] = sscanf(postSelected, "%d", 1);
        postSelected = postSelected(nextIdx : end);
        if ~ismember(methodInfo.numSubRegionPP, [1, 4, 16, 64]) || ~ismember(methodInfo.numNodiSingPP, [1, 3, 7, 12, 19])
            error("Valori nodi non validi")
        end
        methodInfo.numNodiPP = methodInfo.numSubRegionPP * methodInfo.numNodiSingPP; 
        [PPn, PPw] = GaussHammerComposite(methodInfo.numSubRegionPP, methodInfo.numNodiSingPP);

    otherwise
        error("Metodo non esistente")
end

end

