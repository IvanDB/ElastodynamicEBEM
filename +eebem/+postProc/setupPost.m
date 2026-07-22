function [postProcInfo, PPn, PPw] = setupPost()
% Setup post processing quadrature nodes and wheigths for GPU implematation
% using Gauss Hammer Composite quadrature scheme with 64 subregions with 3
% nodes each
import eebem.utility.quadratureRules.*
postProcInfo.typePost = "GPU";

%Setting up space nodes info
postProcInfo.numSubRegionPP = 64;

postProcInfo.numNodeSingPP = 3;

if ~ismember(postProcInfo.numSubRegionPP, [1, 4, 16, 64]) || ~ismember(postProcInfo.numNodeSingPP, [1, 3, 7, 12, 19])
    error("Valori nodi non validi")
end
postProcInfo.numNodesPP = postProcInfo.numSubRegionPP * postProcInfo.numNodeSingPP; 
[PPn, PPw] = GaussHammerComposite(postProcInfo.numSubRegionPP, postProcInfo.numNodeSingPP);

end

