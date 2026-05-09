function subBlockV = calcSingSubBlockV(pbParam, methodSpecs, constValuesCurr, indTemp)
%CALCSINGSUBBLOCKV Summary of this function goes here
%   Detailed explanation goes here

arguments (Input)
    pbParam
    methodSpecs
    constValuesCurr
    indTemp
end

import eebem.utility.quadratureRules.*

subBlockV = zeros(3, 3);
coef = [1, -2, 1];
istTemp = indTemp + [-1, 0, 1];

minDiffTemp = istTemp(1) * pbParam.deltaT;

for indZeta = 1 : 3
    if(istTemp(indZeta) <= 0)
        continue
    end

    intgTot = zeros(3, 3);
    currT = istTemp(indZeta) * pbParam.deltaT;

    for indSRext = 1 : methodSpecs.numSRext
        for indGHext = 1 : methodSpecs.numGHext
            indEXTn = (indSRext-1) * methodSpecs.numGHext + indGHext;
            nodoExt = constValuesCurr.EXTn{indEXTn};
            pesoExt = constValuesCurr.EXTw{indGHext};
        
            intgG2DC = zeros(3, 3);
    
            for indChild = 1 : 3            
                rMin = max(pbParam.velS * minDiffTemp, 0);
                rInt = pbParam.velS * currT;
                rExt = pbParam.velP * currT;
    
                [G2DCnodes, G2DCweights] = generateFinalG2Dnodes(constValuesCurr.childVerts{indEXTn, indChild}, rMin, rInt, rExt, methodSpecs.numSNGLR);
                intgG2DC = intgG2DC + kernelV(length(G2DCweights), G2DCnodes, G2DCweights, nodoExt, currT, pbParam.velP, pbParam.velS);
            end
    
            intgTot = intgTot + pesoExt .* intgG2DC;
        end
    end
    
    subBlockV = subBlockV + coef(indZeta) .* intgTot;
end
subBlockV = subBlockV ./ (4*pi*pbParam.rho);
end