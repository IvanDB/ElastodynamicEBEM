function betaI = calcBetaI(pbParam, domainMesh, constValues, methodSpecs)
%CALCBETAISummary of this function goes here
%   Detailed explanation goes here
arguments
    pbParam     struct
    domainMesh  struct
    constValues cell
    methodSpecs struct
end

import eebem.core.*

numT = domainMesh.numTriangles;

gI = getDatumHandleDirichlet(pbParam);

betaI = cell(pbParam.nT, 1);
parfor indTemp = 1 : pbParam.nT
    betaI{indTemp} = cell(numT, 1);

    tInz = pbParam.deltaT * (indTemp - 1);
    tFin = pbParam.deltaT * indTemp;

    for indT = 1 : numT
        locVect = zeros(3, 1);
    
        for indSRext = 1 : methodSpecs.numSRext
            for indGHext = 1 : methodSpecs.numGHext
                indEXT = (indSRext - 1) * methodSpecs.numGHext + indGHext;

                locVect = locVect + (gI(constValues{indT}.EXTn{indEXT}, tFin) - gI(constValues{indT}.EXTn{indEXT}, tInz)) * constValues{indT}.EXTw{indGHext};
            end
        end
    
        locVect(abs(locVect) < 1.0e-14) = 0;
    
        betaI{indTemp}{indT} = locVect;
    end

    betaI{indTemp} = cell2mat(betaI{indTemp});
end

end