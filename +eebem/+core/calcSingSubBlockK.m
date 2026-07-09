function subBlockK = calcSingSubBlockK(pbParam, domainMesh, methodSpecs, constValuesCurr, G1Dn, G1Dw, indTemp, indM, indV)
%CALCSINGSUBBLOCKK Summary of this function goes here
%   Detailed explanation goes here

arguments (Input)
    pbParam
    domainMesh
    methodSpecs
    constValuesCurr
    G1Dn
    G1Dw
    indTemp
    indM
    indV
end

import eebem.utility.quadratureRules.*

subBlockK = zeros(3, 3);
coef = [-1, 3, -3, 1];
istTemp = indTemp + [-2, -1, 0, 1];

if(indTemp >= ceil(domainMesh.maxL(indM) / (pbParam.velS * pbParam.deltaT)) + 2)
    return;
end

minDiffTemp = istTemp(1) * pbParam.deltaT;
normInt = domainMesh.normal(indM, :);
vettVMS = cross(normInt, constValuesCurr.matCoeff(indV, :));

for indZeta = 1 : 4
    if (istTemp(indZeta) <= 0)
        continue
    end

    intgTot = zeros(3, 3);
    currT = istTemp(indZeta) * pbParam.deltaT;

    for indSRext = 1 : methodSpecs.numSRext
        for indGHext = 1 : methodSpecs.numGHext
            indEXTn = (indSRext-1) * methodSpecs.numGHext + indGHext;
            nodoExt = constValuesCurr.EXTn{indEXTn};
            pesoExt = constValuesCurr.EXTw{indGHext};

            intgIntSing = zeros(3);
        
            for indChild = 1 : 3
                rMin = max(pbParam.velS * minDiffTemp, 0);
                rInt = pbParam.velS * currT;
                rExt = pbParam.velP * currT;

                [G2DCn, G2DCw] = generateFinalG2Dnodes(constValuesCurr.childVerts{indEXTn, indChild}, rMin, rInt, rExt, G1Dn, G1Dw);
                intgIntSing = intgIntSing + kernelK(length(G2DCw), G2DCn, G2DCw, nodoExt, currT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, ...
                                                            vettVMS, constValuesCurr.matCoeff, constValuesCurr.vetCoeff, indV, normInt);
            end

            intgTot = intgTot + pesoExt .* intgIntSing;
        end
    end

    subBlockK = subBlockK + coef(indZeta) .* intgTot;
end

subBlockK = subBlockK ./ (4 * pi * pbParam.deltaT);
end