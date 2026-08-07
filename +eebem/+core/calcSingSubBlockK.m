function subBlockK = calcSingSubBlockK(pbParam, domainMesh, methodSpecs, constValuesCurr, G1Dn, G1Dw, indTemp, indM, indV)
%CALCSINGSUBBLOCKK  Correct the self-triangle (singular) contribution of a double-layer (K) matrix block.
%   SUBBLOCKK = CALCSINGSUBBLOCKK(PBPARAM, DOMAINMESH, METHODSPECS, CONSTVALUESCURR, G1DN, G1DW, INDTEMP, INDM, INDV)
%   evaluates, for triangle INDM, discrete time-lag INDTEMP and trial vertex INDV
%   (1, 2 or 3, local to the triangle), the singular self-interaction 3x3 sub-block
%   of the double-layer operator that CALCMATRIXK's GPU kernel does not handle. 
%   Integration is performed on the light-cone-intersected sub-triangles 
%   via GENERATEFINALG2DNODES, combined with a third-order backward
%   finite-difference in time (coefficients [-1, 3, -3, 1]). 
%   The function returns early (an all-zero block) once the time-lag exceeds the light-cone
%   support implied by the triangle's largest edge length (DOMAINMESH.maxL).
%
%   Input arguments:
%       PBPARAM         - (struct) physical/time-discretization parameters
%                         (deltaT, velP, velS), see READINPUTFILE.
%       DOMAINMESH      - (struct) triangulated boundary mesh
%                         (normal, maxL), see READSPACEMESH.
%       METHODSPECS     - (struct) quadrature scheme sizes, see GENERATEQUADDATA.
%       CONSTVALUESCURR - (struct) precomputed data for triangle INDM,
%                         one entry of CALCCONSTVALUES's output.
%       G1DN            - (double) 1D Gauss-Legendre nodes on [-1, 1].
%       G1DW            - (double) 1D Gauss-Legendre weights.
%       INDTEMP         - (nonnegative integer) discrete time-lag index for this block.
%       INDM            - (positive integer) index of the triangle.
%       INDV            - (1, 2 or 3) local index of the trial vertex.
%
%   Output arguments:
%       SUBBLOCKK - (3x3 double) singular correction to add at the
%                   (triangle INDM, vertex INDV) position of the block.
%
%   See also CALCMATRIXK, CALCCONSTVALUES, GENERATEFINALG2DNODES

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

                [G2DCn, G2DCw] = generateFinalG2Dnodes(constValuesCurr.childVerts{indEXTn, indChild}, rMin, rInt, rExt, G1Dn = G1Dn, G1Dw = G1Dw);
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