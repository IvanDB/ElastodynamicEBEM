function subBlockV = calcSingSubBlockV(pbParam, methodSpecs, constValuesCurr, G1Dn, G1Dw, indTemp)
%CALCSINGSUBBLOCKV  Correct the self-triangle (singular) contribution of a single-layer (V) matrix block.
%   SUBBLOCKV = CALCSINGSUBBLOCKV(PBPARAM, METHODSPECS, CONSTVALUESCURR, G1DN, G1DW, INDTEMP)
%   evaluates, for one triangle and one discrete time-lag INDTEMP, the diagonal self-interaction
%   3x3 sub-block of the single-layer operator that CALCMATRIXV's GPU kernel does not handle.
%   Integration is performed on the light-cone-intersected sub-triangles (childVerts,
%   from CALCCONSTVALUES) via GENERATEFINALG2DNODES, combined with a
%   second-order backward finite-difference in time (coefficients [1, -2, 1])..
%
%   Input arguments:
%       PBPARAM         - (struct) physical/time-discretization parameters
%                         (deltaT, velP, velS, rho), see READINPUTFILE.
%       METHODSPECS     - (struct) quadrature scheme sizes, see GENERATEQUADDATA.
%       CONSTVALUESCURR - (struct) precomputed data for this triangle,
%                         one entry of CALCCONSTVALUES's output.
%       G1DN            - (double) 1D Gauss-Legendre nodes on [-1, 1].
%       G1DW            - (double) 1D Gauss-Legendre weights.
%       INDTEMP         - (nonnegative integer) discrete time-lag index for this block.
%
%   Output arguments:
%       SUBBLOCKV - (3x3 double) singular correction to add on
%                   the diagonal of the triangle's self-block.
%
%   See also CALCMATRIXV, CALCCONSTVALUES, GENERATEFINALG2DNODES

arguments (Input)
    pbParam         (1, 1) struct
    methodSpecs     (1, 1) struct
    constValuesCurr (1, 1) struct
    G1Dn            double
    G1Dw            double
    indTemp         (1, 1) double {mustBeInteger, mustBeNonnegative}
end

arguments (Output)
    subBlockV (3, 3) double
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
    
                [G2DCn, G2DCw] = generateFinalG2Dnodes(constValuesCurr.childVerts{indEXTn, indChild}, rMin, rInt, rExt, G1Dn = G1Dn, G1Dw = G1Dw);
                intgG2DC = intgG2DC + kernelV(length(G2DCw), G2DCn, G2DCw, nodoExt, currT, pbParam.velP, pbParam.velS);
            end
    
            intgTot = intgTot + pesoExt .* intgG2DC;
        end
    end
    
    subBlockV = subBlockV + coef(indZeta) .* intgTot;
end
subBlockV = subBlockV ./ (4*pi*pbParam.rho);
end