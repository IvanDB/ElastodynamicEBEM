function quadData = generateQuadData(methodId, methodSpecs)
%GENERATEQUADDATA  Build the quadrature nodes/weights used throughout the BEM assembly.
%   QUADDATA = GENERATEQUADDATA(METHODID, METHODSPECS) either looks up a predefined
%   quadrature scheme by METHODID (an index into an internal list of preset 
%   "quadType + sizes" strings), or takes the sizes directly from the METHODSPECS
%   name-value arguments (used when METHODID is 0, the default). 
%   From the resulting sizes it builds:
%   a composite Gauss-Hammer rule for the outer/test integrations (QUADDATA.EXTn/EXTw),
%   a composite Gauss-Hammer rule for the inner/trial integrations (QUADDATA.INTn/INTw), and
%   a 1D Gauss-Legendre rule for the singular integrations (QUADDATA.G1Dn/G1Dw).
%
%   Input arguments:
%       METHODID    - (nonnegative integer, default 0) index into the
%                     internal preset list; 0 means "use METHODSPECS".
%       METHODSPECS - (name-value) quadType (string, default "FN"); 
%                     numSRext/numGHext (outer sub-triangles / Gauss- Hammer
%                      nodes per sub-triangle, defaults 16/1); 
%                     numSRint/numGHint (inner counterparts, defaults 64/3);
%                     numSNGLR (radial singular quadrature order, default 256);
%                     numBOUND (boundary/edge integration points, used only by
%                      CALCSINGSUBBLOCKK_C, default 256).
%
%   Output arguments:
%       QUADDATA - (struct) with fields EXTn, EXTw, INTn, INTw, G1Dn,
%                  G1Dw and methodSpecs (the resolved METHODSPECS,
%                  including derived numEXT, numINT and stringID).
%
%   Notes:
%       Asserts that numSRext/numSRint are perfect squares (required by
%       the composite Gauss-Hammer subdivision) and that numGHext/numGHint 
%       are one of the valid base-rule sizes {1, 3, 7, 12, 19}.
%
%   See also eebem.utility.quadratureRules.GaussHammerComposite,
%   eebem.utility.quadratureRules.Gauss1D, eebem.core.calcMatrixSpecs

arguments (Input)
    methodId (1, 1) double {mustBeInteger, mustBeNonnegative} = 0

    methodSpecs.quadType (1, 1) string {mustBeText} = "FN"

    methodSpecs.numSRext (1, 1) double {mustBeInteger, mustBePositive} = 16
    methodSpecs.numGHext (1, 1) double {mustBeInteger, mustBePositive} = 1

    methodSpecs.numSRint (1, 1) double {mustBeInteger, mustBePositive} = 64
    methodSpecs.numGHint (1, 1) double {mustBeInteger, mustBePositive} = 3

    methodSpecs.numSNGLR (1, 1) double {mustBeInteger, mustBePositive} = 256

    methodSpecs.numBOUND (1, 1) double {mustBeInteger, mustBePositive} = 256
end

arguments (Output)
    quadData (1, 1) struct
end

import eebem.utility.quadratureRules.*

% TODO: move list to a .txt file
methodsList = ["SA 01 03", "SA 01 07", "SA 01 12", "SA 01 19", ...                                  % cose
               "MX 01 12 16 03", "MX 01 12 64 03", "MX 01 19 16 03", "MX 01 19 64 03", ...          % cose
               "FN 01 19 64 03 256 256", ...                                                      % old 27
               "FN 16 01 64 03 256 256", ...                                                      % old 40
               "FN 16 03 64 03 256 256", ...                                                      % old 48
              ];

if(methodId > 0)
    assert(methodId <= length(methodsList), sprintf("Method not available. ID must be <= %d.", length(methodsList)))
    methodString = methodsList(methodId);

    [methodSpecs.quadType, ~, ~, nextIdx] = sscanf(methodString, "%s", 1);
    methodString = extractAfter(methodString, nextIdx);

    [methodSpecs.numSRext, ~, ~, nextIdx] = sscanf(methodString, "%d", 1);
    methodString = extractAfter(methodString, nextIdx);
    [methodSpecs.numGHext, ~, ~, nextIdx] = sscanf(methodString, "%d", 1);
    methodString = extractAfter(methodString, nextIdx);

    [methodSpecs.numSRint, ~, ~, nextIdx] = sscanf(methodString, "%d", 1);
    methodString = extractAfter(methodString, nextIdx);
    [methodSpecs.numGHint, ~, ~, nextIdx] = sscanf(methodString, "%d", 1);
    methodString = extractAfter(methodString, nextIdx);


    [methodSpecs.numSNGLR, ~, ~, nextIdx] = sscanf(methodString, "%d", 1);
    methodString = extractAfter(methodString, nextIdx);

    [methodSpecs.numBOUND, ~, ~, ~] = sscanf(methodString, "%d", 1);
end

methodSpecs.numEXT = methodSpecs.numSRext * methodSpecs.numGHext;
methodSpecs.numINT = methodSpecs.numSRint * methodSpecs.numGHint;

methodSpecs.stringID = sprintf("%s-%d-%d-%d-%d-%d-%d", methodSpecs.quadType, methodSpecs.numSRext, methodSpecs.numGHext, ...
                            methodSpecs.numSRint, methodSpecs.numGHint, methodSpecs.numSNGLR, methodSpecs.numBOUND);


assert(isSquare(methodSpecs.numSRext) && isGHvalid(methodSpecs.numGHext), "External quadrature data not valid.")
[quadData.EXTn, quadData.EXTw] = GaussHammerComposite(methodSpecs.numSRext, methodSpecs.numGHext);

assert(isSquare(methodSpecs.numSRint) && isGHvalid(methodSpecs.numGHint), "Internal quadrature data not valid.")
[quadData.INTn, quadData.INTw] = GaussHammerComposite(methodSpecs.numSRint, methodSpecs.numGHint);

assert(isSquare(methodSpecs.numSNGLR), "Singular quadrature data not valid.")
[quadData.G1Dn, quadData.G1Dw] = Gauss1D(1, sqrt(methodSpecs.numSNGLR), 0, 0);

quadData.methodSpecs = methodSpecs;
end

function flag = isSquare(val)
    %True iff VAL is a perfect square (required outer/inner sub-triangle count).
    flag = sqrt(val) == floor(sqrt(val));
    return
end

function flag = isGHvalid(val)
    %True iff VAL is a supported base Gauss-Hammer rule size {1,3,7,12,19}.
    flag = ismember(val, [1, 3, 7, 12, 19]);
    return
end
