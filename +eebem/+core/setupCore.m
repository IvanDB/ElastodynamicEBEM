function quadData = setupCore(methodId, methodSpecs)
arguments
    methodId (1, 1) double {mustBeInteger, mustBeNonnegative} = 0

    methodSpecs.quadType (1, 1) string {mustBeText} = "FN"

    methodSpecs.numSRext (1, 1) double {mustBeInteger, mustBePositive} = 16
    methodSpecs.numGHext (1, 1) double {mustBeInteger, mustBePositive} = 1

    methodSpecs.numSRint (1, 1) double {mustBeInteger, mustBePositive} = 64
    methodSpecs.numGHint (1, 1) double {mustBeInteger, mustBePositive} = 3

    methodSpecs.numSNGLR (1, 1) double {mustBeInteger, mustBePositive} = 256

    methodSpecs.numBOUND (1, 1) double {mustBeInteger, mustBePositive} = 256
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

methodSpecs.stringID = sprintf("%s_%d-%d_%d-%d_%d_%d", methodSpecs.quadType, methodSpecs.numSRext, methodSpecs.numGHext, ...
                            methodSpecs.numSRint, methodSpecs.numGHint, methodSpecs.numSNGLR, methodSpecs.numBOUND);


assert(isSquare(methodSpecs.numSRext) && isGHvalid(methodSpecs.numGHext), "External quadrature data not valid.")
[quadData.EXTn, quadData.EXTw] = GaussHammerComposite(methodSpecs.numSRext, methodSpecs.numGHext);

assert(isSquare(methodSpecs.numSRint) && isGHvalid(methodSpecs.numGHint), "External quadrature data not valid.")
[quadData.INTn, quadData.INTw] = GaussHammerComposite(methodSpecs.numSRint, methodSpecs.numGHint);

%[quadData.G2Dn, quadData.G2Dw] = doppioGauss1D(methodSpecs.numSNGLR);

quadData.methodSpecs = methodSpecs;
end

function flag = isSquare(val)
    flag = sqrt(val) == floor(sqrt(val));
    return
end

function flag = isGHvalid(val)
    flag = ismember(val, [1, 3, 7, 12, 19]);
    return
end
