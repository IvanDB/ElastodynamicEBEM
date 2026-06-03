function [pbParam, domainMesh] = finalizeParameters(pbParam, domainMesh, timeSpecs)
arguments
    pbParam    (1, 1) struct
    domainMesh (1, 1) struct
    timeSpecs.betaVal  (1, 1) double {mustBePositive} = 1
    timeSpecs.TimeMult (1, 1) double {mustBePositive} = 1
end

pbParam.beta = timeSpecs.betaVal;

pbParam.Tfin = timeSpecs.TimeMult * pbParam.defTimeLimit;
pbParam.nT = ceil(pbParam.defNumIntvls * (2 ^ (domainMesh.lev - 1)) * pbParam.beta);
pbParam.deltaT = pbParam.Tfin / pbParam.nT;
end