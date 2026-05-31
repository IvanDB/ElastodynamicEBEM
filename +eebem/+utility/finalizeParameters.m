function [pbParam, domainMesh] = finalizeParameters(pbParam, domainMesh, timeSpecs)
arguments
    pbParam    (1, 1) struct
    domainMesh (1, 1) struct
    timeSpecs.betaVal  (1, 1) double {mustBePositive} = 1
    timeSpecs.TimeMult (1, 1) double {mustBePositive} = 1
end

pbParam.betaVal = timeSpecs.betaVal;

pbParam.Tfin = timeSpecs.TimeMult * pbParam.defTimeLimit;
pbParam.nT = ceil(pbParam.defNumIntvls * domainMesh.lev * timeSpecs.betaVal);
pbParam.deltaT = pbParam.Tfin / pbParam.nT;
end