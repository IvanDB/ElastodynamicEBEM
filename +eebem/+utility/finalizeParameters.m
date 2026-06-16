function [pbParam, domainMesh] = finalizeParameters(pbParam, domainMesh, timeSpecs)
arguments
    pbParam    (1, 1) struct
    domainMesh (1, 1) struct
    timeSpecs.betaMult  (1, 1) double {mustBePositive} = 1
    timeSpecs.timeMult (1, 1) double {mustBePositive} = 1
end

pbParam.beta = timeSpecs.betaMult;
pbParam.tMlt = timeSpecs.timeMult;

pbParam.Tfin = pbParam.tMlt * pbParam.defTimeLimit;

spaceTimeCouplingFactor = 2 ^ (pbParam.STcoupling * (domainMesh.lev - 1));
pbParam.nT = ceil(pbParam.defNumIntvls * spaceTimeCouplingFactor * pbParam.beta * pbParam.tMlt);

pbParam.deltaT = pbParam.Tfin / pbParam.nT;
end