function [pbParam, domainMesh] = finalizeParameters(pbParam, domainMesh, timeSpecs)
%FINALIZEPARAMETERS  Derive the final time-discretization parameters (Tfin, nT, deltaT).
%   [PBPARAM, DOMAINMESH] = FINALIZEPARAMETERS(PBPARAM, DOMAINMESH, TIMESPECS)
%   combines the problem's default time horizon and interval count
%   (PBPARAM.defaultValues.timeLimit/numIntvls) with the user-supplied multipliers and
%   the mesh refinement level to produce the final simulation horizon PBPARAM.Tfin,
%   number of time steps PBPARAM.nT, and time step PBPARAM.deltaT = Tfin / nT.
%   If PBPARAM.STcoupling is true the number of intervals is scaled by 2^(mesh level - 1),
%   implementing the space-time mesh coupling used to keep the CFL-like ratio between
%   mesh size and time step roughly constant when the mesh is refined.
%
%   Input arguments:
%       PBPARAM    - (struct) physical/problem parameters, in particular
%                    defaultValues.timeLimit, defaultValues.numIntvls
%                    and STcoupling, see READINPUTFILE.
%       DOMAINMESH - (struct) triangulated boundary mesh; only DOMAINMESH.lev
%                    (refinement level) is used, see READSPACEMESH.
%       TIMESPECS  - (name-value) betaMult (positive double, default 1)
%                    multiplies the interval count; timeMult (positive
%                    double, default 1) multiplies the time horizon.
%
%   Output arguments:
%       PBPARAM    - (struct) input struct with beta, tMlt, Tfin,
%                    nT and deltaT fields added, see READINPUTFILE.
%       DOMAINMESH - (struct) returned unchanged (kept for a symmetrical,
%                    future-proof call signature), see READSPACEMESH.
%
%   See also GENERATEFILENAMES, eebem.core.calcNumMatrixBlocks

arguments (Input)
    pbParam    (1, 1) struct
    domainMesh (1, 1) struct
    timeSpecs.betaMult (1, 1) double {mustBePositive} = 1
    timeSpecs.timeMult (1, 1) double {mustBePositive} = 1
end

arguments (Output)
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
end

pbParam.beta = timeSpecs.betaMult;
pbParam.tMlt = timeSpecs.timeMult;

pbParam.Tfin = pbParam.tMlt * pbParam.defaultValues.timeLimit;

spaceTimeCouplingFactor = 2 ^ (pbParam.STcoupling * (domainMesh.lev - 1));
pbParam.nT = ceil(pbParam.defaultValues.numIntvls * spaceTimeCouplingFactor * pbParam.beta * pbParam.tMlt);

pbParam.deltaT = pbParam.Tfin / pbParam.nT;
end