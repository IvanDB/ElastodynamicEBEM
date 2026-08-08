function [nDim, tVal, jVal, climCustom] = getPlotProblemParameter(pbParam)
%GETPLOTPROBLEMPARAMETER  Look up the plotting defaults (time instants, components) for a benchmark problem.
%   [NDIM, TVAL, JVAL, CLIMCUSTOM] = GETPLOTPROBLEMPARAMETER(PBPARAM) returns, for each
%   built-in benchmark problem (selected via PBPARAM.domainName), hand-picked defaults
%   describing which solution snapshots PLOTSOLUTIONS should render: the plot
%   dimensionality, which discrete time-step indices to plot, which vector component(s)
%   (or 0 for the magnitude) to plot, and whether a custom color-axis limit should be used.
%
%   Input arguments:
%       PBPARAM - (struct) physical/problem parameters; only PBPARAM.domainName (and,
%                 for several cases, PBPARAM.nT/Tfin) are used, see READINPUTFILE.
%
%   Output arguments:
%       NDIM       - (integer) 2 or 3, the view dimensionality passed to VIEW.
%       TVAL       - (integer vector) discrete time-step indices to plot.
%       JVAL       - (integer vector) vector component(s) to plot (1, 2,
%                    3, or a combination); 0 means "plot the magnitude".
%       CLIMCUSTOM - (logical) whether a custom color-axis limit should
%                    be applied (currently unimplemented, see Notes).
%
%   Notes:
%       Emits a warning and returns without setting any output for PBPARAM.domainName
%       values not in its internal list (e.g. custom user problems): callers must handle
%       that case, as PLOTSOLUTIONS does via its early "plotFigs/saveFigs" flag check. The
%       CLIMCUSTOM branch is marked "still WIP" where it is consumed in PLOTSOLUTIONS.
%
%   See also PLOTSOLUTIONS, PLOTMESH

arguments (Input)
    pbParam (1, 1) struct
end

arguments (Output)
    nDim        (1, 1) double {mustBeInteger}
    tVal        (1, :) double {mustBeInteger}
    jVal        (1, :) double {mustBeInteger}
    climCustom  (1, 1) logical
end

%Assegnazione parametri in base al problema
switch pbParam.domainName
    case 'screenTest'
        nDim = 2;
        tVal = pbParam.nT;
        jVal = 1 : 3;
        climCustom = false;
    case 'screenUniform'
        nDim = 2;
        tVal = pbParam.nT;
        jVal = 3;
        climCustom = false;
    case 'screenGraded'
        nDim = 2;
        tVal = pbParam.nT;
        jVal = 3;
        climCustom = false;
    case 'sphereUniform'
        nDim = 3;
        tVal = [pbParam.nT .* (1/3), pbParam.nT * (2/3), (pbParam.nT * (2/3)) + 1, pbParam.nT];
        jVal = 3;
        climCustom = false;
    case 'sphereNotUniform'
        nDim = 3;
        tVal = [pbParam.nT .* (1/3), pbParam.nT * (2/3), (pbParam.nT * (2/3)) + 1, pbParam.nT];
        jVal = 3;
        climCustom = false;
    case "barH1"
        nDim = 3;
        tVal = pbParam.nT .* (1 : 6) ./ pbParam.Tfin ;
        jVal = 3;
        climCustom = false;
    case "barH3"
        nDim = 3;
        tVal = pbParam.nT .* (1 : pbParam.Tfin) ./ pbParam.Tfin ;
        jVal = 1 : 3;
        climCustom = false;
    case 'sphereWave'
        return
    case "DesCop-cube"
        nDim = 3;
        tVal = floor(pbParam.nT .* [7.50, 8.25, 9.00, 9.75, 10.50, 11.25] ./ pbParam.Tfin);
        jVal = 0;
        climCustom = false;
    case "DesCop-sphere"
        nDim = 3;
        tVal = 1 : 2 : pbParam.nT;
        jVal = [1, 3];
        climCustom = false;
    case "industrialComponent"
        nDim = 3;
        tVal = [1 : 1 : 9, 10 : 10 : 99, 100 : 100 : pbParam.nT];
        jVal = 1 : 3;
        climCustom = false;
    case "ConeBall"
        nDim = 3;
        tVal = pbParam.nT;
        jVal = 1 : 3;
        climCustom = false;
    otherwise
        warning("Problem plots specs not set. Loading from file in not available.")
end
end
