function [nDim, tVal, jVal, climCustom] = getPlotProblemParameter(pbParam)
%GETPLOTPROBLEMPARAMETER Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    pbParam
end

arguments (Output)
    nDim        double {mustBeInteger}
    tVal        double {mustBeInteger}
    jVal        double {mustBeInteger}
    climCustom  logical
end

%Assegnazione parametri in base al problema
switch pbParam.domainType
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
    case {"barH1-symm", "barH1-asym"}
        nDim = 3;
        tVal = [(1 : 1/(pbParam.Tfin/pbParam.nT)), (2 : pbParam.Tfin)./(pbParam.Tfin/pbParam.nT)] ;
        jVal = 3;
        climCustom = false;
    case {"barH3-symm", "barH3-asym"}
        nDim = 3;
        tVal = pbParam.nT .* (1 : pbParam.Tfin) ./ pbParam.Tfin ;
        jVal = 3;
        climCustom = false;
    case 'sphereWave'
        return
    case {"DesCop-cube-symm", "DesCop-cube-asym"}
        nDim = 3;
        tVal = [(1 : 1/(pbParam.Tfin/pbParam.nT)), (2 : pbParam.Tfin)./(pbParam.Tfin/pbParam.nT)]; %floor(pbParam.nT .* [7.50, 8.25, 9.00, 9.75, 10.50, 11.25] ./ pbParam.Tfin);
        jVal = 0;
        climCustom = true;
    case "DesCop-sphere"
        nDim = 3;
        tVal = 1 : 2 : pbParam.nT;
        jVal = 3;
        climCustom = true;
    case 'elementoIndustriale'
        nDim = 3;
        tVal = 1;
        jVal = 3;
        climCustom = false;
end
end