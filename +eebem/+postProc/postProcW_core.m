function solution = postProcW_core(basePath, pbParam, domainMesh, density, numPoints, xVal, tVal, methodInfo, PPn, PPw)
import eebem.postProc.*
% Compute dimentions of solution matrix
xSize = size(xVal, 1);
tSize = length(tVal);
solRaw = cell(xSize * tSize, 1); % "raw" solution initialization

parfor ind = 1 : numPoints
    [indX, indT] = ind2sub([xSize, tSize], ind);
    solRaw{ind} = postProcW_calc(pbParam, domainMesh, density, methodInfo, xVal(indX, :), tVal(indT), PPn, PPw);
end

solution = reshape(solRaw, [xSize, tSize]);
%matrix_max = reshape(maxP_raw, [xSize, tSize]);
%matrix_cond = reshape(condP_raw, [xSize, tSize]);
%matrix_nan = reshape(hasNaN_raw, [xSize, tSize]);

%Salvataggio campo vettoriale calcolato
outPath = fullfile(basePath, "solutionData", "IN_" + pbParam.domainType + pbParam.lev + "FN_16-1_64-3_256_256");
save(outPath + "_solution", 'solution');


return