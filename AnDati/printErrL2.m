clear
clc

function g = detDir(pbParam)
    h = 1 * strcmp(pbParam.domainType, "barH1") + 3*strcmp(pbParam.domainType, "barH3") + 3*strcmp(pbParam.domainType, "barH3sim");
    velP = pbParam.velP;
    hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
    k_val = (0 : hat_k)';
    sol_an = @(x, t) sum((-1).^k_val .* ( (velP.*t - 2.*h.*k_val - (h-x(3))) .* ((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                    - (velP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((velP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
    g = @(x, t) [0; 0; sol_an(x, t)];
end

addpath(genpath("../functions"));
numGH = 7;
[nodes, weights] = GaussHammerComposite(1, numGH);

for beta = 1 %[1, 2, 4]
    listProblems = ["input_barH1_small.txt", "input_barH1_mid.txt", "input_barH1_large.txt", "input_barH1_maxed.txt", ...
                    "input_barH3_small.txt", "input_barH3_mid.txt", "input_barH3_large.txt", "input_barH3_maxed.txt", ...
                    "input_barH3sim_mid.txt", "input_barH3sim_large.txt", "input_barH3sim_maxed.txt", ...
                ];
    
    normL2err = cell(1, length(listProblems));
    for indProblem = 1 : length(listProblems)
        problemFileName = listProblems(indProblem);
        
        pbParam = readInputFile(problemFileName);
        domainMesh = readSpaceMesh(pbParam.domainType, pbParam.lev);
           
        deltaT = pbParam.Tfin / pbParam.nT;
        
        g = detDir(pbParam);
        path = "./outputData/" + pbParam.domainType + pbParam.lev + "_NEW_19_64_3_256_density.mat";
        if ~exist(path, "file")
            normL2err{indProblem} = [NaN; NaN; NaN];
            continue
        end
        density = load(path).density;
    
        tempVal = zeros(3, pbParam.nT+1);

        parfor indT = 1 : domainMesh.numberTriangles
            vertsT = domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
            l21 = (vertsT(2, :) - vertsT(1, :))';
            l31 = (vertsT(3, :) - vertsT(1, :))';
            n = cross(l21, l31);
            matCoeff{indT} = [-1, -1, 0; 1, 0, 0; 0, 1, 0] / [l21, l31, n];
            vetCoeff{indT} = [1; 0; 0] - matCoeff{indT} * vertsT(1, :)';
        end

    
        parfor indTemp = 1 : pbParam.nT
            tCurr = indTemp * deltaT;
            for indT = 1 : domainMesh.numberTriangles
                nodesT = nodes * domainMesh.coordinates(domainMesh.triangles(indT, 1:3), :);
                weightsT = weights * domainMesh.area(indT);
    
                apprVal = reshape(density(9*(indT-1) + (1:9), indTemp), [3 3]);
                apprFun = @(x) apprVal * (matCoeff{indT}*(x') + vetCoeff{indT});
                f = @(x, t) (g(x, t) -  apprFun(x)) .^ 2;
    
                for indGH = 1 : numGH
                    tempVal(:, indTemp+1) = tempVal(:, indTemp+1) + (weightsT(indGH) .* f(nodesT(indGH, :), tCurr));
                end 
            end
        end
    
        normL2err{indProblem} = sqrt(sum(((tempVal(:, 1 : end-1) + tempVal(:, 2 : end)) .* deltaT ./ 2), 2));
    end
    normL2errMAT = cell2mat(normL2err);
    rapL2err = normL2errMAT(:, 1 : (end-1)) ./ normL2errMAT(:, 2 : end);
    
    errL2 = normL2errMAT';
    rap = zeros(11, 3);
    rap(2:4, :) = rapL2err(:, 1:3)';
    rap(6:8, :) = rapL2err(:, 5:7)';
    rap(10:11, :) = rapL2err(:, 9:10)';
    disp("beta = " + beta)
    % disp(table(errL2))
    disp(table(errL2, rap))
end