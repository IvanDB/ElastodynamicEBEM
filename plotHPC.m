
%addpath(genpath("./functions"));

% listProblems = ["input_barH1_largeHPC1.txt", "input_barH1_maxedHPC1.txt", "input_barH1_largeHPC2.txt", "input_barH1_maxedHPC2.txt", ...
%                 "input_barH3_largeHPC1.txt", "input_barH3_maxedHPC1.txt", "input_barH3_largeHPC2.txt", "input_barH3_maxedHPC2.txt"
%                ];

% quadName = "_NEW_37_256_1_256_";
quadName = "_NEW_19_64_3_256_";

for beta = [1 2 4]
    listProblems = "beta=" + beta + "/" + ["input_barH1_small.txt", "input_barH1_mid.txt", "input_barH1_large.txt", "input_barH1_maxed.txt", ...
                                            "input_barH3_small.txt", "input_barH3_mid.txt", "input_barH3_large.txt", "input_barH3_maxed.txt", ...
                                            "input_barH3sim_mid.txt", "input_barH3sim_large.txt", "input_barH3sim_maxed.txt", ...
                                           ];
    
    % listProblems = ["input_barH1_maxedMOD.txt", "input_barH3_maxedMOD.txt"];
    
    for indProblem = 1 : length(listProblems)
        switch indProblem 
            % case {1, 2, 5, 6}
            %     path = "./outputData/HPC - primo round/";
            % case {3, 4, 7, 8}
            %     path = "./outputData/HPC - secondo round/";
    
            case {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
                path = "./outputData/HPC - test beta = " + beta + "/";
            
            % case {1, 2}
            %     path = "./outputData/HPC - testfin beta = 1/";
    
            otherwise
                continue
        end
    
        problemFileName = listProblems(indProblem);
        
        pbParam = readInputFile(problemFileName);
        domainMesh = readSpaceMesh(pbParam.domainType, pbParam.lev);
        
        path = path + pbParam.domainType + pbParam.lev + quadName + "density.mat";

        if ~exist(path, "file")
                continue
        end
        
        density = load(path).density;
    
        deltaT = pbParam.Tfin / pbParam.nT;
        if (beta ~= round(domainMesh.lMin ./ (pbParam.velP*deltaT)))
            error("PROBLEMA");
        end
        
        % switch [pbParam.domainType, beta]
        %     case ['barH1', 2]
        %         tFin = 20;
        %         maxY = 2.25;
        %     case ['barH1', 4]
        %         tFin = 24;
        %         maxY = 2.25;
        %     case ['barH3', 2]
        %         tFin = 20;
        %         maxY = 7.25;
        %     case ['barH3', 4]
        %         tFin = 48;
        %         maxY = 7.25;
        % end
        
        switch pbParam.domainType
            case {'barH3', "barH3sim"}
                h = 3;
                minY = -3.25;
                maxY = 7.25;
            case 'barH1'
                h = 1;
                minY = -0.95;
                maxY = 2.25;
        end
    
        switch [pbParam.domainType, pbParam.lev]
            case {['barH3', 0], ['barH3sim', 0]}
                indTop = 1;
                indBot = 9;
                legendLoc = "northwest";
                levName = "small";
                errY = 1.5;
            case ['barH1', 0]
                indTop = 7;
                indBot = 22;
                legendLoc = "northeast";
                levName = "small";
                errY = 0.5;
            case {['barH3', 1], ['barH3sim', 1]}
                indTop = 2;
                indBot = 38;
                legendLoc = "northwest";
                levName = "mid";
                errY = 1.5;
            case ['barH1', 1]
                indTop = 9;
                indBot = 87;
                legendLoc = "northeast";
                levName = "mid";
                errY = 0.5;
            case {['barH3', 2], ['barH3sim', 2]}
                indTop = 7;
                indBot = 135;
                legendLoc = "northwest";
                levName = "large";
            case ['barH1', 2]
                indTop = 17;
                indBot = 283;
                legendLoc = "northeast";
                levName = "large";
                errY = 0.5;
            case {['barH3', 3], ['barH3sim', 3]}
                indTop = 27;
                indBot = 539;
                legendLoc = "northwest";
                levName = "maxed";
            case ['barH1', 3]
                indTop = 65;
                indBot = 1131;
                legendLoc = "northeast";
                levName = "maxed";
        end
        if ~exist("./outputData/HPC - test beta = 1/" + pbParam.domainType + pbParam.lev + quadName + "density.mat", "file")
                continue
        end

        densityREF = load("./outputData/HPC - test beta = 1/" + pbParam.domainType + pbParam.lev + quadName + "density.mat").density;
        densityTEST = kron(densityREF, ones(1, beta));

        figure(10*beta + indProblem)
        cla
        hold on
        plot(deltaT : deltaT : pbParam.Tfin, density(indTop*3, :), DisplayName = "$\mbox{\boldmath $x$} = [0, 0, h]$");
        plot(deltaT : deltaT : pbParam.Tfin, density(indBot*3, :), DisplayName = "$\mbox{\boldmath $x$} = [0, 0, 0]$");
        
        
        ax = gca;
        legend(ax, ...
                Interpreter = "latex", ...
                Box = "on", ...
                BackgroundAlpha = 1, ...
                Location = legendLoc, ...
                FontSize = 14);
        xlim([0 pbParam.Tfin])
        ylim([minY maxY])
        title("h=" + h +", " + levName + ", beta=" + beta);
        print("plotTestBeta/u3 - beta = " + beta + " - h" + h + " " + levName + ".svg", "-dsvg", "-vector");

        if (beta ~= 1)
            figure(100*beta + indProblem)
            cla
            hold on
            plot(deltaT : deltaT : pbParam.Tfin, density(indTop*3, :) - densityTEST(indTop*3, :), DisplayName = "$\mbox{\boldmath $x$} = [0, 0, h]$");
            plot(deltaT : deltaT : pbParam.Tfin, density(indBot*3, :) - densityTEST(indBot*3, :), DisplayName = "$\mbox{\boldmath $x$} = [0, 0, 0]$");
            
            ax = gca;
            legend(ax, ...
                    Interpreter = "latex", ...
                    Box = "on", ...
                    BackgroundAlpha = 1, ...
                    Location = legendLoc, ...
                    FontSize = 14);
            xlim([0 pbParam.Tfin])
            ylim([-errY errY])
            title("h=" + h +", " + levName + ", beta=" + beta);
            print("plotTestBeta/errVSbeta1 - beta = " + beta + " - h" + h + " " + levName + ".svg", "-dsvg", "-vector");
        end

    end
end

