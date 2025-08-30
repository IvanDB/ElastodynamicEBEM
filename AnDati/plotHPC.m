
addpath(genpath("../functions"));

quadName = "_NEW_12_16_3_144_";

for beta = 1 %[1 2 4]
    listProblems = ["input_barH1_small.txt", "input_barH1_mid.txt", "input_barH1_large.txt", "input_barH1_maxed.txt", ...
                    "input_barH3_small.txt", "input_barH3_mid.txt", "input_barH3_large.txt", "input_barH3_maxed.txt", ...
                    "input_barH3sim_mid.txt", "input_barH3sim_large.txt", "input_barH3sim_maxed.txt", ...
                   ];
    
    for indProblem = 1 : length(listProblems)
        switch indProblem     
            case {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
                path = "./outputData/";    
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
            case 'barH3'
                h = 3;
                hTitle = "3base";
                minY = -3.25;
                maxY = 7.25;
            case 'barH1'
                h = 1;
                hTitle = "1base";
                minY = -0.95;
                maxY = 2.25;
            case 'barH3sim'
                h = 3;
                hTitle = "barH3sim";
                minY = -3.25;
                maxY = 7.25;
        end
    
        switch [pbParam.domainType, pbParam.lev]
            case {['barH3', 0], ['barH3sim', 0]}
                indTop = 1;
                indBot = 9;
                legendLoc = "southwest";
                levName = "small";
                errY = 1.5;
            case ['barH1', 0]
                indTop = 7;
                indBot = 22;
                legendLoc = "southwest";
                levName = "small";
                errY = 0.5;
            case {['barH3', 1], ['barH3sim', 1]}
                indTop = 2;
                indBot = 38;
                legendLoc = "southwest";
                levName = "mid";
                errY = 1.5;
            case ['barH1', 1]
                indTop = 9;
                indBot = 87;
                legendLoc = "southwest";
                levName = "mid";
                errY = 0.5;
            case {['barH3', 2], ['barH3sim', 2]}
                indTop = 7;
                indBot = 135;
                legendLoc = "southwest";
                levName = "large";
            case ['barH1', 2]
                indTop = 17;
                indBot = 283;
                legendLoc = "southwest";
                levName = "large";
                errY = 0.5;
            case {['barH3', 3], ['barH3sim', 3]}
                indTop = 27;
                indBot = 539;
                legendLoc = "southwest";
                levName = "maxed";
            case ['barH1', 3]
                indTop = 65;
                indBot = 1131;
                legendLoc = "southwest";
                levName = "maxed";
        end

        % if ~exist("./outputData/" + pbParam.domainType + pbParam.lev + quadName + "density.mat", "file")
        %         continue
        % end
        % densityREF = load("./outputData/HPC - test beta = 1/" + pbParam.domainType + pbParam.lev + quadName + "density.mat").density;
        % densityTEST = kron(densityREF, ones(1, beta));

        figure(10*indProblem)
        cla
        hold on
        densValueTop = mean(density((indTop-1)*9 + [3 6 9], :), 1);
        densValueBot = mean(density((indBot-1)*9 + [3 6 9], :), 1);
        plot(deltaT : deltaT : pbParam.Tfin, densValueTop, DisplayName = "$\mbox{\boldmath $x$} = [0, 0, h]$");
        plot(deltaT : deltaT : pbParam.Tfin, densValueBot, DisplayName = "$\mbox{\boldmath $x$} = [0, 0, 0]$");


        ax = gca;
        legend(ax, ...
                Interpreter = "latex", ...
                Box = "on", ...
                BackgroundAlpha = 1, ...
                Location = legendLoc, ...
                FontSize = 14);
        xlim([0 pbParam.Tfin])
        ylim([minY maxY])
        title("h=" + hTitle +", " + levName + ", beta=" + beta);
        print("plotTest/u3 - beta = " + beta + " - h" + hTitle + " " + levName + ".svg", "-dsvg", "-vector");
        % 
        % if (beta ~= 1)
        %     figure(100*beta + indProblem)
        %     cla
        %     hold on
        %     plot(deltaT : deltaT : pbParam.Tfin, density(indTop*3, :) - densityTEST(indTop*3, :), DisplayName = "$\mbox{\boldmath $x$} = [0, 0, h]$");
        %     plot(deltaT : deltaT : pbParam.Tfin, density(indBot*3, :) - densityTEST(indBot*3, :), DisplayName = "$\mbox{\boldmath $x$} = [0, 0, 0]$");
        % 
        %     ax = gca;
        %     legend(ax, ...
        %             Interpreter = "latex", ...
        %             Box = "on", ...
        %             BackgroundAlpha = 1, ...
        %             Location = legendLoc, ...
        %             FontSize = 14);
        %     xlim([0 pbParam.Tfin])
        %     ylim([-errY errY])
        %     title("h=" + h +", " + levName + ", beta=" + beta);
        %     print("plotTestBeta/errVSbeta1 - beta = " + beta + " - h" + h + " " + levName + ".svg", "-dsvg", "-vector");
        % end

    end
end

