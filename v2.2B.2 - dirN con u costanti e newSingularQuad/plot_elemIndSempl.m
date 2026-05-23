clc
close all

addpath(genpath("./functions"));

for levAnalisi = [1, 2, 3]

    pbParam = readInputFile("input_elementIndSempl" + levAnalisi + "_materialeBarre.txt");
    
    if ~exist("./outputData/elementoIndustrialeSemplificato" + pbParam.lev + "_" + pbParam.lambda + "_NEW_19_64_3_256_density.mat", "file")
        continue
    end
    density = load("./outputData/elementoIndustrialeSemplificato" + pbParam.lev + "_" + pbParam.lambda + "_NEW_19_64_3_256_density.mat").density;
    
    deltaT = pbParam.Tfin / pbParam.nT;
    
    switch pbParam.lev
        case 1
            indTop = 755;
            indMid = 352;
            indBot = 115;
        case {2, 3}
            indTop = 3018;
            indMid = 1408;
            indBot = 460;
        otherwise
            return
    end
    
    figure(levAnalisi)
    cla
    hold on
    plot(deltaT : deltaT : pbParam.Tfin, density(indTop*3, :), DisplayName = "Faccia superiore");
    plot(deltaT : deltaT : pbParam.Tfin, density(indBot*3, :), DisplayName = "Faccia inferiore");
    plot(deltaT : deltaT : pbParam.Tfin, density(indMid*3, :), DisplayName = "Punto medio bordo laterale");
    ax = gca;
    legend(ax, ...
            Interpreter = "latex", ...
            Box = "on", ...
            BackgroundAlpha = 1, ...
            Location = "northwest", ...
            FontSize = 14);
    
    
    title("Andamento temporale componente $u_3$", Interpreter="latex", FontSize=18);
    xlabel("$t$", FontSize = 18, Interpreter="latex")
    ylabel("$u_3(x, t)$", FontSize = 18, Interpreter="latex")

    %print("u_3(x, t)_" + pbParam.domainType + "_" + "materialeBarre", "-dsvg", "-vector");
    %print("u_3(x, t)_" + pbParam.domainType + "_" + "metallo", "-dsvg", "-vector");
end