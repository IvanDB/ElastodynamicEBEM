function plotNodes(V3, nodi3D, mType, mSize, figMainIndex)
    persistent counter;
    persistent figSubIndex;
    persistent figureMainSavedIndex;
    if isempty(counter)
        counter = 0;
    end
    if isempty(figSubIndex)
        figSubIndex = 0;
    end
    if isempty(figureMainSavedIndex)
        figureMainSavedIndex = figMainIndex;
    end
    if(figureMainSavedIndex ~= figMainIndex)
        figureMainSavedIndex = figMainIndex;
        figSubIndex = 0;
    end
    if(mod(counter, 6) == 0)
        figSubIndex = figSubIndex + 1;
    end
    figure(figMainIndex + figSubIndex);
    title("Nodi su 3D")
    hold on
    view(135, 30);
    fill3(V3(:, 1), V3(:, 2), V3(:, 3), "w", LineStyle="-", FaceColor="w");
    plot3(nodi3D(:, 1), nodi3D(:, 2), nodi3D(:, 3), mType, MarkerSize=mSize);
    axis equal
    grid on
    box off
    counter = counter + 1;
end