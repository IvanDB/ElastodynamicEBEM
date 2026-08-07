function glbIndexFigures = plotMesh(domainMesh, glbIndexFigures, glbFlags, basePath)
%PLOTMESH  Render (and optionally export) a 3D surface plot of the boundary mesh.
%   GLBINDEXFIGURES = PLOTMESH(DOMAINMESH, GLBINDEXFIGURES, GLBFLAGS, BASEPATH) draws
%   DOMAINMESH's triangulation with TRISURF, labels the axes, and titles the figure with the
%   mesh name and refinement level. Returns immediately, without creating a figure, if both
%   GLBFLAGS.plotFigs and GLBFLAGS.saveFigs are false. If GLBFLAGS.saveFigs is true, exports
%   the figure as .jpg (always) and .svg (on Windows only) under BASEPATH/outputPlot/mesh,
%   skipping the export if the target file already exists. If GLBFLAGS.plotFigs is false,
%   the figure is closed right after being (optionally) exported.
%
%   Input arguments:
%       DOMAINMESH      - (struct) triangulated boundary mesh, see READSPACEMESH.
%       GLBINDEXFIGURES - (nonnegative integer) running figure-number
%                         counter, incremented by one if a figure is created.
%       GLBFLAGS        - (struct) must contain plotFigs and
%                         saveFigs logicals, see SETUPWORKSPACE.
%       BASEPATH        - (string, optional, default ".") project root,
%                         used only when GLBFLAGS.saveFigs is true.
%
%   Output arguments:
%       GLBINDEXFIGURES - (nonnegative integer) updated figure counter.
%
%   See also PLOTSOLUTIONS, GETPLOTPROBLEMPARAMETER

arguments
    domainMesh      (1, 1) struct
    glbIndexFigures (1, 1) double {mustBeInteger, mustBeNonnegative}
    glbFlags        (1, 1) struct
    basePath        (1 ,1) string = "."
end

%Early return if both flags are false
if(~glbFlags.plotFigs && ~glbFlags.saveFigs)
    return
end

%Figure initialization
glbIndexFigures = glbIndexFigures + 1;
fig = figure(glbIndexFigures);
fig.Visible = glbFlags.plotFigs;

%Mesh surface plot
trisurf(domainMesh.triangles(:, 1 : 3), domainMesh.coordinates(:, 1), domainMesh.coordinates(:, 2), domainMesh.coordinates(:, 3), ...
            FaceColor = "w", EdgeColor = "b", EdgeAlpha = 0.5, LineWidth = 0.5);

%Set figure parameters
xlabel("$x_1$", Interpreter = 'latex', FontSize = 14, Rotation = 0)
ylabel("$x_2$", Interpreter = 'latex', FontSize = 14, Rotation = 0)
zlabel("$x_3$", Interpreter = 'latex', FontSize = 14, Rotation = 0)

ax = gca;
ax.TickLabelInterpreter = "latex";

title(domainMesh.name + " - lev" + domainMesh.lev, Interpreter = 'latex', FontSize = 12)
axis equal
grid off
box off

%Export image
if(glbFlags.saveFigs)
    folderPath = fullfile(basePath, "outputPlot", "mesh");
    if(~exist(folderPath, 'dir'))
        mkdir(folderPath);
    end
    figName = fullfile(folderPath, domainMesh.name + "-" + domainMesh.lev);
    if(~exist(figName + ".jpg", "file"))
        exportgraphics(fig, figName + ".jpg", ContentType = "image")
    end
    if(ispc && ~exist(figName + ".svg", "file"))
        exportgraphics(fig, figName + ".svg", ContentType = "vector")
    end
end

%Close figure if not visible
if(~glbFlags.plotFigs)
    close(fig)
end
return