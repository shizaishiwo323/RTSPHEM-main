function ExportPhreeqcReactionGridPlot(outputPath, reactionGrid, contourX, contourY, contourPhi, hmaxCm, lengthXAxis, lengthYAxis)
% ExportPhreeqcReactionGridPlot - Save PHREEQC reaction grid over the solid boundary.

if iscell(outputPath) || (isstring(outputPath) && numel(outputPath) > 1)
    outputPaths = cellstr(outputPath);
    for iPath = 1:numel(outputPaths)
        ExportPhreeqcReactionGridPlot(outputPaths{iPath}, reactionGrid, ...
            contourX, contourY, contourPhi, hmaxCm, lengthXAxis, lengthYAxis);
    end
    return;
end

outputPath = char(outputPath);
[outputDir, ~, ~] = fileparts(outputPath);
if ~isempty(outputDir) && ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

coord = reactionGrid.coordV;
triangles = reactionGrid.V0T;
if isfield(reactionGrid, 'numT')
    numTriangles = reactionGrid.numT;
else
    numTriangles = size(triangles, 1);
end

fig = figure('Visible', 'off', 'Position', [100, 100, 1100, 950]);
cleanupFigure = onCleanup(@() closeFigureIfOpen(fig));

ax = axes(fig);
hold(ax, 'on');
patch(ax, 'Faces', triangles, 'Vertices', coord, ...
    'FaceColor', 'none', ...
    'EdgeColor', [0.65 0.72 0.80], ...
    'LineWidth', 0.45);

if ~isempty(contourPhi) && any(isfinite(contourPhi(:)))
    contour(ax, contourX, contourY, contourPhi, [0 0], ...
        'Color', [0.80 0.10 0.08], ...
        'LineWidth', 1.8);
end

plot(ax, [0 lengthXAxis lengthXAxis 0 0], [0 0 lengthYAxis lengthYAxis 0], ...
    'k-', 'LineWidth', 1.0);

axis(ax, 'equal');
xlim(ax, [0 lengthXAxis]);
ylim(ax, [0 lengthYAxis]);
box(ax, 'on');
grid(ax, 'on');
xlabel(ax, 'x [cm]');
ylabel(ax, 'y [cm]');
title(ax, sprintf('PHREEQC reaction grid: Hmax = %.4g cm (%.2f um), %d triangles', ...
    hmaxCm, hmaxCm * 1e4, numTriangles), ...
    'Interpreter', 'none');
legend(ax, {'PHREEQC grid', 'solid boundary (level set = 0)', 'domain'}, ...
    'Location', 'northeastoutside', ...
    'Interpreter', 'none');

print(fig, outputPath, '-dpng', '-r300');
end

function closeFigureIfOpen(fig)
if ~isempty(fig) && ishghandle(fig)
    close(fig);
end
end
