function exported = precip_ExportBenchmarkSnapshots(snapshotDir, targetTimesSeconds, exported, ...
    stepIndex, timeSeconds, grid, currentLevelSet, initialLevelSet, config)
%PRECIP_EXPORTBENCHMARKSNAPSHOTS Export stable benchmark interface snapshots.
%
% The helper writes one PNG per target time once the simulation time reaches
% that target. Filenames are stable for Zhang/Yoon comparison times:
% benchmark_snapshot_013min.png, benchmark_snapshot_018min.png, and
% benchmark_snapshot_118min.png.

if isempty(targetTimesSeconds)
    return;
end

targetTimesSeconds = double(targetTimesSeconds(:));
if nargin < 3 || isempty(exported)
    exported = false(size(targetTimesSeconds));
else
    exported = logical(exported(:));
end
if numel(exported) ~= numel(targetTimesSeconds)
    error('RTSPHEM:Precipitate:InvalidSnapshotExportState', ...
        'The exported flag vector must have the same length as targetTimesSeconds.');
end

snapshotGrid = normalizeSnapshotGrid(grid);
validateSnapshotInputs(snapshotDir, targetTimesSeconds, stepIndex, timeSeconds, ...
    snapshotGrid, currentLevelSet, initialLevelSet);

due = ~exported & timeSeconds >= targetTimesSeconds;
if ~any(due)
    return;
end

if ~exist(snapshotDir, 'dir')
    mkdir(snapshotDir);
end

currentSegments = extractZeroLevelSegments(snapshotGrid, currentLevelSet);
initialSegments = extractZeroLevelSegments(snapshotGrid, initialLevelSet);
for idx = find(due(:))'
    outputPath = fullfile(snapshotDir, precip_BenchmarkSnapshotFilename(targetTimesSeconds(idx)));
    writeSnapshotPng(outputPath, targetTimesSeconds(idx), stepIndex, timeSeconds, ...
        snapshotGrid, currentLevelSet, currentSegments, initialSegments, config);
    exported(idx) = true;
end
end

function snapshotGrid = normalizeSnapshotGrid(grid)
snapshotGrid = struct();
snapshotGrid.coordV = getGridMember(grid, 'coordV');
snapshotGrid.V0T = getGridMember(grid, 'V0T');
end

function value = getGridMember(grid, memberName)
if isstruct(grid) && isfield(grid, memberName)
    value = grid.(memberName);
elseif isobject(grid) && isprop(grid, memberName)
    value = grid.(memberName);
else
    error('RTSPHEM:Precipitate:InvalidSnapshotGrid', ...
        'grid must provide %s as a struct field or object property.', memberName);
end
end

function validateSnapshotInputs(snapshotDir, targetTimesSeconds, stepIndex, timeSeconds, ...
    grid, currentLevelSet, initialLevelSet)
if ~(ischar(snapshotDir) || isstring(snapshotDir)) || strlength(string(snapshotDir)) == 0
    error('RTSPHEM:Precipitate:InvalidSnapshotPath', ...
        'snapshotDir must be a nonempty character vector or string scalar.');
end
if any(~isfinite(targetTimesSeconds)) || any(targetTimesSeconds < 0)
    error('RTSPHEM:Precipitate:InvalidSnapshotTimes', ...
        'targetTimesSeconds must be finite nonnegative values.');
end
if ~isscalar(stepIndex) || ~isnumeric(stepIndex) || ~isfinite(stepIndex) || stepIndex < 0 || stepIndex ~= fix(stepIndex)
    error('RTSPHEM:Precipitate:InvalidSnapshotStep', ...
        'stepIndex must be a finite nonnegative integer.');
end
if ~isscalar(timeSeconds) || ~isnumeric(timeSeconds) || ~isfinite(timeSeconds) || timeSeconds < 0
    error('RTSPHEM:Precipitate:InvalidSnapshotTime', ...
        'timeSeconds must be a finite nonnegative scalar.');
end
coordV = grid.coordV;
V0T = grid.V0T;
if ~isnumeric(coordV) || size(coordV, 2) ~= 2 || isempty(coordV) || any(~isfinite(coordV(:)))
    error('RTSPHEM:Precipitate:InvalidSnapshotGrid', ...
        'grid.coordV must be a finite N x 2 numeric array.');
end
if ~isnumeric(V0T) || size(V0T, 2) ~= 3 || isempty(V0T) || any(V0T(:) < 1) || ...
        any(V0T(:) > size(coordV, 1)) || any(V0T(:) ~= fix(V0T(:)))
    error('RTSPHEM:Precipitate:InvalidSnapshotGrid', ...
        'grid.V0T must contain valid triangle vertex indices.');
end
validateLevelSetVector(currentLevelSet, size(coordV, 1), 'currentLevelSet');
validateLevelSetVector(initialLevelSet, size(coordV, 1), 'initialLevelSet');
end

function validateLevelSetVector(levelSet, numVertices, name)
if ~isnumeric(levelSet) || ~isreal(levelSet) || numel(levelSet) ~= numVertices || any(~isfinite(levelSet(:)))
    error('RTSPHEM:Precipitate:InvalidSnapshotLevelSet', ...
        '%s must be a finite real vector with one value per vertex.', name);
end
end

function writeSnapshotPng(outputPath, targetSeconds, stepIndex, timeSeconds, grid, ...
    currentLevelSet, currentSegments, initialSegments, config)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 650]);
cleanupObj = onCleanup(@() closeFigureIfOpen(fig));
ax = axes(fig);
hold(ax, 'on');

solidFaceData = mean(reshape(currentLevelSet(grid.V0T), size(grid.V0T)), 2) >= 0;
hSolid = patch(ax, 'Faces', grid.V0T, 'Vertices', grid.coordV, ...
    'FaceVertexCData', double(solidFaceData), 'FaceColor', 'flat', ...
    'EdgeColor', [0.88, 0.88, 0.88], 'LineWidth', 0.25, ...
    'DisplayName', 'solid mask');
colormap(ax, [0.93, 0.97, 1.00; 0.84, 0.84, 0.84]);

hInitial = plot(ax, NaN, NaN, 'Color', [0.15, 0.15, 0.15], ...
    'LineStyle', '--', 'LineWidth', 1.0, 'DisplayName', 'initial interface');
hCurrent = plot(ax, NaN, NaN, 'Color', [0.00, 0.32, 0.72], ...
    'LineStyle', '-', 'LineWidth', 1.8, 'DisplayName', 'current interface');
plotSegments(ax, initialSegments, [0.15, 0.15, 0.15], '--', 1.0);
plotSegments(ax, currentSegments, [0.00, 0.32, 0.72], '-', 1.8);

axis(ax, 'equal');
axis(ax, 'tight');
box(ax, 'on');
ax.XGrid = 'on';
ax.YGrid = 'on';
xlabel(ax, 'x [cm]');
ylabel(ax, 'y [cm]');
title(ax, sprintf('CaCO3 precipitation benchmark snapshot: target %.3g min, simulated %.3g min', ...
    targetSeconds / 60, timeSeconds / 60), 'Interpreter', 'none');

annotationText = buildAnnotationText(stepIndex, timeSeconds, config);
text(ax, 0.01, 0.99, annotationText, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'w', 'Margin', 5, 'EdgeColor', [0.7, 0.7, 0.7], ...
    'Interpreter', 'none');

legend(ax, [hSolid, hInitial, hCurrent], ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Interpreter', 'none');

print(fig, outputPath, '-dpng', '-r300');
end

function annotationText = buildAnnotationText(stepIndex, timeSeconds, config)
inletLabel = '';
if isstruct(config) && isfield(config, 'inletA') && isfield(config, 'inletB') && ...
        isstruct(config.inletA) && isstruct(config.inletB)
    nameA = fieldOrDefault(config.inletA, 'name', 'inlet A');
    nameB = fieldOrDefault(config.inletB, 'name', 'inlet B');
    inletLabel = sprintf('\nLower inlet: %s\nUpper inlet: %s', ...
        char(string(nameA)), char(string(nameB)));
end
annotationText = sprintf('step: %d\ntime: %.6g s%s', stepIndex, timeSeconds, inletLabel);
end

function value = fieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function plotSegments(ax, segments, color, lineStyle, lineWidth)
for iSegment = 1:numel(segments)
    segment = segments{iSegment};
    plot(ax, segment(:, 1), segment(:, 2), ...
        'Color', color, 'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
        'HandleVisibility', 'off');
end
end

function segments = extractZeroLevelSegments(grid, levels)
edgePairs = [1, 2; 2, 3; 3, 1];
segments = {};
for iTriangle = 1:size(grid.V0T, 1)
    vertexIds = grid.V0T(iTriangle, :);
    points = grid.coordV(vertexIds, :);
    values = levels(vertexIds);
    crossings = zeros(0, 2);
    for iEdge = 1:size(edgePairs, 1)
        edge = edgePairs(iEdge, :);
        p1 = points(edge(1), :);
        p2 = points(edge(2), :);
        v1 = values(edge(1));
        v2 = values(edge(2));
        if v1 == 0 && v2 == 0
            crossings = [crossings; p1; p2]; %#ok<AGROW>
        elseif v1 == 0
            crossings = [crossings; p1]; %#ok<AGROW>
        elseif v2 == 0
            crossings = [crossings; p2]; %#ok<AGROW>
        elseif v1 * v2 < 0
            fraction = abs(v1) / (abs(v1) + abs(v2));
            crossings = [crossings; p1 + fraction * (p2 - p1)]; %#ok<AGROW>
        end
    end
    if size(crossings, 1) >= 2
        crossings = unique(round(crossings * 1e12) / 1e12, 'rows', 'stable');
        if size(crossings, 1) >= 2
            segments{end + 1} = crossings(1:2, :); %#ok<AGROW>
        end
    end
end
end

function closeFigureIfOpen(fig)
if ishghandle(fig)
    close(fig);
end
end
