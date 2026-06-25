function geometry = BuildCutCellMetrics(mesh, levelSet, options)
%BUILDCUTCELLMETRICS Compute cut-cell water/solid measures from a level set.
%
% Negative level-set values are treated as fluid/pore, matching PNM_beauty3.

if nargin < 3 || isempty(options)
    options = struct();
end

thicknessCm = getOption(options, 'thickness_cm', 1);
if ~(isscalar(thicknessCm) && isfinite(thicknessCm) && thicknessCm > 0)
    error('RTSPHEM:Geometry:InvalidThickness', ...
        'options.thickness_cm must be a positive finite scalar.');
end

[verticesCm, triangles, cellAreaCm2] = unpackMesh(mesh);
levelSet = levelSet(:);
numVertices = size(verticesCm, 1);
if numel(levelSet) ~= numVertices
    error('RTSPHEM:Geometry:LevelSetSizeMismatch', ...
        'levelSet must contain one value per mesh vertex.');
end

numCells = size(triangles, 1);
cellVolumeCm3 = cellAreaCm2(:) .* thicknessCm;
cellCentroidCm = zeros(numCells, 2);
waterVolumeCm3 = zeros(numCells, 1);
solidVolumeCm3 = zeros(numCells, 1);
fluidFraction = zeros(numCells, 1);
interfaceAreaCm2 = zeros(numCells, 1);
interfaceCentroidCm = nan(numCells, 2);
interfaceNormal = nan(numCells, 2);
interfaceHCm = zeros(numCells, 1);
interfaceLengthScaleCm = zeros(numCells, 1);
activeFluidCell = false(numCells, 1);
cutCell = false(numCells, 1);

for iCell = 1:numCells
    vertexIds = triangles(iCell, :);
    points = verticesCm(vertexIds, :);
    values = levelSet(vertexIds);
    totalArea = cellAreaCm2(iCell);
    cellCentroidCm(iCell, :) = mean(points, 1);

    waterArea = clipTriangleAreaBelowZero(points, values);
    waterArea = min(max(waterArea, 0), totalArea);
    solidArea = max(totalArea - waterArea, 0);

    [interfaceLength, segmentPoints] = triangleInterfaceSegment(points, values);
    waterVolumeCm3(iCell) = waterArea * thicknessCm;
    solidVolumeCm3(iCell) = solidArea * thicknessCm;
    interfaceAreaCm2(iCell) = interfaceLength * thicknessCm;
    fluidFraction(iCell) = safeDivide(waterArea, totalArea);
    activeFluidCell(iCell) = waterArea > 0;
    cutCell(iCell) = waterArea > 0 && solidArea > 0 && interfaceLength > 0;

    if size(segmentPoints, 1) >= 2
        interfaceCentroidCm(iCell, :) = mean(segmentPoints, 1);
        interfaceHCm(iCell) = interfaceLength;
        interfaceNormal(iCell, :) = linearLevelSetNormal(points, values);
    end
end

activeInterface = interfaceAreaCm2 > 0;
interfaceLengthScaleCm(activeInterface) = min( ...
    sqrt(waterVolumeCm3(activeInterface) ./ max(interfaceAreaCm2(activeInterface), eps)), ...
    sqrt(cellVolumeCm3(activeInterface)));

geometry = struct();
geometry.cell_volume_cm3 = cellVolumeCm3;
geometry.cell_centroid_cm = cellCentroidCm;
geometry.water_volume_cm3 = waterVolumeCm3;
geometry.solid_volume_cm3 = solidVolumeCm3;
geometry.fluid_fraction = fluidFraction;
geometry.interface_area_cm2 = interfaceAreaCm2;
geometry.interface_centroid_cm = interfaceCentroidCm;
geometry.interface_normal = interfaceNormal;
geometry.interface_h_cm = interfaceHCm;
geometry.interface_length_scale_cm = interfaceLengthScaleCm;
geometry.active_fluid_cell = activeFluidCell;
geometry.cut_cell = cutCell;
geometry.level_set = levelSet;
end

function [verticesCm, triangles, cellAreaCm2] = unpackMesh(mesh)
if hasFieldOrProperty(mesh, 'rectGrid')
    rectGrid = getFieldOrProperty(mesh, 'rectGrid');
    if isstruct(rectGrid) || isobject(rectGrid)
        mesh = rectGrid;
    end
end

if hasFieldOrProperty(mesh, 'vertices_cm')
    verticesCm = getFieldOrProperty(mesh, 'vertices_cm');
elseif hasFieldOrProperty(mesh, 'coordV')
    verticesCm = getFieldOrProperty(mesh, 'coordV');
else
    error('RTSPHEM:Geometry:MissingVertices', ...
        'mesh must contain vertices_cm or coordV.');
end

if hasFieldOrProperty(mesh, 'triangles')
    triangles = getFieldOrProperty(mesh, 'triangles');
elseif hasFieldOrProperty(mesh, 'V0T')
    triangles = getFieldOrProperty(mesh, 'V0T');
else
    error('RTSPHEM:Geometry:MissingTriangles', ...
        'mesh must contain triangles or V0T.');
end

if size(verticesCm, 2) ~= 2 || size(triangles, 2) ~= 3
    error('RTSPHEM:Geometry:InvalidMesh', ...
        'Only two-dimensional triangular meshes are supported.');
end

if hasFieldOrProperty(mesh, 'areaT') && ...
        numel(getFieldOrProperty(mesh, 'areaT')) == size(triangles, 1)
    cellAreaCm2 = getFieldOrProperty(mesh, 'areaT');
    cellAreaCm2 = cellAreaCm2(:);
else
    cellAreaCm2 = triangleAreas(verticesCm, triangles);
end
cellAreaCm2 = max(cellAreaCm2(:), 0);
end

function tf = hasFieldOrProperty(value, name)
tf = (isstruct(value) && isfield(value, name)) || ...
    (isobject(value) && isprop(value, name));
end

function fieldValue = getFieldOrProperty(value, name)
fieldValue = value.(name);
end

function areas = triangleAreas(verticesCm, triangles)
numCells = size(triangles, 1);
areas = zeros(numCells, 1);
for iCell = 1:numCells
    points = verticesCm(triangles(iCell, :), :);
    areas(iCell) = polygonArea(points);
end
end

function area = clipTriangleAreaBelowZero(points, values)
polyPoints = points;
polyValues = values(:);
outputPoints = zeros(0, 2);
outputValues = zeros(0, 1);

for iPoint = 1:size(polyPoints, 1)
    jPoint = mod(iPoint, size(polyPoints, 1)) + 1;
    p1 = polyPoints(iPoint, :);
    p2 = polyPoints(jPoint, :);
    v1 = polyValues(iPoint);
    v2 = polyValues(jPoint);
    inside1 = v1 <= 0;
    inside2 = v2 <= 0;

    if inside1 && inside2
        outputPoints(end + 1, :) = p2; %#ok<AGROW>
        outputValues(end + 1, 1) = v2; %#ok<AGROW>
    elseif inside1 && ~inside2
        outputPoints(end + 1, :) = interpolateLevelZero(p1, p2, v1, v2); %#ok<AGROW>
        outputValues(end + 1, 1) = 0; %#ok<AGROW>
    elseif ~inside1 && inside2
        outputPoints(end + 1, :) = interpolateLevelZero(p1, p2, v1, v2); %#ok<AGROW>
        outputValues(end + 1, 1) = 0; %#ok<AGROW>
        outputPoints(end + 1, :) = p2; %#ok<AGROW>
        outputValues(end + 1, 1) = v2; %#ok<AGROW>
    end
end

if size(outputPoints, 1) < 3
    area = 0;
else
    area = polygonArea(outputPoints);
end
end

function [lengthValue, segmentPoints] = triangleInterfaceSegment(points, values)
if all(values <= 0) || all(values > 0)
    lengthValue = 0;
    segmentPoints = zeros(0, 2);
    return;
end

segmentPoints = zeros(0, 2);
edges = [1 2; 2 3; 3 1];
for iEdge = 1:3
    i1 = edges(iEdge, 1);
    i2 = edges(iEdge, 2);
    v1 = values(i1);
    v2 = values(i2);
    p1 = points(i1, :);
    p2 = points(i2, :);

    if v1 == 0 && v2 == 0
        segmentPoints(end + 1, :) = p1; %#ok<AGROW>
        segmentPoints(end + 1, :) = p2; %#ok<AGROW>
    elseif v1 == 0
        segmentPoints(end + 1, :) = p1; %#ok<AGROW>
    elseif v2 == 0
        segmentPoints(end + 1, :) = p2; %#ok<AGROW>
    elseif v1 * v2 < 0
        segmentPoints(end + 1, :) = interpolateLevelZero(p1, p2, v1, v2); %#ok<AGROW>
    end
end

segmentPoints = uniqueRowsWithTolerance(segmentPoints);
if size(segmentPoints, 1) < 2
    lengthValue = 0;
    return;
end

distances = pairwiseDistances(segmentPoints);
[lengthValue, linearIndex] = max(distances(:));
[rowIndex, colIndex] = ind2sub(size(distances), linearIndex);
segmentPoints = segmentPoints([rowIndex, colIndex], :);
end

function point = interpolateLevelZero(p1, p2, v1, v2)
denominator = v1 - v2;
if abs(denominator) < eps
    lambda = 0.5;
else
    lambda = v1 / denominator;
end
lambda = min(max(lambda, 0), 1);
point = p1 + lambda .* (p2 - p1);
end

function normal = linearLevelSetNormal(points, values)
x = points(:, 1);
y = points(:, 2);
A = [x, y, ones(3, 1)];
coefficients = A \ values(:);
gradient = coefficients(1:2).';
gradientNorm = norm(gradient);
if gradientNorm <= eps
    normal = [NaN, NaN];
else
    normal = -gradient ./ gradientNorm;
end
end

function area = polygonArea(points)
x = points(:, 1);
y = points(:, 2);
area = 0.5 * abs(sum(x .* y([2:end, 1]) - y .* x([2:end, 1])));
end

function rows = uniqueRowsWithTolerance(rows)
if isempty(rows)
    return;
end
scale = 1e14;
[~, uniqueIndices] = unique(round(rows .* scale) ./ scale, 'rows', 'stable');
rows = rows(sort(uniqueIndices), :);
end

function distances = pairwiseDistances(points)
numPoints = size(points, 1);
distances = zeros(numPoints, numPoints);
for iPoint = 1:numPoints
    delta = points - points(iPoint, :);
    distances(iPoint, :) = sqrt(sum(delta .^ 2, 2));
end
end

function value = safeDivide(numerator, denominator)
if denominator <= eps
    value = 0;
else
    value = numerator ./ denominator;
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
