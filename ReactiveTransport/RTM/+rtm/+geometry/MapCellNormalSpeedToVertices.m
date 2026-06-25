function vertexSpeed = MapCellNormalSpeedToVertices(mesh, geometry, cellNormalSpeed, options)
%MAPCELLNORMALSPEEDTOVERTICES Extend interface-normal speeds to mesh vertices.
%
% Interface cell speeds are first averaged onto incident vertices, then any
% remaining vertices are filled from the nearest active interface centroid.

if nargin < 4 || isempty(options)
    options = struct();
end

[verticesCm, triangles] = unpackMesh(mesh);
numVertices = size(verticesCm, 1);
numCells = size(triangles, 1);
cellNormalSpeed = cellNormalSpeed(:);
if numel(cellNormalSpeed) ~= numCells
    error('RTSPHEM:Geometry:SpeedSizeMismatch', ...
        'cellNormalSpeed must contain one value per mesh cell.');
end

interfaceArea = requireGeometryVector(geometry, 'interface_area_cm2', numCells);
active = interfaceArea > 0 & isfinite(cellNormalSpeed) & cellNormalSpeed > 0;
vertexSpeed = zeros(numVertices, 1);
vertexWeight = zeros(numVertices, 1);
for iCell = find(active).'
    vertexIds = triangles(iCell, :);
    weight = interfaceArea(iCell);
    vertexSpeed(vertexIds) = vertexSpeed(vertexIds) + cellNormalSpeed(iCell) .* weight;
    vertexWeight(vertexIds) = vertexWeight(vertexIds) + weight;
end

hasIncidentSpeed = vertexWeight > 0;
vertexSpeed(hasIncidentSpeed) = vertexSpeed(hasIncidentSpeed) ./ vertexWeight(hasIncidentSpeed);
if all(~active)
    return;
end

fillMode = string(getFieldOrDefault(options, 'vertexSpeedExtension', 'nearest_interface'));
switch fillMode
    case "zero"
        return;
    case "nearest_interface"
        vertexSpeed(~hasIncidentSpeed) = nearestInterfaceSpeed( ...
            verticesCm(~hasIncidentSpeed, :), geometry, cellNormalSpeed, active);
    case "mean_interface"
        meanSpeed = sum(cellNormalSpeed(active) .* interfaceArea(active)) ./ ...
            sum(interfaceArea(active));
        vertexSpeed(~hasIncidentSpeed) = meanSpeed;
    otherwise
        error('RTSPHEM:Geometry:UnsupportedVertexSpeedExtension', ...
            'Unsupported vertexSpeedExtension: %s.', fillMode);
end
end

function filled = nearestInterfaceSpeed(queryPoints, geometry, cellNormalSpeed, active)
filled = zeros(size(queryPoints, 1), 1);
if isempty(queryPoints)
    return;
end
centroids = getFieldOrDefault(geometry, 'interface_centroid_cm', []);
if isempty(centroids) || size(centroids, 1) ~= numel(cellNormalSpeed) || ...
        size(centroids, 2) ~= 2
    activeSpeeds = cellNormalSpeed(active);
    filled(:) = mean(activeSpeeds);
    return;
end
activeCentroids = centroids(active, :);
activeSpeeds = cellNormalSpeed(active);
valid = all(isfinite(activeCentroids), 2);
activeCentroids = activeCentroids(valid, :);
activeSpeeds = activeSpeeds(valid);
if isempty(activeSpeeds)
    filled(:) = mean(cellNormalSpeed(active));
    return;
end
for iPoint = 1:size(queryPoints, 1)
    delta = activeCentroids - queryPoints(iPoint, :);
    [~, nearestIndex] = min(sum(delta .^ 2, 2));
    filled(iPoint) = activeSpeeds(nearestIndex);
end
end

function values = requireGeometryVector(geometry, fieldName, expectedLength)
if ~isstruct(geometry) || ~isfield(geometry, fieldName)
    error('RTSPHEM:Geometry:MissingGeometryField', ...
        'geometry.%s is required.', fieldName);
end
values = geometry.(fieldName)(:);
if numel(values) ~= expectedLength
    error('RTSPHEM:Geometry:GeometrySizeMismatch', ...
        'geometry.%s must contain one value per mesh cell.', fieldName);
end
values(~isfinite(values)) = 0;
values = max(values, 0);
end

function [verticesCm, triangles] = unpackMesh(mesh)
if isstruct(mesh) && isfield(mesh, 'rectGrid')
    mesh = mesh.rectGrid;
end
if isstruct(mesh) && isfield(mesh, 'vertices_cm')
    verticesCm = mesh.vertices_cm;
elseif isstruct(mesh) && isfield(mesh, 'coordV')
    verticesCm = mesh.coordV;
elseif isobject(mesh) && isprop(mesh, 'coordV')
    verticesCm = mesh.coordV;
else
    error('RTSPHEM:Geometry:MissingVertices', ...
        'mesh must contain vertices_cm or coordV.');
end

if isstruct(mesh) && isfield(mesh, 'triangles')
    triangles = mesh.triangles;
elseif isstruct(mesh) && isfield(mesh, 'V0T')
    triangles = mesh.V0T;
elseif isobject(mesh) && isprop(mesh, 'V0T')
    triangles = mesh.V0T;
else
    error('RTSPHEM:Geometry:MissingTriangles', ...
        'mesh must contain triangles or V0T.');
end

if size(verticesCm, 2) ~= 2 || size(triangles, 2) ~= 3
    error('RTSPHEM:Geometry:InvalidMesh', ...
        'Only two-dimensional triangular meshes are supported.');
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
