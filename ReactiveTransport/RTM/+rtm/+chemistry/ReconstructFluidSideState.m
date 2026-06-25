function interfaceState = ReconstructFluidSideState(state, geometry, connectivity, options)
%RECONSTRUCTFLUIDSIDESTATE Reconstruct fluid-side component concentrations.
%
% Concentrations are derived from conserved component moles. For each
% interface cell, a one-sided linear least-squares reconstruction is limited
% to the source-cell/neighbour stencil bounds. Degenerate stencils fall back
% to the source cell average.

if nargin < 4 || isempty(options)
    options = struct();
end
rtm.state.ValidateState(state);
[waterVolume, cellCentroid, interfaceCentroid, interfaceArea, interfaceNormal] = ...
    validateGeometry(geometry, state);
neighbors = getNeighbors(connectivity, size(state.component_moles, 1));
useLimiter = logical(getFieldOrDefault(options, 'useLimiter', true));

concentration = rtm.state.ConvertMolesConcentrations(state, geometry, 'moles_to_concentration');
numCells = size(state.component_moles, 1);
numComponents = size(state.component_moles, 2);

values = nan(numCells, numComponents);
sourceCell = nan(numCells, 1);
qualityFlag = strings(numCells, 1);
qualityFlag(:) = "inactive";
fallbackUsed = false(numCells, 1);

for iCell = 1:numCells
    if waterVolume(iCell) <= 0 || interfaceArea(iCell) <= 0 || ...
            any(~isfinite(interfaceCentroid(iCell, :)))
        continue;
    end

    sourceCell(iCell) = iCell;
    fluidSideNeighbors = filterFluidSideNeighbors( ...
        neighbors{iCell}(:), waterVolume, cellCentroid, interfaceCentroid(iCell, :), ...
        interfaceNormal(iCell, :));
    stencil = unique([iCell; fluidSideNeighbors(:)]);
    [reconstructed, usedFallback] = reconstructCell( ...
        iCell, stencil, concentration, cellCentroid, interfaceCentroid(iCell, :), useLimiter);
    values(iCell, :) = reconstructed;
    fallbackUsed(iCell) = usedFallback;
    if usedFallback
        qualityFlag(iCell) = "cell_average";
    else
        qualityFlag(iCell) = "linear_limited";
    end
end

interfaceState = struct();
interfaceState.component_names = state.component_names;
interfaceState.component_concentration_mol_cm3 = values;
interfaceState.source_cell = sourceCell;
interfaceState.quality_flag = qualityFlag;
interfaceState.fallback_used = fallbackUsed;
end

function [reconstructed, usedFallback] = reconstructCell(sourceCell, stencil, ...
    concentration, cellCentroid, interfacePoint, useLimiter)
sourceValue = concentration(sourceCell, :);
usedFallback = true;
reconstructed = sourceValue;

if numel(stencil) < 3
    return;
end

delta = cellCentroid(stencil, :) - cellCentroid(sourceCell, :);
if rank(delta, 1e-12) < 2
    return;
end

targetDelta = interfacePoint - cellCentroid(sourceCell, :);
numComponents = size(concentration, 2);
candidate = zeros(1, numComponents);
for iComponent = 1:numComponents
    rhs = concentration(stencil, iComponent) - sourceValue(iComponent);
    gradient = delta \ rhs;
    candidate(iComponent) = sourceValue(iComponent) + targetDelta * gradient;
end

candidate(~isfinite(candidate)) = sourceValue(~isfinite(candidate));
candidate = max(candidate, 0);
if useLimiter
    stencilMin = min(concentration(stencil, :), [], 1);
    stencilMax = max(concentration(stencil, :), [], 1);
    candidate = min(max(candidate, stencilMin), stencilMax);
end
reconstructed = candidate;
usedFallback = false;
end

function fluidSideNeighbors = filterFluidSideNeighbors( ...
        neighborIds, waterVolume, cellCentroid, interfacePoint, normal)
neighborIds = neighborIds(waterVolume(neighborIds) > 0);
fluidSideNeighbors = neighborIds;
if isempty(neighborIds) || any(~isfinite(normal)) || norm(normal) <= eps
    return;
end
normal = normal ./ norm(normal);
relativePosition = cellCentroid(neighborIds, :) - interfacePoint;
fluidSide = relativePosition * normal(:) > 0;
fluidSideNeighbors = neighborIds(fluidSide);
end

function [waterVolume, cellCentroid, interfaceCentroid, interfaceArea, interfaceNormal] = ...
    validateGeometry(geometry, state)
numCells = size(state.component_moles, 1);
requiredFields = {'water_volume_cm3', 'cell_centroid_cm', ...
    'interface_centroid_cm', 'interface_area_cm2'};
for iField = 1:numel(requiredFields)
    if ~isfield(geometry, requiredFields{iField})
        error('RTSPHEM:Chemistry:MissingGeometryField', ...
            'geometry.%s is required.', requiredFields{iField});
    end
end
waterVolume = geometry.water_volume_cm3(:);
cellCentroid = geometry.cell_centroid_cm;
interfaceCentroid = geometry.interface_centroid_cm;
interfaceArea = geometry.interface_area_cm2(:);
interfaceNormal = nan(numCells, 2);
if isfield(geometry, 'interface_normal') && ~isempty(geometry.interface_normal)
    interfaceNormal = geometry.interface_normal;
end
if numel(waterVolume) ~= numCells || numel(interfaceArea) ~= numCells || ...
        size(cellCentroid, 1) ~= numCells || size(interfaceCentroid, 1) ~= numCells || ...
        size(cellCentroid, 2) ~= 2 || size(interfaceCentroid, 2) ~= 2 || ...
        size(interfaceNormal, 1) ~= numCells || size(interfaceNormal, 2) ~= 2
    error('RTSPHEM:Chemistry:GeometrySizeMismatch', ...
        'Geometry fields must match state cell count and use 2D centroids.');
end
if any(~isfinite(waterVolume)) || any(~isfinite(interfaceArea)) || ...
        any(~isfinite(cellCentroid(:)))
    error('RTSPHEM:Chemistry:InvalidGeometryMeasure', ...
        'Geometry volume, area, and cell centroid fields must be finite.');
end
if any(waterVolume < 0) || any(interfaceArea < 0)
    error('RTSPHEM:Chemistry:NegativeGeometryMeasure', ...
        'Geometry water volume and interface area fields must be nonnegative.');
end
end

function neighbors = getNeighbors(connectivity, numCells)
if ~isstruct(connectivity) || ~isfield(connectivity, 'cell_neighbors')
    neighbors = repmat({zeros(0, 1)}, numCells, 1);
    return;
end
neighbors = connectivity.cell_neighbors;
if numel(neighbors) ~= numCells
    error('RTSPHEM:Chemistry:ConnectivitySizeMismatch', ...
        'connectivity.cell_neighbors must contain one entry per cell.');
end
neighbors = neighbors(:);
for iCell = 1:numCells
    ids = neighbors{iCell}(:);
    if any(ids < 1) || any(ids > numCells) || any(ids ~= round(ids))
        error('RTSPHEM:Chemistry:InvalidConnectivity', ...
            'cell_neighbors contains invalid cell indices.');
    end
    neighbors{iCell} = ids;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
