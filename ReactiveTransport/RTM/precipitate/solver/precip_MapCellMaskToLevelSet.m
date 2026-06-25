function mapped = precip_MapCellMaskToLevelSet(state, spec, options)
% precip_MapCellMaskToLevelSet - Map Yoon substrate/Vm masks to a level set.
%
% The sign convention matches the legacy HyPHM level-set path: positive values
% are mobile fluid cells and negative values are solid/blocked cells.

if nargin < 3
    options = struct();
end

substrateMask = logicalMask(state, 'substrateMask');
blockedMask = logicalMask(state, 'blockedMask');
vmThreshold = getFieldOrDefault(spec, 'blockedVmThreshold', 0.6);
if isfield(options, 'blockedVmThreshold') && ~isempty(options.blockedVmThreshold)
    vmThreshold = options.blockedVmThreshold;
end
vmBlockedMask = false(size(substrateMask));
if isfield(state, 'Vm') && ~isempty(state.Vm)
    vmBlockedMask = state.Vm >= vmThreshold;
end

solidMask = substrateMask | blockedMask | vmBlockedMask;
[xCenters, yCenters] = cellCenters(state, spec);
[xGrid, yGrid] = meshgrid(xCenters, yCenters);
levelSet = signedDistanceLikeField(xGrid, yGrid, solidMask, spec);

mapped = struct();
mapped.levelSetAtCellCenters_cm = levelSet;
mapped.solidMask = solidMask;
mapped.fluidMask = ~solidMask;
mapped.substrateMask = substrateMask;
mapped.blockedMask = blockedMask | vmBlockedMask;
mapped.xCenters_cm = xCenters;
mapped.yCenters_cm = yCenters;
mapped.cellCenters_cm = [xGrid(:), yGrid(:)];
mapped.blockedVmThreshold = vmThreshold;
mapped.levelSetSignConvention = 'positive_fluid_negative_solid';
mapped.source = 'yoon_cell_mask_to_level_set';
end

function levelSet = signedDistanceLikeField(xGrid, yGrid, solidMask, spec)
spacing = min(getFieldOrDefault(spec, 'dx_cm', 1), ...
    getFieldOrDefault(spec, 'dy_cm', 1));
levelSet = ones(size(solidMask)) .* (0.5 * spacing);
if ~any(solidMask(:))
    return;
end
if all(solidMask(:))
    levelSet(:) = -0.5 * spacing;
    return;
end

solidPoints = [xGrid(solidMask), yGrid(solidMask)];
fluidPoints = [xGrid(~solidMask), yGrid(~solidMask)];
for iCell = 1:numel(solidMask)
    point = [xGrid(iCell), yGrid(iCell)];
    if solidMask(iCell)
        levelSet(iCell) = -nearestDistance(point, fluidPoints, spacing);
    else
        levelSet(iCell) = nearestDistance(point, solidPoints, spacing);
    end
end
end

function distance = nearestDistance(point, points, fallback)
if isempty(points)
    distance = 0.5 * fallback;
    return;
end
d = hypot(points(:, 1) - point(1), points(:, 2) - point(2));
distance = max(min(d), 0.5 * fallback);
end

function [xCenters, yCenters] = cellCenters(state, spec)
if isfield(state, 'grid') && isfield(state.grid, 'xCenters_cm') && ...
        isfield(state.grid, 'yCenters_cm')
    xCenters = state.grid.xCenters_cm(:)';
    yCenters = state.grid.yCenters_cm(:);
else
    dx = getFieldOrDefault(spec, 'dx_cm', ...
        getFieldOrDefault(spec, 'lengthXAxis_cm', 1) / spec.numX);
    dy = getFieldOrDefault(spec, 'dy_cm', ...
        getFieldOrDefault(spec, 'lengthYAxis_cm', 1) / spec.numY);
    xCenters = ((1:spec.numX) - 0.5) .* dx;
    yCenters = ((1:spec.numY)' - 0.5) .* dy;
end
end

function mask = logicalMask(state, fieldName)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    mask = logical(state.(fieldName));
elseif isfield(state, 'Vm')
    mask = false(size(state.Vm));
else
    error('RTSPHEM:Precipitate:MissingVm', ...
        'state.Vm is required when %s is absent.', fieldName);
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
