function [newState, ledger] = ConservativeRemap(oldState, oldGeometry, newGeometry, options)
%CONSERVATIVEREMAP Remap conserved component moles between cut-cell geometries.
%
% options.overlap_volume_cm3 is an nOld-by-nNew matrix of old water volume
% portions that overlap new cells. Component moles are distributed according
% to each old cell's overlap fractions.

if nargin < 4 || isempty(options)
    options = struct();
end

rtm.state.ValidateState(oldState);
oldWater = requireWaterVolume(oldGeometry, 'oldGeometry');
newWater = requireWaterVolume(newGeometry, 'newGeometry');
numOld = numel(oldWater);
numNew = numel(newWater);
if size(oldState.component_moles, 1) ~= numOld
    error('RTSPHEM:Geometry:StateGeometrySizeMismatch', ...
        'oldState rows must match oldGeometry.water_volume_cm3.');
end

if ~(isfield(options, 'overlap_volume_cm3') && ~isempty(options.overlap_volume_cm3)) && ...
        numOld == numNew
    [newState, ledger] = sameGridDefaultRemap( ...
        oldState, oldWater, newWater, options);
    return;
end

overlapVolume = resolveOverlapVolume(options, oldWater, newWater);
transferWeights = zeros(size(overlapVolume));
activeOld = oldWater > 0;
transferWeights(activeOld, :) = overlapVolume(activeOld, :) ./ oldWater(activeOld);
transferWeights(~isfinite(transferWeights)) = 0;
transferWeights = max(transferWeights, 0);

newComponentMoles = transferWeights.' * oldState.component_moles;
newState = oldState;
newState.component_moles = newComponentMoles;
newState.temperature_C = remapCellVector(oldState.temperature_C(:), transferWeights, numNew, 25);
newState.pressure_atm = remapCellVector(oldState.pressure_atm(:), transferWeights, numNew, 1);
newState.mineral_moles = resolveNewMineralMoles(oldState, transferWeights, numNew, options);

rtm.state.ValidateState(newState);
ledger = buildLedger(oldState, newState, oldWater, newWater, overlapVolume, transferWeights);
end

function waterVolume = requireWaterVolume(geometry, label)
if ~isstruct(geometry) || ~isfield(geometry, 'water_volume_cm3')
    error('RTSPHEM:Geometry:MissingWaterVolume', ...
        '%s.water_volume_cm3 is required.', label);
end
waterVolume = geometry.water_volume_cm3(:);
if any(~isfinite(waterVolume))
    error('RTSPHEM:Geometry:InvalidWaterVolume', ...
        '%s.water_volume_cm3 must contain finite values.', label);
end
if any(waterVolume < 0)
    error('RTSPHEM:Geometry:NegativeWaterVolume', ...
        '%s.water_volume_cm3 must be nonnegative.', label);
end
end

function overlapVolume = resolveOverlapVolume(options, oldWater, newWater)
numOld = numel(oldWater);
numNew = numel(newWater);
if isfield(options, 'overlap_volume_cm3') && ~isempty(options.overlap_volume_cm3)
    overlapVolume = options.overlap_volume_cm3;
else
    if numOld ~= numNew
        error('RTSPHEM:Geometry:MissingOverlapVolume', ...
            'options.overlap_volume_cm3 is required when cell counts differ.');
    end
    overlapVolume = diag(min(oldWater, newWater));
end
overlapVolume = overlapVolume(:,:);
if ~isequal(size(overlapVolume), [numOld, numNew])
    error('RTSPHEM:Geometry:OverlapSizeMismatch', ...
        'overlap_volume_cm3 must have size nOld-by-nNew.');
end
if any(~isfinite(overlapVolume(:))) || any(overlapVolume(:) < 0)
    error('RTSPHEM:Geometry:InvalidOverlapVolume', ...
        'overlap_volume_cm3 must contain nonnegative finite values.');
end
rowOverlap = sum(overlapVolume, 2);
if any(rowOverlap > oldWater + max(1e-14, 1e-12 .* max(oldWater, 1)))
    error('RTSPHEM:Geometry:InvalidOverlapVolume', ...
        'overlap_volume_cm3 cannot exceed old water volume by row.');
end
end

function values = remapCellVector(oldValues, transferWeights, numNew, fallbackValue)
weightedSum = transferWeights.' * oldValues;
weightTotal = sum(transferWeights, 1).';
values = fallbackValue .* ones(numNew, 1);
active = weightTotal > 0;
values(active) = weightedSum(active) ./ weightTotal(active);
values(~isfinite(values)) = fallbackValue;
end

function mineralMoles = resolveNewMineralMoles(oldState, transferWeights, numNew, options)
if isfield(options, 'new_mineral_moles') && ~isempty(options.new_mineral_moles)
    mineralMoles = options.new_mineral_moles(:,:);
else
    mineralMoles = transferWeights.' * oldState.mineral_moles;
end
if size(mineralMoles, 1) ~= numNew || size(mineralMoles, 2) ~= size(oldState.mineral_moles, 2)
    error('RTSPHEM:Geometry:MineralSizeMismatch', ...
        'new mineral moles must match the new cell count and mineral count.');
end
if any(~isfinite(mineralMoles(:)))
    error('RTSPHEM:Geometry:NonfiniteMineralMoles', ...
        'new mineral moles must be finite.');
end
tolerance = max(1e-14, 1e-12 .* max(abs(mineralMoles), 1));
if any(mineralMoles(:) < -tolerance(:))
    error('RTSPHEM:Geometry:NegativeMineralMoles', ...
        'new mineral moles must be nonnegative.');
end
mineralMoles(abs(mineralMoles) <= tolerance) = 0;
end

function ledger = buildLedger(oldState, newState, oldWater, newWater, overlapVolume, transferWeights)
oldTotals = sum(oldState.component_moles, 1);
newTotals = sum(newState.component_moles, 1);
assignedOldMoles = sum((sum(transferWeights, 2) .* oldState.component_moles), 1);
unassignedComponentMoles = oldTotals - assignedOldMoles;
componentResidual = newTotals - oldTotals;

newWaterWithoutOverlap = find(newWater > 0 & sum(overlapVolume, 1).' <= 0);
oldWaterUnmapped = find(oldWater > 0 & sum(overlapVolume, 2) <= 0);

ledger = struct();
ledger.initial_component_moles_total = oldTotals;
ledger.final_component_moles_total = newTotals;
ledger.assigned_component_moles = assignedOldMoles;
ledger.unassigned_component_moles = unassignedComponentMoles;
ledger.component_residual_moles = componentResidual;
ledger.max_abs_component_residual_moles = max(abs(componentResidual));
ledger.new_water_without_overlap_cells = newWaterWithoutOverlap(:).';
ledger.old_water_unmapped_cells = oldWaterUnmapped(:).';
ledger.overlap_volume_cm3 = overlapVolume;
end

function [newState, ledger] = sameGridDefaultRemap(oldState, oldWater, newWater, options)
numCells = numel(oldWater);
overlapDiagonal = min(oldWater, newWater);
fraction = zeros(numCells, 1);
activeOld = oldWater > 0;
fraction(activeOld) = overlapDiagonal(activeOld) ./ oldWater(activeOld);
fraction(~isfinite(fraction)) = 0;
fraction = max(fraction, 0);

newState = oldState;
newState.component_moles = oldState.component_moles .* fraction;
newState.temperature_C = oldState.temperature_C(:);
newState.pressure_atm = oldState.pressure_atm(:);
newState.mineral_moles = resolveNewMineralMolesForSameGrid( ...
    oldState, fraction, numCells, options);

rtm.state.ValidateState(newState);
overlapVolume = sparse(1:numCells, 1:numCells, overlapDiagonal, numCells, numCells);
ledger = buildSameGridLedger(oldState, newState, oldWater, newWater, ...
    overlapVolume, fraction);
end

function mineralMoles = resolveNewMineralMolesForSameGrid( ...
        oldState, fraction, numCells, options)
if isfield(options, 'new_mineral_moles') && ~isempty(options.new_mineral_moles)
    mineralMoles = options.new_mineral_moles(:,:);
else
    mineralMoles = oldState.mineral_moles .* fraction;
end
if size(mineralMoles, 1) ~= numCells || size(mineralMoles, 2) ~= size(oldState.mineral_moles, 2)
    error('RTSPHEM:Geometry:MineralSizeMismatch', ...
        'new mineral moles must match the new cell count and mineral count.');
end
if any(~isfinite(mineralMoles(:)))
    error('RTSPHEM:Geometry:NonfiniteMineralMoles', ...
        'new mineral moles must be finite.');
end
tolerance = max(1e-14, 1e-12 .* max(abs(mineralMoles), 1));
if any(mineralMoles(:) < -tolerance(:))
    error('RTSPHEM:Geometry:NegativeMineralMoles', ...
        'new mineral moles must be nonnegative.');
end
mineralMoles(abs(mineralMoles) <= tolerance) = 0;
end

function ledger = buildSameGridLedger( ...
        oldState, newState, oldWater, newWater, overlapVolume, fraction)
oldTotals = sum(oldState.component_moles, 1);
newTotals = sum(newState.component_moles, 1);
assignedOldMoles = sum(fraction .* oldState.component_moles, 1);
unassignedComponentMoles = oldTotals - assignedOldMoles;
componentResidual = newTotals - oldTotals;

newWaterWithoutOverlap = find(newWater > 0 & diag(overlapVolume) <= 0);
oldWaterUnmapped = find(oldWater > 0 & fraction <= 0);

ledger = struct();
ledger.initial_component_moles_total = oldTotals;
ledger.final_component_moles_total = newTotals;
ledger.assigned_component_moles = assignedOldMoles;
ledger.unassigned_component_moles = unassignedComponentMoles;
ledger.component_residual_moles = componentResidual;
ledger.max_abs_component_residual_moles = max(abs(componentResidual));
ledger.new_water_without_overlap_cells = newWaterWithoutOverlap(:).';
ledger.old_water_unmapped_cells = oldWaterUnmapped(:).';
ledger.overlap_volume_cm3 = overlapVolume;
end
