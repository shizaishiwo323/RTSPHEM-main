function result = PhreeqcKineticsBackend(state, geometry, dtSeconds, options)
%PHREEQCKINETICSBACKEND PHREEQC database kinetics with conserved deltas.

if nargin < 4 || isempty(options)
    options = struct();
end
rtm.state.ValidateState(state);
validateGeometry(geometry, size(state.component_moles, 1));
if ~(isscalar(dtSeconds) && isfinite(dtSeconds) && dtSeconds >= 0)
    error('RTSPHEM:Chemistry:InvalidTimeStep', ...
        'dtSeconds must be a nonnegative finite scalar.');
end

numCells = size(state.component_moles, 1);
waterVolume = geometry.water_volume_cm3(:);
interfaceArea = geometry.interface_area_cm2(:);
[clusterIds, clusterInfo] = rtm.chemistry.ReactionClusterDiagnostics( ...
    options, numCells);
calciteIndex = find(strcmp(state.mineral_names, 'Calcite'), 1);
if isempty(calciteIndex)
    error('RTSPHEM:Chemistry:MissingCalcite', ...
        'phreeqc_kinetics requires mineral_names to include Calcite.');
end

hMolCm3 = requiredVectorOption(options, 'h_mol_cm3', numCells, ...
    'RTSPHEM:Chemistry:MissingPhreeqcHydrogenInput');
batchState = buildBatchState(state, geometry, hMolCm3);
batchOptions = options;
batchOptions.timeStepSize = dtSeconds;
batchOptions.rateLaw = char(getOption(options, 'rateLaw', 'database_calcite'));
runner = getRunBatchFunction(options);
batchResult = runner(batchState, batchOptions);

rawDissolved = realizedFromBatch(batchResult, numCells);
actualDissolved = min(max(rawDissolved, 0), max(state.mineral_moles(:, calciteIndex), 0));
componentDelta = componentDeltaFromBatch(state, batchState, batchResult, ...
    waterVolume, rawDissolved, actualDissolved);
[componentDelta, roundoffInfo] = suppressTinyNegativeComponentRoundoff( ...
    componentDelta, state.component_moles, options);
mineralDelta = zeros(size(state.mineral_moles));
mineralDelta(:, calciteIndex) = -actualDissolved;

result = struct();
result.component_delta_moles = componentDelta;
result.mineral_delta_moles = mineralDelta;
result.realized_interface_moles = actualDissolved;
result.candidate_interface_moles = rawDissolved;
result.interface_rate_mol_cm2_s = perAreaRate(actualDissolved, interfaceArea, dtSeconds);
result.candidate_interface_rate_mol_cm2_s = perAreaRate(rawDissolved, interfaceArea, dtSeconds);
result.inventory_limited = rawDissolved > state.mineral_moles(:, calciteIndex) + ...
    (1e-14 + 1e-12 .* max(rawDissolved, 1));
result.reactant_limited = false(numCells, 1);
result.pH = optionalResultVector(batchResult, 'pH', numCells, NaN);
result.saturation_index = optionalResultVector(batchResult, 'calciteSI', numCells, NaN);
result.charge_balance_residual_eq = optionalResultVector(batchResult, ...
    'chargeBalance', numCells, 0);
result.converged = true;
result.failed_cells = failedCellsFromBatch(batchResult);
result.error_message = string(getFieldOrDefault(batchResult, ...
    'errorMessage', ""));
result.cluster_ids = clusterIds;
result.aux = struct('chemistry_mode', "phreeqc_kinetics", ...
    'mineral_delta_source', "phreeqc", ...
    'reaction_cluster_count', clusterInfo.count, ...
    'reaction_cluster_max_membership', clusterInfo.max_membership, ...
    'reaction_cluster_overlapping_cell_count', clusterInfo.overlapping_cell_count, ...
    'component_delta_roundoff_suppressed_moles', ...
        roundoffInfo.suppressed_moles_total, ...
    'component_delta_roundoff_suppressed_entries', ...
        roundoffInfo.suppressed_entries);
end

function [componentDelta, info] = suppressTinyNegativeComponentRoundoff( ...
        componentDelta, componentMoles, options)
info = struct('suppressed_moles_total', 0, 'suppressed_entries', 0);
tolerance = getOption(options, 'componentNegativeRoundoffTolerance_mol', 1e-24);
if ~(isscalar(tolerance) && isfinite(tolerance) && tolerance >= 0)
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.componentNegativeRoundoffTolerance_mol must be a nonnegative finite scalar.');
end
updatedMoles = componentMoles + componentDelta;
roundoffMask = updatedMoles < 0 & abs(updatedMoles) <= tolerance;
if ~any(roundoffMask(:))
    return;
end
suppressed = -updatedMoles(roundoffMask);
componentDelta(roundoffMask) = -componentMoles(roundoffMask);
info.suppressed_moles_total = sum(suppressed, 'omitnan');
info.suppressed_entries = nnz(roundoffMask);
end

function failedCells = failedCellsFromBatch(batchResult)
failedCells = getFieldOrDefault(batchResult, 'failedCells', []);
failedCells = failedCells(:);
failedCells = failedCells(isfinite(failedCells));
end

function batchState = buildBatchState(state, geometry, hMolCm3)
waterVolume = geometry.water_volume_cm3(:);
batchState = struct();
batchState.h_mol_cm3 = hMolCm3(:);
batchState.ca_total_mol_cm3 = componentConcentration(state, 'Ca', waterVolume);
batchState.c_total_mol_cm3 = componentConcentration(state, 'C', waterVolume);
batchState.na_total_mol_cm3 = componentConcentration(state, 'Na', waterVolume);
batchState.cl_total_mol_cm3 = componentConcentration(state, 'Cl', waterVolume);
batchState.ca_mol_cm3 = batchState.ca_total_mol_cm3;
batchState.c_mol_cm3 = batchState.c_total_mol_cm3;
batchState.na_mol_cm3 = batchState.na_total_mol_cm3;
batchState.cl_mol_cm3 = batchState.cl_total_mol_cm3;
batchState.water_volume_cm3 = waterVolume;
batchState.reaction_water_volume_cm3 = waterVolume;
batchState.interface_area_cm2 = geometry.interface_area_cm2(:);
batchState.calcite_moles = mineralMoles(state, 'Calcite', size(state.component_moles, 1));
end

function values = componentConcentration(state, componentName, waterVolume)
idx = find(strcmp(state.component_names, componentName), 1);
values = zeros(size(waterVolume));
if isempty(idx)
    return;
end
activeWater = waterVolume > 0;
values(activeWater) = state.component_moles(activeWater, idx) ./ waterVolume(activeWater);
end

function values = mineralMoles(state, mineralName, numCells)
idx = find(strcmp(state.mineral_names, mineralName), 1);
if isempty(idx)
    values = zeros(numCells, 1);
else
    values = state.mineral_moles(:, idx);
end
end

function componentDelta = componentDeltaFromBatch(state, batchState, batchResult, ...
        waterVolume, rawDissolved, actualDissolved)
componentDelta = zeros(size(state.component_moles));
scale = ones(size(rawDissolved));
activeRaw = rawDissolved > 0;
scale(activeRaw) = actualDissolved(activeRaw) ./ rawDissolved(activeRaw);
scale(~isfinite(scale)) = 0;
for iComponent = 1:numel(state.component_names)
    componentName = char(state.component_names{iComponent});
    before = batchConcentration(batchState, componentName, numel(waterVolume));
    after = batchConcentration(batchResult, componentName, numel(waterVolume));
    delta = (after - before) .* waterVolume;
    if strcmp(componentName, 'Ca') || strcmp(componentName, 'C')
        delta = delta .* scale;
    end
    componentDelta(:, iComponent) = delta;
end
if any(~isfinite(componentDelta(:)))
    error('RTSPHEM:Chemistry:InvalidPhreeqcComponentDelta', ...
        'PHREEQC component deltas must be finite.');
end
end

function values = batchConcentration(batchStruct, componentName, numCells)
switch componentName
    case 'Ca'
        fieldNames = {'ca_total_mol_cm3', 'ca_mol_cm3'};
    case 'C'
        fieldNames = {'c_total_mol_cm3', 'c_mol_cm3'};
    case 'Na'
        fieldNames = {'na_total_mol_cm3', 'na_mol_cm3'};
    case 'Cl'
        fieldNames = {'cl_total_mol_cm3', 'cl_mol_cm3'};
    otherwise
        fieldNames = {};
end
values = zeros(numCells, 1);
for iField = 1:numel(fieldNames)
    if isfield(batchStruct, fieldNames{iField}) && ~isempty(batchStruct.(fieldNames{iField}))
        values = batchStruct.(fieldNames{iField})(:);
        return;
    end
end
end

function values = realizedFromBatch(batchResult, numCells)
if ~isfield(batchResult, 'calciteDissolvedMoles') || isempty(batchResult.calciteDissolvedMoles)
    error('RTSPHEM:Chemistry:MissingPhreeqcKineticDelta', ...
        'PHREEQC kinetics result must include calciteDissolvedMoles.');
end
values = batchResult.calciteDissolvedMoles(:);
if numel(values) ~= numCells
    error('RTSPHEM:Chemistry:PhreeqcResultSizeMismatch', ...
        'PHREEQC calciteDissolvedMoles must match the state cell count.');
end
end

function rate = perAreaRate(moles, interfaceArea, dtSeconds)
rate = zeros(size(moles));
active = interfaceArea > 0 & dtSeconds > 0;
rate(active) = moles(active) ./ interfaceArea(active) ./ dtSeconds;
rate(~isfinite(rate)) = 0;
end

function runner = getRunBatchFunction(options)
if isfield(options, 'runBatchFunction') && ~isempty(options.runBatchFunction)
    runner = options.runBatchFunction;
else
    runner = @RunPhreeqcCalciteBatch;
end
end

function values = requiredVectorOption(options, fieldName, numCells, errorId)
if ~isfield(options, fieldName) || isempty(options.(fieldName))
    error(errorId, 'options.%s is required.', fieldName);
end
values = options.(fieldName)(:);
if isscalar(values) && numCells > 1
    values = repmat(values, numCells, 1);
end
if numel(values) ~= numCells || any(~isfinite(values))
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.%s must be finite and match the state cell count.', fieldName);
end
if any(values < 0)
    error('RTSPHEM:Chemistry:NegativePhreeqcChemistryInput', ...
        'options.%s must be nonnegative.', fieldName);
end
end

function values = optionalResultVector(result, fieldName, numCells, defaultValue)
if isfield(result, fieldName) && ~isempty(result.(fieldName))
    values = result.(fieldName)(:);
else
    values = repmat(defaultValue, numCells, 1);
end
if isscalar(values) && numCells > 1
    values = repmat(values, numCells, 1);
end
if numel(values) ~= numCells
    values = repmat(defaultValue, numCells, 1);
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function validateGeometry(geometry, numCells)
requiredFields = {'water_volume_cm3', 'interface_area_cm2'};
for iField = 1:numel(requiredFields)
    fieldName = requiredFields{iField};
    if ~isstruct(geometry) || ~isfield(geometry, fieldName)
        error('RTSPHEM:Chemistry:MissingGeometryField', ...
            'geometry.%s is required.', fieldName);
    end
    if numel(geometry.(fieldName)) ~= numCells
        error('RTSPHEM:Chemistry:GeometrySizeMismatch', ...
            'geometry.%s must match the state cell count.', fieldName);
    end
    values = geometry.(fieldName)(:);
    if any(~isfinite(values))
        error('RTSPHEM:Chemistry:InvalidGeometryMeasure', ...
            'geometry.%s must contain finite values.', fieldName);
    end
    if any(values < 0)
        error('RTSPHEM:Chemistry:NegativeGeometryMeasure', ...
            'geometry.%s must be nonnegative.', fieldName);
    end
end
end
