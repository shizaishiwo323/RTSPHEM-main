function result = ExternalTstPhreeqcBackend(state, geometry, dtSeconds, options)
%EXTERNALTSTPHREEQCBACKEND External TST rate with PHREEQC aqueous closure.
%
% The transported state remains conserved component moles. Hydrogen activity
% used by the external TST rate is supplied in the chemistry input, so free H+
% is not required as a transported component.

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
        'external_tst_phreeqc requires mineral_names to include Calcite.');
end

hMolCm3 = requiredVectorOption(options, 'h_mol_cm3', numCells, ...
    'RTSPHEM:Chemistry:MissingPhreeqcHydrogenInput');
hActivityMolCm3 = optionalVectorOption(options, 'h_activity_mol_cm3', ...
    hMolCm3, numCells);
rateConstant = scalarOption(options, 'rate_constant_cm_s', 0);
if rateConstant < 0
    error('RTSPHEM:Chemistry:InvalidRateConstant', ...
        'options.rate_constant_cm_s must be nonnegative.');
end
maxMineralFraction = fractionOption(options, 'maxMineralFraction', Inf);

candidateMoles = hActivityMolCm3 .* 1000 .* rateConstant .* ...
    interfaceArea .* dtSeconds;
if any(~isfinite(candidateMoles))
    error('RTSPHEM:Chemistry:InvalidCandidateMoles', ...
        'external TST candidate moles must be finite.');
end
mineralMoles = state.mineral_moles(:, calciteIndex);
mineralCapacity = mineralMoles .* maxMineralFraction;
realizedMoles = min(candidateMoles, mineralMoles);
realizedMoles = min(realizedMoles, mineralCapacity);

batchState = buildBatchState(state, geometry, hMolCm3, realizedMoles, options);
batchOptions = options;
batchOptions.timeStepSize = dtSeconds;
batchOptions.rateLaw = 'external_tst_phreeqc';
runner = getRunBatchFunction(options);
batchResult = runner(batchState, batchOptions);

componentDelta = componentDeltaFromBatch(state, batchState, batchResult, waterVolume);
[componentDelta, roundoffInfo] = suppressTinyNegativeComponentRoundoff( ...
    componentDelta, state.component_moles, options);
mineralDelta = zeros(size(state.mineral_moles));
actualMoles = realizedFromBatch(batchResult, realizedMoles, numCells);
actualMoles = min(max(actualMoles, 0), max(mineralMoles, 0));
actualMoles = min(actualMoles, max(mineralCapacity, 0));
mineralDelta(:, calciteIndex) = -actualMoles;

result = struct();
result.component_delta_moles = componentDelta;
result.mineral_delta_moles = mineralDelta;
result.realized_interface_moles = actualMoles;
result.candidate_interface_moles = candidateMoles;
result.interface_rate_mol_cm2_s = perAreaRate(actualMoles, interfaceArea, dtSeconds);
result.candidate_interface_rate_mol_cm2_s = perAreaRate(candidateMoles, interfaceArea, dtSeconds);
result.inventory_limited = candidateMoles > mineralCapacity + ...
    (1e-14 + 1e-12 .* max(candidateMoles, 1));
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
result.aux = struct('chemistry_mode', "external_tst_phreeqc", ...
    'hydrogen_source', "reactionInput", ...
    'reaction_cluster_count', clusterInfo.count, ...
    'reaction_cluster_max_membership', clusterInfo.max_membership, ...
    'reaction_cluster_overlapping_cell_count', clusterInfo.overlapping_cell_count, ...
    'component_delta_roundoff_suppressed_moles', ...
        roundoffInfo.suppressed_moles_total, ...
    'component_delta_roundoff_suppressed_entries', ...
        roundoffInfo.suppressed_entries, ...
    'phreeqc_session_reused', logical(getFieldOrDefault(batchResult, ...
        'phreeqcSessionReused', false)), ...
    'phreeqc_input_written', logical(getFieldOrDefault(batchResult, ...
        'inputWritten', true)), ...
    'phreeqc_run_method', string(getFieldOrDefault(batchResult, ...
        'phreeqcRunMethod', "")), ...
    'phreeqc_database_path', string(getFieldOrDefault(batchResult, ...
        'databasePath', "")));
end

function [componentDelta, info] = suppressTinyNegativeComponentRoundoff( ...
        componentDelta, componentMoles, options)
info = struct('suppressed_moles_total', 0, 'suppressed_entries', 0);
tolerance = scalarOption(options, ...
    'componentNegativeRoundoffTolerance_mol', 1e-24);
if tolerance < 0
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.componentNegativeRoundoffTolerance_mol must be nonnegative.');
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

function batchState = buildBatchState(state, geometry, hMolCm3, prescribedMoles, options)
numCells = size(state.component_moles, 1);
waterVolume = geometry.water_volume_cm3(:);
reactionWaterVolume = reactionWaterVolumeForPhreeqc(waterVolume, options);
batchState = struct();
batchState.h_mol_cm3 = hMolCm3(:);
batchState.ca_total_mol_cm3 = componentConcentration(state, 'Ca', reactionWaterVolume);
batchState.c_total_mol_cm3 = componentConcentration(state, 'C', reactionWaterVolume);
batchState.na_total_mol_cm3 = componentConcentration(state, 'Na', reactionWaterVolume);
batchState.cl_total_mol_cm3 = componentConcentration(state, 'Cl', reactionWaterVolume);
batchState.ca_mol_cm3 = batchState.ca_total_mol_cm3;
batchState.c_mol_cm3 = batchState.c_total_mol_cm3;
batchState.na_mol_cm3 = batchState.na_total_mol_cm3;
batchState.cl_mol_cm3 = batchState.cl_total_mol_cm3;
batchState.water_volume_cm3 = waterVolume;
batchState.reaction_water_volume_cm3 = reactionWaterVolume;
batchState.interface_area_cm2 = geometry.interface_area_cm2(:);
batchState.calcite_moles = mineralMoles(state, 'Calcite', numCells);
batchState.prescribed_calcite_dissolved_moles = prescribedMoles(:);
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

function values = reactionWaterVolumeForPhreeqc(waterVolume, options)
values = waterVolume(:);
if isfield(options, 'reactionWaterVolumeCm3') && ~isempty(options.reactionWaterVolumeCm3)
    values = vectorOption(options.reactionWaterVolumeCm3, ...
        'reactionWaterVolumeCm3', numel(waterVolume));
end
minReactionWaterVolume = getFieldOrDefault(options, ...
    'minReactionWaterVolumeCm3', 0);
if ~(isscalar(minReactionWaterVolume) && isfinite(minReactionWaterVolume) && ...
        minReactionWaterVolume >= 0)
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.minReactionWaterVolumeCm3 must be a nonnegative finite scalar.');
end
values = max(values, minReactionWaterVolume);
end

function componentDelta = componentDeltaFromBatch(state, batchState, batchResult, waterVolume)
componentDelta = zeros(size(state.component_moles));
reactionWaterVolume = batchConservedWaterVolume(batchState, waterVolume);
for iComponent = 1:numel(state.component_names)
    componentName = char(state.component_names{iComponent});
    before = batchConcentration(batchState, componentName, numel(waterVolume));
    after = batchConcentration(batchResult, componentName, numel(waterVolume));
    componentDelta(:, iComponent) = (after - before) .* reactionWaterVolume;
end
if any(~isfinite(componentDelta(:)))
    error('RTSPHEM:Chemistry:InvalidPhreeqcComponentDelta', ...
        'PHREEQC component deltas must be finite.');
end
end

function values = batchConservedWaterVolume(batchState, waterVolume)
if isfield(batchState, 'reaction_water_volume_cm3') && ...
        ~isempty(batchState.reaction_water_volume_cm3)
    values = batchState.reaction_water_volume_cm3(:);
else
    values = waterVolume(:);
end
if numel(values) ~= numel(waterVolume)
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'batchState.reaction_water_volume_cm3 must match the state cell count.');
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

function values = realizedFromBatch(batchResult, fallbackMoles, numCells)
if isfield(batchResult, 'calciteDissolvedMoles') && ~isempty(batchResult.calciteDissolvedMoles)
    values = batchResult.calciteDissolvedMoles(:);
else
    values = fallbackMoles(:);
end
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
values = vectorOption(options.(fieldName), fieldName, numCells);
end

function values = optionalVectorOption(options, fieldName, defaultValues, numCells)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    values = vectorOption(options.(fieldName), fieldName, numCells);
else
    values = defaultValues(:);
end
end

function values = vectorOption(values, fieldName, numCells)
values = values(:);
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

function value = scalarOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && isfinite(value))
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.%s must be a finite scalar.', fieldName);
end
end

function value = fractionOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && (isfinite(value) || isequal(value, Inf)))
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.%s must be a finite scalar or Inf.', fieldName);
end
if value < 0
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.%s must be nonnegative.', fieldName);
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
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
