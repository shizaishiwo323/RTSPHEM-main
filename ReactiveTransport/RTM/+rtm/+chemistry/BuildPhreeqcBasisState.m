function basisState = BuildPhreeqcBasisState(state, geometry, options)
%BUILDPHREEQCBASISSTATE Convert conserved RTM state to PHREEQC batch fields.
%
% Transported quantities remain conserved component moles. This helper
% exposes concentration-like totals for the PHREEQC text builder while keeping
% pH/free H/carbonate species out of state.component_names.

if nargin < 3 || isempty(options)
    options = struct();
end
rtm.state.ValidateState(state);
rtm.chemistry.ValidatePhreeqcTransportBasis(state, options);
numCells = size(state.component_moles, 1);
waterVolume = requiredGeometryVector(geometry, 'water_volume_cm3', numCells);
componentWaterVolume = componentVolume(waterVolume, options, numCells);
reactionWaterVolume = reactionVolume(waterVolume, options, numCells);

basisState = struct();
basisState.ca_total_mol_cm3 = componentConcentration(state, 'Ca', componentWaterVolume);
basisState.c_total_mol_cm3 = componentConcentration(state, 'C', componentWaterVolume);
basisState.na_total_mol_cm3 = componentConcentration(state, 'Na', componentWaterVolume);
basisState.cl_total_mol_cm3 = componentConcentration(state, 'Cl', componentWaterVolume);
basisState.alkalinity_mol_cm3 = componentConcentration(state, ...
    'Alkalinity', componentWaterVolume);
basisState.ca_mol_cm3 = basisState.ca_total_mol_cm3;
basisState.c_mol_cm3 = basisState.c_total_mol_cm3;
basisState.na_mol_cm3 = basisState.na_total_mol_cm3;
basisState.cl_mol_cm3 = basisState.cl_total_mol_cm3;
basisState.water_volume_cm3 = waterVolume;
basisState.reaction_water_volume_cm3 = reactionWaterVolume;
basisState.interface_area_cm2 = optionalGeometryVector( ...
    geometry, 'interface_area_cm2', numCells, zeros(numCells, 1));
basisState.calcite_moles = mineralMoles(state, 'Calcite', numCells);
basisState.initial_calcite_moles = initialMineralMoles(state, 'Calcite', ...
    numCells);
if isfield(options, 'h_mol_cm3') && ~isempty(options.h_mol_cm3)
    basisState.h_mol_cm3 = requiredOptionVector(options.h_mol_cm3, ...
        numCells, 'options.h_mol_cm3');
end
end

function values = componentConcentration(state, componentName, waterVolume)
idx = find(strcmp(state.component_names, componentName), 1);
values = zeros(size(waterVolume));
if isempty(idx)
    return;
end
active = waterVolume > 0;
values(active) = state.component_moles(active, idx) ./ waterVolume(active);
end

function values = mineralMoles(state, mineralName, numCells)
idx = find(strcmp(state.mineral_names, mineralName), 1);
if isempty(idx)
    values = zeros(numCells, 1);
else
    values = state.mineral_moles(:, idx);
end
end

function values = initialMineralMoles(state, mineralName, numCells)
values = mineralMoles(state, mineralName, numCells);
if ~isfield(state, 'chemistry_aux') || ~isstruct(state.chemistry_aux) || ...
        ~isfield(state.chemistry_aux, 'initial_mineral_moles') || ...
        isempty(state.chemistry_aux.initial_mineral_moles)
    return;
end
idx = find(strcmp(state.mineral_names, mineralName), 1);
candidate = state.chemistry_aux.initial_mineral_moles;
if isempty(idx) || size(candidate, 1) ~= numCells || size(candidate, 2) < idx
    return;
end
valuesCandidate = candidate(:, idx);
if isnumeric(valuesCandidate) && all(isfinite(valuesCandidate(:))) && ...
        all(valuesCandidate(:) >= 0)
    values = valuesCandidate(:);
end
end

function values = componentVolume(waterVolume, options, numCells)
if isfield(options, 'componentWaterVolumeCm3') && ~isempty(options.componentWaterVolumeCm3)
    values = requiredOptionVector(options.componentWaterVolumeCm3, numCells, ...
        'options.componentWaterVolumeCm3');
else
    values = waterVolume;
end
if any(values < 0)
    error('RTSPHEM:Chemistry:InvalidPhreeqcBasisState', ...
        'options.componentWaterVolumeCm3 must be nonnegative.');
end
end

function values = reactionVolume(waterVolume, options, numCells)
if isfield(options, 'reactionWaterVolumeCm3') && ~isempty(options.reactionWaterVolumeCm3)
    values = requiredOptionVector(options.reactionWaterVolumeCm3, numCells, ...
        'options.reactionWaterVolumeCm3');
else
    values = waterVolume;
end
if isfield(options, 'minReactionWaterVolumeCm3') && ...
        ~isempty(options.minReactionWaterVolumeCm3)
    floorValue = options.minReactionWaterVolumeCm3;
    if ~(isscalar(floorValue) && isfinite(floorValue) && floorValue >= 0)
        error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
            'options.minReactionWaterVolumeCm3 must be nonnegative and finite.');
    end
    values = max(values, floorValue);
end
end

function values = requiredGeometryVector(geometry, fieldName, numCells)
if ~isstruct(geometry) || ~isfield(geometry, fieldName) || ...
        isempty(geometry.(fieldName))
    error('RTSPHEM:Chemistry:InvalidPhreeqcBasisState', ...
        'geometry.%s is required.', fieldName);
end
values = requiredOptionVector(geometry.(fieldName), numCells, ...
    ['geometry.', fieldName]);
if any(values < 0)
    error('RTSPHEM:Chemistry:InvalidPhreeqcBasisState', ...
        'geometry.%s must be nonnegative.', fieldName);
end
end

function values = optionalGeometryVector(geometry, fieldName, numCells, defaultValue)
if isstruct(geometry) && isfield(geometry, fieldName) && ~isempty(geometry.(fieldName))
    values = requiredOptionVector(geometry.(fieldName), numCells, ...
        ['geometry.', fieldName]);
else
    values = defaultValue;
end
if any(values < 0)
    error('RTSPHEM:Chemistry:InvalidPhreeqcBasisState', ...
        'geometry.%s must be nonnegative.', fieldName);
end
end

function values = requiredOptionVector(values, numCells, fieldName)
values = values(:);
if isscalar(values) && numCells > 1
    values = repmat(values, numCells, 1);
end
if numel(values) ~= numCells || any(~isfinite(values))
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        '%s must be finite and match the state cell count.', fieldName);
end
end
