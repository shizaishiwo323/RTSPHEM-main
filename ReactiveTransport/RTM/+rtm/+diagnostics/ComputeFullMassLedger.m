function ledger = ComputeFullMassLedger(stepInfo, options)
%COMPUTEFULLMASSLEDGER Build a full coupled RTM mass/geometry ledger.

if nargin < 2 || isempty(options)
    options = struct(); %#ok<NASGU>
end
mass = getStructField(stepInfo, 'mass');
transport = firstStructField(stepInfo, {'transport', 'transport_ledger'});
reaction = firstStructField(stepInfo, {'reaction', 'chemistry_ledger'});
remap = getStructField(stepInfo, 'remap');
geometry = getStructField(stepInfo, 'geometry');

componentNames = string(getFieldOrDefault(mass, 'component_names', ...
    getFieldOrDefault(transport, 'component_names', {})));
initialComponents = rowVector(getRequiredField(mass, ...
    'initial_component_moles_total'));
finalComponents = rowVector(getRequiredField(mass, ...
    'final_component_moles_total'));
numComponents = numel(initialComponents);
if isempty(componentNames)
    componentNames = "component_" + string(1:numComponents);
end
componentNames = componentNames(:).';

transportSource = vectorField(transport, 'source_delta_moles_total', ...
    numComponents, 0);
transportInternal = vectorField(transport, ...
    'internal_flux_delta_moles_total', numComponents, 0);
transportBoundary = vectorField(transport, ...
    'boundary_delta_moles_total', numComponents, 0);
reactionDelta = vectorField(reaction, 'component_delta_moles_total', ...
    numComponents, 0);
remapDelta = vectorField(remap, 'component_delta_moles_total', ...
    numComponents, 0);
componentResidual = finalComponents - initialComponents - ...
    transportSource - transportInternal - transportBoundary - ...
    reactionDelta - remapDelta;
componentScale = max(abs([initialComponents(:); finalComponents(:); 1e-300]));

ledger = struct();
ledger.component_names = componentNames;
ledger.initial_component_moles_total = initialComponents;
ledger.final_component_moles_total = finalComponents;
ledger.transport_source_delta_moles_total = transportSource;
ledger.transport_internal_flux_delta_moles_total = transportInternal;
ledger.transport_boundary_delta_moles_total = transportBoundary;
ledger.reaction_delta_moles_total = reactionDelta;
ledger.remap_delta_moles_total = remapDelta;
ledger.transport_roundoff_suppressed_moles_total = ...
    scalarField(transport, 'roundoff_suppressed_moles_total', 0);
ledger.component_residual_moles = componentResidual;
ledger.max_abs_component_residual_moles = max(abs(componentResidual));
ledger.component_relative_residual = ...
    ledger.max_abs_component_residual_moles ./ max(componentScale, eps);

ledger.initial_mineral_moles_total = vectorField(reaction, ...
    'initial_mineral_moles_total', [], []);
ledger.final_mineral_moles_total = vectorField(reaction, ...
    'final_mineral_moles_total', [], []);
ledger.reaction_mineral_delta_moles_total = vectorField(reaction, ...
    'mineral_delta_moles_total', [], []);
ledger.mineral_residual_moles = mineralResidual(ledger);
ledger.max_abs_mineral_residual_moles = max(abs(ledger.mineral_residual_moles));
mineralScale = max(abs([ledger.initial_mineral_moles_total(:); ...
    ledger.final_mineral_moles_total(:); ...
    ledger.reaction_mineral_delta_moles_total(:); 1e-300]));
ledger.mineral_relative_residual = ...
    ledger.max_abs_mineral_residual_moles ./ max(mineralScale, eps);

ledger.initial_solid_volume_cm3 = scalarField(geometry, ...
    'solid_volume_before_cm3', NaN);
ledger.final_solid_volume_cm3 = scalarField(geometry, ...
    'solid_volume_after_cm3', NaN);
ledger.expected_solid_volume_change_cm3 = expectedSolidChange(geometry);
ledger.actual_solid_volume_change_cm3 = scalarField(geometry, ...
    'actual_solid_volume_change_cm3', NaN);
ledger.solid_volume_residual_cm3 = ledger.actual_solid_volume_change_cm3 - ...
    ledger.expected_solid_volume_change_cm3;
solidScale = max(abs([ledger.initial_solid_volume_cm3, ...
    ledger.final_solid_volume_cm3, ...
    ledger.expected_solid_volume_change_cm3, 1e-300]));
ledger.solid_volume_relative_residual = ...
    abs(ledger.solid_volume_residual_cm3) ./ max(solidScale, eps);
ledger.max_displacement_over_h = scalarField(geometry, ...
    'max_displacement_over_h', 0);
end

function value = expectedSolidChange(geometry)
if isfield(geometry, 'expected_solid_volume_change_cm3') && ...
        ~isempty(geometry.expected_solid_volume_change_cm3)
    value = geometry.expected_solid_volume_change_cm3;
else
    dissolvedMoles = scalarField(geometry, 'mineral_dissolved_moles', 0);
    molarVolume = scalarField(geometry, 'molar_volume_cm3_mol', 1);
    value = -dissolvedMoles .* molarVolume;
end
end

function residual = mineralResidual(ledger)
initialMinerals = ledger.initial_mineral_moles_total;
finalMinerals = ledger.final_mineral_moles_total;
deltaMinerals = ledger.reaction_mineral_delta_moles_total;
if isempty(initialMinerals) || isempty(finalMinerals) || isempty(deltaMinerals)
    residual = [];
    return;
end
numMinerals = max([numel(initialMinerals), numel(finalMinerals), numel(deltaMinerals)]);
initialMinerals = expandVector(initialMinerals, numMinerals);
finalMinerals = expandVector(finalMinerals, numMinerals);
deltaMinerals = expandVector(deltaMinerals, numMinerals);
residual = finalMinerals - initialMinerals - deltaMinerals;
end

function values = expandVector(values, expectedLength)
values = values(:).';
if isscalar(values) && expectedLength > 1
    values = repmat(values, 1, expectedLength);
elseif numel(values) ~= expectedLength
    values = nan(1, expectedLength);
end
end

function value = firstStructField(structValue, fieldNames)
value = struct();
for iField = 1:numel(fieldNames)
    candidate = getStructField(structValue, fieldNames{iField});
    if ~isempty(fieldnames(candidate))
        value = candidate;
        return;
    end
end
end

function value = getStructField(structValue, fieldName)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        isstruct(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = struct();
end
end

function value = getRequiredField(structValue, fieldName)
if ~isstruct(structValue) || ~isfield(structValue, fieldName) || ...
        isempty(structValue.(fieldName))
    error('RTSPHEM:Diagnostics:MissingField', ...
        'Missing required field: %s.', fieldName);
end
value = structValue.(fieldName);
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function values = vectorField(structValue, fieldName, expectedLength, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    values = structValue.(fieldName);
else
    values = defaultValue;
end
values = rowVector(values);
if isempty(expectedLength)
    return;
end
if isscalar(values) && expectedLength > 1
    values = repmat(values, 1, expectedLength);
elseif numel(values) ~= expectedLength
    values = repmat(defaultValue, 1, expectedLength);
end
end

function value = scalarField(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    values = structValue.(fieldName);
    value = values(1);
else
    value = defaultValue;
end
end

function values = rowVector(values)
values = values(:).';
end
