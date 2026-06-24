function result = precip_ScaleSignedCalciteDeltaToCellInventory(result, state, options)
% precip_ScaleSignedCalciteDeltaToCellInventory - Scale signed PHREEQC deltas.

if nargin < 3
    options = struct();
end
numCells = numel(result.calciteDeltaMoles);
if numCells == 0
    return;
end

waterVolume = optionalColumn(state, 'water_volume_cm3', numCells, 0);
calciteMoles = optionalColumn(state, 'calcite_moles', numCells, Inf);
solutionWaterKg = getOption(options, 'solutionWaterKg', 1);
timeStepSize = getOption(options, 'timeStepSize', 1);
capPrecipitationByInventory = logical(getOption(options, ...
    'capPrecipitationByInitialAqueousInventory', false));

scale = max(waterVolume(:), 0) * 1e-3 ./ max(solutionWaterKg, eps);
rawDeltaCell = result.calciteDeltaMoles(:) .* scale;
deltaCell = rawDeltaCell;

negative = deltaCell < 0;
deltaCell(negative) = -min(abs(deltaCell(negative)), max(calciteMoles(negative), 0));

positive = deltaCell > 0;
if capPrecipitationByInventory && any(positive)
    availableCa = availableMoles(state, result, numCells, waterVolume, ...
        {'ca_total_mol_cm3', 'ca_mol_cm3'});
    availableC = availableCarbonMoles(state, result, numCells, waterVolume);
    maxPrecipitation = min(availableCa, availableC);
    deltaCell(positive) = min(deltaCell(positive), max(maxPrecipitation(positive), 0));
end

deltaCell(~isfinite(deltaCell)) = 0;
result = scaleAqueousFieldsIfDeltaWasLimited(result, state, rawDeltaCell, deltaCell);

result.calciteDeltaMoles = deltaCell;
result.calcitePrecipitatedMoles = max(deltaCell, 0);
result.calciteDissolvedMoles = max(-deltaCell, 0);
result.calciteSignedRate_mol_s = deltaCell ./ max(timeStepSize, eps);
result.calciteDissolutionRate_mol_s = result.calciteDissolvedMoles ./ max(timeStepSize, eps);
result.calciteRate_mol_s = result.calciteDissolvedMoles ./ max(timeStepSize, eps);
result.calciteRate_mol_dm2_s = result.calciteRate_mol_s;
end

function result = scaleAqueousFieldsIfDeltaWasLimited(result, state, rawDeltaCell, deltaCell)
limited = abs(rawDeltaCell) > eps & abs(deltaCell - rawDeltaCell) > eps;
if ~any(limited)
    return;
end

ratio = ones(size(deltaCell));
ratio(limited) = deltaCell(limited) ./ rawDeltaCell(limited);
fieldsToBlend = {'h_mol_cm3', 'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'na_total_mol_cm3', 'cl_total_mol_cm3', 'ca_mol_cm3', 'hco3_mol_cm3', ...
    'co3_mol_cm3', 'cl_mol_cm3', 'na_mol_cm3'};
for iField = 1:numel(fieldsToBlend)
    fieldName = fieldsToBlend{iField};
    if ~isfield(result, fieldName) || isempty(result.(fieldName))
        continue;
    end
    postValue = result.(fieldName)(:);
    [preValue, foundPreValue] = concentrationFromFields( ...
        state, struct(), numel(deltaCell), preFieldCandidates(fieldName));
    if ~foundPreValue
        continue;
    end
    blended = preValue + ratio .* (postValue - preValue);
    result.(fieldName) = blended;
end
end

function candidates = preFieldCandidates(fieldName)
switch fieldName
    case 'ca_total_mol_cm3'
        candidates = {'ca_total_mol_cm3', 'ca_mol_cm3'};
    case 'c_total_mol_cm3'
        candidates = {'c_total_mol_cm3', 'c_mol_cm3'};
    case 'na_total_mol_cm3'
        candidates = {'na_total_mol_cm3', 'na_mol_cm3'};
    case 'cl_total_mol_cm3'
        candidates = {'cl_total_mol_cm3', 'cl_mol_cm3'};
    otherwise
        candidates = {fieldName};
end
end

function moles = availableCarbonMoles(state, result, numCells, waterVolume)
cTotal = concentrationFromFields(state, result, numCells, ...
    {'c_total_mol_cm3', 'c_mol_cm3'});
if all(cTotal == 0)
    hco3 = concentrationFromFields(state, result, numCells, {'hco3_mol_cm3'});
    co3 = concentrationFromFields(state, result, numCells, {'co3_mol_cm3'});
    cTotal = hco3 + co3;
end
moles = max(cTotal(:), 0) .* max(waterVolume(:), 0);
end

function moles = availableMoles(state, result, numCells, waterVolume, fieldNames)
concentration = concentrationFromFields(state, result, numCells, fieldNames);
moles = max(concentration(:), 0) .* max(waterVolume(:), 0);
end

function [concentration, found] = concentrationFromFields(state, result, numCells, fieldNames)
concentration = zeros(numCells, 1);
found = false;
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    if isfield(state, fieldName) && ~isempty(state.(fieldName))
        concentration = state.(fieldName)(:);
        found = true;
        return;
    end
    if isfield(result, fieldName) && ~isempty(result.(fieldName))
        concentration = result.(fieldName)(:);
        found = true;
        return;
    end
end
end

function values = optionalColumn(state, fieldName, numCells, defaultValue)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    values = state.(fieldName)(:);
else
    values = repmat(defaultValue, numCells, 1);
end
if numel(values) ~= numCells
    error('RTSPHEM:Precipitate:StateSizeMismatch', ...
        'State field %s has %d values, expected %d.', fieldName, numel(values), numCells);
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
