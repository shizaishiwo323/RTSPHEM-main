function coarseState = AggregatePhreeqcStateToReactionGrid(fineState, projection)
% Aggregate fine transport/geometry fields onto a PHREEQC reaction grid.

fineToReaction = projection.fineToReactionCell(:);
numFine = numel(fineToReaction);
numReaction = projection.numReactionCells;
valid = fineToReaction > 0 & fineToReaction <= numReaction;

coarseState = struct();
waterVolume = optionalVector(fineState, 'water_volume_cm3', numFine, 0);
reactionWaterVolume = optionalVector(fineState, 'reaction_water_volume_cm3', numFine, waterVolume);
interfaceArea = optionalVector(fineState, 'interface_area_cm2', numFine, 0);
calciteMoles = optionalVector(fineState, 'calcite_moles', numFine, 0);

coarseState.water_volume_cm3 = accumByReaction(waterVolume, fineToReaction, numReaction, valid);
coarseState.reaction_water_volume_cm3 = accumByReaction(reactionWaterVolume, fineToReaction, numReaction, valid);
coarseState.interface_area_cm2 = accumByReaction(interfaceArea, fineToReaction, numReaction, valid);
coarseState.calcite_moles = accumByReaction(calciteMoles, fineToReaction, numReaction, valid);

concentrationFields = { ...
    'h_mol_cm3', 'ca_mol_cm3', 'c_mol_cm3', 'na_mol_cm3', 'cl_mol_cm3', ...
    'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'na_total_mol_cm3', 'cl_total_mol_cm3', ...
    'hco3_mol_cm3', 'co3_mol_cm3', 'alkalinity_mol_cm3'};

for iField = 1:numel(concentrationFields)
    fieldName = concentrationFields{iField};
    if ~isfield(fineState, fieldName) || isempty(fineState.(fieldName))
        continue;
    end
    values = fineState.(fieldName)(:);
    assert(numel(values) == numFine, ...
        'RTSPHEM:Phreeqc:ProjectionFieldSizeMismatch', ...
        'Field %s has %d entries, expected %d.', fieldName, numel(values), numFine);
    weighted = accumByReaction(values .* waterVolume, fineToReaction, numReaction, valid);
    coarseState.(fieldName) = safeDivide(weighted, coarseState.water_volume_cm3);
end

if isfield(fineState, 'initial_calcite_moles') && ~isempty(fineState.initial_calcite_moles)
    values = fineState.initial_calcite_moles(:);
    coarseState.initial_calcite_moles = accumByReaction(values, fineToReaction, numReaction, valid);
end
end

function values = optionalVector(s, fieldName, n, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    values = s.(fieldName)(:);
else
    if isscalar(defaultValue)
        values = repmat(defaultValue, n, 1);
    else
        values = defaultValue(:);
    end
end
assert(numel(values) == n, ...
    'RTSPHEM:Phreeqc:ProjectionFieldSizeMismatch', ...
    'Field %s has %d entries, expected %d.', fieldName, numel(values), n);
values(~isfinite(values)) = 0;
values = max(values, 0);
end

function out = accumByReaction(values, fineToReaction, numReaction, valid)
out = accumarray(fineToReaction(valid), values(valid), [numReaction, 1], @sum, 0);
end

function out = safeDivide(numerator, denominator)
out = zeros(size(numerator));
active = denominator > 0;
out(active) = numerator(active) ./ denominator(active);
out(~isfinite(out)) = 0;
end
