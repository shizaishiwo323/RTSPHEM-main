function mixedState = BuildPhreeqcAgglomeratedState(state, waterVolumeCm3, agglomerateWeightMatrix)
% Build mixed source-cell solutions for agglomerated PHREEQC reactions.

numCells = numel(waterVolumeCm3);
mixedState = state;
waterVolumeCm3 = max(waterVolumeCm3(:), 0);
if isempty(agglomerateWeightMatrix) || nnz(agglomerateWeightMatrix) == 0
    return;
end
if size(agglomerateWeightMatrix, 1) ~= numCells
    error('RTSPHEM:Phreeqc:AgglomerateSizeMismatch', ...
        'Agglomerate matrix has %d rows, expected %d.', ...
        size(agglomerateWeightMatrix, 1), numCells);
end

sourceWaterCm3 = full(sum(spones(agglomerateWeightMatrix) .* waterVolumeCm3, 1))';
sourceCells = find(sourceWaterCm3 > 0);
columnSums = full(sum(agglomerateWeightMatrix, 1))';
badColumns = sourceCells(abs(columnSums(sourceCells) - 1) > 1e-10);
if ~isempty(badColumns)
    error('RTSPHEM:Phreeqc:AgglomerateWeightsNotNormalized', ...
        'Agglomerate weights must sum to one for each active source column.');
end
concentrationFields = { ...
    'h_mol_cm3', 'ca_mol_cm3', 'c_mol_cm3', ...
    'na_mol_cm3', 'cl_mol_cm3', ...
    'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'na_total_mol_cm3', 'cl_total_mol_cm3', ...
    'hco3_mol_cm3', 'co3_mol_cm3'};

for iField = 1:numel(concentrationFields)
    fieldName = concentrationFields{iField};
    if ~isfield(state, fieldName)
        continue;
    end
    values = state.(fieldName)(:);
    if numel(values) ~= numCells
        error('RTSPHEM:Phreeqc:AgglomerateFieldSizeMismatch', ...
            'Field %s has incompatible size for agglomerated PHREEQC input.', ...
            fieldName);
    end
    mixedValues = full(agglomerateWeightMatrix' * values);
    values(sourceCells) = mixedValues(sourceCells);
    mixedState.(fieldName) = values;
end

if isfield(state, 'water_volume_cm3')
    values = state.water_volume_cm3(:);
    values(sourceCells) = max(values(sourceCells), sourceWaterCm3(sourceCells));
    mixedState.water_volume_cm3 = values;
end
if isfield(state, 'reaction_water_volume_cm3')
    values = state.reaction_water_volume_cm3(:);
    values(sourceCells) = max(values(sourceCells), sourceWaterCm3(sourceCells));
    mixedState.reaction_water_volume_cm3 = values;
end
end
