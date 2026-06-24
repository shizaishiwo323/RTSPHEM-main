function updated = ApplyPhreeqcAgglomeratedDelta(baseState, reactionInput, reactionOutput, ...
    waterVolumeCm3, agglomerateWeightMatrix)
% Apply PHREEQC concentration deltas from agglomerated cut-cell reactions.
%
% reactionInput/reactionOutput are cell-centered PHREEQC states for reaction
% source cells. agglomerateWeightMatrix(i,j) gives the fraction of source
% reaction j assigned to water cell i. The update converts each source
% concentration delta to a mole delta over the agglomerated water volume and
% distributes that mole delta back to member cells.

numCells = numel(waterVolumeCm3);
updated = baseState;
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
safeWaterCm3 = max(waterVolumeCm3, eps);
concentrationFields = { ...
    'h_mol_cm3', 'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'na_total_mol_cm3', 'cl_total_mol_cm3', ...
    'ca_mol_cm3', 'hco3_mol_cm3', 'co3_mol_cm3', ...
    'na_mol_cm3', 'cl_mol_cm3'};

for iField = 1:numel(concentrationFields)
    fieldName = concentrationFields{iField};
    if ~isfield(baseState, fieldName) || ~isfield(reactionInput, fieldName) || ...
            ~isfield(reactionOutput, fieldName)
        continue;
    end
    baseValues = baseState.(fieldName)(:);
    inputValues = reactionInput.(fieldName)(:);
    outputValues = reactionOutput.(fieldName)(:);
    if numel(baseValues) ~= numCells || ...
            numel(inputValues) ~= size(agglomerateWeightMatrix, 2) || ...
            numel(outputValues) ~= size(agglomerateWeightMatrix, 2)
        error('RTSPHEM:Phreeqc:AgglomerateFieldSizeMismatch', ...
            'Field %s has incompatible size for agglomerated PHREEQC update.', ...
            fieldName);
    end
    sourceMoleDelta = (outputValues - inputValues) .* sourceWaterCm3;
    cellMoleDelta = agglomerateWeightMatrix * sourceMoleDelta;
    updated.(fieldName) = max(baseValues + cellMoleDelta ./ safeWaterCm3, 0);
end

if isfield(updated, 'h_mol_cm3')
    updated.pH = -log10(max(updated.h_mol_cm3(:) * 1000, 1e-14));
end
end
