function fineResult = ExpandPhreeqcReactionGridResult(coarseResult, fineTemplate, projection)
% Expand PHREEQC reaction-grid output back to fine transport cells.

fineToReaction = projection.fineToReactionCell(:);
numFine = numel(fineToReaction);
numReaction = projection.numReactionCells;
valid = fineToReaction > 0 & fineToReaction <= numReaction;

fineResult = fineTemplate;
fields = fieldnames(coarseResult);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    coarseValues = coarseResult.(fieldName);
    if ~(isnumeric(coarseValues) && isvector(coarseValues) && numel(coarseValues) == numReaction)
        fineResult.(fieldName) = coarseValues;
        continue;
    end

    fineValues = zeros(numFine, 1);
    if isfield(fineTemplate, fieldName) && isnumeric(fineTemplate.(fieldName)) && ...
            isvector(fineTemplate.(fieldName)) && numel(fineTemplate.(fieldName)) == numFine
        fineValues = fineTemplate.(fieldName)(:);
    end
    fineValues(valid) = coarseValues(fineToReaction(valid));
    fineValues(~isfinite(fineValues)) = 0;
    fineResult.(fieldName) = fineValues;
end

fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'calciteDissolvedMoles', 'interface_area_cm2');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'calciteDeltaMoles', 'interface_area_cm2');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'calciteRate_mol_s', 'interface_area_cm2');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'calciteKinDeltaRate_mol_s', 'interface_area_cm2');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'calciteRawKinDeltaRate_mol_s', 'interface_area_cm2');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'water_phase_h_delta_moles', 'water_volume_cm3');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'water_phase_ca_delta_moles', 'water_volume_cm3');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'water_phase_c_delta_moles', 'water_volume_cm3');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'water_phase_na_delta_moles', 'water_volume_cm3');
fineResult = distributeAmountField(fineResult, coarseResult, fineTemplate, projection, ...
    'water_phase_cl_delta_moles', 'water_volume_cm3');

if isfield(fineTemplate, 'solutionNumber')
    fineResult.solutionNumber = fineTemplate.solutionNumber(:);
else
    fineResult.solutionNumber = (1:numFine)';
end
end

function result = distributeAmountField(result, coarseResult, fineTemplate, projection, fieldName, weightField)
if ~isfield(coarseResult, fieldName) || isempty(coarseResult.(fieldName)) || ...
        ~isfield(fineTemplate, weightField) || isempty(fineTemplate.(weightField))
    return;
end

fineToReaction = projection.fineToReactionCell(:);
numFine = numel(fineToReaction);
numReaction = projection.numReactionCells;
valid = fineToReaction > 0 & fineToReaction <= numReaction;
coarseValues = coarseResult.(fieldName)(:);
if numel(coarseValues) ~= numReaction
    return;
end

weights = max(fineTemplate.(weightField)(:), 0);
if numel(weights) ~= numFine
    return;
end
weightSums = accumarray(fineToReaction(valid), weights(valid), [numReaction, 1], @sum, 0);
fineValues = zeros(numFine, 1);
active = valid & weightSums(fineToReaction) > 0;
fineValues(active) = coarseValues(fineToReaction(active)) .* ...
    weights(active) ./ weightSums(fineToReaction(active));
fineValues(~isfinite(fineValues)) = 0;
result.(fieldName) = fineValues;
end
