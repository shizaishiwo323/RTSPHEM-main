function result = InferCalciteDissolutionFromTotals(result, state, timeStepSize)
% InferCalciteDissolutionFromTotals - Backfill calcite dissolution from Ca/C totals.
%
% PHREEQC database RATES can update solution totals even when KIN_DELTA is not
% available in selected output for a given run. This fallback preserves the
% preferred KIN_DELTA value when present, and otherwise infers CaCO3 moles from
% the stoichiometric increase in aqueous Ca and C.

if nargin < 3 || isempty(timeStepSize)
    timeStepSize = 1;
end

if ~needsFallback(result)
    return;
end

numCells = numel(result.ca_total_mol_cm3);
waterVolume = optionalColumn(state, 'water_volume_cm3', numCells, 1000);
caBefore = optionalColumn(state, 'ca_mol_cm3', numCells, 0);
cBefore = optionalColumn(state, 'c_mol_cm3', numCells, 0);

deltaCaMoles = (result.ca_total_mol_cm3(:) - caBefore(:)) .* waterVolume(:);
deltaCMoles = (result.c_total_mol_cm3(:) - cBefore(:)) .* waterVolume(:);
dissolvedMoles = max(min(deltaCaMoles, deltaCMoles), 0);
dissolvedMoles(~isfinite(dissolvedMoles)) = 0;

result.calciteDissolvedMoles = dissolvedMoles;
result.calciteDeltaMoles = -dissolvedMoles;
result.calciteRate_mol_s = dissolvedMoles ./ max(timeStepSize, eps);
result.calciteRate_mol_dm2_s = result.calciteRate_mol_s;
end

function tf = needsFallback(result)
if ~isfield(result, 'calciteDissolvedMoles') || isempty(result.calciteDissolvedMoles)
    tf = true;
    return;
end
if any(result.calciteDissolvedMoles > 0 & isfinite(result.calciteDissolvedMoles))
    tf = false;
    return;
end
if ~isfield(result, 'calciteRate_mol_s') || all(~isfinite(result.calciteRate_mol_s))
    tf = true;
    return;
end
tf = all(result.calciteRate_mol_s == 0 | ~isfinite(result.calciteRate_mol_s));
end

function values = optionalColumn(state, fieldName, numCells, defaultValue)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    values = state.(fieldName)(:);
else
    values = repmat(defaultValue, numCells, 1);
end
if numel(values) ~= numCells
    error('RTSPHEM:Phreeqc:StateSizeMismatch', ...
        'State field %s has %d values, expected %d.', fieldName, numel(values), numCells);
end
end
