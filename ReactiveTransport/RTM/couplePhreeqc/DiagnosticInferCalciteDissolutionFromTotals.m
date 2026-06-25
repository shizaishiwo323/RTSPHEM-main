function result = DiagnosticInferCalciteDissolutionFromTotals(result, state, timeStepSize)
%DIAGNOSTICINFERCALCITEDISSOLUTIONFROMTOTALS Infer CaCO3 only for diagnostics.
%
% This helper estimates a stoichiometric calcite amount from aqueous Ca/C
% total changes, but it deliberately writes only diagnostic fields. Production
% mineral deltas must come from PHREEQC KIN_DELTA or explicitly prescribed
% external reaction amounts.

if nargin < 3 || isempty(timeStepSize)
    timeStepSize = 1;
end

numCells = numel(result.ca_total_mol_cm3);
waterVolume = optionalColumn(state, 'water_volume_cm3', numCells, 1000);
caBefore = optionalColumn(state, 'ca_mol_cm3', numCells, 0);
cBefore = optionalColumn(state, 'c_mol_cm3', numCells, 0);

deltaCaMoles = (result.ca_total_mol_cm3(:) - caBefore(:)) .* waterVolume(:);
deltaCMoles = (result.c_total_mol_cm3(:) - cBefore(:)) .* waterVolume(:);
inferredDissolved = max(min(deltaCaMoles, deltaCMoles), 0);
inferredDissolved(~isfinite(inferredDissolved)) = 0;

productionDissolved = optionalColumn(result, 'calciteDissolvedMoles', ...
    numCells, 0);

result.diagnostic_inferred_calcite_dissolved_moles = inferredDissolved;
result.diagnostic_inferred_calcite_delta_moles = -inferredDissolved;
result.diagnostic_inferred_calcite_cell_rate_mol_s = ...
    inferredDissolved ./ max(timeStepSize, eps);
result.diagnostic_inferred_calcite_rate_mol_s = ...
    result.diagnostic_inferred_calcite_cell_rate_mol_s;
result.diagnostic_calcite_inference_minus_reported_moles = ...
    inferredDissolved - productionDissolved(:);
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
