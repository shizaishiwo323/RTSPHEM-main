function result = precip_ApplyPrescribedCalciteDissolutionSigned(result, state, options)
% precip_ApplyPrescribedCalciteDissolutionSigned - Preserve TST-match moles.

if nargin < 3
    options = struct();
end
if ~isfield(state, 'prescribed_calcite_dissolved_moles') || ...
        isempty(state.prescribed_calcite_dissolved_moles)
    return;
end

numCells = numel(result.calciteDeltaMoles);
timeStepSize = getOption(options, 'timeStepSize', 1);
calciteMoles = optionalColumn(state, 'calcite_moles', numCells, Inf);
prescribed = optionalColumn(state, 'prescribed_calcite_dissolved_moles', numCells, 0);
dissolvedMoles = min(max(prescribed(:), 0), max(calciteMoles(:), 0));
dissolvedMoles(~isfinite(dissolvedMoles)) = 0;

result.calciteDeltaMoles = -dissolvedMoles;
result.calcitePrecipitatedMoles = zeros(numCells, 1);
result.calciteDissolvedMoles = dissolvedMoles;
result.calciteSignedRate_mol_s = -dissolvedMoles ./ max(timeStepSize, eps);
result.calciteDissolutionRate_mol_s = dissolvedMoles ./ max(timeStepSize, eps);
result.calciteRate_mol_s = result.calciteDissolutionRate_mol_s;
result.calciteRate_mol_dm2_s = result.calciteRate_mol_s;
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
