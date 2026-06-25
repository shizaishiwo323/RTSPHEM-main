function diagnostics = ValidateCalciteComponentStoichiometry( ...
        state, componentDelta, calciteDissolvedMoles, clusterIds, options)
%VALIDATECALCITECOMPONENTSTOICHIOMETRY Check PHREEQC Ca/C deltas.
%
% Calcite dissolution adds one mole of Ca and one mole of C to the aqueous
% conserved totals for every mole of calcite dissolved. For clustered
% PHREEQC calls, the check is applied to each disjoint cluster.

if nargin < 5 || isempty(options)
    options = struct();
end
rtm.state.ValidateState(state);
numCells = size(state.component_moles, 1);
componentDelta = validateComponentDelta(componentDelta, state);
calciteDissolvedMoles = validateVector(calciteDissolvedMoles, numCells, ...
    'calciteDissolvedMoles');
clusterIds = validateClusterIds(clusterIds, numCells);

caIndex = find(strcmp(state.component_names, 'Ca'), 1);
cIndex = find(strcmp(state.component_names, 'C'), 1);
if isempty(caIndex) || isempty(cIndex)
    error('RTSPHEM:Chemistry:MissingCalciteStoichiometryComponents', ...
        'PHREEQC calcite chemistry requires conserved Ca and C components.');
end

absoluteTolerance = scalarOption(options, ...
    'calciteStoichiometryAbsoluteTolerance_mol', 1e-14);
relativeTolerance = scalarOption(options, ...
    'calciteStoichiometryRelativeTolerance', 1e-8);

groups = stoichiometryGroups(clusterIds);
residuals = zeros(numel(groups), 2);
relativeResiduals = zeros(numel(groups), 2);
labels = strings(numel(groups), 1);
for iGroup = 1:numel(groups)
    members = groups{iGroup};
    expected = sum(calciteDissolvedMoles(members), 'omitnan');
    observed = [sum(componentDelta(members, caIndex), 'omitnan'), ...
        sum(componentDelta(members, cIndex), 'omitnan')];
    residuals(iGroup, :) = observed - expected;
    scale = max([abs(expected), abs(observed), 1e-300]);
    relativeResiduals(iGroup, :) = abs(residuals(iGroup, :)) ./ max(scale, eps);
    labels(iGroup) = groupLabel(clusterIds(members));
end

maxAbsResidual = max(abs(residuals(:)));
maxRelResidual = max(relativeResiduals(:));
diagnostics = struct();
diagnostics.accepted = maxAbsResidual <= absoluteTolerance || ...
    maxRelResidual <= relativeTolerance;
diagnostics.group_labels = labels;
diagnostics.component_names = ["Ca", "C"];
diagnostics.residual_moles = residuals;
diagnostics.relative_residual = relativeResiduals;
diagnostics.max_absolute_residual_moles = maxAbsResidual;
diagnostics.max_relative_residual = maxRelResidual;
diagnostics.absolute_tolerance_mol = absoluteTolerance;
diagnostics.relative_tolerance = relativeTolerance;

if ~diagnostics.accepted
    error('RTSPHEM:Chemistry:CalciteStoichiometryMismatch', ...
        ['PHREEQC calcite component deltas are inconsistent with ', ...
        'calciteDissolvedMoles: max residual %.6g mol, relative %.6g.'], ...
        maxAbsResidual, maxRelResidual);
end
end

function values = validateComponentDelta(values, state)
values = values(:, :);
if ~isequal(size(values), size(state.component_moles))
    error('RTSPHEM:Chemistry:InvalidCalciteStoichiometryInput', ...
        'componentDelta must match state.component_moles size.');
end
if any(~isfinite(values(:)))
    error('RTSPHEM:Chemistry:InvalidCalciteStoichiometryInput', ...
        'componentDelta must contain finite values.');
end
end

function values = validateVector(values, numCells, fieldName)
values = values(:);
if numel(values) ~= numCells || any(~isfinite(values))
    error('RTSPHEM:Chemistry:InvalidCalciteStoichiometryInput', ...
        '%s must be finite and match the state cell count.', fieldName);
end
end

function values = validateClusterIds(values, numCells)
if nargin < 1 || isempty(values)
    values = nan(numCells, 1);
else
    values = values(:);
end
if numel(values) ~= numCells
    error('RTSPHEM:Chemistry:InvalidCalciteStoichiometryInput', ...
        'clusterIds must match the state cell count.');
end
if any(~isfinite(values) & ~isnan(values))
    error('RTSPHEM:Chemistry:InvalidCalciteStoichiometryInput', ...
        'clusterIds must contain finite values or NaN.');
end
end

function groups = stoichiometryGroups(clusterIds)
if all(isnan(clusterIds))
    groups = {find(true(size(clusterIds)))};
    return;
end
finiteLabels = unique(clusterIds(isfinite(clusterIds)), 'stable');
groups = cell(numel(finiteLabels), 1);
for iLabel = 1:numel(finiteLabels)
    groups{iLabel} = find(clusterIds == finiteLabels(iLabel));
end
unclustered = find(isnan(clusterIds));
if ~isempty(unclustered)
    groups{end + 1, 1} = unclustered;
end
end

function label = groupLabel(values)
finiteValues = unique(values(isfinite(values)), 'stable');
if isempty(finiteValues)
    label = "unclustered";
else
    label = "cluster_" + string(finiteValues(1));
end
end

function value = scalarOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && isfinite(value) && value >= 0)
    error('RTSPHEM:Chemistry:InvalidCalciteStoichiometryInput', ...
        'options.%s must be a nonnegative finite scalar.', fieldName);
end
end
