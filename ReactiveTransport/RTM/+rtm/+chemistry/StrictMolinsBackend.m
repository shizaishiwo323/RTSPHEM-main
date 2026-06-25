function result = StrictMolinsBackend(state, geometry, dtSeconds, options)
%STRICTMOLINSBACKEND Irreversible Molins-style calcite dissolution backend.
%
% The strict benchmark mode transports an abstract H_reactant component only.
% It does not invoke PHREEQC or carbonate speciation.

if nargin < 4 || isempty(options)
    options = struct();
end
rtm.state.ValidateState(state);
validateGeometry(geometry, size(state.component_moles, 1));
if ~(isscalar(dtSeconds) && isfinite(dtSeconds) && dtSeconds >= 0)
    error('RTSPHEM:Chemistry:InvalidTimeStep', ...
        'dtSeconds must be a nonnegative finite scalar.');
end

hIndex = find(strcmp(state.component_names, 'H_reactant'), 1);
if isempty(hIndex)
    error('RTSPHEM:Chemistry:MissingHReactant', ...
        'strict_molins requires component_names to include H_reactant.');
end
calciteIndex = find(strcmp(state.mineral_names, 'Calcite'), 1);
if isempty(calciteIndex)
    error('RTSPHEM:Chemistry:MissingCalcite', ...
        'strict_molins requires mineral_names to include Calcite.');
end

rateConstantCmS = getOption(options, 'rate_constant_cm_s', 0);
if ~(isscalar(rateConstantCmS) && isfinite(rateConstantCmS) && rateConstantCmS >= 0)
    error('RTSPHEM:Chemistry:InvalidRateConstant', ...
        'options.rate_constant_cm_s must be a nonnegative finite scalar.');
end

waterVolumeCm3 = geometry.water_volume_cm3(:);
interfaceAreaCm2 = geometry.interface_area_cm2(:);
hMoles = state.component_moles(:, hIndex);
calciteMoles = state.mineral_moles(:, calciteIndex);
numCells = size(state.component_moles, 1);
numComponents = size(state.component_moles, 2);
numMinerals = size(state.mineral_moles, 2);
maxReactantFraction = fractionOption(options, 'maxReactantFraction', Inf);
maxMineralFraction = fractionOption(options, 'maxMineralFraction', Inf);

[candidateMoles, candidateRatePerArea, interfaceStateSource] = computeCandidateMoles( ...
    state, geometry, options, hIndex, rateConstantCmS, dtSeconds);
candidateMoles(~isfinite(candidateMoles)) = 0;
candidateMoles = max(candidateMoles, 0);
explicitCandidateMoles = candidateMoles;
clusters = getReactionClusters(options, numCells);
if isempty(clusters) && strcmp(reactionTimeIntegration(options), 'exact_first_order')
    candidateMoles = exactFirstOrderCandidate(candidateMoles, hMoles);
end

componentDelta = zeros(numCells, numComponents);
mineralDelta = zeros(numCells, numMinerals);
if isempty(clusters)
    reactantCapacity = hMoles .* maxReactantFraction;
    mineralCapacity = calciteMoles .* maxMineralFraction;
    realizedMoles = min(candidateMoles, hMoles);
    realizedMoles = min(realizedMoles, reactantCapacity);
    realizedMoles = min(realizedMoles, calciteMoles);
    realizedMoles = min(realizedMoles, mineralCapacity);
    realizedMoles = max(realizedMoles, 0);
    componentDelta(:, hIndex) = -realizedMoles;
    mineralDelta(:, calciteIndex) = -realizedMoles;
else
    [candidateMoles, realizedMoles, componentDelta(:, hIndex), ...
        mineralDelta(:, calciteIndex)] = realizeClusteredReaction( ...
        explicitCandidateMoles, hMoles, calciteMoles, clusters, ...
        reactionTimeIntegration(options), maxReactantFraction, maxMineralFraction);
end

interfaceRate = zeros(numCells, 1);
activeInterface = interfaceAreaCm2 > 0 & dtSeconds > 0;
interfaceRate(activeInterface) = realizedMoles(activeInterface) ./ ...
    interfaceAreaCm2(activeInterface) ./ dtSeconds;
interfaceRate(~isfinite(interfaceRate)) = 0;

result = struct();
result.component_delta_moles = componentDelta;
result.mineral_delta_moles = mineralDelta;
result.realized_interface_moles = realizedMoles;
result.candidate_interface_moles = candidateMoles;
result.explicit_candidate_interface_moles = explicitCandidateMoles;
result.interface_rate_mol_cm2_s = interfaceRate;
result.candidate_interface_rate_mol_cm2_s = candidateRatePerArea;
limitTolerance = 1e-14 + 1e-12 .* max(candidateMoles, 1);
result.reactant_limited = candidateMoles > hMoles + limitTolerance;
result.inventory_limited = candidateMoles > calciteMoles + limitTolerance;
result.pH = nan(numCells, 1);
result.saturation_index = nan(numCells, numMinerals);
result.converged = true;
result.failed_cells = [];
result.error_message = "";
result.aux = struct('chemistry_mode', 'strict_molins', ...
    'interface_state_source', interfaceStateSource, ...
    'clustered_reaction', ~isempty(clusters));
end

function [candidateMoles, realizedMoles, componentDelta, mineralDelta] = ...
    realizeClusteredReaction(explicitCandidateMoles, hMoles, calciteMoles, ...
        clusters, timeIntegration, maxReactantFraction, maxMineralFraction)
numCells = numel(hMoles);
candidateMoles = zeros(numCells, 1);
realizedMoles = zeros(numCells, 1);
componentDelta = zeros(numCells, 1);
mineralDelta = zeros(numCells, 1);

for iCluster = 1:numel(clusters)
    sourceCells = clusters(iCluster).source_cells(:);
    memberCells = clusters(iCluster).member_cells(:);
    sourceDemand = explicitCandidateMoles(sourceCells);
    sourceDemand(~isfinite(sourceDemand)) = 0;
    sourceDemand = max(sourceDemand, 0);
    if ~any(sourceDemand > 0)
        continue;
    end

    clusterH = sum(hMoles(memberCells));
    if strcmp(timeIntegration, 'exact_first_order')
        effectiveTotalCandidate = exactFirstOrderScalar(sum(sourceDemand), clusterH);
        sourceCandidate = distributeByWeights(effectiveTotalCandidate, sourceDemand);
    else
        sourceCandidate = sourceDemand;
    end
    candidateMoles(sourceCells) = candidateMoles(sourceCells) + sourceCandidate;

    sourceMineralCapacity = min(calciteMoles(sourceCells), ...
        calciteMoles(sourceCells) .* maxMineralFraction);
    clusterReactantCapacity = clusterH .* maxReactantFraction;
    totalRealized = min(sum(sourceCandidate), clusterH);
    totalRealized = min(totalRealized, clusterReactantCapacity);
    totalRealized = min(totalRealized, sum(sourceMineralCapacity));
    totalRealized = max(totalRealized, 0);
    sourceRealized = allocateByCapacity(sourceCandidate, ...
        sourceMineralCapacity, totalRealized);

    realizedMoles(sourceCells) = realizedMoles(sourceCells) + sourceRealized;
    mineralDelta(sourceCells) = mineralDelta(sourceCells) - sourceRealized;
    memberConsumption = allocateByCapacity(hMoles(memberCells), ...
        hMoles(memberCells), sum(sourceRealized));
    componentDelta(memberCells) = componentDelta(memberCells) - memberConsumption;
end
end

function value = exactFirstOrderScalar(explicitCandidateMoles, hMoles)
if hMoles <= 0 || explicitCandidateMoles <= 0
    value = 0;
else
    value = hMoles .* (1 - exp(-explicitCandidateMoles ./ hMoles));
end
value(~isfinite(value)) = 0;
value = max(value, 0);
end

function values = distributeByWeights(totalValue, weights)
weights = max(weights(:), 0);
if totalValue <= 0 || sum(weights) <= 0
    values = zeros(size(weights));
else
    values = totalValue .* weights ./ sum(weights);
end
end

function allocation = allocateByCapacity(weights, capacities, totalValue)
weights = max(weights(:), 0);
capacities = max(capacities(:), 0);
allocation = zeros(size(weights));
remaining = max(totalValue, 0);
active = capacities > 0 & weights > 0;
while remaining > 0 && any(active)
    proposal = remaining .* weights(active) ./ sum(weights(active));
    activeIndices = find(active);
    cappedLocal = proposal >= capacities(active) - allocation(active);
    increment = proposal;
    increment(cappedLocal) = capacities(activeIndices(cappedLocal)) - ...
        allocation(activeIndices(cappedLocal));
    increment = max(increment, 0);
    allocation(activeIndices) = allocation(activeIndices) + increment;
    remaining = totalValue - sum(allocation);
    active = capacities - allocation > 1e-18 & weights > 0;
    if ~any(cappedLocal)
        break;
    end
end
if remaining > 0 && sum(capacities - allocation) > 0
    spare = max(capacities - allocation, 0);
    allocation = allocation + min(remaining, sum(spare)) .* spare ./ sum(spare);
end
allocation = min(allocation, capacities);
end

function clusters = getReactionClusters(options, numCells)
clusters = struct([]);
if ~isfield(options, 'reactionClusters') || isempty(options.reactionClusters)
    return;
end
clusters = options.reactionClusters;
requiredFields = {'source_cells', 'member_cells'};
for iCluster = 1:numel(clusters)
    for iField = 1:numel(requiredFields)
        if ~isfield(clusters(iCluster), requiredFields{iField})
            error('RTSPHEM:Chemistry:InvalidReactionCluster', ...
                'reactionClusters.%s is required.', requiredFields{iField});
        end
    end
    validateClusterCells(clusters(iCluster).source_cells, numCells);
    validateClusterCells(clusters(iCluster).member_cells, numCells);
end
end

function validateClusterCells(cellIds, numCells)
cellIds = cellIds(:);
if any(~isfinite(cellIds)) || any(cellIds < 1) || any(cellIds > numCells) || ...
        any(cellIds ~= round(cellIds))
    error('RTSPHEM:Chemistry:InvalidReactionCluster', ...
        'reactionClusters contain invalid cell indices.');
end
end

function candidateMoles = exactFirstOrderCandidate(explicitCandidateMoles, hMoles)
candidateMoles = zeros(size(explicitCandidateMoles));
active = hMoles > 0 & explicitCandidateMoles > 0;
exponent = -explicitCandidateMoles(active) ./ hMoles(active);
candidateMoles(active) = hMoles(active) .* (1 - exp(exponent));
candidateMoles(~isfinite(candidateMoles)) = 0;
candidateMoles = max(candidateMoles, 0);
end

function value = reactionTimeIntegration(options)
value = getOption(options, 'reaction_time_integration', 'explicit_euler');
value = lower(strrep(strtrim(char(value)), '-', '_'));
if ~any(strcmp(value, {'explicit_euler', 'exact_first_order'}))
    error('RTSPHEM:Chemistry:UnknownReactionTimeIntegration', ...
        'Unsupported strict_molins reaction time integration: %s.', value);
end
end

function [candidateMoles, candidateRatePerArea, source] = computeCandidateMoles( ...
    state, geometry, options, hIndex, rateConstantCmS, dtSeconds)
numCells = size(state.component_moles, 1);
interfaceAreaCm2 = geometry.interface_area_cm2(:);
if isfield(options, 'interfaceState') && ~isempty(options.interfaceState)
    [candidateMoles, candidateRatePerArea] = candidateFromInterfaceState( ...
        options.interfaceState, getOption(options, 'interfaceQuadrature', []), ...
        state, hIndex, rateConstantCmS, dtSeconds, interfaceAreaCm2);
    source = "interfaceState";
    return;
end

waterVolumeCm3 = geometry.water_volume_cm3(:);
hMoles = state.component_moles(:, hIndex);
hConcentration = zeros(numCells, 1);
activeWater = waterVolumeCm3 > 0;
hConcentration(activeWater) = hMoles(activeWater) ./ waterVolumeCm3(activeWater);
candidateRatePerArea = rateConstantCmS .* hConcentration;
candidateMoles = candidateRatePerArea .* interfaceAreaCm2 .* dtSeconds;
source = "cellAverage";
end

function [candidateMoles, candidateRatePerArea] = candidateFromInterfaceState( ...
    interfaceState, quadrature, state, hIndex, rateConstantCmS, dtSeconds, interfaceAreaCm2)
numCells = size(state.component_moles, 1);
componentIndex = find(strcmp(interfaceState.component_names, state.component_names{hIndex}), 1);
if isempty(componentIndex)
    error('RTSPHEM:Chemistry:MissingHReactant', ...
        'interfaceState must include H_reactant concentration.');
end
values = interfaceState.component_concentration_mol_cm3(:, componentIndex);
if numel(values) ~= numCells
    error('RTSPHEM:Chemistry:InvalidInterfaceConcentration', ...
        'interfaceState concentrations must match the state cell count.');
end
candidateMoles = zeros(numCells, 1);

if isempty(quadrature)
    active = interfaceAreaCm2 > 0;
    validateInterfaceConcentrations(values(active));
    candidateMoles(active) = rateConstantCmS .* values(active) .* ...
        interfaceAreaCm2(active) .* dtSeconds;
else
    validateQuadrature(quadrature, numCells);
    validateInterfaceConcentrations(values(unique(quadrature.cell_id(:))));
    for iPoint = 1:numel(quadrature.cell_id)
        cellId = quadrature.cell_id(iPoint);
        concentration = values(cellId);
        candidateMoles(cellId) = candidateMoles(cellId) + ...
            rateConstantCmS .* concentration .* ...
            quadrature.weight_cm2(iPoint) .* dtSeconds;
    end
end

candidateRatePerArea = zeros(numCells, 1);
activeInterface = interfaceAreaCm2 > 0 & dtSeconds > 0;
candidateRatePerArea(activeInterface) = candidateMoles(activeInterface) ./ ...
    interfaceAreaCm2(activeInterface) ./ dtSeconds;
candidateRatePerArea(~isfinite(candidateRatePerArea)) = 0;
end

function validateInterfaceConcentrations(values)
if any(~isfinite(values(:)))
    error('RTSPHEM:Chemistry:InvalidInterfaceConcentration', ...
        'active interfaceState concentrations must be finite.');
end
if any(values(:) < 0)
    error('RTSPHEM:Chemistry:NegativeInterfaceConcentration', ...
        'active interfaceState concentrations must be nonnegative.');
end
end

function validateQuadrature(quadrature, numCells)
requiredFields = {'cell_id', 'weight_cm2'};
for iField = 1:numel(requiredFields)
    if ~isfield(quadrature, requiredFields{iField})
        error('RTSPHEM:Chemistry:InvalidInterfaceQuadrature', ...
            'interfaceQuadrature.%s is required.', requiredFields{iField});
    end
end
cellId = quadrature.cell_id(:);
if any(cellId < 1) || any(cellId > numCells) || any(cellId ~= round(cellId))
    error('RTSPHEM:Chemistry:InvalidInterfaceQuadrature', ...
        'interfaceQuadrature.cell_id contains invalid cell indices.');
end
if numel(quadrature.weight_cm2) ~= numel(cellId)
    error('RTSPHEM:Chemistry:InvalidInterfaceQuadrature', ...
        'interfaceQuadrature weights must match cell_id length.');
end
weights = quadrature.weight_cm2(:);
if any(~isfinite(weights)) || any(weights < 0)
    error('RTSPHEM:Chemistry:InvalidInterfaceQuadrature', ...
        'interfaceQuadrature.weight_cm2 must contain nonnegative finite values.');
end
end

function validateGeometry(geometry, numCells)
requiredFields = {'water_volume_cm3', 'interface_area_cm2'};
for iField = 1:numel(requiredFields)
    fieldName = requiredFields{iField};
    if ~isstruct(geometry) || ~isfield(geometry, fieldName)
        error('RTSPHEM:Chemistry:MissingGeometryField', ...
            'geometry.%s is required.', fieldName);
    end
    if numel(geometry.(fieldName)) ~= numCells
        error('RTSPHEM:Chemistry:GeometrySizeMismatch', ...
            'geometry.%s must match the state cell count.', fieldName);
    end
    values = geometry.(fieldName)(:);
    if any(~isfinite(values))
        error('RTSPHEM:Chemistry:InvalidGeometryMeasure', ...
            'geometry.%s must contain finite values.', fieldName);
    end
    if any(values < 0)
        error('RTSPHEM:Chemistry:NegativeGeometryMeasure', ...
            'geometry.%s must be nonnegative.', fieldName);
    end
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function value = fractionOption(options, fieldName, defaultValue)
value = getOption(options, fieldName, defaultValue);
if ~(isscalar(value) && isfinite(value) || isequal(value, Inf))
    error('RTSPHEM:Chemistry:InvalidFractionLimit', ...
        'options.%s must be a nonnegative finite scalar or Inf.', fieldName);
end
if value < 0
    error('RTSPHEM:Chemistry:InvalidFractionLimit', ...
        'options.%s must be nonnegative.', fieldName);
end
end
