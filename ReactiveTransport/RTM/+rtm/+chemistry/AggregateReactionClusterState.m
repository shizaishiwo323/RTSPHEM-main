function [clusterState, map] = AggregateReactionClusterState( ...
        state, geometry, hMolCm3, prescribedMoles, clusters, options)
%AGGREGATEREACTIONCLUSTERSTATE Build PHREEQC batch state for reaction clusters.

if nargin < 6 || isempty(options)
    options = struct();
end
numCells = size(state.component_moles, 1);
numClusters = numel(clusters);
waterVolume = requiredVector(geometry.water_volume_cm3, numCells, ...
    'geometry.water_volume_cm3');
reactionWaterVolume = reactionWaterVolumeForPhreeqc(waterVolume, options);
hMolCm3 = requiredVector(hMolCm3, numCells, 'hMolCm3');
prescribedMoles = requiredVector(prescribedMoles, numCells, 'prescribedMoles');
interfaceArea = requiredVector(geometry.interface_area_cm2, numCells, ...
    'geometry.interface_area_cm2');

map = struct();
map.clusters = clusters;
map.cluster_id_by_cell = zeros(numCells, 1);
map.membership_count = zeros(numCells, 1);
map.source_cells = cell(numClusters, 1);
map.member_cells = cell(numClusters, 1);
map.reaction_water_volume_cm3 = reactionWaterVolume;
map.water_volume_cm3 = waterVolume;
map.prescribed_calcite_dissolved_moles = prescribedMoles;
map.calcite_moles = mineralMoles(state, 'Calcite', numCells);
map.initial_calcite_moles = initialMineralMoles(state, 'Calcite', numCells);
map.before = buildCellBatchState(state, geometry, hMolCm3, prescribedMoles, ...
    reactionWaterVolume);

clusterState = emptyBatchState(numClusters);
for iCluster = 1:numClusters
    [sourceCells, memberCells] = clusterCells(clusters(iCluster), numCells);
    map.source_cells{iCluster} = sourceCells;
    map.member_cells{iCluster} = memberCells;
    map.membership_count(memberCells) = map.membership_count(memberCells) + 1;
    map.cluster_id_by_cell(memberCells) = iCluster;

    clusterWater = sum(waterVolume(memberCells), 'omitnan');
    clusterReactionWater = sum(reactionWaterVolume(memberCells), 'omitnan');
    if clusterReactionWater <= 0
        componentWater = 0;
    else
        componentWater = clusterReactionWater;
    end

    clusterState.h_mol_cm3(iCluster, 1) = weightedMean( ...
        hMolCm3(memberCells), reactionWaterVolume(memberCells));
    clusterState.ca_total_mol_cm3(iCluster, 1) = clusterComponentConcentration( ...
        state, 'Ca', memberCells, componentWater);
    clusterState.c_total_mol_cm3(iCluster, 1) = clusterComponentConcentration( ...
        state, 'C', memberCells, componentWater);
    clusterState.na_total_mol_cm3(iCluster, 1) = clusterComponentConcentration( ...
        state, 'Na', memberCells, componentWater);
    clusterState.cl_total_mol_cm3(iCluster, 1) = clusterComponentConcentration( ...
        state, 'Cl', memberCells, componentWater);
    clusterState.alkalinity_mol_cm3(iCluster, 1) = clusterComponentConcentration( ...
        state, 'Alkalinity', memberCells, componentWater);
    clusterState.water_volume_cm3(iCluster, 1) = clusterWater;
    clusterState.reaction_water_volume_cm3(iCluster, 1) = clusterReactionWater;
    clusterState.interface_area_cm2(iCluster, 1) = sum(interfaceArea(sourceCells), ...
        'omitnan');
    clusterState.calcite_moles(iCluster, 1) = sum(map.calcite_moles(sourceCells), ...
        'omitnan');
    clusterState.initial_calcite_moles(iCluster, 1) = ...
        sum(map.initial_calcite_moles(sourceCells), 'omitnan');
    clusterState.prescribed_calcite_dissolved_moles(iCluster, 1) = ...
        sum(prescribedMoles(sourceCells), 'omitnan');
end

if any(map.membership_count > 1)
    error('RTSPHEM:Chemistry:OverlappingReactionClusters', ...
        'Reaction cluster member cells must be disjoint.');
end
clusterState.ca_mol_cm3 = clusterState.ca_total_mol_cm3;
clusterState.c_mol_cm3 = clusterState.c_total_mol_cm3;
clusterState.na_mol_cm3 = clusterState.na_total_mol_cm3;
clusterState.cl_mol_cm3 = clusterState.cl_total_mol_cm3;
end

function batchState = buildCellBatchState(state, geometry, hMolCm3, ...
        prescribedMoles, reactionWaterVolume)
numCells = size(state.component_moles, 1);
batchState = struct();
batchState.h_mol_cm3 = hMolCm3(:);
batchState.ca_total_mol_cm3 = componentConcentration(state, 'Ca', ...
    reactionWaterVolume);
batchState.c_total_mol_cm3 = componentConcentration(state, 'C', ...
    reactionWaterVolume);
batchState.na_total_mol_cm3 = componentConcentration(state, 'Na', ...
    reactionWaterVolume);
batchState.cl_total_mol_cm3 = componentConcentration(state, 'Cl', ...
    reactionWaterVolume);
batchState.alkalinity_mol_cm3 = componentConcentration(state, 'Alkalinity', ...
    reactionWaterVolume);
batchState.ca_mol_cm3 = batchState.ca_total_mol_cm3;
batchState.c_mol_cm3 = batchState.c_total_mol_cm3;
batchState.na_mol_cm3 = batchState.na_total_mol_cm3;
batchState.cl_mol_cm3 = batchState.cl_total_mol_cm3;
batchState.water_volume_cm3 = geometry.water_volume_cm3(:);
batchState.reaction_water_volume_cm3 = reactionWaterVolume(:);
batchState.interface_area_cm2 = geometry.interface_area_cm2(:);
batchState.calcite_moles = mineralMoles(state, 'Calcite', numCells);
batchState.initial_calcite_moles = initialMineralMoles(state, 'Calcite', ...
    numCells);
batchState.prescribed_calcite_dissolved_moles = prescribedMoles(:);
end

function batchState = emptyBatchState(numRows)
zerosColumn = zeros(numRows, 1);
batchState = struct();
batchState.h_mol_cm3 = zerosColumn;
batchState.ca_total_mol_cm3 = zerosColumn;
batchState.c_total_mol_cm3 = zerosColumn;
batchState.na_total_mol_cm3 = zerosColumn;
batchState.cl_total_mol_cm3 = zerosColumn;
batchState.alkalinity_mol_cm3 = zerosColumn;
batchState.water_volume_cm3 = zerosColumn;
batchState.reaction_water_volume_cm3 = zerosColumn;
batchState.interface_area_cm2 = zerosColumn;
batchState.calcite_moles = zerosColumn;
batchState.initial_calcite_moles = zerosColumn;
batchState.prescribed_calcite_dissolved_moles = zerosColumn;
end

function values = componentConcentration(state, componentName, waterVolume)
idx = find(strcmp(state.component_names, componentName), 1);
values = zeros(size(waterVolume));
if isempty(idx)
    return;
end
active = waterVolume > 0;
values(active) = state.component_moles(active, idx) ./ waterVolume(active);
end

function value = clusterComponentConcentration(state, componentName, cells, waterVolume)
idx = find(strcmp(state.component_names, componentName), 1);
value = 0;
if isempty(idx) || waterVolume <= 0
    return;
end
value = sum(state.component_moles(cells, idx), 'omitnan') ./ waterVolume;
end

function values = mineralMoles(state, mineralName, numCells)
idx = find(strcmp(state.mineral_names, mineralName), 1);
if isempty(idx)
    values = zeros(numCells, 1);
else
    values = state.mineral_moles(:, idx);
end
end

function values = initialMineralMoles(state, mineralName, numCells)
values = mineralMoles(state, mineralName, numCells);
if ~isfield(state, 'chemistry_aux') || ~isstruct(state.chemistry_aux) || ...
        ~isfield(state.chemistry_aux, 'initial_mineral_moles') || ...
        isempty(state.chemistry_aux.initial_mineral_moles)
    return;
end
idx = find(strcmp(state.mineral_names, mineralName), 1);
initialMineralMoles = state.chemistry_aux.initial_mineral_moles;
if isempty(idx) || size(initialMineralMoles, 1) ~= numCells || ...
        size(initialMineralMoles, 2) < idx
    return;
end
candidate = initialMineralMoles(:, idx);
if isnumeric(candidate) && all(isfinite(candidate(:))) && all(candidate(:) >= 0)
    values = candidate(:);
end
end

function [sourceCells, memberCells] = clusterCells(cluster, numCells)
sourceCells = requiredCellIndices(cluster.source_cells, numCells, ...
    'source_cells');
memberCells = requiredCellIndices(cluster.member_cells, numCells, ...
    'member_cells');
if isempty(memberCells)
    memberCells = sourceCells;
end
end

function values = requiredCellIndices(values, numCells, fieldName)
values = unique(values(:), 'stable');
if any(values < 1 | values > numCells | values ~= fix(values))
    error('RTSPHEM:Chemistry:InvalidReactionCluster', ...
        'Reaction cluster %s must contain valid cell indices.', fieldName);
end
end

function values = requiredVector(values, numCells, fieldName)
values = values(:);
if isscalar(values) && numCells > 1
    values = repmat(values, numCells, 1);
end
if numel(values) ~= numCells || any(~isfinite(values))
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        '%s must be finite and match the state cell count.', fieldName);
end
end

function values = reactionWaterVolumeForPhreeqc(waterVolume, options)
values = waterVolume(:);
if isfield(options, 'reactionWaterVolumeCm3') && ~isempty(options.reactionWaterVolumeCm3)
    values = requiredVector(options.reactionWaterVolumeCm3, numel(waterVolume), ...
        'options.reactionWaterVolumeCm3');
end
minReactionWaterVolume = getFieldOrDefault(options, ...
    'minReactionWaterVolumeCm3', 0);
if ~(isscalar(minReactionWaterVolume) && isfinite(minReactionWaterVolume) && ...
        minReactionWaterVolume >= 0)
    error('RTSPHEM:Chemistry:InvalidPhreeqcChemistryInput', ...
        'options.minReactionWaterVolumeCm3 must be a nonnegative finite scalar.');
end
values = max(values, minReactionWaterVolume);
end

function value = weightedMean(values, weights)
totalWeight = sum(weights, 'omitnan');
if totalWeight <= 0
    value = 0;
else
    value = sum(values(:) .* weights(:), 'omitnan') ./ totalWeight;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
