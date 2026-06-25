function tests = test_PhreeqcReactionClusters
tests = functiontests(localfunctions);
end

function setupOnce(~)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rtmDir);
end

function testAggregatesClusterStateExtensiveQuantities(testCase)
state = carbonateState();
state.chemistry_aux.initial_mineral_moles = [2e-6; 0; 5e-6];
geometry = clusterGeometry();
clusters = sharedClusters();
prescribed = [1e-10; 0; 3e-10];
hMolCm3 = [1e-7; 2e-7; 4e-7];

[clusterState, map] = rtm.chemistry.AggregateReactionClusterState( ...
    state, geometry, hMolCm3, prescribed, clusters);

verifyEqual(testCase, numel(clusterState.h_mol_cm3), 1);
verifyEqual(testCase, clusterState.water_volume_cm3, 6, 'AbsTol', 1e-18);
verifyEqual(testCase, clusterState.reaction_water_volume_cm3, 6, 'AbsTol', 1e-18);
verifyEqual(testCase, clusterState.interface_area_cm2, 4, 'AbsTol', 1e-18);
verifyEqual(testCase, clusterState.prescribed_calcite_dissolved_moles, 4e-10, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, clusterState.calcite_moles, 3e-6, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, clusterState.initial_calcite_moles, 7e-6, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, clusterState.ca_total_mol_cm3, ...
    sum(state.component_moles(:, 1)) ./ 6, 'RelTol', 1e-12);
verifyEqual(testCase, clusterState.c_total_mol_cm3, ...
    sum(state.component_moles(:, 2)) ./ 6, 'RelTol', 1e-12);
verifyEqual(testCase, clusterState.alkalinity_mol_cm3, ...
    sum(state.component_moles(:, 5)) ./ 6, 'RelTol', 1e-12);
verifyEqual(testCase, clusterState.h_mol_cm3, ...
    sum(hMolCm3(:) .* geometry.water_volume_cm3(:)) ./ 6, ...
    'RelTol', 1e-12);
verifyEqual(testCase, map.cluster_id_by_cell, [1; 1; 1]);
verifyEqual(testCase, map.membership_count, [1; 1; 1]);
verifyEqual(testCase, map.source_cells{1}, [1; 3]);
verifyEqual(testCase, map.member_cells{1}, [1; 2; 3]);
end

function testScatterClusterResultConservesComponentAndMineralDeltas(testCase)
state = carbonateState();
geometry = clusterGeometry();
clusters = sharedClusters();
prescribed = [1e-10; 0; 3e-10];
hMolCm3 = [1e-7; 2e-7; 4e-7];
[clusterState, map] = rtm.chemistry.AggregateReactionClusterState( ...
    state, geometry, hMolCm3, prescribed, clusters);

clusterResult = struct();
clusterResult.ca_total_mol_cm3 = clusterState.ca_total_mol_cm3 + 4e-10 ./ 6;
clusterResult.c_total_mol_cm3 = clusterState.c_total_mol_cm3 + 4e-10 ./ 6;
clusterResult.na_total_mol_cm3 = clusterState.na_total_mol_cm3;
clusterResult.cl_total_mol_cm3 = clusterState.cl_total_mol_cm3;
clusterResult.alkalinity_mol_cm3 = clusterState.alkalinity_mol_cm3 + 5e-10 ./ 6;
clusterResult.calciteDissolvedMoles = 4e-10;
clusterResult.pH = 6.5;
clusterResult.calciteSI = -1;
clusterResult.chargeBalance = 2e-12;

cellResult = rtm.chemistry.ScatterClusterReactionResult( ...
    clusterResult, map, size(state.component_moles, 1));

expectedCaDelta = 4e-10 .* geometry.water_volume_cm3(:) ./ ...
    sum(geometry.water_volume_cm3(:));
expectedMineral = [1e-10; 0; 3e-10];
verifyEqual(testCase, ...
    (cellResult.ca_total_mol_cm3(:) - map.before.ca_total_mol_cm3(:)) .* ...
    map.reaction_water_volume_cm3(:), expectedCaDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ...
    (cellResult.c_total_mol_cm3(:) - map.before.c_total_mol_cm3(:)) .* ...
    map.reaction_water_volume_cm3(:), expectedCaDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
expectedAlkalinityDelta = 5e-10 .* geometry.water_volume_cm3(:) ./ ...
    sum(geometry.water_volume_cm3(:));
verifyEqual(testCase, ...
    (cellResult.alkalinity_mol_cm3(:) - map.before.alkalinity_mol_cm3(:)) .* ...
    map.reaction_water_volume_cm3(:), expectedAlkalinityDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, cellResult.calciteDissolvedMoles, expectedMineral, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, sum(cellResult.calciteDissolvedMoles), ...
    clusterResult.calciteDissolvedMoles, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, cellResult.pH, repmat(6.5, 3, 1));
verifyEqual(testCase, cellResult.chargeBalance, repmat(2e-12, 3, 1));
end

function testRunPhreeqcReactionClustersAggregatesRunsAndScatters(testCase)
captured = struct('call_count', 0);
state = carbonateState();
geometry = clusterGeometry();
clusters = sharedClusters();
prescribed = [1e-10; 0; 3e-10];
hMolCm3 = [1e-7; 2e-7; 4e-7];
options = struct('runBatchFunction', @mockRunBatch);

[cellResult, info] = rtm.chemistry.RunPhreeqcReactionClusters( ...
    state, geometry, hMolCm3, prescribed, clusters, options);

verifyEqual(testCase, captured.call_count, 1);
verifyEqual(testCase, numel(captured.batch_state.h_mol_cm3), 1);
verifyEqual(testCase, captured.batch_state.prescribed_calcite_dissolved_moles, ...
    sum(prescribed), 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, info.cluster_count, 1);
verifyEqual(testCase, info.membership_count, [1; 1; 1]);
verifyEqual(testCase, info.cluster_state.water_volume_cm3, 6, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, cellResult.calciteDissolvedMoles, prescribed, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, sum(cellResult.calciteDissolvedMoles), sum(prescribed), ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ...
    sum((cellResult.ca_total_mol_cm3 - info.map.before.ca_total_mol_cm3) .* ...
    info.map.reaction_water_volume_cm3), sum(prescribed), ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);

    function batchResult = mockRunBatch(batchState, batchOptions)
        captured.call_count = captured.call_count + 1;
        captured.batch_state = batchState;
        captured.batch_options = batchOptions;
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        water = batchState.reaction_water_volume_cm3(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ water;
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ water;
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.alkalinity_mol_cm3 = batchState.alkalinity_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
        batchResult.pH = 6.4;
        batchResult.calciteSI = -0.5;
        batchResult.chargeBalance = 1e-12;
    end
end

function testRunPhreeqcReactionClustersCanOmitPrescribedReactionForKinetics(testCase)
captured = struct();
state = carbonateState();
geometry = clusterGeometry();
clusters = sharedClusters();
hMolCm3 = [1e-7; 2e-7; 4e-7];
options = struct('runBatchFunction', @mockRunBatch, ...
    'omitPrescribedCalciteReaction', true);

[cellResult, info] = rtm.chemistry.RunPhreeqcReactionClusters( ...
    state, geometry, hMolCm3, zeros(3, 1), clusters, options);

verifyFalse(testCase, isfield(captured.batch_state, ...
    'prescribed_calcite_dissolved_moles'));
verifyEqual(testCase, info.cluster_count, 1);
verifyEqual(testCase, sum(cellResult.calciteDissolvedMoles), 6e-10, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);

    function batchResult = mockRunBatch(batchState, ~)
        captured.batch_state = batchState;
        dissolved = 6e-10;
        water = batchState.reaction_water_volume_cm3(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ water;
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ water;
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.alkalinity_mol_cm3 = batchState.alkalinity_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function testScatterMapsFailedClusterIdsToMemberCells(testCase)
state = carbonateState();
geometry = clusterGeometry();
clusters = struct( ...
    'source_cells', {1, 3}, ...
    'member_cells', {[1; 2], 3});
prescribed = [1e-10; 0; 3e-10];
hMolCm3 = [1e-7; 2e-7; 4e-7];
[clusterState, map] = rtm.chemistry.AggregateReactionClusterState( ...
    state, geometry, hMolCm3, prescribed, clusters);

clusterResult = struct();
clusterResult.ca_total_mol_cm3 = clusterState.ca_total_mol_cm3;
clusterResult.c_total_mol_cm3 = clusterState.c_total_mol_cm3;
clusterResult.na_total_mol_cm3 = clusterState.na_total_mol_cm3;
clusterResult.cl_total_mol_cm3 = clusterState.cl_total_mol_cm3;
clusterResult.alkalinity_mol_cm3 = clusterState.alkalinity_mol_cm3;
clusterResult.calciteDissolvedMoles = [0; 0];
clusterResult.failedCells = 2;
clusterResult.errorMessage = "cluster 2 PHREEQC failure";

cellResult = rtm.chemistry.ScatterClusterReactionResult( ...
    clusterResult, map, size(state.component_moles, 1));

verifyEqual(testCase, cellResult.failedCells, 3);
verifyEqual(testCase, cellResult.errorMessage, ...
    "cluster 2 PHREEQC failure");
end

function state = carbonateState()
state = struct();
state.component_names = {'Ca', 'C', 'Na', 'Cl', 'Alkalinity'};
state.component_moles = [
    1e-9, 2e-9, 1e-8, 1e-8, 4e-9
    2e-9, 3e-9, 2e-8, 2e-8, 5e-9
    3e-9, 4e-9, 3e-8, 3e-8, 6e-9
    ];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1e-6; 0; 2e-6];
state.temperature_C = [25; 25; 25];
state.pressure_atm = [1; 1; 1];
state.time_s = 0;
end

function geometry = clusterGeometry()
geometry = struct();
geometry.water_volume_cm3 = [1; 2; 3];
geometry.interface_area_cm2 = [1; 0; 3];
geometry.solid_volume_cm3 = [1; 0; 1];
geometry.interface_h_cm = [1; 1; 1];
end

function clusters = sharedClusters()
clusters = struct('source_cells', [1; 3], 'member_cells', [1; 2; 3]);
end
