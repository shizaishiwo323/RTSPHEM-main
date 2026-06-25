function tests = test_BuildStrictMolinsReactionInput
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testBuildsInterfaceStateAndQuadratureForStrictBackend(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1; 1; 1];
geometry.cell_centroid_cm = [0 0; 1 0; 0 1];
geometry.interface_centroid_cm = [0.5 0; NaN NaN; NaN NaN];
geometry.interface_area_cm2 = [2; 0; 0];
geometry.interface_normal = [1 0; NaN NaN; NaN NaN];
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = [1e-6; 2e-6; 1e-6];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1e-5; 1e-5; 1e-5];
state.temperature_C = [25; 25; 25];
state.pressure_atm = [1; 1; 1];
state.time_s = 0;
connectivity.cell_neighbors = {[2 3], 1, 1};
options.rate_constant_cm_s = 0.2;

reactionInput = rtm.chemistry.BuildStrictMolinsReactionInput( ...
    state, geometry, connectivity, options);

verifyEqual(testCase, reactionInput.rate_constant_cm_s, 0.2);
verifyEqual(testCase, reactionInput.interfaceQuadrature.cell_id, 1);
verifyEqual(testCase, reactionInput.interfaceQuadrature.weight_cm2, 2);
verifyEqual(testCase, reactionInput.interfaceState.source_cell(1), 1);

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 1, reactionInput);

verifyGreaterThan(testCase, result.candidate_interface_moles(1), 0);
verifyEqual(testCase, result.aux.interface_state_source, "interfaceState");
end

function testDefaultRateConstantIsZero(testCase)
geometry = struct();
geometry.water_volume_cm3 = 1;
geometry.cell_centroid_cm = [0 0];
geometry.interface_centroid_cm = [0 0];
geometry.interface_area_cm2 = 1;
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = 1e-6;
state.mineral_names = {'Calcite'};
state.mineral_moles = 1e-5;
state.temperature_C = 25;
state.pressure_atm = 1;
state.time_s = 0;

reactionInput = rtm.chemistry.BuildStrictMolinsReactionInput(state, geometry, struct());

verifyEqual(testCase, reactionInput.rate_constant_cm_s, 0);
end

function testBuildsReactionClustersFromInterfaceNeighborhood(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1e-6; 1; 1];
geometry.cell_centroid_cm = [0 0; 1 0; 0 1];
geometry.interface_centroid_cm = [0.1 0; NaN NaN; NaN NaN];
geometry.interface_area_cm2 = [1; 0; 0];
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = [1e-9; 1e-6; 2e-6];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1e-3; 1e-3; 1e-3];
state.temperature_C = [25; 25; 25];
state.pressure_atm = [1; 1; 1];
state.time_s = 0;
connectivity.cell_neighbors = {[2 3], 1, 1};

reactionInput = rtm.chemistry.BuildStrictMolinsReactionInput( ...
    state, geometry, connectivity, struct());

verifyTrue(testCase, isfield(reactionInput, 'reactionClusters'));
verifyEqual(testCase, reactionInput.reactionClusters.source_cells, 1);
verifyEqual(testCase, reactionInput.reactionClusters.member_cells, [1; 2; 3]);
end

function testReactionClusterDepthCanIncludeSecondRingNeighbors(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1e-6; 1; 1];
geometry.cell_centroid_cm = [0 0; 1 0; 2 0];
geometry.interface_centroid_cm = [0.1 0; NaN NaN; NaN NaN];
geometry.interface_area_cm2 = [1; 0; 0];
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = [1e-9; 1e-6; 2e-6];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1e-3; 1e-3; 1e-3];
state.temperature_C = [25; 25; 25];
state.pressure_atm = [1; 1; 1];
state.time_s = 0;
connectivity.cell_neighbors = {2, [1 3], 2};
options.reactionClusterDepthCells = 2;

reactionInput = rtm.chemistry.BuildStrictMolinsReactionInput( ...
    state, geometry, connectivity, options);

verifyEqual(testCase, reactionInput.reactionClusters.member_cells, [1; 2; 3]);
end
