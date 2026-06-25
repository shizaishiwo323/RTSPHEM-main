function tests = test_StrictMolinsBackend
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

function testStrictMolinsConsumesHReactantAndDissolvesCalcite(testCase)
state = validStrictState([3e-6; 4e-6], [1e-5; 1e-5]);
geometry = validGeometry();
options.rate_constant_cm_s = 0.1;

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 5, options);

verifyEqual(testCase, result.component_delta_moles, [-3e-6; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.mineral_delta_moles, [-3e-6; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.realized_interface_moles, [3e-6; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.interface_rate_mol_cm2_s, [3e-7; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, result.converged);
verifyFalse(testCase, result.reactant_limited(1));
verifyFalse(testCase, result.inventory_limited(1));
end

function testStrictMolinsLimitsReactionByReactantMoles(testCase)
state = validStrictState([2e-7; 4e-6], [1e-3; 1e-5]);
geometry = validGeometry();
geometry.water_volume_cm3(1) = 0.01;
options.rate_constant_cm_s = 0.1;

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 5, options);

verifyEqual(testCase, result.candidate_interface_moles, [2e-5; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.realized_interface_moles, [2e-7; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles, [-2e-7; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, result.reactant_limited(1));
verifyFalse(testCase, result.inventory_limited(1));
end

function testStrictMolinsLimitsReactionByReactantFraction(testCase)
state = validStrictState([2e-7; 4e-6], [1e-3; 1e-5]);
geometry = validGeometry();
geometry.water_volume_cm3(1) = 0.01;
options.rate_constant_cm_s = 0.1;
options.maxReactantFraction = 0.25;

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 5, options);

verifyEqual(testCase, result.candidate_interface_moles, [2e-5; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.realized_interface_moles, [5e-8; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles, [-5e-8; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, result.reactant_limited(1));
end

function testStrictMolinsClusterSharesReactantCapacity(testCase)
state = validStrictState([1e-9; 1e-5], [1e-3; 1e-3]);
geometry = validGeometry();
geometry.water_volume_cm3 = [1e-6; 1];
geometry.interface_area_cm2 = [1; 0];
options.rate_constant_cm_s = 0.1;
options.reaction_time_integration = 'explicit_euler';
options.reactionClusters = struct('source_cells', 1, ...
    'member_cells', [1; 2]);

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 1, options);

verifyGreaterThan(testCase, result.realized_interface_moles(1), state.component_moles(1));
verifyEqual(testCase, sum(result.component_delta_moles, 'all'), ...
    -result.realized_interface_moles(1), 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyLessThan(testCase, result.component_delta_moles(2), 0);
verifyTrue(testCase, result.aux.clustered_reaction);
end

function testStrictMolinsExactFirstOrderReactionAvoidsFullExplicitDepletion(testCase)
state = validStrictState([10; 0], [100; 100]);
geometry = validGeometry();
geometry.water_volume_cm3 = [2; 2];
geometry.interface_area_cm2 = [4; 0];
options.rate_constant_cm_s = 0.5;
options.maxReactantFraction = Inf;
options.reaction_time_integration = 'exact_first_order';

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 3, options);

expected = 10 .* (1 - exp(-0.5 .* 4 ./ 2 .* 3));
verifyEqual(testCase, result.realized_interface_moles(1), expected, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyLessThan(testCase, result.realized_interface_moles(1), 10);
verifyEqual(testCase, result.component_delta_moles(1), -expected, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testStrictMolinsLimitsReactionByCalciteInventory(testCase)
state = validStrictState([3e-6; 4e-6], [4e-7; 1e-5]);
geometry = validGeometry();
options.rate_constant_cm_s = 0.1;

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 5, options);

verifyEqual(testCase, result.realized_interface_moles, [4e-7; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.mineral_delta_moles, [-4e-7; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyFalse(testCase, result.reactant_limited(1));
verifyTrue(testCase, result.inventory_limited(1));
end

function testStrictMolinsUsesInterfaceStateWhenProvided(testCase)
state = validStrictState([3e-6; 4e-6], [1e-5; 1e-5]);
geometry = validGeometry();
interfaceState = struct();
interfaceState.component_names = {'H_reactant'};
interfaceState.component_concentration_mol_cm3 = [1e-7; 0];
quadrature = struct();
quadrature.cell_id = [1; 1];
quadrature.weight_cm2 = [0.5; 1.5];
options.rate_constant_cm_s = 0.1;
options.interfaceState = interfaceState;
options.interfaceQuadrature = quadrature;

result = rtm.chemistry.StrictMolinsBackend(state, geometry, 5, options);

verifyEqual(testCase, result.candidate_interface_moles, [1e-7; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.interface_rate_mol_cm2_s, [1e-8; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.aux.interface_state_source, "interfaceState");
end

function testStrictMolinsRejectsMissingHReactantBasis(testCase)
state = validStrictState([3e-6; 4e-6], [1e-5; 1e-5]);
state.component_names = {'Ca_total'};
geometry = validGeometry();

verifyError(testCase, @() rtm.chemistry.StrictMolinsBackend(state, geometry, 5, struct()), ...
    'RTSPHEM:Chemistry:MissingHReactant');
end

function testStrictMolinsRejectsNegativeGeometryMeasure(testCase)
state = validStrictState([3e-6; 4e-6], [1e-5; 1e-5]);
geometry = validGeometry();
geometry.interface_area_cm2(1) = -2;

verifyError(testCase, @() rtm.chemistry.StrictMolinsBackend(state, geometry, 5, struct()), ...
    'RTSPHEM:Chemistry:NegativeGeometryMeasure');
end

function testStrictMolinsRejectsNegativeInterfaceStateConcentration(testCase)
state = validStrictState([3e-6; 4e-6], [1e-5; 1e-5]);
geometry = validGeometry();
interfaceState = struct();
interfaceState.component_names = {'H_reactant'};
interfaceState.component_concentration_mol_cm3 = [-1e-7; 2e-7];
options = struct('rate_constant_cm_s', 0.1, 'interfaceState', interfaceState);

verifyError(testCase, @() rtm.chemistry.StrictMolinsBackend(state, geometry, 5, options), ...
    'RTSPHEM:Chemistry:NegativeInterfaceConcentration');
end

function state = validStrictState(hMoles, calciteMoles)
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = hMoles(:);
state.mineral_names = {'Calcite'};
state.mineral_moles = calciteMoles(:);
state.temperature_C = 25 * ones(numel(hMoles), 1);
state.pressure_atm = ones(numel(hMoles), 1);
state.time_s = 0;
end

function geometry = validGeometry()
geometry = struct();
geometry.water_volume_cm3 = [1; 2];
geometry.interface_area_cm2 = [2; 0];
geometry.cell_volume_cm3 = [1; 2];
geometry.solid_volume_cm3 = [0.5; 1];
geometry.fluid_fraction = [1; 1];
end
