function tests = test_CutCellTransportFV
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

function testZeroFluxAndZeroSourceKeepsMolesAndAdvancesTime(testCase)
state = validState();
geometry = validGeometry();

[updated, ledger] = rtm.transport.CutCellTransportFV(state, geometry, 2, struct());

verifyEqual(testCase, updated.component_moles, state.component_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, updated.time_s, 7);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], 'AbsTol', 1e-18);
end

function testCellSourceAddsComponentMoles(testCase)
state = validState();
geometry = validGeometry();
options.component_source_mol_s = [1e-9, 0; 0, -2e-9];

[updated, ledger] = rtm.transport.CutCellTransportFV(state, geometry, 3, options);

verifyEqual(testCase, updated.component_moles, ...
    [1.003e-6, 2e-6; 3e-6, 3.994e-6], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.source_delta_moles_total, [3e-9, -6e-9], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], 'AbsTol', 1e-18);
end

function testInternalFaceFluxTransfersMolesConservatively(testCase)
state = validState();
geometry = validGeometry();
options.internal_face_cells = [1 2];
options.internal_face_flux_mol_s = [2e-9, -1e-9];

[updated, ledger] = rtm.transport.CutCellTransportFV(state, geometry, 4, options);

verifyEqual(testCase, updated.component_moles, ...
    [0.992e-6, 2.004e-6; 3.008e-6, 3.996e-6], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, sum(updated.component_moles, 1), ...
    sum(state.component_moles, 1), 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.internal_flux_delta_moles_total, [0 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], 'AbsTol', 1e-18);
end

function testInternalFaceAdvectionDiffusionUsesCellConcentrations(testCase)
state = validState();
geometry = validGeometry();
options.internal_face_cells = [1 2];
options.internal_face_area_cm2 = 2;
options.internal_face_distance_cm = 0.5;
options.internal_face_velocity_cm_s = 0.25;
options.diffusion_coefficient_cm2_s = [0.01, 0.02];

[updated, ledger] = rtm.transport.CutCellTransportFV(state, geometry, 0.1, options);

leftConcentration = state.component_moles(1, :) ./ geometry.water_volume_cm3(1);
rightConcentration = state.component_moles(2, :) ./ geometry.water_volume_cm3(2);
advectiveFlux = options.internal_face_velocity_cm_s .* ...
    options.internal_face_area_cm2 .* leftConcentration;
diffusiveFlux = options.diffusion_coefficient_cm2_s .* ...
    options.internal_face_area_cm2 ./ options.internal_face_distance_cm .* ...
    (leftConcentration - rightConcentration);
expectedFlux = advectiveFlux + diffusiveFlux;
expectedDelta = expectedFlux .* 0.1;

verifyEqual(testCase, updated.component_moles(1, :), ...
    state.component_moles(1, :) - expectedDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.component_moles(2, :), ...
    state.component_moles(2, :) + expectedDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, sum(updated.component_moles, 1), ...
    sum(state.component_moles, 1), 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.internal_flux_delta_moles_total, [0 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], 'AbsTol', 1e-18);
end

function testImplicitEulerInternalTransportAllowsLargeStableStep(testCase)
state = validState();
geometry = validGeometry();
options.time_integration = 'implicit_euler';
options.internal_face_cells = [1 2];
options.internal_face_area_cm2 = 2;
options.internal_face_distance_cm = 0.5;
options.internal_face_velocity_cm_s = 0.25;
options.diffusion_coefficient_cm2_s = [0.01, 0.02];

[updated, ledger] = rtm.transport.CutCellTransportFV(state, geometry, 4, options);

verifyGreaterThanOrEqual(testCase, updated.component_moles, ...
    zeros(size(updated.component_moles)));
verifyEqual(testCase, sum(updated.component_moles, 1), ...
    sum(state.component_moles, 1), 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.internal_flux_delta_moles_total, [0 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], 'AbsTol', 1e-18);
verifyGreaterThan(testCase, updated.component_moles(2, 1), ...
    state.component_moles(2, 1));
end

function testDirichletBoundaryAddsAdvectiveAndDiffusiveFlux(testCase)
state = validState();
geometry = validGeometry();
options.boundary_face_cells = 1;
options.boundary_face_area_cm2 = 2;
options.boundary_face_distance_cm = 0.5;
options.boundary_face_velocity_cm_s = 0.25;
options.boundary_concentration_mol_cm3 = [3e-6, 1e-6];
options.diffusion_coefficient_cm2_s = [0.1, 0.2];

[updated, ledger] = rtm.transport.CutCellTransportFV(state, geometry, 4, options);

cellConcentration = state.component_moles(1, :) ./ geometry.water_volume_cm3(1);
advectiveFlux = options.boundary_face_velocity_cm_s .* ...
    options.boundary_face_area_cm2 .* options.boundary_concentration_mol_cm3;
diffusiveFlux = options.diffusion_coefficient_cm2_s .* ...
    options.boundary_face_area_cm2 ./ options.boundary_face_distance_cm .* ...
    (options.boundary_concentration_mol_cm3 - cellConcentration);
expectedDelta = (advectiveFlux + diffusiveFlux) .* 4;

verifyEqual(testCase, updated.component_moles(1, :), ...
    state.component_moles(1, :) + expectedDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.boundary_delta_moles_total, expectedDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.boundary_advective_delta_moles_total, ...
    advectiveFlux .* 4, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.boundary_diffusive_delta_moles_total, ...
    diffusiveFlux .* 4, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], 'AbsTol', 1e-18);
end

function testNegativeTransportStateIsRejected(testCase)
state = validState();
geometry = validGeometry();
options.component_source_mol_s = [-2e-6, 0; 0, 0];

verifyError(testCase, @() rtm.transport.CutCellTransportFV(state, geometry, 1, options), ...
    'RTSPHEM:Transport:NegativeComponentMoles');
end

function testTinyNegativeTransportRoundoffIsRecorded(testCase)
state = validState();
geometry = validGeometry();
options.component_source_mol_s = [-(1e-6 + 5e-16), 0; 0, 0];
expectedSuppressed = -(state.component_moles(1, 1) + ...
    options.component_source_mol_s(1, 1));

[updated, ledger] = rtm.transport.CutCellTransportFV(state, geometry, 1, options);

verifyEqual(testCase, updated.component_moles(1, 1), 0, 'AbsTol', 0);
verifyEqual(testCase, ledger.roundoff_suppressed_moles_total, expectedSuppressed, ...
    'RelTol', 1e-12, 'AbsTol', 1e-24);
verifyEqual(testCase, ledger.roundoff_suppressed_entries, 1);
verifyEqual(testCase, ledger.component_residual_moles(1), expectedSuppressed, ...
    'RelTol', 1e-12, 'AbsTol', 1e-24);
end

function state = validState()
state = struct();
state.component_names = {'Ca_total', 'C_total'};
state.component_moles = [1e-6, 2e-6; 3e-6, 4e-6];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1; 1];
state.temperature_C = [25; 25];
state.pressure_atm = [1; 1];
state.time_s = 5;
end

function geometry = validGeometry()
geometry = struct();
geometry.water_volume_cm3 = [1; 2];
geometry.cell_volume_cm3 = [1; 2];
geometry.solid_volume_cm3 = [0; 0];
geometry.fluid_fraction = [1; 1];
end
