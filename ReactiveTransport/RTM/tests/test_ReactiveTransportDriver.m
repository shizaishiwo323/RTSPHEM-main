function tests = test_ReactiveTransportDriver
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc', 'tests'));
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testStrictMolinsAcceptedStepCommitsState(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0.1;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 0.25;
state = strictState(1e-6, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyTrue(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transaction_status, 'committed');
verifyLessThan(testCase, result.state.component_moles(1), state.component_moles(1));
verifyLessThan(testCase, result.state.mineral_moles(1), state.mineral_moles(1));
verifyLessThan(testCase, result.geometry.solid_volume_cm3, ...
    geometry.solid_volume_cm3);
verifyGreaterThan(testCase, result.geometry.water_volume_cm3, ...
    geometry.water_volume_cm3);
verifyEqual(testCase, result.geometry.cell_volume_cm3, ...
    geometry.cell_volume_cm3, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.solver_state.dt_s, 1);
verifyEqual(testCase, result.geometry_info.expected_solid_volume_change_cm3, ...
    -result.reaction_result.realized_interface_moles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, isfield(result.geometry_info, ...
    'cell_actual_solid_volume_change_cm3'));
verifyEqual(testCase, result.geometry_info.cell_actual_solid_volume_change_cm3, ...
    result.geometry.solid_volume_cm3 - geometry.solid_volume_cm3, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.geometry_info.actual_solid_volume_change_cm3, ...
    result.geometry.solid_volume_cm3 - geometry.solid_volume_cm3, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testStrictMolinsRejectedStepRollsBackAndShrinksDt(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 10;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 0.01;
cfg.failure.shrinkFactor = 0.5;
cfg.failure.minDt_s = 1e-8;
cfg.failure.maxRetries = 12;
state = strictState(1e-3, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyFalse(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transaction_status, 'rolled_back');
verifyEqual(testCase, result.state.component_moles, state.component_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.state.mineral_moles, state.mineral_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.solver_state.dt_s, 0.5);
verifyEqual(testCase, result.solver_state.rejected_steps, 1);
verifyTrue(testCase, any(result.diagnostics.reasons == "geometry displacement exceeds tolerance"));
end

function testStrictMolinsDriverPropagatesReactantFractionLimit(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 10;
cfg.time.rt.maxReactantFraction = 0.2;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1e6;
state = strictState(1e-3, 1);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyTrue(testCase, result.diagnostics.accepted);
verifyGreaterThan(testCase, result.reaction_result.candidate_interface_moles, ...
    state.component_moles);
verifyEqual(testCase, result.reaction_result.realized_interface_moles, ...
    2e-4, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.component_moles, 8e-4, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testSolidVolumeOvershootRejectsAndRollsBackGeometry(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0.1;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1e6;
cfg.failure.shrinkFactor = 0.5;
state = strictState(10, 10);
geometry = singleInterfaceGeometry();
geometry.solid_volume_cm3 = 1e-6;
geometry.cell_volume_cm3 = geometry.water_volume_cm3 + geometry.solid_volume_cm3;
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyFalse(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transaction_status, 'rolled_back');
verifyEqual(testCase, result.geometry.solid_volume_cm3, geometry.solid_volume_cm3, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.geometry.water_volume_cm3, geometry.water_volume_cm3, ...
    'AbsTol', 1e-18);
verifyTrue(testCase, any(result.diagnostics.reasons == ...
    "solid volume would become negative"));
end

function testTransportStepRunsBeforeStrictChemistry(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0.1;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 0.25;
cfg.transport.options.component_source_mol_s = 1e-7;
state = strictState(1e-6, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyTrue(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transport_ledger.source_delta_moles_total, 1e-7, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyGreaterThan(testCase, result.reaction_result.candidate_interface_moles, 1e-7);
verifyEqual(testCase, result.diagnostics.component_absolute_residual_moles, 0, ...
    'AbsTol', 1e-18);
end

function testBoundaryDirichletFluxIsIncludedInDriverMassLedger(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0;
cfg.transport.options.boundary_face_cells = 1;
cfg.transport.options.boundary_face_area_cm2 = 2;
cfg.transport.options.boundary_face_distance_cm = 0.5;
cfg.transport.options.boundary_face_velocity_cm_s = 0.25;
cfg.transport.options.boundary_concentration_mol_cm3 = 3e-6;
cfg.transport.options.diffusion_coefficient_cm2_s = 0.1;
state = strictState(1e-6, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(4);

cellConcentration = state.component_moles ./ geometry.water_volume_cm3;
advectiveDelta = 0.25 .* 2 .* 3e-6 .* 4;
diffusiveDelta = 0.1 .* 2 ./ 0.5 .* (3e-6 - cellConcentration) .* 4;
expectedBoundaryDelta = advectiveDelta + diffusiveDelta;

verifyTrue(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transport_ledger.boundary_delta_moles_total, ...
    expectedBoundaryDelta, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.component_moles, ...
    state.component_moles + expectedBoundaryDelta, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.diagnostics.component_absolute_residual_moles, ...
    0, 'AbsTol', 1e-18);
end

function testPhreeqcKineticsModeRunsThroughSharedDriver(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.benchmark.enabled = false;
cfg.chemistry.mode = 'phreeqc_kinetics';
cfg.chemistry.options = struct();
cfg.chemistry.options.h_mol_cm3 = 1e-7;
cfg.chemistry.options.runBatchFunction = @mockRunBatch;
cfg.phreeqc.engine = 'mock';
cfg.phreeqc.databasePolicy = 'not_used';
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;
state = carbonateState(0, 0, 1e-8, 1e-8, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

expectedDissolved = 2e-8;
verifyTrue(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.reaction_result.realized_interface_moles, ...
    expectedDissolved, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.component_moles(1, 1), ...
    expectedDissolved, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.component_moles(1, 2), ...
    expectedDissolved, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.mineral_moles, ...
    state.mineral_moles - expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);

    function batchResult = mockRunBatch(batchState, batchOptions)
        verifyEqual(testCase, batchOptions.rateLaw, 'database_calcite');
        dissolved = 2e-8;
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function testExternalTstPhreeqcModeRunsThroughSharedDriver(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.benchmark.enabled = false;
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.chemistry.rate_constant_cm_s = 1e-4;
cfg.chemistry.options = struct();
cfg.chemistry.options.h_mol_cm3 = 1e-7;
cfg.chemistry.options.h_activity_mol_cm3 = 2e-7;
cfg.chemistry.options.runBatchFunction = @mockRunBatch;
cfg.phreeqc.engine = 'mock';
cfg.phreeqc.databasePolicy = 'not_used';
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;
state = carbonateState(0, 0, 1e-8, 1e-8, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

expectedDissolved = 2e-7 * 1000 * 1e-4;
verifyTrue(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transaction_status, 'committed');
verifyEqual(testCase, result.reaction_result.realized_interface_moles, ...
    expectedDissolved, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.component_moles(1, 1), ...
    expectedDissolved, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.component_moles(1, 2), ...
    expectedDissolved, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.state.mineral_moles, ...
    state.mineral_moles - expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function testExternalTstDriverBuildsDisjointReactionClusters(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.benchmark.enabled = false;
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.chemistry.rate_constant_cm_s = 1e-4;
cfg.chemistry.reaction_cluster_depth_cells = 1;
cfg.chemistry.options = struct();
cfg.chemistry.options.h_mol_cm3 = [1e-7; 1e-7; 1e-7];
cfg.chemistry.options.h_activity_mol_cm3 = [1e-7; 1e-7; 1e-7];
cfg.chemistry.options.runBatchFunction = @mockRunBatch;
cfg.phreeqc.engine = 'mock';
cfg.phreeqc.databasePolicy = 'not_used';
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;
state = carbonateState(0, 0, 3e-8, 3e-8, 3e-3);
state.component_moles = repmat(state.component_moles ./ 3, 3, 1);
state.mineral_moles = repmat(state.mineral_moles ./ 3, 3, 1);
state.temperature_C = repmat(25, 3, 1);
state.pressure_atm = repmat(1, 3, 1);
geometry = threeCellSharedClusterGeometry();
connectivity = struct('cell_neighbors', {{2; [1; 3]; 2}});
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, connectivity);

result = driver.runOneStep(1);

verifyTrue(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.reaction_result.cluster_ids, [1; 1; 1]);
verifyEqual(testCase, result.reaction_result.aux.reaction_cluster_count, 1);
verifyEqual(testCase, ...
    result.reaction_result.aux.reaction_cluster_max_membership, 1);

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function testPhreeqcChargeResidualRejectsAndRollsBack(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.benchmark.enabled = false;
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.chemistry.rate_constant_cm_s = 1e-4;
cfg.chemistry.chargeAbsoluteTolerance_eq = 1e-8;
cfg.chemistry.options = struct();
cfg.chemistry.options.h_mol_cm3 = 1e-7;
cfg.chemistry.options.h_activity_mol_cm3 = 2e-7;
cfg.chemistry.options.runBatchFunction = @mockRunBatch;
cfg.phreeqc.engine = 'mock';
cfg.phreeqc.databasePolicy = 'not_used';
cfg.failure.shrinkFactor = 0.5;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;
state = carbonateState(0, 0, 1e-8, 1e-8, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyFalse(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transaction_status, 'rolled_back');
verifyEqual(testCase, result.state.component_moles, state.component_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.geometry.solid_volume_cm3, geometry.solid_volume_cm3, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.diagnostics.charge_absolute_residual_eq, 2e-6, ...
    'AbsTol', 1e-18);
verifyTrue(testCase, any(result.diagnostics.reasons == ...
    "chemistry charge balance residual exceeds tolerance"));

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
        batchResult.chargeBalance = 2e-6;
    end
end

function testPhreeqcFailedCellsRejectAndRollBack(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.benchmark.enabled = false;
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.chemistry.rate_constant_cm_s = 1e-4;
cfg.chemistry.options = struct();
cfg.chemistry.options.h_mol_cm3 = 1e-7;
cfg.chemistry.options.h_activity_mol_cm3 = 2e-7;
cfg.chemistry.options.runBatchFunction = @mockRunBatch;
cfg.phreeqc.engine = 'mock';
cfg.phreeqc.databasePolicy = 'not_used';
cfg.failure.shrinkFactor = 0.5;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;
state = carbonateState(0, 0, 1e-8, 1e-8, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyFalse(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transaction_status, 'rolled_back');
verifyEqual(testCase, result.state.component_moles, state.component_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.diagnostics.failed_cells, [4; 9]);
verifyTrue(testCase, any(result.diagnostics.reasons == ...
    "chemistry has failed cells"));

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
        batchResult.chargeBalance = 0;
        batchResult.failedCells = [4; 9];
        batchResult.errorMessage = "mock failed cells";
    end
end

function testNegativeStateUpdateRejectsAndRollsBack(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.benchmark.enabled = false;
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.chemistry.rate_constant_cm_s = 1e-4;
cfg.chemistry.options = struct();
cfg.chemistry.options.h_mol_cm3 = 1e-7;
cfg.chemistry.options.h_activity_mol_cm3 = 2e-7;
cfg.chemistry.options.runBatchFunction = @mockRunBatch;
cfg.phreeqc.engine = 'mock';
cfg.phreeqc.databasePolicy = 'not_used';
cfg.failure.shrinkFactor = 0.5;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;
state = carbonateState(1e-8, 0, 1e-8, 1e-8, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

result = driver.runOneStep(1);

verifyFalse(testCase, result.diagnostics.accepted);
verifyEqual(testCase, result.transaction_status, 'rolled_back');
verifyEqual(testCase, result.state.component_moles, state.component_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.geometry.solid_volume_cm3, geometry.solid_volume_cm3, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.solver_state.dt_s, 0.5, 'AbsTol', 1e-18);
verifyTrue(testCase, any(contains(result.diagnostics.reasons, ...
    "state update failed")));
verifyEqual(testCase, result.diagnostics.error_message, ...
    "component_moles must be nonnegative.");

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = -1e-8;
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
        batchResult.chargeBalance = 0;
    end
end

function testExternalTstPhreeqcDriverReusesPersistentSession(testCase)
databasePath = createTempDatabase(testCase);
mockEngine = MockIPhreeqcEngine(minimalSelectedOutput());

cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.benchmark.enabled = false;
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.chemistry.rate_constant_cm_s = 1e-4;
cfg.chemistry.options = struct();
cfg.chemistry.options.h_mol_cm3 = 1e-7;
cfg.chemistry.options.h_activity_mol_cm3 = 2e-7;
cfg.phreeqc.engine = 'mock';
cfg.phreeqc.databasePolicy = 'exact_local';
cfg.phreeqc.databasePath = databasePath;
cfg.phreeqc.persistSession = true;
cfg.phreeqc.useRunString = true;
cfg.phreeqc.engineFactory = @() mockEngine;
cfg.time.rt.initialDt_s = 0.5;
cfg.time.rt.maxDt_s = 0.5;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;

state = carbonateState(0, 0, 1e-8, 1e-8, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());
cleanupObj = onCleanup(@() delete(driver));

summary = driver.runRtSubcycles(1.0);

verifyEqual(testCase, summary.accepted_steps, 2);
verifyEqual(testCase, mockEngine.LoadDatabaseCount, 1);
verifyEqual(testCase, mockEngine.RunStringCount, 2);
verifyEqual(testCase, mockEngine.LoadedDatabasePaths, string(databasePath));
verifyTrue(testCase, all(arrayfun(@(result) ...
    result.reaction_result.aux.phreeqc_session_reused, summary.step_results)));
verifyTrue(testCase, all(arrayfun(@(result) ...
    result.reaction_result.aux.phreeqc_input_written == false, summary.step_results)));
clear cleanupObj;
end

function testRunRtSubcyclesAdvancesToRequestedTime(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0;
cfg.time.rt.initialDt_s = 0.4;
cfg.time.rt.maxDt_s = 0.4;
state = strictState(1e-6, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

summary = driver.runRtSubcycles(1.0);

verifyEqual(testCase, summary.accepted_steps, 3);
verifyEqual(testCase, summary.rejected_steps, 0);
verifyEqual(testCase, summary.time_s, 1.0, 'AbsTol', 1e-14);
verifyEqual(testCase, summary.state.time_s, 1.0, 'AbsTol', 1e-14);
accepted = arrayfun(@(result) result.diagnostics.accepted, summary.step_results);
verifyEqual(testCase, accepted, [true true true]);
end

function testRunRtSubcyclesRetriesRejectedStepWithSmallerDt(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0.002;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 0.015;
cfg.failure.shrinkFactor = 0.5;
cfg.failure.minDt_s = 1e-8;
cfg.failure.maxRetries = 4;
cfg.time.rt.initialDt_s = 1.0;
cfg.time.rt.maxDt_s = 1.0;
state = strictState(1e-3, 1);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

summary = driver.runRtSubcycles(1.0);

verifyEqual(testCase, summary.rejected_steps, 1);
verifyEqual(testCase, summary.accepted_steps, 2);
verifyEqual(testCase, summary.time_s, 1.0, 'AbsTol', 1e-14);
verifyEqual(testCase, summary.step_results(1).transaction_status, 'rolled_back');
verifyEqual(testCase, summary.step_results(2).transaction_status, 'committed');
verifyEqual(testCase, summary.step_results(3).transaction_status, 'committed');
end

function testRunRtSubcyclesWritesDiagnosticCsvWhenConfigured(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupDiagnostics(outputDir));

cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0;
cfg.time.rt.initialDt_s = 0.5;
cfg.time.rt.maxDt_s = 0.5;
cfg.output.diagnosticsDir = outputDir;
state = strictState(1e-6, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

summary = driver.runRtSubcycles(1.0);

verifyTrue(testCase, isfield(summary, 'diagnostic_paths'));
verifyTrue(testCase, exist(summary.diagnostic_paths.mass_balance_components, 'file') == 2);
verifyTrue(testCase, exist(summary.diagnostic_paths.solid_geometry_balance, 'file') == 2);
verifyTrue(testCase, exist(summary.diagnostic_paths.chemistry_step_status, 'file') == 2);
verifyTrue(testCase, exist(summary.diagnostic_paths.step_rejection_log, 'file') == 2);
verifyTrue(testCase, exist(summary.diagnostic_paths.transport_flux_balance, 'file') == 2);
verifyTrue(testCase, exist(summary.diagnostic_paths.reaction_cluster_diagnostics, 'file') == 2);
verifyTrue(testCase, exist(summary.diagnostic_paths.runtime_manifest, 'file') == 2);
massTable = readtable(summary.diagnostic_paths.mass_balance_components, ...
    'TextType', 'string');
verifyEqual(testCase, height(massTable), 2);
verifyEqual(testCase, massTable.step_index, [1; 2]);
manifest = jsondecode(fileread(summary.diagnostic_paths.runtime_manifest));
verifyEqual(testCase, string(manifest.solver_architecture), "conservative_v2");
verifyEqual(testCase, string(manifest.chemistry_mode), "strict_molins");
end

function testRunRtSubcyclesWritesDiagnosticsBeforeAbort(testCase)
outputDir = tempname;
cleanup = onCleanup(@() cleanupDiagnostics(outputDir));

cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 10;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 0.01;
cfg.failure.maxRetries = 0;
cfg.output.diagnosticsDir = outputDir;
state = strictState(1e-3, 1e-3);
geometry = singleInterfaceGeometry();
driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, struct());

verifyError(testCase, @() driver.runRtSubcycles(1.0), ...
    'RTSPHEM:Driver:SubcycleRetryAborted');

rejectionPath = fullfile(outputDir, 'step_rejection_log.csv');
manifestPath = fullfile(outputDir, 'runtime_manifest.json');
verifyTrue(testCase, exist(rejectionPath, 'file') == 2);
verifyTrue(testCase, exist(manifestPath, 'file') == 2);
rejectionTable = readtable(rejectionPath, 'TextType', 'string');
verifyEqual(testCase, height(rejectionTable), 1);
verifyTrue(testCase, contains(rejectionTable.reason(1), ...
    "geometry displacement exceeds tolerance"));
end

function state = strictState(hMoles, calciteMoles)
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = hMoles;
state.mineral_names = {'Calcite'};
state.mineral_moles = calciteMoles;
state.temperature_C = 25;
state.pressure_atm = 1;
state.time_s = 0;
end

function state = carbonateState(caMoles, cMoles, naMoles, clMoles, calciteMoles)
state = struct();
state.component_names = {'Ca', 'C', 'Na', 'Cl'};
state.component_moles = [caMoles, cMoles, naMoles, clMoles];
state.mineral_names = {'Calcite'};
state.mineral_moles = calciteMoles;
state.temperature_C = 25;
state.pressure_atm = 1;
state.time_s = 0;
end

function geometry = singleInterfaceGeometry()
geometry = struct();
geometry.water_volume_cm3 = 1;
geometry.solid_volume_cm3 = 1e-3;
geometry.cell_volume_cm3 = 1.001;
geometry.fluid_fraction = 1;
geometry.cell_centroid_cm = [0 0];
geometry.interface_centroid_cm = [0 0];
geometry.interface_area_cm2 = 1;
geometry.interface_h_cm = 1e-4;
geometry.interface_normal = [1 0];
geometry.active_fluid_cell = true;
end

function geometry = threeCellSharedClusterGeometry()
geometry = struct();
geometry.water_volume_cm3 = [1; 1; 1];
geometry.solid_volume_cm3 = [1e-3; 0; 1e-3];
geometry.cell_volume_cm3 = geometry.water_volume_cm3 + geometry.solid_volume_cm3;
geometry.fluid_fraction = geometry.water_volume_cm3 ./ geometry.cell_volume_cm3;
geometry.cell_centroid_cm = [0 0; 1 0; 2 0];
geometry.interface_centroid_cm = [0 0; NaN NaN; 2 0];
geometry.interface_area_cm2 = [1; 0; 1];
geometry.interface_h_cm = [1e-4; 1e-4; 1e-4];
geometry.interface_normal = [1 0; NaN NaN; -1 0];
geometry.active_fluid_cell = true(3, 1);
end

function cleanupDiagnostics(outputDir)
massPath = fullfile(outputDir, 'mass_balance_components.csv');
solidPath = fullfile(outputDir, 'solid_geometry_balance.csv');
chemistryPath = fullfile(outputDir, 'chemistry_step_status.csv');
rejectionPath = fullfile(outputDir, 'step_rejection_log.csv');
fluxPath = fullfile(outputDir, 'transport_flux_balance.csv');
clusterPath = fullfile(outputDir, 'reaction_cluster_diagnostics.csv');
manifestPath = fullfile(outputDir, 'runtime_manifest.json');
if exist(massPath, 'file') == 2
    delete(massPath);
end
if exist(solidPath, 'file') == 2
    delete(solidPath);
end
if exist(chemistryPath, 'file') == 2
    delete(chemistryPath);
end
if exist(rejectionPath, 'file') == 2
    delete(rejectionPath);
end
if exist(fluxPath, 'file') == 2
    delete(fluxPath);
end
if exist(clusterPath, 'file') == 2
    delete(clusterPath);
end
if exist(manifestPath, 'file') == 2
    delete(manifestPath);
end
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end

function databasePath = createTempDatabase(testCase)
databasePath = [tempname, '.dat'];
fid = fopen(databasePath, 'w');
if fid == -1
    error('test_ReactiveTransportDriver:TempFileOpenFailed', ...
        'Cannot create temp database.');
end
fprintf(fid, 'TITLE mock database\n');
fclose(fid);
testCase.addTeardown(@() deleteExistingFile(databasePath));
end

function deleteExistingFile(pathValue)
if exist(pathValue, 'file') == 2
    delete(pathValue);
end
end

function rawOutput = minimalSelectedOutput()
rawOutput = {
    'sim', 'state', 'soln', 'pH', 'charge(eq)', 'Ca(mol/kgw)', 'C(mol/kgw)', 'Na(mol/kgw)', 'Cl(mol/kgw)', 'm_H+(mol/kgw)', 'm_Ca+2(mol/kgw)', 'm_HCO3-(mol/kgw)', 'm_CO3-2(mol/kgw)', 'm_Cl-(mol/kgw)', 'm_Na+(mol/kgw)', 'si_Calcite', 'KIN_DELTA_Calcite', 'RATE_Calcite';
    1, 'react', 1, 2.1, 0, 1e-4, 1e-4, 0.01, 0.11, 0.08, 1e-4, 8e-5, 1e-8, 0.11, 0.01, -3, 0, 0
    };
end
