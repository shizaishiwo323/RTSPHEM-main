function tests = test_RunDriverBenchmarkCase
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

function testRunsDriverCaseAndReportsBenchmarkObservables(testCase)
caseOptions = struct();
caseOptions.totalTime_s = 1.0;
caseOptions.configFactory = @configFactory;
caseOptions.stateFactory = @stateFactory;
caseOptions.geometryFactory = @geometryFactory;
caseOptions.connectivityFactory = @(~, ~) struct();

summary = rtm.benchmark.RunDriverBenchmarkCase(0.5, ...
    struct('name', "dt0p5", 'index', 2), caseOptions);

verifyTrue(testCase, summary.accepted);
verifyEqual(testCase, summary.run_name, "dt0p5");
verifyEqual(testCase, summary.refinement_scale, 0.5);
verifyEqual(testCase, summary.time_s, 1.0, 'AbsTol', 1e-14);
verifyEqual(testCase, summary.accepted_steps, 2);
verifyEqual(testCase, summary.rejected_steps, 0);
verifyGreaterThan(testCase, summary.final_porosity, summary.initial_porosity);
verifyGreaterThan(testCase, summary.mineral_dissolved_moles, 0);
verifyTrue(testCase, isfield(summary, 'final_solid_volume_cm3'));
verifyTrue(testCase, isfield(summary, 'solid_volume_change_cm3'));
verifyTrue(testCase, isfield(summary, 'final_surface_area_cm2'));
verifyEqual(testCase, summary.solid_volume_change_cm3, ...
    -summary.mineral_dissolved_moles, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, summary.final_solid_volume_cm3, ...
    summary.initial_solid_volume_cm3 + summary.solid_volume_change_cm3, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, summary.final_surface_area_cm2, ...
    summary.initial_surface_area_cm2, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, isfield(summary, 'mean_effective_rate_mol_cm2_s'));
verifyEqual(testCase, summary.mean_effective_rate_mol_cm2_s, ...
    summary.mineral_dissolved_moles ./ ...
    summary.initial_surface_area_cm2 ./ summary.time_s, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, summary.final_mineral_moles, ...
    summary.initial_mineral_moles - summary.mineral_dissolved_moles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, isfield(summary, 'runtime_manifest'));
verifyEqual(testCase, summary.runtime_manifest.solver_architecture, ...
    "conservative_v2");
verifyEqual(testCase, summary.runtime_manifest.chemistry_mode, ...
    "strict_molins");
verifyEqual(testCase, summary.runtime_manifest.units.component_state, "mol");

    function cfg = configFactory(refinementScale, runInfo)
        verifyEqual(testCase, refinementScale, 0.5);
        verifyEqual(testCase, runInfo.name, "dt0p5");
        cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
        cfg.chemistry.rate_constant_cm_s = 0.1;
        cfg.geometry.molarVolume_cm3_mol = 1;
        cfg.geometry.maxDisplacementOverH = 1;
        cfg.time.rt.initialDt_s = refinementScale;
        cfg.time.rt.maxDt_s = refinementScale;
    end
end

function testRejectedSubcycleMarksCaseRejected(testCase)
caseOptions = struct();
caseOptions.totalTime_s = 1.0;
caseOptions.configFactory = @rejectingConfigFactory;
caseOptions.stateFactory = @stateFactory;
caseOptions.geometryFactory = @rejectingGeometryFactory;
caseOptions.connectivityFactory = @(~, ~) struct();

summary = rtm.benchmark.RunDriverBenchmarkCase(1.0, ...
    struct('name', "rejecting"), caseOptions);

verifyFalse(testCase, summary.accepted);
verifyGreaterThan(testCase, summary.rejected_steps, 0);
verifyEqual(testCase, summary.run_name, "rejecting");
verifyTrue(testCase, isfield(summary, 'runtime_manifest'));
verifyEqual(testCase, summary.runtime_manifest.chemistry_mode, "strict_molins");
end

function testAdaptiveRejectedSubcycleCanStillAcceptCase(testCase)
caseOptions = struct();
caseOptions.totalTime_s = 1.0;
caseOptions.configFactory = @adaptiveRetryConfigFactory;
caseOptions.stateFactory = @stateFactory;
caseOptions.geometryFactory = @adaptiveRetryGeometryFactory;
caseOptions.connectivityFactory = @(~, ~) struct();

summary = rtm.benchmark.RunDriverBenchmarkCase(1.0, ...
    struct('name', "adaptive_retry"), caseOptions);

verifyTrue(testCase, summary.accepted, summary.failure_message);
verifyGreaterThan(testCase, summary.rejected_steps, 0);
verifyEqual(testCase, summary.time_s, 1.0, 'AbsTol', 1e-14);
verifyGreaterThan(testCase, summary.accepted_steps, 0);
verifyLessThanOrEqual(testCase, summary.max_displacement_over_h, 0.2);
verifyFalse(testCase, isfield(summary.solver_state, 'abort') && ...
    logical(summary.solver_state.abort));
end

function testRequiresFactoryFunctions(testCase)
verifyError(testCase, ...
    @() rtm.benchmark.RunDriverBenchmarkCase(1.0, struct(), struct()), ...
    'RTSPHEM:Benchmark:MissingDriverCaseFactory');
end

function testReportsObservationRecordsAtRequestedTimes(testCase)
caseOptions = struct();
caseOptions.totalTime_s = 1.0;
caseOptions.observationTimes_s = [0.5; 1.0];
caseOptions.configFactory = @configFactoryForObservation;
caseOptions.stateFactory = @stateFactory;
caseOptions.geometryFactory = @geometryFactory;
caseOptions.connectivityFactory = @(~, ~) struct();

summary = rtm.benchmark.RunDriverBenchmarkCase(0.5, ...
    struct('name', "dt0p5_observed"), caseOptions);

verifyEqual(testCase, summary.observation_times_s, [0.5; 1.0]);
verifyEqual(testCase, numel(summary.observation_records), 2);
verifyEqual(testCase, [summary.observation_records.time_s]', [0.5; 1.0], ...
    'AbsTol', 1e-14);
verifyLessThan(testCase, summary.observation_records(2).final_solid_volume_cm3, ...
    summary.observation_records(1).final_solid_volume_cm3);
verifyGreaterThan(testCase, summary.observation_records(1).mean_effective_rate_mol_cm2_s, 0);
verifyGreaterThan(testCase, summary.observation_records(2).mean_effective_rate_mol_cm2_s, 0);
end

function testDoesNotReportUnreachedObservationTimes(testCase)
caseOptions = struct();
caseOptions.totalTime_s = 1.0;
caseOptions.observationTimes_s = [0.5; 900];
caseOptions.configFactory = @configFactoryForObservation;
caseOptions.stateFactory = @stateFactory;
caseOptions.geometryFactory = @geometryFactory;
caseOptions.connectivityFactory = @(~, ~) struct();

summary = rtm.benchmark.RunDriverBenchmarkCase(0.5, ...
    struct('name', "dt0p5_observed"), caseOptions);

verifyEqual(testCase, summary.observation_times_s, [0.5; 900]);
verifyEqual(testCase, numel(summary.observation_records), 1);
verifyEqual(testCase, summary.observation_records.requested_time_s, 0.5, ...
    'AbsTol', 1e-14);
verifyEqual(testCase, summary.observation_records.time_s, 0.5, ...
    'AbsTol', 1e-14);
end

function testLightweightHistoryRecordsObservationsWithoutKeepingAllSteps(testCase)
caseOptions = struct();
caseOptions.totalTime_s = 2.0;
caseOptions.observationTimes_s = [0.5; 1.5; 2.0];
caseOptions.lightweightStepHistory = true;
caseOptions.configFactory = @configFactoryForObservation;
caseOptions.stateFactory = @stateFactory;
caseOptions.geometryFactory = @geometryFactory;
caseOptions.connectivityFactory = @(~, ~) struct();

summary = rtm.benchmark.RunDriverBenchmarkCase(0.5, ...
    struct('name', "dt0p5_lightweight"), caseOptions);

verifyTrue(testCase, summary.accepted);
verifyTrue(testCase, summary.lightweight_step_history);
verifyEmpty(testCase, summary.step_results);
verifyEqual(testCase, summary.accepted_steps, 4);
verifyEqual(testCase, summary.rejected_steps, 0);
verifyEqual(testCase, numel(summary.observation_records), 3);
verifyEqual(testCase, [summary.observation_records.requested_time_s]', ...
    [0.5; 1.5; 2.0], 'AbsTol', 1e-14);
verifyEqual(testCase, [summary.observation_records.time_s]', ...
    [0.5; 1.5; 2.0], 'AbsTol', 1e-14);
verifyGreaterThan(testCase, summary.max_displacement_over_h, 0);
verifyEqual(testCase, summary.max_component_mass_residual_moles, 0, ...
    'AbsTol', 1e-18);
end

function testPartIIInitializerUsesPartISteadyStateWithoutMovingGeometry(testCase)
caseOptions = struct();
caseOptions.totalTime_s = 0;
caseOptions.steadyOptions = struct('window_s', 1.0, 'maxWindows', 1, ...
    'rateAbsoluteTolerance_mol_s', 0, 'rateRelativeTolerance', 0, ...
    'requiredConsecutivePasses', 1);
caseOptions.configFactory = @partIIConfigFactory;
caseOptions.stateFactory = @largeReactantStateFactory;
caseOptions.geometryFactory = @largeSolidGeometryFactory;
caseOptions.connectivityFactory = @(~, ~) struct();

summary = rtm.benchmark.RunDriverBenchmarkCase(1.0, ...
    struct('name', "partII_initialized"), caseOptions);

verifyTrue(testCase, isfield(summary, 'steady_initializer'));
verifyEqual(testCase, summary.steady_initializer.windows, 1);
verifyLessThan(testCase, summary.state.component_moles, 10);
verifyEqual(testCase, summary.state.mineral_moles, 100, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, summary.geometry.solid_volume_cm3, 100, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, summary.initial_mineral_moles, 100, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, summary.final_mineral_moles, 100, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testIntegrationInitializerDoesNotLeaveMolesInDryCells(testCase)
matrix = struct();
matrix.time_steps_s = 0.054;
matrix.grid_resolutions = "128x64";
suiteOptions = rtm.benchmark.CreateMolinsDriverCaseOptions( ...
    'integration_phreeqc', struct( ...
    'acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true, ...
    'phreeqcRunBatchFunction', @rtm.benchmark.MockPhreeqcBatch, ...
    'totalTime_s', 0.054, ...
    'circle_cut_subsamples', 64, ...
    'interface_arc_samples', 8192));

summary = rtm.benchmark.RunDriverBenchmarkCase( ...
    suiteOptions.refinementScales(1), ...
    struct('name', suiteOptions.runNames(1), 'index', 1), ...
    suiteOptions.caseOptions);

dryCells = summary.geometry.water_volume_cm3(:) <= 0;
verifyTrue(testCase, summary.accepted, summary.failure_message);
verifyEqual(testCase, summary.rejected_steps, 0);
verifyEqual(testCase, summary.state.component_moles(dryCells, :), ...
    zeros(nnz(dryCells), numel(summary.state.component_names)), ...
    'AbsTol', 1e-18);
end

function testReportsAcceptanceMeshMetadata(testCase)
matrix = struct();
matrix.time_steps_s = [54; 5.4; 0.54; 0.054];
matrix.grid_resolutions = ["128x64"; "256x128"; "512x256"];
suiteOptions = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('refinementScales', matrix.time_steps_s, ...
    'acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true, ...
    'totalTime_s', 0));

summary = rtm.benchmark.RunDriverBenchmarkCase(54, ...
    struct('name', "dt_54"), suiteOptions.caseOptions);

verifyTrue(testCase, isfield(summary, 'benchmark_mesh'));
verifyEqual(testCase, summary.benchmark_mesh.nx, 128);
verifyEqual(testCase, summary.benchmark_mesh.ny, 64);
verifyEqual(testCase, summary.benchmark_mesh.resolution_label, "128x64");
verifyEqual(testCase, summary.benchmark_mesh.cell_count, 128 * 64);
verifyEqual(testCase, summary.initial_surface_area_cm2, ...
    2 * pi * 0.01, 'RelTol', 1e-12);
verifyGreaterThan(testCase, summary.initial_solid_volume_cm3, 0);
verifyTrue(testCase, isfield(summary, 'final_solid_volume_cm3'));
verifyTrue(testCase, isfield(summary, 'final_surface_area_cm2'));
verifyEqual(testCase, summary.final_solid_volume_cm3, ...
    summary.initial_solid_volume_cm3, 'RelTol', 1e-12);
verifyEqual(testCase, summary.final_surface_area_cm2, ...
    summary.initial_surface_area_cm2, 'RelTol', 1e-12);
verifyTrue(testCase, isnan(summary.mean_effective_rate_mol_cm2_s));
verifyTrue(testCase, isfield(summary, 'acceptance_matrix'));
verifyEqual(testCase, summary.acceptance_matrix.time_steps_s, ...
    [54; 5.4; 0.54; 0.054]);
end

function cfg = partIIConfigFactory(refinementScale, ~)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partII_strict');
cfg.chemistry.rate_constant_cm_s = 0.1;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1e6;
cfg.time.rt.initialDt_s = refinementScale;
cfg.time.rt.maxDt_s = refinementScale;
end

function cfg = configFactoryForObservation(refinementScale, ~)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partII_strict');
cfg.chemistry.rate_constant_cm_s = 0.1;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 1;
cfg.time.rt.initialDt_s = refinementScale;
cfg.time.rt.maxDt_s = refinementScale;
end

function cfg = rejectingConfigFactory(~, ~)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 10;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 0.01;
cfg.failure.maxRetries = 1;
cfg.time.rt.initialDt_s = 1.0;
cfg.time.rt.maxDt_s = 1.0;
end

function cfg = adaptiveRetryConfigFactory(~, ~)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0.5;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 0.2;
cfg.failure.maxRetries = 6;
cfg.time.rt.initialDt_s = 1.0;
cfg.time.rt.maxDt_s = 1.0;
cfg.time.rt.maxReactantFraction = Inf;
cfg.time.geometry.maxMineralFraction = Inf;
end

function state = stateFactory(~, ~)
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = 1e-3;
state.mineral_names = {'Calcite'};
state.mineral_moles = 1;
state.temperature_C = 25;
state.pressure_atm = 1;
state.time_s = 0;
end

function state = largeReactantStateFactory(~, ~)
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = 10;
state.mineral_names = {'Calcite'};
state.mineral_moles = 100;
state.temperature_C = 25;
state.pressure_atm = 1;
state.time_s = 0;
end

function geometry = geometryFactory(~, ~)
geometry = struct();
geometry.water_volume_cm3 = 1;
geometry.solid_volume_cm3 = 1;
geometry.cell_volume_cm3 = 2;
geometry.fluid_fraction = 0.5;
geometry.cell_centroid_cm = [0 0];
geometry.interface_centroid_cm = [0 0];
geometry.interface_area_cm2 = 1;
geometry.interface_h_cm = 1;
geometry.interface_normal = [1 0];
geometry.active_fluid_cell = true;
end

function geometry = largeSolidGeometryFactory(~, ~)
geometry = geometryFactory([], []);
geometry.solid_volume_cm3 = 100;
geometry.cell_volume_cm3 = geometry.water_volume_cm3 + geometry.solid_volume_cm3;
geometry.fluid_fraction = geometry.water_volume_cm3 ./ geometry.cell_volume_cm3;
geometry.interface_h_cm = 1;
end

function geometry = rejectingGeometryFactory(refinementScale, runInfo)
geometry = geometryFactory(refinementScale, runInfo);
geometry.interface_h_cm = 1e-4;
end

function geometry = adaptiveRetryGeometryFactory(refinementScale, runInfo)
geometry = geometryFactory(refinementScale, runInfo);
geometry.interface_h_cm = 1e-3;
end
