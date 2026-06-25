function tests = test_CreateMolinsDriverCaseOptions
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

function testPartIStrictOptionsRunDriverConvergenceSmoke(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupOutput(outputDir));

options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('refinementScales', [0.5, 0.25], ...
    'totalTime_s', 0.5, ...
    'outputDir', outputDir));

verifyEqual(testCase, options.suiteName, "molins_partI_strict_driver");
verifyEqual(testCase, options.observableName, "final_porosity");
verifyEqual(testCase, options.refinementScales, [0.5; 0.25]);
verifyEqual(testCase, options.runNames, ["dt_0p5"; "dt_0p25"]);
verifyTrue(testCase, isa(options.runFunction, 'function_handle'));
verifyTrue(testCase, isa(options.observableFunction, 'function_handle'));

suite = rtm.benchmark.RunConvergenceSuite(options);

verifyEqual(testCase, numel(suite.runs), 2);
verifyTrue(testCase, all([suite.runs.accepted]));
verifyTrue(testCase, exist(suite.report_path, 'file') == 2);
verifyGreaterThan(testCase, suite.runs(1).summary.final_porosity, ...
    suite.runs(1).summary.initial_porosity);
end

function testSupportsMeanEffectiveRateObservable(testCase)
options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('refinementScales', [0.5, 0.25], ...
    'observableName', 'mean_effective_rate_mol_cm2_s'));

summary = struct();
summary.mean_effective_rate_mol_cm2_s = 4.32e-8;

verifyEqual(testCase, options.observableName, "mean_effective_rate_mol_cm2_s");
verifyEqual(testCase, options.observableFunction(summary), 4.32e-8, ...
    'AbsTol', 1e-18);
end

function testPropagatesReactionFractionLimitOptions(testCase)
options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('refinementScales', 0.1, ...
    'maxReactantFraction', Inf, ...
    'maxMineralFraction', 0.05, ...
    'reactionClusterDepthCells', 2));

cfg = options.caseOptions.configFactory(0.1, struct('name', "dt_0p1"));

verifyEqual(testCase, cfg.time.rt.maxReactantFraction, Inf);
verifyEqual(testCase, cfg.time.geometry.maxMineralFraction, 0.05);
verifyEqual(testCase, cfg.chemistry.reaction_cluster_depth_cells, 2);
end

function testIntegrationPhreeqcOptionsUseExternalTstModeWithChemistryInput(testCase)
options = rtm.benchmark.CreateMolinsDriverCaseOptions('integration_phreeqc', ...
    struct('refinementScales', 0.1, ...
    'h_mol_cm3', 1e-7, ...
    'h_activity_mol_cm3', 2e-7, ...
    'phreeqcRunBatchFunction', @mockRunBatch));

cfg = options.caseOptions.configFactory(0.1, struct('name', "dt_0p1"));
state = options.caseOptions.stateFactory(0.1, struct());
geometry = options.caseOptions.geometryFactory(0.1, struct());

verifyEqual(testCase, cfg.chemistry.mode, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.chemistry.options.h_mol_cm3, 1e-7);
verifyEqual(testCase, cfg.chemistry.options.h_activity_mol_cm3, 2e-7);
verifyEqual(testCase, cfg.phreeqc.engine, 'mock');
verifyEqual(testCase, state.component_names, {'Ca', 'C', 'Na', 'Cl'});
verifyTrue(testCase, isfield(geometry, 'interface_area_cm2'));

    function result = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        result = struct();
        result.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        result.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        result.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        result.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        result.calciteDissolvedMoles = dissolved;
    end
end

function testIntegrationPhreeqcOptionsPassExactLocalDatabasePath(testCase)
options = rtm.benchmark.CreateMolinsDriverCaseOptions('integration_phreeqc', ...
    struct('refinementScales', [0.5, 0.25]));

cfg = options.caseOptions.configFactory(0.5, struct('name', "dt_0p5"));

verifyEqual(testCase, cfg.phreeqc.engine, 'iphreeqc_com');
verifyEqual(testCase, cfg.phreeqc.databasePolicy, 'exact_local');
verifyEqual(testCase, cfg.phreeqc.databaseName, 'phreeqc-m.dat');
verifyTrue(testCase, isfield(cfg.phreeqc, 'databasePath'));
verifyTrue(testCase, isfield(cfg.chemistry.options, 'databasePath'));
verifyEqual(testCase, string(cfg.phreeqc.databasePath), ...
    string(cfg.chemistry.options.databasePath));
verifyTrue(testCase, exist(cfg.chemistry.options.databasePath, 'file') == 2);
verifyTrue(testCase, endsWith(string(cfg.chemistry.options.databasePath), ...
    fullfile('phreeqc', 'database', 'phreeqc-m.dat')));
end

function testIntegrationAcceptanceGridDoesNotPlaceSoluteInDryCells(testCase)
matrix = struct();
matrix.time_steps_s = 0.054;
matrix.grid_resolutions = "128x64";

options = rtm.benchmark.CreateMolinsDriverCaseOptions('integration_phreeqc', ...
    struct('acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true, ...
    'phreeqcRunBatchFunction', @rtm.benchmark.MockPhreeqcBatch));

state = options.caseOptions.stateFactory(0.054, ...
    struct('name', "grid_128x64_dt_0p054", 'index', 1));
geometry = options.caseOptions.geometryFactory(0.054, ...
    struct('name', "grid_128x64_dt_0p054", 'index', 1));

dryCells = geometry.water_volume_cm3(:) <= 0;
verifyGreaterThan(testCase, nnz(dryCells), 0);
verifyEqual(testCase, state.component_moles(dryCells, :), ...
    zeros(nnz(dryCells), numel(state.component_names)), 'AbsTol', 0);
end

function testIntegrationAcceptanceGridBoundaryConcentrationMatchesComponents(testCase)
matrix = struct();
matrix.time_steps_s = 0.054;
matrix.grid_resolutions = "128x64";

options = rtm.benchmark.CreateMolinsDriverCaseOptions('integration_phreeqc', ...
    struct('acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true, ...
    'phreeqcRunBatchFunction', @rtm.benchmark.MockPhreeqcBatch));

cfg = options.caseOptions.configFactory(0.054, ...
    struct('name', "grid_128x64_dt_0p054", 'index', 1));

verifyEqual(testCase, size(cfg.transport.options.boundary_concentration_mol_cm3, 2), 4);
verifyEqual(testCase, cfg.transport.options.boundary_concentration_mol_cm3, ...
    zeros(size(cfg.transport.options.boundary_concentration_mol_cm3)), ...
    'AbsTol', 0);
end

function testLightweightStepHistoryOptionPropagatesToCaseOptions(testCase)
options = rtm.benchmark.CreateMolinsDriverCaseOptions('partII_strict', ...
    struct('lightweightStepHistory', true));

verifyTrue(testCase, options.caseOptions.lightweightStepHistory);
end

function testAcceptanceMatrixBuildsMolinsGridFactories(testCase)
matrix = struct();
matrix.time_steps_s = [54; 5.4; 0.54; 0.054];
matrix.grid_resolutions = ["128x64"; "256x128"; "512x256"];
matrix.partII_observation_times_min = [15; 30; 45];

options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('refinementScales', matrix.time_steps_s, ...
    'acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true));

cfg = options.caseOptions.configFactory(54, struct('name', "dt_54"));
state = options.caseOptions.stateFactory(54, struct('name', "dt_54"));
geometry = options.caseOptions.geometryFactory(54, struct('name', "dt_54"));
connectivity = options.caseOptions.connectivityFactory(54, struct('name', "dt_54"));

verifyEqual(testCase, cfg.benchmark.mesh.nx, 128);
verifyEqual(testCase, cfg.benchmark.mesh.ny, 64);
verifyEqual(testCase, cfg.benchmark.mesh.resolution_label, "128x64");
verifyEqual(testCase, cfg.benchmark.mesh.domain_length_cm, 0.1);
verifyEqual(testCase, cfg.benchmark.mesh.domain_height_cm, 0.05);
verifyEqual(testCase, cfg.benchmark.mesh.circle_center_cm, [0.05 0.025]);
verifyEqual(testCase, cfg.chemistry.reaction_time_integration, 'exact_first_order');
verifyEqual(testCase, cfg.chemistry.reaction_cluster_depth_cells, 0);
verifyEqual(testCase, cfg.time.rt.maxReactantFraction, Inf);
verifyEqual(testCase, numel(geometry.water_volume_cm3), 128 * 64);
verifyEqual(testCase, size(state.component_moles, 1), 128 * 64);
verifyEqual(testCase, size(state.mineral_moles, 1), 128 * 64);
verifyTrue(testCase, all(geometry.water_volume_cm3 >= 0));
verifyTrue(testCase, all(geometry.solid_volume_cm3 >= 0));
verifyEqual(testCase, geometry.mesh_resolution, [128 64]);
verifyEqual(testCase, string(geometry.mesh_resolution_label), "128x64");
verifyEqual(testCase, sum(state.mineral_moles), ...
    sum(geometry.solid_volume_cm3) ./ cfg.geometry.molarVolume_cm3_mol, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, state.mineral_moles(geometry.solid_volume_cm3 == 0), ...
    zeros(nnz(geometry.solid_volume_cm3 == 0), 1), 'AbsTol', 0);
verifyTrue(testCase, isfield(connectivity, 'cell_neighbors'));
verifyEqual(testCase, numel(connectivity.cell_neighbors), 128 * 64);
verifyTrue(testCase, isfield(cfg.transport, 'options'));
verifyTrue(testCase, isfield(cfg.transport.options, 'internal_face_cells'));
verifyTrue(testCase, isfield(cfg.transport.options, 'boundary_face_cells'));
verifyEqual(testCase, cfg.transport.options.time_integration, 'implicit_euler');
verifyEqual(testCase, unique(cfg.transport.options.boundary_concentration_mol_cm3), ...
    1.255e-6, 'RelTol', 1e-12);
verifyGreaterThan(testCase, size(cfg.transport.options.internal_face_cells, 1), 0);
verifyGreaterThan(testCase, numel(cfg.transport.options.boundary_face_cells), 0);
verifyEqual(testCase, size(cfg.transport.options.boundary_concentration_mol_cm3, 2), 1);

partIIOptions = rtm.benchmark.CreateMolinsDriverCaseOptions('partII_strict', ...
    struct('refinementScales', matrix.time_steps_s, ...
    'acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true));
verifyEqual(testCase, partIIOptions.caseOptions.observationTimes_s, ...
    [900; 1800; 2700]);
end

function testAcceptanceGridInterfaceCentroidsLieOnCircle(testCase)
matrix = struct();
matrix.time_steps_s = 0.054;
matrix.grid_resolutions = "128x64";

options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true));

geometry = options.caseOptions.geometryFactory(0.054, ...
    struct('name', "grid_128x64_dt_0p054", 'index', 1));
activeInterface = geometry.interface_area_cm2(:) > 0;
radius = 0.01;
center = [0.05 0.025];
distanceToCenter = hypot( ...
    geometry.interface_centroid_cm(activeInterface, 1) - center(1), ...
    geometry.interface_centroid_cm(activeInterface, 2) - center(2));

verifyGreaterThan(testCase, nnz(activeInterface), 0);
verifyEqual(testCase, max(abs(distanceToCenter - radius)), 0, ...
    'AbsTol', 1e-12);
end

function testAcceptanceGridActiveInterfaceCellsHaveWater(testCase)
matrix = struct();
matrix.time_steps_s = 0.054;
matrix.grid_resolutions = ["128x64"; "512x256"];

for iGrid = 1:numel(matrix.grid_resolutions)
    options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
        struct('acceptanceMatrix', matrix, ...
        'useAcceptanceGrid', true));

    geometry = options.caseOptions.geometryFactory(iGrid, ...
        struct('name', "grid_" + matrix.grid_resolutions(iGrid) + "_dt_0p054", ...
        'index', iGrid));
    activeInterface = geometry.interface_area_cm2(:) > 0;

    verifyGreaterThan(testCase, nnz(activeInterface), 0);
    verifyGreaterThan(testCase, min(geometry.water_volume_cm3(activeInterface)), 0);
    verifyTrue(testCase, all(isfinite(geometry.interface_centroid_cm(activeInterface, :)), 'all'));
end
end

function testAcceptanceGridVerticalFaceVelocitiesPreserveInletFlow(testCase)
matrix = struct();
matrix.time_steps_s = 0.054;
matrix.grid_resolutions = "128x64";

options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true));

cfg = options.caseOptions.configFactory(0.054, ...
    struct('name', "grid_128x64_dt_0p054", 'index', 1));
mesh = cfg.benchmark.mesh;
transport = cfg.transport.options;
verticalFace = transport.internal_face_cells(:, 2) - ...
    transport.internal_face_cells(:, 1) == 1 & ...
    transport.internal_face_velocity_cm_s(:) > 0;
leftCells = transport.internal_face_cells(verticalFace, 1);
ix = mod(leftCells - 1, mesh.nx) + 1;
faceFluxCapacity = accumarray(ix, ...
    transport.internal_face_area_cm2(verticalFace) .* ...
    transport.internal_face_velocity_cm_s(verticalFace), ...
    [mesh.nx, 1], @sum, 0);
activeColumns = faceFluxCapacity > 0;
expectedFluxCapacity = 0.12 .* mesh.domain_height_cm;

verifyGreaterThan(testCase, nnz(activeColumns), 0);
verifyEqual(testCase, faceFluxCapacity(activeColumns), ...
    repmat(expectedFluxCapacity, nnz(activeColumns), 1), ...
    'RelTol', 1e-12, 'AbsTol', 1e-14);
end

function testAcceptanceGridCanEnablePotentialFlowVerticalVelocity(testCase)
matrix = struct();
matrix.time_steps_s = 0.054;
matrix.grid_resolutions = "128x64";

options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true, ...
    'usePotentialFlowVelocity', true));

cfg = options.caseOptions.configFactory(0.054, ...
    struct('name', "grid_128x64_dt_0p054", 'index', 1));
mesh = cfg.benchmark.mesh;
transport = cfg.transport.options;
horizontalFace = transport.internal_face_cells(:, 2) - ...
    transport.internal_face_cells(:, 1) == mesh.nx;
verticalVelocity = transport.internal_face_velocity_cm_s(horizontalFace);

verifyTrue(testCase, any(verticalVelocity > 0));
verifyTrue(testCase, any(verticalVelocity < 0));
end

function testAcceptanceMatrixExpandsGridTimeStepCrossProduct(testCase)
matrix = struct();
matrix.time_steps_s = [54; 5.4; 0.54; 0.054];
matrix.grid_resolutions = ["128x64"; "256x128"; "512x256"];

options = rtm.benchmark.CreateMolinsDriverCaseOptions('partI_strict', ...
    struct('acceptanceMatrix', matrix, ...
    'useAcceptanceGrid', true));

verifyEqual(testCase, numel(options.refinementScales), 12);
verifyEqual(testCase, numel(options.runNames), 12);
verifyTrue(testCase, any(contains(options.runNames, "128x64") & ...
    contains(options.runNames, "dt_54")));
verifyTrue(testCase, any(contains(options.runNames, "512x256") & ...
    contains(options.runNames, "dt_0p054")));

runIndex = find(contains(options.runNames, "512x256") & ...
    contains(options.runNames, "dt_0p054"), 1);
runInfo = struct('name', options.runNames(runIndex), ...
    'index', runIndex, ...
    'refinement_scale', options.refinementScales(runIndex));
cfg = options.caseOptions.configFactory(options.refinementScales(runIndex), runInfo);
geometry = options.caseOptions.geometryFactory(options.refinementScales(runIndex), runInfo);

verifyEqual(testCase, cfg.time.rt.requestedDt_s, 0.054, 'AbsTol', 1e-14);
verifyEqual(testCase, cfg.time.rt.initialDt_s, 0.054, 'AbsTol', 1e-14);
verifyEqual(testCase, cfg.benchmark.mesh.resolution_label, "512x256");
verifyEqual(testCase, cfg.chemistry.reaction_cluster_depth_cells, 6);
verifyEqual(testCase, numel(geometry.water_volume_cm3), 512 * 256);
end

function testUnknownMolinsDriverCaseKindIsRejected(testCase)
verifyError(testCase, ...
    @() rtm.benchmark.CreateMolinsDriverCaseOptions('mystery', struct()), ...
    'RTSPHEM:Benchmark:UnknownMolinsDriverCase');
end

function cleanupOutput(outputDir)
reportPath = fullfile(outputDir, 'benchmark_convergence_report.json');
if exist(reportPath, 'file') == 2
    delete(reportPath);
end
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end
