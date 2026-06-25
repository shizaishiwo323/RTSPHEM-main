function tests = test_ExternalTstPhreeqcBackend
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testExternalTstUsesChemistryInputHWithoutTransportingH(testCase)
captured = struct();
state = carbonateState();
geometry = twoCellGeometry();
dtSeconds = 2;
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1e-4;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, ...
    dtSeconds, options);

expectedCandidate = options.h_activity_mol_cm3(:) .* ...
    options.rate_constant_cm_s .* geometry.interface_area_cm2(:) .* dtSeconds;
expectedRealized = min(expectedCandidate, state.mineral_moles(:));
verifyFalse(testCase, any(strcmp(state.component_names, 'H')));
verifyFalse(testCase, any(strcmp(state.component_names, 'H+')));
verifyEqual(testCase, captured.batchState.h_mol_cm3, options.h_mol_cm3(:), ...
    'AbsTol', 1e-18);
verifyEqual(testCase, captured.batchState.prescribed_calcite_dissolved_moles, ...
    expectedRealized, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.realized_interface_moles, expectedRealized, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.candidate_interface_moles, expectedCandidate, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 1), expectedRealized, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 2), expectedRealized, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 3:4), zeros(2, 2), ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.mineral_delta_moles, -expectedRealized, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.charge_balance_residual_eq, [3e-9; -4e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.aux.chemistry_mode, "external_tst_phreeqc");

    function batchResult = mockRunBatch(batchState, batchOptions)
        captured.batchState = batchState;
        captured.batchOptions = batchOptions;
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        waterVolume = batchState.water_volume_cm3(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ waterVolume;
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ waterVolume;
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
        batchResult.pH = [6.5; 6.8];
        batchResult.calciteSI = [-1; -0.5];
        batchResult.chargeBalance = [3e-9; -4e-9];
    end
end

function testRequiresExplicitHydrogenChemistryInput(testCase)
state = carbonateState();
geometry = twoCellGeometry();

verifyError(testCase, ...
    @() rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, struct()), ...
    'RTSPHEM:Chemistry:MissingPhreeqcHydrogenInput');
end

function testExternalTstCarriesAlkalinityComponentDelta(testCase)
state = carbonateState();
state.component_names{end + 1} = 'Alkalinity';
state.component_moles(:, end + 1) = [6e-9; 8e-9];
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1e-4;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

verifyEqual(testCase, result.component_delta_moles(:, 5), [3e-9; 4e-9], ...
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
        batchResult.alkalinity_mol_cm3 = batchState.alkalinity_mol_cm3(:) + ...
            [3e-9; 4e-9] ./ batchState.reaction_water_volume_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function testExternalTstLimitsPrescribedReactionByMineralFraction(testCase)
captured = struct();
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1;
options.maxMineralFraction = 0.2;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, ...
    10, options);

expectedRealized = state.mineral_moles(:) .* options.maxMineralFraction;
verifyEqual(testCase, captured.batchState.prescribed_calcite_dissolved_moles, ...
    expectedRealized, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.realized_interface_moles, expectedRealized, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, all(result.inventory_limited));

    function batchResult = mockRunBatch(batchState, ~)
        captured.batchState = batchState;
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

function testExternalTstScalesComponentDeltaWithReactionWaterVolume(testCase)
captured = struct();
state = carbonateState();
geometry = twoCellGeometry();
geometry.water_volume_cm3 = [1e-6; 2e-6];
state.mineral_moles = [1e-9; 2e-9];
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1;
options.minReactionWaterVolumeCm3 = 1000;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, ...
    10, options);

expectedDissolved = state.mineral_moles(:);
verifyEqual(testCase, captured.batchState.reaction_water_volume_cm3, ...
    [1000; 1000], 'AbsTol', 0);
verifyEqual(testCase, captured.batchState.ca_total_mol_cm3, ...
    state.component_moles(:, 1) ./ captured.batchState.reaction_water_volume_cm3, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, captured.batchState.c_total_mol_cm3, ...
    state.component_moles(:, 2) ./ captured.batchState.reaction_water_volume_cm3, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 1), expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-16);
verifyEqual(testCase, result.component_delta_moles(:, 2), expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-16);
verifyEqual(testCase, result.realized_interface_moles, expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);

    function batchResult = mockRunBatch(batchState, ~)
        captured.batchState = batchState;
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        reactionWaterVolume = batchState.reaction_water_volume_cm3(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ reactionWaterVolume;
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ reactionWaterVolume;
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function testExternalTstDoesNotReactDryInterfaceCells(testCase)
state = carbonateState();
state.component_moles(2, :) = 0;
geometry = twoCellGeometry();
geometry.water_volume_cm3(2) = 0;
options = struct();
options.h_mol_cm3 = [1e-7; 1e-7];
options.h_activity_mol_cm3 = [1e-7; 1e-7];
options.rate_constant_cm_s = 0.1;
options.minReactionWaterVolumeCm3 = 1000;
options.runBatchFunction = @rtm.benchmark.MockPhreeqcBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

verifyGreaterThan(testCase, result.realized_interface_moles(1), 0);
verifyEqual(testCase, result.realized_interface_moles(2), 0, 'AbsTol', 0);
verifyEqual(testCase, result.component_delta_moles(2, :), zeros(1, 4), ...
    'AbsTol', 0);
end

function testBenchmarkMockUsesReactionWaterVolumeScaling(testCase)
state = carbonateState();
geometry = twoCellGeometry();
geometry.water_volume_cm3 = [1e-6; 2e-6];
state.mineral_moles = [1e-9; 2e-9];
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1;
options.minReactionWaterVolumeCm3 = 1000;
options.runBatchFunction = @rtm.benchmark.MockPhreeqcBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, ...
    10, options);

verifyEqual(testCase, result.component_delta_moles(:, 1), state.mineral_moles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 2), state.mineral_moles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 3:4), zeros(2, 2), ...
    'AbsTol', 1e-18);
end

function testExternalTstSuppressesTinyNegativeComponentRoundoff(testCase)
state = carbonateState();
state.component_moles(:, 3:4) = 0;
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1e-4;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);
[updated, ledger] = rtm.chemistry.ApplyReactionResult(state, result);

verifyEqual(testCase, result.component_delta_moles(:, 3:4), zeros(2, 2), ...
    'AbsTol', 0);
verifyEqual(testCase, updated.component_moles(:, 3:4), zeros(2, 2), ...
    'AbsTol', 0);
verifyGreaterThan(testCase, ...
    result.aux.component_delta_roundoff_suppressed_moles, 0);
verifyEqual(testCase, ledger.component_delta_moles_total(3:4), [0 0], ...
    'AbsTol', 0);

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = [-5e-26; 0];
        batchResult.cl_total_mol_cm3 = [0; -5e-26 ./ 2];
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function testExternalTstPropagatesFailedCellsFromBatchResult(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1e-4;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

verifyEqual(testCase, result.failed_cells, [2; 5]);
verifyEqual(testCase, result.error_message, "mock PHREEQC cell failure");
verifyFalse(testCase, result.converged);

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
        batchResult.failedCells = [2; 5];
        batchResult.errorMessage = "mock PHREEQC cell failure";
    end
end

function testExternalTstMarksErrorMessageAsNotConverged(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1e-4;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

verifyEqual(testCase, result.failed_cells, zeros(0, 1));
verifyEqual(testCase, result.error_message, "mock PHREEQC warning");
verifyFalse(testCase, result.converged);

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
        batchResult.errorMessage = "mock PHREEQC warning";
    end
end

function testExternalTstRecordsNonzeroPhreeqcRunStatus(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [2e-7; 4e-8];
options.rate_constant_cm_s = 1e-4;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

verifyFalse(testCase, result.converged);
verifyEqual(testCase, result.aux.phreeqc_run_status, 7);

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
        batchResult.runStatus = 7;
    end
end

function testExternalTstReportsReactionClusterDiagnostics(testCase)
state = carbonateState();
state.component_moles = [state.component_moles; state.component_moles(2, :)];
state.mineral_moles = [state.mineral_moles; 2e-9];
state.temperature_C = [state.temperature_C; 25];
state.pressure_atm = [state.pressure_atm; 1];
geometry = twoCellGeometry();
geometry.water_volume_cm3 = [geometry.water_volume_cm3; 3];
geometry.interface_area_cm2 = [geometry.interface_area_cm2; 4];
geometry.solid_volume_cm3 = [geometry.solid_volume_cm3; 1];
geometry.interface_h_cm = [geometry.interface_h_cm; 1];
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8; 2e-8];
options.h_activity_mol_cm3 = [1e-7; 5e-8; 2e-8];
options.rate_constant_cm_s = 1e-4;
options.reactionClusters = struct( ...
    'source_cells', {1, 3}, ...
    'member_cells', {[1; 2], 3});
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

verifyEqual(testCase, result.cluster_ids, [1; 1; 2]);
verifyEqual(testCase, result.aux.reaction_cluster_count, 2);
verifyEqual(testCase, result.aux.reaction_cluster_max_membership, 1);
verifyEqual(testCase, result.aux.reaction_cluster_overlapping_cell_count, 0);

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

function testExternalTstSolvesPhreeqcOnReactionClusters(testCase)
captured = struct();
state = carbonateState();
state.component_moles = [state.component_moles; 3e-9, 4e-9, 1e-8, 1e-8];
state.mineral_moles = [state.mineral_moles; 2e-6];
state.temperature_C = [state.temperature_C; 25];
state.pressure_atm = [state.pressure_atm; 1];
geometry = twoCellGeometry();
geometry.water_volume_cm3 = [geometry.water_volume_cm3; 3];
geometry.interface_area_cm2 = [geometry.interface_area_cm2; 4];
geometry.solid_volume_cm3 = [geometry.solid_volume_cm3; 1];
geometry.interface_h_cm = [geometry.interface_h_cm; 1];
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8; 2e-8];
options.h_activity_mol_cm3 = [1e-7; 5e-8; 2e-8];
options.rate_constant_cm_s = 1e-4;
options.reactionClusters = struct('source_cells', [1; 3], ...
    'member_cells', [1; 2; 3]);
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

expectedCandidate = options.h_activity_mol_cm3(:) .* ...
    options.rate_constant_cm_s .* geometry.interface_area_cm2(:);
expectedClusterDissolved = expectedCandidate(1) + expectedCandidate(3);
expectedComponentDelta = expectedClusterDissolved .* ...
    geometry.water_volume_cm3(:) ./ sum(geometry.water_volume_cm3(:));
verifyEqual(testCase, numel(captured.batchState.h_mol_cm3), 1);
verifyEqual(testCase, captured.batchState.water_volume_cm3, 6, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, captured.batchState.interface_area_cm2, 6, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ...
    captured.batchState.prescribed_calcite_dissolved_moles, ...
    expectedClusterDissolved, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 1), ...
    expectedComponentDelta, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 2), ...
    expectedComponentDelta, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.realized_interface_moles, ...
    [expectedCandidate(1); 0; expectedCandidate(3)], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.pH, repmat(6.4, 3, 1));
verifyTrue(testCase, result.aux.phreeqc_clustered_reaction);

    function batchResult = mockRunBatch(batchState, ~)
        captured.batchState = batchState;
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.reaction_water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.reaction_water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
        batchResult.pH = 6.4;
        batchResult.calciteSI = -0.7;
        batchResult.chargeBalance = 1e-12;
    end
end

function testExternalTstRejectsOverlappingReactionClusters(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [1e-7; 5e-8];
options.rate_constant_cm_s = 1e-4;
options.reactionClusters = struct( ...
    'source_cells', {1, 2}, ...
    'member_cells', {[1; 2], 2});
options.runBatchFunction = @rtm.benchmark.MockPhreeqcBatch;

verifyError(testCase, ...
    @() rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:OverlappingReactionClusters');
end

function testExternalTstRejectsNegativeGeometryMeasure(testCase)
state = carbonateState();
geometry = twoCellGeometry();
geometry.water_volume_cm3(2) = -2;
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];

verifyError(testCase, ...
    @() rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:NegativeGeometryMeasure');
end

function testExternalTstRejectsNegativeHydrogenActivity(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [1e-7; -5e-8];

verifyError(testCase, ...
    @() rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:NegativePhreeqcChemistryInput');
end

function testExternalTstRejectsNonfinitePhreeqcComponentOutput(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [1e-7; 5e-8];
options.runBatchFunction = @mockRunBatch;

verifyError(testCase, ...
    @() rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:InvalidPhreeqcComponentDelta');

    function batchResult = mockRunBatch(batchState, ~)
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = [NaN; batchState.ca_total_mol_cm3(2)];
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = batchState.prescribed_calcite_dissolved_moles(:);
    end
end

function testExternalTstRejectsCalciteStoichiometryMismatch(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.h_activity_mol_cm3 = [1e-7; 5e-8];
options.rate_constant_cm_s = 1e-4;
options.runBatchFunction = @mockRunBatch;
options.calciteStoichiometryAbsoluteTolerance_mol = 1e-18;

verifyError(testCase, ...
    @() rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:CalciteStoichiometryMismatch');

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            2 .* dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function state = carbonateState()
state = struct();
state.component_names = {'Ca', 'C', 'Na', 'Cl'};
state.component_moles = [
    0, 0, 1e-8, 1e-8
    1e-9, 2e-9, 1e-8, 1e-8
    ];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1e-6; 5e-10];
state.temperature_C = [25; 25];
state.pressure_atm = [1; 1];
state.time_s = 0;
end

function geometry = twoCellGeometry()
geometry = struct();
geometry.water_volume_cm3 = [1; 2];
geometry.interface_area_cm2 = [2; 3];
geometry.solid_volume_cm3 = [1; 1];
geometry.interface_h_cm = [1; 1];
end
