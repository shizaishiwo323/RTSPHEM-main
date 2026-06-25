function tests = test_PhreeqcKineticsBackend
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

function testKineticsUsesPhreeqcActualCalciteDelta(testCase)
captured = struct();
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 3, options);

expectedDissolved = [2e-8; 1e-6];
expectedDissolved(2) = state.mineral_moles(2);
verifyFalse(testCase, isfield(captured.batchState, ...
    'prescribed_calcite_dissolved_moles'));
verifyEqual(testCase, captured.batchOptions.rateLaw, 'database_calcite');
verifyEqual(testCase, result.realized_interface_moles, expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 1), expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.component_delta_moles(:, 2), expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.mineral_delta_moles, -expectedDissolved, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.charge_balance_residual_eq, [5e-9; -6e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, result.aux.chemistry_mode, "phreeqc_kinetics");

    function batchResult = mockRunBatch(batchState, batchOptions)
        captured.batchState = batchState;
        captured.batchOptions = batchOptions;
        rawDissolved = [2e-8; 1e-6];
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = rawDissolved(:);
        batchResult.pH = [6.9; 7.1];
        batchResult.calciteSI = [-0.8; -0.2];
        batchResult.chargeBalance = [5e-9; -6e-9];
    end
end

function testRequiresExplicitInitialHydrogenForSpeciation(testCase)
state = carbonateState();
geometry = twoCellGeometry();

verifyError(testCase, ...
    @() rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, struct()), ...
    'RTSPHEM:Chemistry:MissingPhreeqcHydrogenInput');
end

function testKineticsPropagatesFailedCellsFromBatchResult(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options);

verifyEqual(testCase, result.failed_cells, [1; 4]);
verifyEqual(testCase, result.error_message, "mock kinetics failure");

    function batchResult = mockRunBatch(batchState, ~)
        rawDissolved = [2e-8; 1e-10];
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = rawDissolved(:);
        batchResult.failedCells = [1; 4];
        batchResult.errorMessage = "mock kinetics failure";
    end
end

function testKineticsSuppressesTinyNegativeComponentRoundoff(testCase)
state = carbonateState();
state.component_moles(:, 3:4) = 0;
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options);
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
        rawDissolved = [2e-8; 1e-10];
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = [-5e-26; 0];
        batchResult.cl_total_mol_cm3 = [0; -5e-26 ./ 2];
        batchResult.calciteDissolvedMoles = rawDissolved(:);
    end
end

function testKineticsReportsReactionClusterDiagnostics(testCase)
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
options.reactionClusters = struct( ...
    'source_cells', {1, 3}, ...
    'member_cells', {[1; 2], 3});
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options);

verifyEqual(testCase, result.cluster_ids, [1; 1; 2]);
verifyEqual(testCase, result.aux.reaction_cluster_count, 2);
verifyEqual(testCase, result.aux.reaction_cluster_max_membership, 1);
verifyEqual(testCase, result.aux.reaction_cluster_overlapping_cell_count, 0);

    function batchResult = mockRunBatch(batchState, ~)
        rawDissolved = [2e-8; 1e-10; 5e-10];
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            rawDissolved(:) ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = rawDissolved(:);
    end
end

function testKineticsRejectsOverlappingReactionClusters(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.reactionClusters = struct( ...
    'source_cells', {1, 2}, ...
    'member_cells', {[1; 2], 2});
options.runBatchFunction = @rtm.benchmark.MockPhreeqcBatch;

verifyError(testCase, ...
    @() rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:OverlappingReactionClusters');
end

function testKineticsRejectsNegativeGeometryMeasure(testCase)
state = carbonateState();
geometry = twoCellGeometry();
geometry.interface_area_cm2(1) = -2;
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];

verifyError(testCase, ...
    @() rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:NegativeGeometryMeasure');
end

function testKineticsRejectsNegativeHydrogenInput(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; -5e-8];

verifyError(testCase, ...
    @() rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:NegativePhreeqcChemistryInput');
end

function testKineticsRejectsNonfinitePhreeqcComponentOutput(testCase)
state = carbonateState();
geometry = twoCellGeometry();
options = struct();
options.h_mol_cm3 = [1e-7; 5e-8];
options.runBatchFunction = @mockRunBatch;

verifyError(testCase, ...
    @() rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:InvalidPhreeqcComponentDelta');

    function batchResult = mockRunBatch(batchState, ~)
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = [NaN; batchState.ca_total_mol_cm3(2)];
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = [1e-10; 1e-10];
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
