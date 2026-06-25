function tests = test_ApplyReactionResult
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

function testAppliesComponentAndMineralDeltas(testCase)
state = validState();
reactionResult = struct();
reactionResult.component_delta_moles = [-1e-7, 2e-7; 3e-7, 4e-7];
reactionResult.mineral_delta_moles = [-2e-7; -1e-7];
reactionResult.realized_interface_moles = [2e-7; 1e-7];
reactionResult.converged = true;

[updated, ledger] = rtm.chemistry.ApplyReactionResult(state, reactionResult);

verifyEqual(testCase, updated.component_moles, ...
    [9e-7, 2.2e-6; 3.3e-6, 4.4e-6], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.mineral_moles, [9.998e-4; 1.9999e-3], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_delta_moles_total, [2e-7, 6e-7], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.mineral_delta_moles_total, -3e-7, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.realized_interface_moles_total, 3e-7, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, ledger.converged);
end

function testRejectsNegativeComponentResult(testCase)
state = validState();
reactionResult = struct();
reactionResult.component_delta_moles = [-2e-6, 0; 0, 0];
reactionResult.mineral_delta_moles = [0; 0];

verifyError(testCase, @() rtm.chemistry.ApplyReactionResult(state, reactionResult), ...
    'RTSPHEM:State:NegativeComponentMoles');
end

function testRejectsNegativeMineralResult(testCase)
state = validState();
reactionResult = struct();
reactionResult.component_delta_moles = zeros(2, 2);
reactionResult.mineral_delta_moles = [-2e-3; 0];

verifyError(testCase, @() rtm.chemistry.ApplyReactionResult(state, reactionResult), ...
    'RTSPHEM:State:NegativeMineralMoles');
end

function testRejectsDeltaSizeMismatch(testCase)
state = validState();
reactionResult = struct();
reactionResult.component_delta_moles = zeros(1, 2);
reactionResult.mineral_delta_moles = zeros(2, 1);

verifyError(testCase, @() rtm.chemistry.ApplyReactionResult(state, reactionResult), ...
    'RTSPHEM:Chemistry:ReactionResultSizeMismatch');
end

function state = validState()
state = struct();
state.component_names = {'H_reactant', 'Ca_total'};
state.component_moles = [1e-6, 2e-6; 3e-6, 4e-6];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1e-3; 2e-3];
state.temperature_C = [25; 25];
state.pressure_atm = [1; 1];
state.time_s = 0;
end
