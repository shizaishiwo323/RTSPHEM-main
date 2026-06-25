function tests = test_RtmConservedState
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

function testCreateConservedStateStoresComponentMoles(testCase)
geometry = struct();
geometry.water_volume_cm3 = [2; 0.5; 0];
componentNames = {'Ca_total', 'C_total'};
concentration = [1e-6, 2e-6; 3e-6, 4e-6; 7e-6, 8e-6];
mineralMoles = [1e-4; 2e-4; 3e-4];

state = rtm.state.CreateConservedState( ...
    concentration, geometry, componentNames, ...
    struct('mineralNames', {{'Calcite'}}, 'mineralMoles', mineralMoles));

verifyEqual(testCase, state.component_names, componentNames);
verifyEqual(testCase, state.component_moles, ...
    [2e-6, 4e-6; 1.5e-6, 2e-6; 0, 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, state.mineral_names, {'Calcite'});
verifyEqual(testCase, state.mineral_moles, mineralMoles);
verifyEqual(testCase, state.time_s, 0);
end

function testConvertMolesToConcentrationsUsesWaterVolume(testCase)
geometry = struct();
geometry.water_volume_cm3 = [2; 0.5; 0];
state = struct();
state.component_names = {'Ca_total', 'C_total'};
state.component_moles = [2e-6, 4e-6; 1.5e-6, 2e-6; 3e-9, 4e-9];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1; 1; 1];
state.temperature_C = [25; 25; 25];
state.pressure_atm = [1; 1; 1];
state.time_s = 0;

concentration = rtm.state.ConvertMolesConcentrations(state, geometry, 'moles_to_concentration');

verifyEqual(testCase, concentration, ...
    [1e-6, 2e-6; 3e-6, 4e-6; 0, 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testConvertConcentrationsToMolesRoundTrips(testCase)
geometry = struct();
geometry.water_volume_cm3 = [2; 0.5; 0];
componentNames = {'Ca_total', 'C_total'};
concentration = [1e-6, 2e-6; 3e-6, 4e-6; 0, 0];

state = rtm.state.CreateConservedState(concentration, geometry, componentNames);
roundTrip = rtm.state.ConvertMolesConcentrations(state, geometry, 'moles_to_concentration');

verifyEqual(testCase, roundTrip, concentration, 'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testCreateConservedStateRejectsNegativeWaterVolume(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1; -0.5];
componentNames = {'Ca_total'};
concentration = [1e-6; 2e-6];

verifyError(testCase, ...
    @() rtm.state.CreateConservedState(concentration, geometry, componentNames), ...
    'RTSPHEM:State:NegativeWaterVolume');
end

function testConvertConcentrationsRejectsNegativeWaterVolume(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1; -0.5];
concentration = [1e-6; 2e-6];

verifyError(testCase, ...
    @() rtm.state.ConvertMolesConcentrations( ...
        concentration, geometry, 'concentration_to_moles'), ...
    'RTSPHEM:State:NegativeWaterVolume');
end

function testConvertConcentrationsRejectsNegativeConcentration(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1; 0.5];
concentration = [1e-6; -2e-6];

verifyError(testCase, ...
    @() rtm.state.ConvertMolesConcentrations( ...
        concentration, geometry, 'concentration_to_moles'), ...
    'RTSPHEM:State:NegativeConcentration');
end

function testCreateConservedStateRejectsNegativeConcentrationInDryCell(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1; 0];
componentNames = {'Ca_total'};
concentration = [1e-6; -2e-6];

verifyError(testCase, ...
    @() rtm.state.CreateConservedState(concentration, geometry, componentNames), ...
    'RTSPHEM:State:NegativeConcentration');
end

function testConvertMolesToConcentrationsRejectsNegativeWaterVolume(testCase)
geometry = struct();
geometry.water_volume_cm3 = [1; -0.5];
state = validState();

verifyError(testCase, ...
    @() rtm.state.ConvertMolesConcentrations( ...
        state, geometry, 'moles_to_concentration'), ...
    'RTSPHEM:State:NegativeWaterVolume');
end

function testValidateStateRejectsNegativeMoles(testCase)
state = validState();
state.component_moles(1, 1) = -1e-12;

verifyError(testCase, @() rtm.state.ValidateState(state), ...
    'RTSPHEM:State:NegativeComponentMoles');
end

function testValidateStateRejectsMismatchedComponentNames(testCase)
state = validState();
state.component_names = {'Ca_total'};

verifyError(testCase, @() rtm.state.ValidateState(state), ...
    'RTSPHEM:State:ComponentNameMismatch');
end

function testApplyComponentDeltaUpdatesNamedComponents(testCase)
state = validState();
deltaNames = {'C_total', 'Ca_total'};
deltaMoles = [1e-7, -2e-7; 3e-7, 4e-7];

updated = rtm.state.ApplyComponentDelta(state, deltaNames, deltaMoles);

verifyEqual(testCase, updated.component_moles, ...
    [8e-7, 2.1e-6; 3.4e-6, 4.3e-6], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.time_s, state.time_s);
end

function testApplyComponentDeltaRejectsNegativeResult(testCase)
state = validState();
deltaNames = {'Ca_total'};
deltaMoles = [-2e-6; 0];

verifyError(testCase, @() rtm.state.ApplyComponentDelta(state, deltaNames, deltaMoles), ...
    'RTSPHEM:State:NegativeComponentMoles');
end

function state = validState()
state = struct();
state.component_names = {'Ca_total', 'C_total'};
state.component_moles = [1e-6, 2e-6; 3e-6, 4e-6];
state.mineral_names = {'Calcite'};
state.mineral_moles = [1; 1];
state.temperature_C = [25; 25];
state.pressure_atm = [1; 1];
state.time_s = 0;
end
