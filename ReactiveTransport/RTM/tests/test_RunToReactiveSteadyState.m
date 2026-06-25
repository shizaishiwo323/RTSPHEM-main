function tests = test_RunToReactiveSteadyState
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

function testZeroReactionConvergesAfterOneWindow(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0;
cfg.time.rt.initialDt_s = 0.25;
cfg.time.rt.maxDt_s = 0.25;
driver = rtm.driver.ReactiveTransportDriver(cfg, strictState(1e-6, 1e-3), ...
    singleInterfaceGeometry(), struct());

options = struct('window_s', 1.0, 'maxWindows', 3, ...
    'rateAbsoluteTolerance_mol_s', 1e-15, ...
    'rateRelativeTolerance', 1e-6, ...
    'requiredConsecutivePasses', 1);
info = rtm.driver.RunToReactiveSteadyState(driver, options);

verifyTrue(testCase, info.converged);
verifyEqual(testCase, info.reason, "absolute_rate_tolerance");
verifyEqual(testCase, info.windows, 1);
verifyEqual(testCase, info.time_s, 1.0, 'AbsTol', 1e-14);
verifyEqual(testCase, info.steady_rate_mol_s, 0, 'AbsTol', 1e-18);
verifyEqual(testCase, info.state.time_s, 1.0, 'AbsTol', 1e-14);
verifyEqual(testCase, info.window_rates_mol_s, 0, 'AbsTol', 1e-18);
end

function testRelativeRateConvergenceRequiresConsecutiveWindows(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0.1;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 2e4;
cfg.time.rt.initialDt_s = 1.0;
cfg.time.rt.maxDt_s = 1.0;
driver = rtm.driver.ReactiveTransportDriver(cfg, strictState(10, 10), ...
    largeSolidInterfaceGeometry(), struct());

options = struct('window_s', 1.0, 'maxWindows', 5, ...
    'rateAbsoluteTolerance_mol_s', 0, ...
    'rateRelativeTolerance', 1.0, ...
    'requiredConsecutivePasses', 2);
info = rtm.driver.RunToReactiveSteadyState(driver, options);

verifyTrue(testCase, info.converged);
verifyEqual(testCase, info.reason, "relative_rate_tolerance");
verifyEqual(testCase, info.windows, 3);
verifyEqual(testCase, numel(info.window_rates_mol_s), 3);
verifyEqual(testCase, info.window_rates_mol_s(1), 1.0, 'RelTol', 1e-12);
expectedSecondRate = 0.1 .* 9 ./ 2;
expectedThirdRate = 0.1 .* (9 - expectedSecondRate) ./ (2 + expectedSecondRate);
verifyEqual(testCase, info.window_rates_mol_s(2), expectedSecondRate, ...
    'RelTol', 1e-12);
verifyEqual(testCase, info.window_rates_mol_s(3), expectedThirdRate, ...
    'RelTol', 1e-12);
end

function testReportsMaxWindowsWithoutConvergence(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.chemistry.rate_constant_cm_s = 0.1;
cfg.geometry.molarVolume_cm3_mol = 1;
cfg.geometry.maxDisplacementOverH = 2e4;
cfg.time.rt.initialDt_s = 1.0;
cfg.time.rt.maxDt_s = 1.0;
driver = rtm.driver.ReactiveTransportDriver(cfg, strictState(10, 10), ...
    largeSolidInterfaceGeometry(), struct());

options = struct('window_s', 1.0, 'maxWindows', 2, ...
    'rateAbsoluteTolerance_mol_s', 0, ...
    'rateRelativeTolerance', 1e-12, ...
    'requiredConsecutivePasses', 1);
info = rtm.driver.RunToReactiveSteadyState(driver, options);

verifyFalse(testCase, info.converged);
verifyEqual(testCase, info.reason, "max_windows");
verifyEqual(testCase, info.windows, 2);
verifyEqual(testCase, info.time_s, 2.0, 'AbsTol', 1e-14);
verifyEqual(testCase, info.steady_rate_mol_s, 0.1 .* 9 ./ 2, ...
    'RelTol', 1e-12);
verifyEqual(testCase, numel(info.summaries), 2);
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

function geometry = largeSolidInterfaceGeometry()
geometry = singleInterfaceGeometry();
geometry.solid_volume_cm3 = 20;
geometry.cell_volume_cm3 = geometry.water_volume_cm3 + geometry.solid_volume_cm3;
geometry.fluid_fraction = geometry.water_volume_cm3 ./ geometry.cell_volume_cm3;
end
