function tests = test_CreateMolinsBenchmarkConfig
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

function testStrictPartIConfigUsesStrictMolinsDefaults(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');

verifyEqual(testCase, cfg.benchmark.name, 'molins_partI_strict');
verifyEqual(testCase, cfg.benchmark.part, 'I');
verifyEqual(testCase, cfg.chemistry.mode, 'strict_molins');
verifyEqual(testCase, cfg.phreeqc.engine, 'none');
verifyEqual(testCase, cfg.transport.backend, 'cut_cell_fv');
verifyTrue(testCase, cfg.strictMassConservation);
verifyEqual(testCase, cfg.geometry.solidRelativeTolerance, 1e-6);
rtm.config.ValidateReactiveTransportConfig(cfg);
end

function testStrictPartIIConfigRequiresPartISteadyInitializer(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partII_strict');

verifyEqual(testCase, cfg.benchmark.name, 'molins_partII_strict');
verifyEqual(testCase, cfg.benchmark.part, 'II');
verifyTrue(testCase, cfg.benchmark.initializeFromPartI);
verifyEqual(testCase, cfg.chemistry.mode, 'strict_molins');
rtm.config.ValidateReactiveTransportConfig(cfg);
end

function testPhreeqcIntegrationConfigIsNotStrictBenchmark(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('integration_phreeqc');

verifyEqual(testCase, cfg.benchmark.name, 'molins_geometry_phreeqc_integration');
verifyFalse(testCase, cfg.benchmark.strictBenchmark);
verifyEqual(testCase, cfg.chemistry.mode, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.phreeqc.engine, 'iphreeqc_com');
verifyEqual(testCase, cfg.phreeqc.databasePolicy, 'exact_local');
rtm.config.ValidateReactiveTransportConfig(cfg);
end

function testUnknownBenchmarkConfigIsRejected(testCase)
verifyError(testCase, @() rtm.config.CreateMolinsBenchmarkConfig('mystery'), ...
    'RTSPHEM:Config:UnknownMolinsBenchmark');
end
