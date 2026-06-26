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

function testStrictPartIConfigIsRemoved(testCase)
verifyError(testCase, @() rtm.config.CreateMolinsBenchmarkConfig('partI_strict'), ...
    'RTSPHEM:Config:UnknownMolinsBenchmark');
end

function testStrictPartIIConfigIsRemoved(testCase)
verifyError(testCase, @() rtm.config.CreateMolinsBenchmarkConfig('partII_strict'), ...
    'RTSPHEM:Config:UnknownMolinsBenchmark');
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
