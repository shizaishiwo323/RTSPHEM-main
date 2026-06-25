function tests = test_ValidateReactiveTransportConfig
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

function testMinimalConfigReceivesConservativeDefaults(testCase)
cfg = struct();

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.solverArchitecture, 'conservative_v2');
verifyEqual(testCase, validated.chemistry.mode, 'strict_molins');
verifyEqual(testCase, validated.chemistry.chargeAbsoluteTolerance_eq, 1e-8);
verifyEqual(testCase, validated.transport.backend, 'cut_cell_fv');
verifyEqual(testCase, validated.time.mode, 'quasi_steady_geometry');
verifyEqual(testCase, validated.cutCell.minFluidFraction, 0.1);
verifyEqual(testCase, validated.failure.maxRetries, 12);
verifyEqual(testCase, validated.failure.shrinkFactor, 0.5);
verifyEqual(testCase, validated.failure.minDt_s, 1e-8);
verifyEqual(testCase, validated.strictMassConservation, true);
end

function testBenchmarkRejectsFallbackDatabasePolicy(testCase)
cfg = conservativeConfig();
cfg.benchmark.enabled = true;
cfg.phreeqc.databasePolicy = 'allow_fallback';
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.phreeqc.engine = 'iphreeqc_com';

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:BenchmarkRequiresExactLocalDatabase');
end

function testStrictMolinsRejectsPhreeqcEngine(testCase)
cfg = conservativeConfig();
cfg.chemistry.mode = 'strict_molins';
cfg.phreeqc.engine = 'iphreeqc_com';

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:StrictMolinsUsesNoPhreeqc');
end

function testLegacyTransportCannotClaimStrictMassConservation(testCase)
cfg = conservativeConfig();
cfg.transport.backend = 'legacy_hyphm';
cfg.strictMassConservation = true;

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:LegacyTransportNotStrict');
end

function testLegacyArchitectureAllowsLegacyTransportWithoutStrictClaim(testCase)
cfg = conservativeConfig();
cfg.solverArchitecture = 'legacy';
cfg.transport.backend = 'legacy_hyphm';
cfg.strictMassConservation = false;
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.phreeqc.engine = 'iphreeqc_com';
cfg.phreeqc.databasePolicy = 'allow_fallback';

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.solverArchitecture, 'legacy');
verifyEqual(testCase, validated.transport.backend, 'legacy_hyphm');
verifyFalse(testCase, validated.strictMassConservation);
end

function cfg = conservativeConfig()
cfg = struct();
cfg.solverArchitecture = 'conservative_v2';
cfg.chemistry = struct('mode', 'external_tst_phreeqc');
cfg.transport = struct('backend', 'cut_cell_fv');
cfg.time = struct('mode', 'quasi_steady_geometry');
cfg.phreeqc = struct('engine', 'iphreeqc_com', ...
    'databasePolicy', 'exact_local', 'databaseName', 'phreeqc-m.dat');
cfg.strictMassConservation = true;
cfg.benchmark = struct('enabled', false);
end
