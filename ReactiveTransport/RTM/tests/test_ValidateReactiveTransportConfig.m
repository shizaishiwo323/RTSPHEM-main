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
verifyEqual(testCase, ...
    validated.chemistry.calciteStoichiometryAbsoluteTolerance_mol, 1e-14);
verifyEqual(testCase, ...
    validated.chemistry.calciteStoichiometryRelativeTolerance, 1e-8);
verifyEqual(testCase, validated.transport.backend, 'cut_cell_fv');
verifyEqual(testCase, validated.time.mode, 'quasi_steady_geometry');
verifyEqual(testCase, validated.time.rt.initialDt_s, 0.01);
verifyEqual(testCase, validated.time.rt.maxDt_s, 0.1);
verifyEqual(testCase, validated.time.rt.advectiveCfl, 0.5);
verifyEqual(testCase, validated.time.rt.diffusiveNumber, 0.5);
verifyEqual(testCase, validated.time.rt.maxReactantFraction, 0.1);
verifyEqual(testCase, validated.time.geometry.boundaryCfl, 0.25);
verifyEqual(testCase, validated.time.geometry.maxMineralFraction, 0.2);
verifyEqual(testCase, validated.time.geometry.maxDt_s, 60);
verifyEqual(testCase, validated.cutCell.minFluidFraction, 0.1);
verifyEqual(testCase, validated.failure.maxRetries, 12);
verifyEqual(testCase, validated.failure.shrinkFactor, 0.5);
verifyEqual(testCase, validated.failure.minDt_s, 1e-8);
verifyEqual(testCase, validated.mass.absoluteTolerance_mol, 1e-14);
verifyEqual(testCase, validated.mass.relativeTolerance, 1e-8);
verifyEqual(testCase, validated.mass.globalRelativeTolerance, 1e-6);
verifyEqual(testCase, validated.geometry.maxDisplacementOverH, 0.25);
verifyEqual(testCase, validated.geometry.solidAbsoluteTolerance_cm3, 1e-14);
verifyEqual(testCase, validated.geometry.solidRelativeTolerance, 1e-8);
verifyEqual(testCase, validated.flow.absoluteTolerance_cm3_s, 1e-12);
verifyEqual(testCase, validated.flow.relativeTolerance, 1e-8);
verifyEqual(testCase, validated.strictMassConservation, true);
end

function testRejectsNegativeMassTolerance(testCase)
cfg = conservativeConfig();
cfg.mass.relativeTolerance = -1;

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidMassTolerance');
end

function testCutCellOptionsCanBeConfigured(testCase)
cfg = conservativeConfig();
cfg.cutCell.minFluidFraction = 0.2;
cfg.cutCell.smallCellMethod = 'disjoint_agglomeration';

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.cutCell.minFluidFraction, 0.2);
verifyEqual(testCase, validated.cutCell.smallCellMethod, ...
    'disjoint_agglomeration');
end

function testRejectsInvalidCutCellOptions(testCase)
cfg = conservativeConfig();
cfg.cutCell.minFluidFraction = -0.1;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidCutCellConfig');

cfg = conservativeConfig();
cfg.cutCell.minFluidFraction = 1.1;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidCutCellConfig');

cfg = conservativeConfig();
cfg.cutCell.smallCellMethod = 'overlapping_agglomeration';
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:UnknownSmallCellMethod');
end

function testTimeStepOptionsCanBeConfigured(testCase)
cfg = conservativeConfig();
cfg.time.rt.initialDt_s = 0.2;
cfg.time.rt.maxDt_s = 1.5;
cfg.time.rt.advectiveCfl = 0.4;
cfg.time.rt.diffusiveNumber = 0.3;
cfg.time.rt.maxReactantFraction = Inf;
cfg.time.geometry.boundaryCfl = 0.2;
cfg.time.geometry.maxMineralFraction = Inf;
cfg.time.geometry.maxDt_s = Inf;

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.time.rt.initialDt_s, 0.2);
verifyEqual(testCase, validated.time.rt.maxDt_s, 1.5);
verifyEqual(testCase, validated.time.rt.advectiveCfl, 0.4);
verifyEqual(testCase, validated.time.rt.diffusiveNumber, 0.3);
verifyEqual(testCase, validated.time.rt.maxReactantFraction, Inf);
verifyEqual(testCase, validated.time.geometry.boundaryCfl, 0.2);
verifyEqual(testCase, validated.time.geometry.maxMineralFraction, Inf);
verifyEqual(testCase, validated.time.geometry.maxDt_s, Inf);
end

function testRejectsInvalidTimeStepOptions(testCase)
cfg = conservativeConfig();
cfg.time.rt.initialDt_s = 0;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidTimeConfig');

cfg = conservativeConfig();
cfg.time.rt.maxDt_s = -1;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidTimeConfig');

cfg = conservativeConfig();
cfg.time.rt.maxReactantFraction = -1;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidTimeConfig');

cfg = conservativeConfig();
cfg.time.geometry.boundaryCfl = -1;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidTimeConfig');
end

function testFailureRetryOptionsCanBeConfigured(testCase)
cfg = conservativeConfig();
cfg.failure.maxRetries = 0;
cfg.failure.shrinkFactor = 0.25;
cfg.failure.minDt_s = 1e-9;

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.failure.maxRetries, 0);
verifyEqual(testCase, validated.failure.shrinkFactor, 0.25);
verifyEqual(testCase, validated.failure.minDt_s, 1e-9);
end

function testRejectsInvalidFailureRetryOptions(testCase)
cfg = conservativeConfig();
cfg.failure.maxRetries = 1.5;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidFailureRetry');

cfg = conservativeConfig();
cfg.failure.shrinkFactor = 1;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidFailureRetry');

cfg = conservativeConfig();
cfg.failure.minDt_s = 0;
verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidFailureRetry');
end

function testRejectsNegativeGeometryTolerance(testCase)
cfg = conservativeConfig();
cfg.geometry.solidAbsoluteTolerance_cm3 = -1;

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidGeometryTolerance');
end

function testGeometryDisplacementCanBeUnlimited(testCase)
cfg = conservativeConfig();
cfg.geometry.maxDisplacementOverH = Inf;

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.geometry.maxDisplacementOverH, Inf);
end

function testRejectsNegativeChargeTolerance(testCase)
cfg = conservativeConfig();
cfg.chemistry.chargeAbsoluteTolerance_eq = -1;

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidChemistryTolerance');
end

function testCalciteStoichiometryTolerancesCanBeConfigured(testCase)
cfg = conservativeConfig();
cfg.chemistry.calciteStoichiometryAbsoluteTolerance_mol = 2e-11;
cfg.chemistry.calciteStoichiometryRelativeTolerance = 3e-6;

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, ...
    validated.chemistry.calciteStoichiometryAbsoluteTolerance_mol, 2e-11);
verifyEqual(testCase, ...
    validated.chemistry.calciteStoichiometryRelativeTolerance, 3e-6);
end

function testRejectsNegativeCalciteStoichiometryTolerance(testCase)
cfg = conservativeConfig();
cfg.chemistry.calciteStoichiometryRelativeTolerance = -1;

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidChemistryTolerance');
end

function testFlowTolerancesCanBeConfigured(testCase)
cfg = conservativeConfig();
cfg.flow.absoluteTolerance_cm3_s = 2e-10;
cfg.flow.relativeTolerance = 3e-7;

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.flow.absoluteTolerance_cm3_s, 2e-10);
verifyEqual(testCase, validated.flow.relativeTolerance, 3e-7);
end

function testRejectsNegativeFlowTolerance(testCase)
cfg = conservativeConfig();
cfg.flow.absoluteTolerance_cm3_s = -1;

verifyError(testCase, @() rtm.config.ValidateReactiveTransportConfig(cfg), ...
    'RTSPHEM:Config:InvalidFlowTolerance');
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

function testSupportsFixedGeometrySteadyRtMode(testCase)
cfg = conservativeConfig();
cfg.chemistry.mode = 'strict_molins';
cfg.phreeqc.engine = 'none';
cfg.phreeqc.databasePolicy = 'not_used';
cfg.time.mode = 'fixed_geometry_steady_rt';

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.time.mode, 'fixed_geometry_steady_rt');
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
