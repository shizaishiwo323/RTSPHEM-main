function tests = test_BenchmarkConvergenceReport
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

function testReportsSecondOrderGridConvergence(testCase)
runs = [
    benchmarkRun("grid32", 1/32, 1.0 + 4/32^2, true)
    benchmarkRun("grid64", 1/64, 1.0 + 4/64^2, true)
    benchmarkRun("grid128", 1/128, 1.0 + 4/128^2, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'porosity', ...
    'errorTolerance', 0.005, ...
    'minObservedOrder', 1.8));

verifyTrue(testCase, report.accepted);
verifyEqual(testCase, report.observable_name, "porosity");
verifyEqual(testCase, report.reference_value, 1.0 + 4/128^2, ...
    'RelTol', 1e-12);
verifyEqual(testCase, report.observed_orders, 2, 'RelTol', 1e-12);
verifyLessThan(testCase, report.max_error_to_reference, 0.005);
verifyEqual(testCase, report.failure_reasons, strings(0, 1));
end

function testRejectsWeakTimeStepConvergence(testCase)
runs = [
    benchmarkRun("dt5p4", 5.4, 1.20, true)
    benchmarkRun("dt0p54", 0.54, 1.15, true)
    benchmarkRun("dt0p054", 0.054, 1.10, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'front_position_cm', ...
    'errorTolerance', 0.01, ...
    'minObservedOrder', 0.5));

verifyFalse(testCase, report.accepted);
verifyTrue(testCase, any(report.failure_reasons == "error exceeds tolerance"));
verifyTrue(testCase, any(report.failure_reasons == "observed order below tolerance"));
end

function testRejectsRejectedBenchmarkRuns(testCase)
runs = [
    benchmarkRun("dt5p4", 5.4, 1.20, true)
    benchmarkRun("dt0p54", 0.54, 1.10, false)
    benchmarkRun("dt0p054", 0.054, 1.01, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'k_k0'));

verifyFalse(testCase, report.accepted);
verifyTrue(testCase, any(report.failure_reasons == "benchmark run rejected"));
end

function testRejectsValuesOutsideReferenceTargetTolerance(testCase)
runs = [
    benchmarkRun("grid128", 1/128, 4.80e-8, true)
    benchmarkRun("grid256", 1/256, 4.70e-8, true)
    benchmarkRun("grid512", 1/512, 4.60e-8, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'mean_effective_rate_mol_cm2_s', ...
    'referenceTargetValue', 4.32e-8, ...
    'referenceRelativeTolerance', 0.05));

verifyFalse(testCase, report.accepted);
verifyEqual(testCase, report.reference_target_value, 4.32e-8, ...
    'RelTol', 1e-12);
verifyEqual(testCase, report.reference_relative_tolerance, 0.05);
verifyGreaterThan(testCase, report.max_relative_error_to_target, 0.05);
verifyTrue(testCase, any(report.failure_reasons == ...
    "reference target error exceeds tolerance"));
end

function testReferenceTargetCanApplyOnlyToSelectedRuns(testCase)
runs = [
    benchmarkRun("grid128_dt_54", 12, 2.0e-10, true)
    benchmarkRun("grid128_dt_0p054", 11, 4.12e-8, true)
    benchmarkRun("grid256_dt_54", 10, 2.0e-10, true)
    benchmarkRun("grid256_dt_0p054", 9, 4.44e-8, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'mean_effective_rate_mol_cm2_s', ...
    'referenceTargetValue', 4.32e-8, ...
    'referenceRelativeTolerance', 0.05, ...
    'referenceTargetRunNamePattern', 'dt_0p054'));

verifyTrue(testCase, report.accepted);
verifyEqual(testCase, report.reference_target_applies, ...
    [false; true; false; true]);
verifyLessThan(testCase, report.max_relative_error_to_target, 0.05);
end

function testSelectedReferenceTargetStillRejectsFailingSelectedRuns(testCase)
runs = [
    benchmarkRun("grid128_dt_54", 12, 2.0e-10, true)
    benchmarkRun("grid128_dt_0p054", 11, 5.00e-8, true)
    benchmarkRun("grid256_dt_54", 10, 2.0e-10, true)
    benchmarkRun("grid256_dt_0p054", 9, 4.44e-8, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'mean_effective_rate_mol_cm2_s', ...
    'referenceTargetValue', 4.32e-8, ...
    'referenceRelativeTolerance', 0.05, ...
    'referenceTargetRunNamePattern', 'dt_0p054'));

verifyFalse(testCase, report.accepted);
verifyTrue(testCase, any(report.failure_reasons == ...
    "reference target error exceeds tolerance"));
end

function testRejectsObservableValuesAboveMaximumTolerance(testCase)
runs = [
    benchmarkRun("dt5p4", 5.4, 5e-9, true)
    benchmarkRun("dt0p54", 0.54, 2e-8, true)
    benchmarkRun("dt0p054", 0.054, 4e-9, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'max_component_mass_residual_moles', ...
    'observableMaximumTolerance', 1e-8));

verifyFalse(testCase, report.accepted);
verifyEqual(testCase, report.observable_maximum_tolerance, 1e-8, ...
    'RelTol', 1e-12);
verifyEqual(testCase, report.max_observable_value, 2e-8, ...
    'RelTol', 1e-12);
verifyTrue(testCase, any(report.failure_reasons == ...
    "observable exceeds maximum tolerance"));
end

function testRequiresStrictlyDecreasingRefinementScale(testCase)
runs = [
    benchmarkRun("coarse", 0.1, 1.2, true)
    benchmarkRun("duplicate", 0.1, 1.1, true)
    benchmarkRun("fine", 0.01, 1.0, true)
    ];

verifyError(testCase, ...
    @() rtm.benchmark.BuildConvergenceReport(runs, struct()), ...
    'RTSPHEM:Benchmark:InvalidRefinementScale');
end

function run = benchmarkRun(name, refinementScale, value, accepted)
run = struct();
run.name = name;
run.refinement_scale = refinementScale;
run.observable_value = value;
run.accepted = accepted;
run.summary = struct('accepted_steps', 1, 'rejected_steps', double(~accepted));
end
