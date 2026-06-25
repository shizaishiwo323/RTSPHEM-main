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

function testSeparatesTimeAndGridRefinementDimensions(testCase)
runs = [
    benchmarkRunWithDimensions("grid128_dt5p4", 1/128, 5.4, 1.0 + 5.4^2, true)
    benchmarkRunWithDimensions("grid128_dt0p54", 1/128, 0.54, 1.0 + 0.54^2, true)
    benchmarkRunWithDimensions("grid128_dt0p054", 1/128, 0.054, 1.0 + 0.054^2, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'front_position_cm', ...
    'refinementDimension', 'time', ...
    'errorTolerance', Inf, ...
    'minObservedOrder', 1.8));

verifyTrue(testCase, report.accepted);
verifyEqual(testCase, report.refinement_dimension, "time");
verifyEqual(testCase, report.refinement_field, "time_step_s");
verifyEqual(testCase, report.refinement_scales, [5.4; 0.54; 0.054], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, report.grid_spacings_cm, repmat(1/128, 3, 1), ...
    'AbsTol', 1e-18);
verifyEqual(testCase, report.time_steps_s, [5.4; 0.54; 0.054], ...
    'AbsTol', 1e-18);
verifyGreaterThan(testCase, min(report.observed_orders), 1.8);
end

function testUsesGridSpacingForGridRefinement(testCase)
runs = [
    benchmarkRunWithDimensions("grid128_dt0p054", 1/128, 0.054, 1.0 + 4/128^2, true)
    benchmarkRunWithDimensions("grid256_dt0p054", 1/256, 0.054, 1.0 + 4/256^2, true)
    benchmarkRunWithDimensions("grid512_dt0p054", 1/512, 0.054, 1.0 + 4/512^2, true)
    ];

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'solid_volume_cm3', ...
    'refinementDimension', 'grid', ...
    'errorTolerance', Inf, ...
    'minObservedOrder', 1.8));

verifyTrue(testCase, report.accepted);
verifyEqual(testCase, report.refinement_dimension, "grid");
verifyEqual(testCase, report.refinement_field, "grid_spacing_cm");
verifyEqual(testCase, report.refinement_scales, [1/128; 1/256; 1/512], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, report.time_steps_s, repmat(0.054, 3, 1), ...
    'AbsTol', 1e-18);
verifyGreaterThan(testCase, min(report.observed_orders), 1.8);
end

function testReportsTwoDimensionalGridTimeMatrix(testCase)
runs = [
    benchmarkRunWithDimensions("grid128_dt5p4", 1/128, 5.4, 10.0, true)
    benchmarkRunWithDimensions("grid128_dt0p54", 1/128, 0.54, 1.0, true)
    benchmarkRunWithDimensions("grid256_dt5p4", 1/256, 5.4, 8.0, false)
    benchmarkRunWithDimensions("grid256_dt0p54", 1/256, 0.54, 0.8, true)
    ];
runs(1).grid_resolution = "128x64";
runs(2).grid_resolution = "128x64";
runs(3).grid_resolution = "256x128";
runs(4).grid_resolution = "256x128";

report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'mean_effective_rate_mol_cm2_s', ...
    'refinementDimension', 'time', ...
    'errorTolerance', Inf));

verifyTrue(testCase, isfield(report, 'refinement_matrix'));
verifyEqual(testCase, report.refinement_matrix.grid_resolutions, ...
    ["128x64"; "256x128"]);
verifyEqual(testCase, report.refinement_matrix.grid_spacings_cm, ...
    [1/128; 1/256], 'AbsTol', 1e-18);
verifyEqual(testCase, report.refinement_matrix.time_steps_s, [5.4, 0.54], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, report.refinement_matrix.observable_values, ...
    [10.0, 1.0; 8.0, 0.8], 'AbsTol', 1e-18);
verifyEqual(testCase, report.refinement_matrix.accepted_runs, ...
    [true, true; false, true]);
verifyTrue(testCase, isfield(report, 'joint_acceptance_matrix'));
verifyEqual(testCase, report.joint_acceptance_matrix.accepted_runs, ...
    [true, true; false, true]);
verifyTrue(testCase, isfield(report, 'time_convergence_report'));
verifyEqual(testCase, string({report.time_convergence_report.grid_resolution}'), ...
    ["128x64"; "256x128"]);
verifyEqual(testCase, report.time_convergence_report(1).time_steps_s, ...
    [5.4, 0.54], 'AbsTol', 1e-18);
verifyEqual(testCase, report.time_convergence_report(1).observable_values, ...
    [10.0, 1.0], 'AbsTol', 1e-18);
verifyTrue(testCase, isfield(report, 'grid_convergence_report'));
verifyEqual(testCase, [report.grid_convergence_report.time_step_s]', ...
    [5.4; 0.54], 'AbsTol', 1e-18);
verifyEqual(testCase, report.grid_convergence_report(2).grid_resolutions, ...
    ["128x64"; "256x128"]);
verifyEqual(testCase, report.grid_convergence_report(2).observable_values, ...
    [1.0; 0.8], 'AbsTol', 1e-18);
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

function run = benchmarkRunWithDimensions(name, gridSpacingCm, timeStepS, ...
        value, accepted)
run = benchmarkRun(name, timeStepS, value, accepted);
run.grid_spacing_cm = gridSpacingCm;
run.time_step_s = timeStepS;
end
