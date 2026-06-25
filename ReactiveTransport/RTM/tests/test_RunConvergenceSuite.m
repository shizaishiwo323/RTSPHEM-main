function tests = test_RunConvergenceSuite
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

function testRunsConvergenceSuiteAndWritesReport(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupOutput(outputDir));
calls = [];

options = struct();
options.suiteName = 'molins_partII_time_step';
options.refinementScales = [5.4, 0.54, 0.054];
options.runNames = ["dt5p4", "dt0p54", "dt0p054"];
options.observableName = 'porosity';
options.errorTolerance = 0.2;
options.minObservedOrder = 0.5;
options.outputDir = outputDir;
options.runFunction = @mockRun;
options.observableFunction = @(summary) summary.final_porosity;

suite = rtm.benchmark.RunConvergenceSuite(options);

verifyEqual(testCase, calls, options.refinementScales);
verifyEqual(testCase, suite.suite_name, "molins_partII_time_step");
verifyEqual(testCase, numel(suite.runs), 3);
verifyEqual(testCase, string({suite.runs.name})', options.runNames(:));
verifyEqual(testCase, [suite.runs.refinement_scale]', options.refinementScales(:));
verifyEqual(testCase, [suite.runs.observable_value]', [1.20; 1.02; 1.002], ...
    'AbsTol', 1e-14);
verifyTrue(testCase, suite.report.accepted);
verifyTrue(testCase, exist(suite.report_path, 'file') == 2);
decoded = jsondecode(fileread(suite.report_path));
verifyEqual(testCase, string(decoded.observable_name), "porosity");

    function summary = mockRun(scale, runInfo)
        calls(end + 1) = scale;
        summary = struct();
        summary.final_porosity = 1 + scale / 27;
        summary.accepted = true;
        summary.run_name = runInfo.name;
    end
end

function testMarksFailedRunAsRejectedRecord(testCase)
options = struct();
options.refinementScales = [1, 0.5, 0.25];
options.runFunction = @mockRunWithFailure;
options.observableFunction = @(summary) summary.value;

suite = rtm.benchmark.RunConvergenceSuite(options);

verifyFalse(testCase, suite.runs(2).accepted);
verifyFalse(testCase, suite.report.accepted);
verifyTrue(testCase, any(suite.report.failure_reasons == "benchmark run rejected"));
verifyEqual(testCase, suite.runs(2).failure_message, "intentional failure");

    function summary = mockRunWithFailure(scale, ~)
        if scale == 0.5
            error('test_RunConvergenceSuite:IntentionalFailure', ...
                'intentional failure');
        end
        summary = struct('value', scale, 'accepted', true);
    end
end

function testPreservesSummaryWhenObservableIsInvalid(testCase)
options = struct();
options.refinementScales = [1, 0.5];
options.runFunction = @mockRunWithInvalidObservable;
options.observableFunction = @(summary) summary.value;

suite = rtm.benchmark.RunConvergenceSuite(options);

verifyFalse(testCase, suite.runs(1).accepted);
verifyEqual(testCase, suite.runs(1).summary.diagnostic, "kept");
verifyEqual(testCase, suite.runs(1).failure_message, ...
    "Observable function must return a finite scalar.");
verifyTrue(testCase, suite.runs(2).accepted);

    function summary = mockRunWithInvalidObservable(scale, ~)
        summary = struct();
        summary.value = scale;
        summary.accepted = true;
        if scale == 1
            summary.value = NaN;
            summary.diagnostic = "kept";
        end
    end
end

function testWritesPartialCheckpointAfterEachRun(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupOutput(outputDir));

options = struct();
options.refinementScales = [1, 0.5];
options.outputDir = outputDir;
options.writePartialCheckpoint = true;
options.runFunction = @mockRunWithFailure;
options.observableFunction = @(summary) summary.value;

suite = rtm.benchmark.RunConvergenceSuite(options);

verifyFalse(testCase, suite.runs(2).accepted);
checkpointPath = fullfile(outputDir, 'benchmark_convergence_partial.json');
verifyTrue(testCase, exist(checkpointPath, 'file') == 2);
decoded = jsondecode(fileread(checkpointPath));
verifyEqual(testCase, numel(decoded.runs), 2);
verifyEqual(testCase, string(decoded.runs(1).name), "run1");
verifyTrue(testCase, logical(decoded.runs(1).accepted));
verifyFalse(testCase, logical(decoded.runs(2).accepted));
verifyEqual(testCase, string(decoded.runs(2).failure_message), ...
    "intentional failure");

    function summary = mockRunWithFailure(scale, ~)
        if scale == 0.5
            error('test_RunConvergenceSuite:IntentionalFailure', ...
                'intentional failure');
        end
        summary = struct('value', scale, 'accepted', true);
    end
end

function testRejectsMissingRunFunction(testCase)
options = struct('refinementScales', [1, 0.5]);
verifyError(testCase, @() rtm.benchmark.RunConvergenceSuite(options), ...
    'RTSPHEM:Benchmark:MissingRunFunction');
end

function cleanupOutput(outputDir)
reportPath = fullfile(outputDir, 'benchmark_convergence_report.json');
if exist(reportPath, 'file') == 2
    delete(reportPath);
end
partialPath = fullfile(outputDir, 'benchmark_convergence_partial.json');
if exist(partialPath, 'file') == 2
    delete(partialPath);
end
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end
