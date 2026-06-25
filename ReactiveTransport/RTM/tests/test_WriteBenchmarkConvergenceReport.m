function tests = test_WriteBenchmarkConvergenceReport
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

function testWritesBenchmarkConvergenceReportJson(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupOutput(outputDir));

runs = [
    benchmarkRun("coarse", 0.2, 1.04, true)
    benchmarkRun("medium", 0.1, 1.01, true)
    benchmarkRun("fine", 0.05, 1.0025, true)
    ];
report = rtm.benchmark.BuildConvergenceReport(runs, ...
    struct('observableName', 'porosity', ...
    'errorTolerance', 0.05, ...
    'minObservedOrder', 1.8));

path = rtm.benchmark.WriteConvergenceReport(outputDir, report);

verifyEqual(testCase, path, ...
    string(fullfile(outputDir, 'benchmark_convergence_report.json')));
verifyTrue(testCase, exist(path, 'file') == 2);
decoded = jsondecode(fileread(path));
verifyEqual(testCase, string(decoded.observable_name), "porosity");
verifyTrue(testCase, logical(decoded.accepted));
verifyEmpty(testCase, string(decoded.failure_reasons));
end

function testRejectsMissingOutputDirectory(testCase)
report = struct('observable_name', "porosity", 'accepted', true);
verifyError(testCase, ...
    @() rtm.benchmark.WriteConvergenceReport(fullfile(tempname, 'missing'), report), ...
    'RTSPHEM:Benchmark:MissingOutputDirectory');
end

function run = benchmarkRun(name, refinementScale, value, accepted)
run = struct();
run.name = name;
run.refinement_scale = refinementScale;
run.observable_value = value;
run.accepted = accepted;
end

function cleanupOutput(outputDir)
reportPath = fullfile(outputDir, 'benchmark_convergence_report.json');
if exist(reportPath, 'file') == 2
    delete(reportPath);
end
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end
