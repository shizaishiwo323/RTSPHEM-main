function tests = test_MolinsBenchmarkEntryPoints
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
end

function testDiagnosticConvergenceRunnerUsesStrictSuites(testCase)
outputRoot = tempname;
mkdir(outputRoot);
cleanup = onCleanup(@() cleanupOutput(outputRoot));

run = run_convergence_molins('diagnostic', quickOptions(outputRoot));

verifyEqual(testCase, run.mode, "diagnostic");
verifyEqual(testCase, run.kinds, ["partI_strict"; "partII_strict"]);
verifyTrue(testCase, exist(run.index_path, 'file') == 2);
decoded = jsondecode(fileread(run.index_path));
verifyEqual(testCase, string(decoded.mode), "diagnostic");
verifyEqual(testCase, string(decoded.kinds), ["partI_strict"; "partII_strict"]);
end

function testPartIStrictEntryPointRunsOnlyPartI(testCase)
outputRoot = tempname;
mkdir(outputRoot);
cleanup = onCleanup(@() cleanupOutput(outputRoot));

run = run_benchmark_molins_partI_strict('diagnostic', quickOptions(outputRoot));

verifyEqual(testCase, run.mode, "diagnostic");
verifyEqual(testCase, run.kinds, "partI_strict");
verifyEqual(testCase, numel(run.suites), 1);
verifyTrue(testCase, contains(run.suites(1).suite_name, "partI"));
end

function testPartIIStrictEntryPointRunsOnlyPartII(testCase)
outputRoot = tempname;
mkdir(outputRoot);
cleanup = onCleanup(@() cleanupOutput(outputRoot));

run = run_benchmark_molins_partII_strict('diagnostic', quickOptions(outputRoot));

verifyEqual(testCase, run.mode, "diagnostic");
verifyEqual(testCase, run.kinds, "partII_strict");
verifyEqual(testCase, numel(run.suites), 1);
verifyTrue(testCase, contains(run.suites(1).suite_name, "partII"));
end

function testIntegrationEntryPointCanUseMockPhreeqcRunner(testCase)
outputRoot = tempname;
mkdir(outputRoot);
cleanup = onCleanup(@() cleanupOutput(outputRoot));

options = quickOptions(outputRoot);
options.phreeqcRunBatchFunction = @rtm.benchmark.MockPhreeqcBatch;

run = run_integration_molins_geometry_phreeqc('diagnostic', options);

verifyEqual(testCase, run.mode, "diagnostic");
verifyEqual(testCase, run.kinds, "integration_phreeqc");
verifyEqual(testCase, numel(run.suites), 1);
verifyTrue(testCase, contains(run.suites(1).suite_name, "phreeqc"));
verifyTrue(testCase, run.accepted);
end

function testConvergenceRunnerRejectsUnknownMode(testCase)
verifyError(testCase, ...
    @() run_convergence_molins('not_a_mode', struct('outputRoot', tempname)), ...
    'RTSPHEM:Benchmark:UnknownMolinsConvergenceMode');
end

function options = quickOptions(outputRoot)
options = struct();
options.outputRoot = outputRoot;
options.refinementScales = [0.5; 0.25];
options.totalTime_s = 0.5;
options.rate_constant_cm_s = 0.1;
options.maxDisplacementOverH = 1;
options.errorTolerance = Inf;
options.minObservedOrder = -Inf;
end

function cleanupOutput(outputRoot)
suiteDirs = {'partI_strict', 'partII_strict', 'integration_phreeqc'};
for iDir = 1:numel(suiteDirs)
    reportPath = fullfile(outputRoot, suiteDirs{iDir}, ...
        'benchmark_convergence_report.json');
    if exist(reportPath, 'file') == 2
        delete(reportPath);
    end
    suiteDir = fullfile(outputRoot, suiteDirs{iDir});
    if exist(suiteDir, 'dir') == 7
        rmdir(suiteDir);
    end
end
indexPath = fullfile(outputRoot, 'molins_convergence_index.json');
if exist(indexPath, 'file') == 2
    delete(indexPath);
end
if exist(outputRoot, 'dir') == 7
    rmdir(outputRoot);
end
end
