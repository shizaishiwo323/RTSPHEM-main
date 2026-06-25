function tests = test_RunMolinsConvergenceSuites
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

function testRunsSelectedMolinsSuitesAndWritesIndex(testCase)
outputRoot = tempname;
mkdir(outputRoot);
cleanup = onCleanup(@() cleanupOutput(outputRoot));

options = struct();
options.outputRoot = outputRoot;
options.kinds = ["partI_strict", "partII_strict"];
options.refinementScales = [0.5, 0.25];
options.totalTime_s = 0.5;
options.rate_constant_cm_s = 0.1;
options.maxDisplacementOverH = 1;
options.errorTolerance = Inf;
options.minObservedOrder = -Inf;

run = rtm.benchmark.RunMolinsConvergenceSuites(options);

verifyEqual(testCase, numel(run.suites), 2);
verifyEqual(testCase, run.kinds, options.kinds(:));
verifyTrue(testCase, exist(run.index_path, 'file') == 2);
for iSuite = 1:numel(run.suites)
    verifyTrue(testCase, exist(run.suites(iSuite).report_path, 'file') == 2);
    verifyEqual(testCase, numel(run.suites(iSuite).runs), 2);
end
decoded = jsondecode(fileread(run.index_path));
verifyEqual(testCase, string(decoded.kinds), options.kinds(:));
verifyEqual(testCase, numel(decoded.suites), 2);
verifyTrue(testCase, isfield(decoded.suites, 'report_path'));
end

function testRejectsMissingOutputRoot(testCase)
verifyError(testCase, ...
    @() rtm.benchmark.RunMolinsConvergenceSuites(struct()), ...
    'RTSPHEM:Benchmark:MissingOutputRoot');
end

function testFullAcceptancePresetRecordsPlanAcceptanceMatrix(testCase)
outputRoot = tempname;
mkdir(outputRoot);
cleanup = onCleanup(@() cleanupOutput(outputRoot));

options = struct();
options.outputRoot = outputRoot;
options.totalTime_s = 0.1;
options.maxDisplacementOverH = 1;
options.errorTolerance = Inf;
options.minObservedOrder = -Inf;
options.phreeqcRunBatchFunction = @mockRunBatch;
options.useAcceptanceGrid = false;

run = run_molins_convergence_suites('full_acceptance', options);

verifyEqual(testCase, run.mode, "full_acceptance");
verifyEqual(testCase, run.kinds, ...
    ["partI_strict"; "partII_strict"; "integration_phreeqc"]);
verifyEqual(testCase, run.acceptance_matrix.time_steps_s, ...
    [54; 5.4; 0.54; 0.054]);
verifyEqual(testCase, run.acceptance_matrix.grid_resolutions, ...
    ["128x64"; "256x128"; "512x256"]);
verifyEqual(testCase, numel(run.suites), 3);
verifyEqual(testCase, run.suites(1).report.observable_name, ...
    "mean_effective_rate_mol_cm2_s");
verifyEqual(testCase, run.suites(1).report.reference_target_value, ...
    4.32e-8, 'RelTol', 1e-12);
verifyEqual(testCase, run.suites(1).report.reference_relative_tolerance, ...
    0.05);
verifyEqual(testCase, run.suites(1).report.reference_target_run_name_pattern, ...
    "dt_0p054");
verifyEqual(testCase, run.suites(2).report.observable_name, ...
    "final_solid_volume_cm3");
verifyEqual(testCase, run.suites(3).report.observable_name, ...
    "max_component_mass_residual_moles");
verifyEqual(testCase, run.suites(3).report.observable_maximum_tolerance, ...
    1e-8, 'RelTol', 1e-12);
for iSuite = 1:numel(run.suites)
    verifyEqual(testCase, run.suites(iSuite).acceptance_matrix.time_steps_s, ...
        [54; 5.4; 0.54; 0.054]);
    verifyEqual(testCase, run.suites(iSuite).acceptance_matrix.grid_resolutions, ...
        ["128x64"; "256x128"; "512x256"]);
    verifyEqual(testCase, run.suites(iSuite).report.refinement_dimension, "time");
    verifyEqual(testCase, run.suites(iSuite).report.refinement_field, ...
        "time_step_s");
    verifyTrue(testCase, all(isfinite(run.suites(iSuite).report.grid_spacings_cm)));
    verifyTrue(testCase, all(isfinite(run.suites(iSuite).report.time_steps_s)));
    verifyEqual(testCase, numel(run.suites(iSuite).runs), 12);
    runNames = string({run.suites(iSuite).runs.name}');
    verifyTrue(testCase, any(contains(runNames, "128x64") & contains(runNames, "dt_54")));
    verifyTrue(testCase, any(contains(runNames, "512x256") & contains(runNames, "dt_0p054")));
    partialPath = fullfile(outputRoot, char(run.kinds(iSuite)), ...
        'benchmark_convergence_partial.json');
    verifyTrue(testCase, exist(partialPath, 'file') == 2);
    partial = jsondecode(fileread(partialPath));
    verifyTrue(testCase, isfield(partial.runs(1).summary, ...
        'reaction_realized_moles'));
    verifyGreaterThanOrEqual(testCase, ...
        partial.runs(1).summary.reaction_realized_moles, 0);
end

decoded = jsondecode(fileread(run.index_path));
verifyEqual(testCase, string(decoded.mode), "full_acceptance");
verifyEqual(testCase, decoded.acceptance_matrix.time_steps_s, ...
    [54; 5.4; 0.54; 0.054]);
verifyEqual(testCase, string(decoded.acceptance_matrix.grid_resolutions), ...
    ["128x64"; "256x128"; "512x256"]);
end

function batchResult = mockRunBatch(batchState, ~)
dissolved = batchState.prescribed_calcite_dissolved_moles(:);
waterVolume = batchState.water_volume_cm3(:);
batchResult = struct();
batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
    dissolved ./ waterVolume;
batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
    dissolved ./ waterVolume;
batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
batchResult.calciteDissolvedMoles = dissolved;
batchResult.pH = 7 .* ones(size(dissolved));
batchResult.calciteSI = zeros(size(dissolved));
end

function cleanupOutput(outputRoot)
suiteDirs = {'partI_strict', 'partII_strict', 'integration_phreeqc'};
for iDir = 1:numel(suiteDirs)
    reportPath = fullfile(outputRoot, suiteDirs{iDir}, ...
        'benchmark_convergence_report.json');
    if exist(reportPath, 'file') == 2
        delete(reportPath);
    end
    partialPath = fullfile(outputRoot, suiteDirs{iDir}, ...
        'benchmark_convergence_partial.json');
    if exist(partialPath, 'file') == 2
        delete(partialPath);
    end
    if exist(fullfile(outputRoot, suiteDirs{iDir}), 'dir') == 7
        rmdir(fullfile(outputRoot, suiteDirs{iDir}));
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
