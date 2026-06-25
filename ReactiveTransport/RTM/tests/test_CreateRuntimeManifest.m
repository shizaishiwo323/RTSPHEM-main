function tests = test_CreateRuntimeManifest
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

function testRuntimeManifestRecordsCoreProvenance(testCase)
fixturePath = tempname;
cleanup = onCleanup(@() deleteIfExists(fixturePath));
fid = fopen(fixturePath, 'w');
fprintf(fid, 'database fixture\n');
fclose(fid);

cfg = struct();
cfg.solverArchitecture = 'conservative_v2';
cfg.chemistry = struct('mode', 'external_tst_phreeqc', ...
    'chargeAbsoluteTolerance_eq', 1e-8);
cfg.transport = struct('backend', 'cut_cell_fv');
cfg.phreeqc = struct('engine', 'iphreeqc_com', ...
    'databasePath', fixturePath, 'databasePolicy', 'exact_local', ...
    'databaseName', 'phreeqc-m.dat', ...
    'comProgId', 'IPhreeqcCOM.Object', ...
    'persistSession', true, 'useRunString', true);
cfg.mass = struct('relativeTolerance', 1e-8);
cfg.geometry = struct('maxDisplacementOverH', 0.25, ...
    'solidRelativeTolerance', 1e-8);
cfg.time = struct('mode', 'quasi_steady_geometry', ...
    'rt', struct('initialDt_s', 0.01, 'maxDt_s', 0.1));
cfg.failure = struct('maxRetries', 12, 'shrinkFactor', 0.5, ...
    'minDt_s', 1e-8);

manifest = rtm.diagnostics.CreateRuntimeManifest(cfg, ...
    struct('gitCommit', 'abc123', 'gitBranch', 'codex/test'));

verifyEqual(testCase, manifest.solver_architecture, "conservative_v2");
verifyEqual(testCase, manifest.chemistry_mode, "external_tst_phreeqc");
verifyEqual(testCase, manifest.chemistry.chargeAbsoluteTolerance_eq, 1e-8);
verifyEqual(testCase, manifest.transport_backend, "cut_cell_fv");
verifyEqual(testCase, manifest.operator_order, "quasi_steady_geometry");
verifyEqual(testCase, manifest.phreeqc.database_path, string(fixturePath));
verifyEqual(testCase, manifest.phreeqc.database_policy, "exact_local");
verifyEqual(testCase, manifest.phreeqc.engine, "iphreeqc_com");
verifyEqual(testCase, manifest.phreeqc.database_name, "phreeqc-m.dat");
verifyEqual(testCase, manifest.phreeqc.com_progid, "IPhreeqcCOM.Object");
verifyTrue(testCase, manifest.phreeqc.persist_session);
verifyTrue(testCase, manifest.phreeqc.use_run_string);
verifyTrue(testCase, isfield(manifest.phreeqc, 'engine_version'));
verifyEqual(testCase, strlength(manifest.phreeqc.database_sha256), 64);
verifyEqual(testCase, manifest.git.commit, "abc123");
verifyEqual(testCase, manifest.git.branch, "codex/test");
verifyEqual(testCase, manifest.time.mode, "quasi_steady_geometry");
verifyEqual(testCase, manifest.time.rt.initialDt_s, 0.01);
verifyEqual(testCase, manifest.time.rt.maxDt_s, 0.1);
verifyEqual(testCase, manifest.geometry.maxDisplacementOverH, 0.25);
verifyEqual(testCase, manifest.geometry.solidRelativeTolerance, 1e-8);
verifyEqual(testCase, manifest.failure.maxRetries, 12);
verifyEqual(testCase, manifest.failure.shrinkFactor, 0.5);
verifyEqual(testCase, manifest.failure.minDt_s, 1e-8);
verifyTrue(testCase, isfield(manifest, 'matlab_version'));
verifyTrue(testCase, isfield(manifest, 'platform'));
verifyEqual(testCase, manifest.units.length, "cm");
verifyEqual(testCase, manifest.units.component_state, "mol");
end

function testRuntimeManifestHandlesMissingDatabaseWithoutFallback(testCase)
cfg = struct();
cfg.solverArchitecture = 'conservative_v2';
cfg.chemistry = struct('mode', 'strict_molins');
cfg.transport = struct('backend', 'cut_cell_fv');
cfg.phreeqc = struct('engine', 'none', 'databasePath', '', ...
    'databasePolicy', 'not_used');

manifest = rtm.diagnostics.CreateRuntimeManifest(cfg);

verifyEqual(testCase, manifest.phreeqc.database_path, "");
verifyEqual(testCase, manifest.phreeqc.database_sha256, "");
verifyEqual(testCase, manifest.phreeqc.database_size_bytes, 0);
end

function testRuntimeManifestSerializesFunctionHandles(testCase)
cfg = struct();
cfg.solverArchitecture = 'conservative_v2';
cfg.chemistry = struct('mode', 'external_tst_phreeqc', ...
    'options', struct('runBatchFunction', @sin));
cfg.transport = struct('backend', 'cut_cell_fv');
cfg.phreeqc = struct('engine', 'mock', 'databasePath', '', ...
    'databasePolicy', 'not_used');

manifest = rtm.diagnostics.CreateRuntimeManifest(cfg);
encoded = jsonencode(manifest);
decoded = jsondecode(encoded);

verifyEqual(testCase, string(decoded.chemistry.options.runBatchFunction), ...
    "<function_handle:sin>");
end

function deleteIfExists(pathValue)
if exist(pathValue, 'file') == 2
    delete(pathValue);
end
end
