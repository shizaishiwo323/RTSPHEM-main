function tests = test_CreateRtmRuntimeConfig
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

function testCreatesValidatedRuntimeConfigWithPhreeqcManifest(testCase)
databasePath = tempname;
cleanup = onCleanup(@() deleteIfExists(databasePath));
writeTextFile(databasePath, 'TITLE runtime config fixture');
cfg = struct();
cfg.solverArchitecture = 'conservative_v2';
cfg.chemistry = struct('mode', 'external_tst_phreeqc');
cfg.transport = struct('backend', 'cut_cell_fv');
cfg.time = struct('mode', 'quasi_steady_geometry');
cfg.phreeqc = struct('engine', 'mock', ...
    'databasePolicy', 'exact_local', ...
    'databaseName', 'phreeqc_rates.dat', ...
    'databasePath', databasePath, ...
    'persistSession', true, ...
    'useRunString', true);

runtimeConfig = rtm.config.CreateRtmRuntimeConfig(cfg, ...
    struct('gitCommit', 'abc123', 'gitBranch', 'codex/runtime'));

verifyEqual(testCase, runtimeConfig.config.chemistry.mode, ...
    'external_tst_phreeqc');
verifyEqual(testCase, runtimeConfig.config.transport.backend, 'cut_cell_fv');
verifyEqual(testCase, runtimeConfig.phreeqc_runtime.engineType, "mock");
verifyEqual(testCase, runtimeConfig.phreeqc_runtime.databasePath, ...
    string(databasePath));
verifyEqual(testCase, strlength(runtimeConfig.phreeqc_runtime.databaseSha256), 64);
verifyEqual(testCase, runtimeConfig.runtime_manifest.phreeqc.database_sha256, ...
    runtimeConfig.phreeqc_runtime.databaseSha256);
verifyEqual(testCase, runtimeConfig.runtime_manifest.git.commit, "abc123");
verifyEqual(testCase, runtimeConfig.runtime_manifest.git.branch, "codex/runtime");
verifyEqual(testCase, runtimeConfig.operator_order, "quasi_steady_geometry");
verifyEqual(testCase, runtimeConfig.chemistry_mode, "external_tst_phreeqc");
end

function testStrictMolinsRuntimeConfigDoesNotRequirePhreeqcDatabase(testCase)
runtimeConfig = rtm.config.CreateRtmRuntimeConfig(struct());

verifyEqual(testCase, runtimeConfig.config.chemistry.mode, 'strict_molins');
verifyEqual(testCase, runtimeConfig.phreeqc_runtime.engineType, "none");
verifyEqual(testCase, runtimeConfig.phreeqc_runtime.databasePath, "");
verifyEqual(testCase, runtimeConfig.runtime_manifest.phreeqc.database_path, "");
end

function writeTextFile(pathValue, textValue)
fid = fopen(pathValue, 'w');
if fid == -1
    error('test_CreateRtmRuntimeConfig:WriteFailed', ...
        'Cannot write test file: %s', pathValue);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', textValue);
clear cleanupObj;
end

function deleteIfExists(pathValue)
if exist(pathValue, 'file') == 2
    delete(pathValue);
end
end
