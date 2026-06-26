function tests = test_ResolvePhreeqcRuntime
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

function testExactLocalRuntimeRecordsDatabaseManifest(testCase)
rtmDir = createTempRtmTree(testCase);
databasePath = fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc_rates.dat');
writeTextFile(databasePath, 'TITLE exact local runtime');
cfg = struct();
cfg.phreeqc = struct('engine', 'mock', ...
    'databasePolicy', 'exact_local', ...
    'databaseName', 'phreeqc_rates.dat', ...
    'comProgId', 'MockIPhreeqc.Object');

runtime = rtm.phreeqc.ResolveRuntime(rtmDir, cfg);

verifyEqual(testCase, runtime.databasePath, string(databasePath));
verifyEqual(testCase, runtime.databaseName, "phreeqc_rates.dat");
verifyEqual(testCase, runtime.databasePolicy, "exact_local");
verifyEqual(testCase, strlength(runtime.databaseSha256), 64);
verifyGreaterThan(testCase, runtime.databaseSizeBytes, 0);
verifyEqual(testCase, runtime.engineType, "mock");
verifyEqual(testCase, runtime.comProgId, "MockIPhreeqc.Object");
verifyTrue(testCase, runtime.isAvailable);
end

function testExactLocalRuntimeRejectsMissingDatabase(testCase)
rtmDir = createTempRtmTree(testCase);
cfg = struct();
cfg.phreeqc = struct('engine', 'mock', ...
    'databasePolicy', 'exact_local', ...
    'databaseName', 'phreeqc_rates.dat');

verifyError(testCase, @() rtm.phreeqc.ResolveRuntime(rtmDir, cfg), ...
    'RTSPHEM:Phreeqc:MissingExactLocalDatabase');
end

function testNotUsedRuntimeDoesNotRequireDatabase(testCase)
runtime = rtm.phreeqc.ResolveRuntime("", struct('phreeqc', struct( ...
    'engine', 'none', 'databasePolicy', 'not_used')));

verifyEqual(testCase, runtime.databasePath, "");
verifyEqual(testCase, runtime.databaseSha256, "");
verifyEqual(testCase, runtime.databaseSizeBytes, 0);
verifyEqual(testCase, runtime.engineType, "none");
verifyTrue(testCase, runtime.isAvailable);
end

function rtmDir = createTempRtmTree(testCase)
rtmDir = tempname;
databaseDir = fullfile(rtmDir, 'phreeqc', 'database');
mkdir(databaseDir);
testCase.addTeardown(@() removeTempTreeFiles(rtmDir));
end

function writeTextFile(pathValue, textValue)
fid = fopen(pathValue, 'w');
if fid == -1
    error('test_ResolvePhreeqcRuntime:WriteFailed', ...
        'Cannot write test file: %s', pathValue);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', textValue);
clear cleanupObj;
end

function removeTempTreeFiles(rtmDir)
databaseDir = fullfile(rtmDir, 'phreeqc', 'database');
deleteIfFile(fullfile(databaseDir, 'phreeqc_rates.dat'));
deleteIfFile(fullfile(databaseDir, 'phreeqc-m.dat'));
deleteIfFile(fullfile(databaseDir, 'phreeqc.dat'));
if exist(databaseDir, 'dir') == 7
    rmdir(databaseDir);
end
phreeqcDir = fullfile(rtmDir, 'phreeqc');
if exist(phreeqcDir, 'dir') == 7
    rmdir(phreeqcDir);
end
if exist(rtmDir, 'dir') == 7
    rmdir(rtmDir);
end
end

function deleteIfFile(pathValue)
if exist(pathValue, 'file') == 2
    delete(pathValue);
end
end
