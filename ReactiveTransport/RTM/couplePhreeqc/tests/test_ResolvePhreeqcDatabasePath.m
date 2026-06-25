function tests = test_ResolvePhreeqcDatabasePath
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
helperDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.helperDir = helperDir;
addpath(helperDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testExactLocalReturnsOnlyBundledPreferredDatabase(testCase)
rtmDir = createTempRtmTree(testCase);
databasePath = fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc-m.dat');
writeTextFile(databasePath, 'TITLE exact local');

resolved = ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc-m.dat', ...
    struct('databasePolicy', 'exact_local'));

verifyEqual(testCase, string(resolved), string(databasePath));
end

function testExactLocalRejectsFallbackDatabase(testCase)
rtmDir = createTempRtmTree(testCase);
writeTextFile(fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc.dat'), ...
    'TITLE fallback should not be used');

verifyError(testCase, ...
    @() ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc-m.dat', ...
    struct('databasePolicy', 'exact_local')), ...
    'RTSPHEM:Phreeqc:MissingExactLocalDatabase');
end

function testAllowFallbackKeepsLegacyFallbackBehavior(testCase)
rtmDir = createTempRtmTree(testCase);
fallbackPath = fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc.dat');
writeTextFile(fallbackPath, 'TITLE fallback allowed');

resolved = ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc-m.dat', ...
    struct('databasePolicy', 'allow_fallback'));

verifyEqual(testCase, string(resolved), string(fallbackPath));
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
    error('test_ResolvePhreeqcDatabasePath:WriteFailed', ...
        'Cannot write test file: %s', pathValue);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', textValue);
clear cleanupObj;
end

function removeTempTreeFiles(rtmDir)
databaseDir = fullfile(rtmDir, 'phreeqc', 'database');
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
