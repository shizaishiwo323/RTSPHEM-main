function tests = test_WriteRuntimeManifest
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

function testWritesRuntimeManifestJson(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupManifest(outputDir));

manifest = struct();
manifest.solver_architecture = "conservative_v2";
manifest.chemistry_mode = "external_tst_phreeqc";
manifest.transport_backend = "cut_cell_fv";
manifest.phreeqc = struct('engine', "iphreeqc_com", ...
    'database_sha256', "abc123");

path = rtm.diagnostics.WriteRuntimeManifest(outputDir, manifest);

verifyEqual(testCase, path, string(fullfile(outputDir, 'runtime_manifest.json')));
verifyTrue(testCase, exist(path, 'file') == 2);
decoded = jsondecode(fileread(path));
verifyEqual(testCase, string(decoded.solver_architecture), "conservative_v2");
verifyEqual(testCase, string(decoded.chemistry_mode), "external_tst_phreeqc");
verifyEqual(testCase, string(decoded.phreeqc.engine), "iphreeqc_com");
verifyEqual(testCase, string(decoded.phreeqc.database_sha256), "abc123");
end

function testRejectsMissingOutputDirectory(testCase)
missingDir = fullfile(tempname, 'missing');
manifest = struct('solver_architecture', "conservative_v2");

verifyError(testCase, ...
    @() rtm.diagnostics.WriteRuntimeManifest(missingDir, manifest), ...
    'RTSPHEM:Diagnostics:MissingOutputDirectory');
end

function cleanupManifest(outputDir)
manifestPath = fullfile(outputDir, 'runtime_manifest.json');
if exist(manifestPath, 'file') == 2
    delete(manifestPath);
end
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end
