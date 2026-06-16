function tests = test_LoadRandomGeometryConfig
tests = functiontests(localfunctions);
end

function testLoadedGeometryOverridesDefaults(testCase)
tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() cleanupTmpDir(tmpDir));

configPath = fullfile(tmpDir, 'random_geometry_config.mat');
circleCenters = [0.01 0.01; 0.03 0.02];
circleRadii = [0.002; 0.004];
circleRadius = 0.003;
circleSpacing = 0.0007;
targetAvgSpacing = 0.0012;
lengthXAxis = 0.075;
lengthYAxis = 0.045;
save(configPath, 'circleCenters', 'circleRadii', 'circleRadius', ...
    'circleSpacing', 'targetAvgSpacing', 'lengthXAxis', 'lengthYAxis');

defaults = struct( ...
    'circleRadius', 0.009, ...
    'circleSpacing', 0.002, ...
    'targetAvgSpacing', 0.003, ...
    'targetLengthYAxis', 0.04, ...
    'targetAspectRatio', 1.5, ...
    'useRandomParticleRadii', false, ...
    'randomParticleRadiusMin', 0.009, ...
    'randomParticleRadiusMax', 0.009);

[loadedData, params] = LoadRandomGeometryConfig(configPath, defaults);

verifyEqual(testCase, loadedData.circleCenters, circleCenters);
verifyEqual(testCase, params.circleRadius, circleRadius);
verifyEqual(testCase, params.circleSpacing, circleSpacing);
verifyEqual(testCase, params.targetAvgSpacing, targetAvgSpacing);
verifyEqual(testCase, params.targetLengthYAxis, lengthYAxis);
verifyEqual(testCase, params.targetAspectRatio, lengthXAxis / lengthYAxis);
verifyTrue(testCase, params.useRandomParticleRadii);
verifyEqual(testCase, params.randomParticleRadiusMin, min(circleRadii));
verifyEqual(testCase, params.randomParticleRadiusMax, max(circleRadii));
end

function testGitLfsPointerGetsClearError(testCase)
tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() cleanupTmpDir(tmpDir));

configPath = fullfile(tmpDir, 'random_geometry_config.mat');
fid = fopen(configPath, 'w');
fprintf(fid, 'version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 4720\n');
fclose(fid);

defaults = struct();
verifyError(testCase, @() LoadRandomGeometryConfig(configPath, defaults), ...
    'RTM:RandomGeometryConfig:GitLfsPointer');
end

function cleanupTmpDir(tmpDir)
if exist(tmpDir, 'dir')
    files = dir(tmpDir);
    for i = 1:numel(files)
        if ~files(i).isdir
            delete(fullfile(tmpDir, files(i).name));
        end
    end
    rmdir(tmpDir);
end
end
