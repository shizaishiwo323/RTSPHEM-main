function tests = test_ShouldExportTargetDissolutionSlice
tests = functiontests(localfunctions);
end

function testFirstStepCountsAsFirstTargetSlice(testCase)
[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    1, 0.55, 0.55, 100, 0, false);

verifyTrue(testCase, doExport);
verifyEqual(testCase, nextIndex, 1);
end

function testDoesNotExportUntilProgressCrossesNextSlice(testCase)
[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    5, 0.552, 0.55, 100, 1, false);

verifyFalse(testCase, doExport);
verifyEqual(testCase, nextIndex, 1);
end

function testExportsWhenProgressCrossesNextSlice(testCase)
[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    6, 0.5546, 0.55, 100, 1, false);

verifyTrue(testCase, doExport);
verifyEqual(testCase, nextIndex, 2);
end

function testExportsWhenVeryCloseToNextSliceBoundary(testCase)
initialPorosity = 0.545;
targetPorosity = initialPorosity + (1 - initialPorosity) / 99;
nearlyAtBoundary = targetPorosity - (1 - initialPorosity) * 0.04 / 99;

[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    8, nearlyAtBoundary, initialPorosity, 100, 1, false);

verifyTrue(testCase, doExport);
verifyEqual(testCase, nextIndex, 2);
end

function testFinalForceExportsLastSlice(testCase)
[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    200, 0.76, 0.55, 100, 46, true);

verifyTrue(testCase, doExport);
verifyEqual(testCase, nextIndex, 100);
end

function testPartialDissolutionEndpointMapsToFinalSlice(testCase)
[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    20, 0.78, 0.55, 10, 9, false, 0.78);

verifyTrue(testCase, doExport);
verifyEqual(testCase, nextIndex, 10);
end

function testFullDissolutionReservesLastSliceForGrainVolumeStop(testCase)
[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    200, 1.0, 0.55, 100, 99, false, 1.0);

verifyFalse(testCase, doExport);
verifyEqual(testCase, nextIndex, 99);

[doExport, nextIndex] = ShouldExportTargetDissolutionSlice( ...
    201, 1.0, 0.55, 100, 99, true, 1.0);
verifyTrue(testCase, doExport);
verifyEqual(testCase, nextIndex, 100);
end
