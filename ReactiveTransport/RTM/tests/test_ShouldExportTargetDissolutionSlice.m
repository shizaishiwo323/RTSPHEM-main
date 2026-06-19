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
    6, 0.5591, 0.55, 100, 1, false);

verifyTrue(testCase, doExport);
verifyEqual(testCase, nextIndex, 2);
end

function testExportsWhenVeryCloseToNextSliceBoundary(testCase)
initialPorosity = 0.545;
targetPorosity = initialPorosity + (1 - initialPorosity) * 2 / 100;
nearlyAtBoundary = targetPorosity - (1 - initialPorosity) * 0.04 / 100;

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
