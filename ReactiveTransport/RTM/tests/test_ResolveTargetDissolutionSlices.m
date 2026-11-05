function tests = test_ResolveTargetDissolutionSlices
tests = functiontests(localfunctions);
end

function testDisabledWhenTargetIsEmpty(testCase)
settings = ResolveTargetDissolutionSlices([], 0.35, [], 0.01);

verifyFalse(testCase, settings.enabled);
verifyEqual(testCase, settings.porosityStepTarget, 0.01);
verifyEqual(testCase, settings.maxTotalTimeSteps, []);
end

function testComputesPorosityIncrementFromInitialPorosity(testCase)
settings = ResolveTargetDissolutionSlices(100, 0.35, [], 0.01);

verifyTrue(testCase, settings.enabled);
verifyEqual(testCase, settings.targetSliceCount, 100);
verifyEqual(testCase, settings.porosityStepTarget, 0.65 / 99, 'AbsTol', 1e-12);
verifyGreaterThanOrEqual(testCase, settings.maxTotalTimeSteps, 150);
end

function testKeepsExplicitSafetyStepLimit(testCase)
settings = ResolveTargetDissolutionSlices(100, 0.35, 120, 0.01);

verifyTrue(testCase, settings.enabled);
verifyEqual(testCase, settings.maxTotalTimeSteps, 120);
end

function testUsesConfiguredPartialDissolutionEndpoint(testCase)
settings = ResolveTargetDissolutionSlices(10, 0.55, [], 0.01, 2, 20, 0.78);

verifyEqual(testCase, settings.finalPorosity, 0.78);
verifyEqual(testCase, settings.porosityStepTarget, 0.23 / 9, 'AbsTol', 1e-12);
end

function testRejectsInvalidTarget(testCase)
verifyError(testCase, @() ResolveTargetDissolutionSlices(1, 0.35, [], 0.01), ...
    'RTM:InvalidTargetDissolutionSlices');
end
