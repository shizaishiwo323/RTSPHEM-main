function tests = test_ResolveAdaptivePorosityStepSettings
tests = functiontests(localfunctions);
end

function testDefaultUpperFactorKeepsTargetSlicesTight(testCase)
settings = ResolveAdaptivePorosityStepSettings(struct());

verifyEqual(testCase, settings.porosityStepUpperFactor, 1.0);
end

function testDefaultGrowthFactorPreservesLegacyAdaptiveMode(testCase)
settings = ResolveAdaptivePorosityStepSettings(struct());

verifyEqual(testCase, settings.adaptiveGrowthFactor, 2.0);
end

function testTargetSliceModeUsesConservativeGrowth(testCase)
settings = ResolveAdaptivePorosityStepSettings(struct('targetDissolutionSlices', 100));

verifyEqual(testCase, settings.adaptiveGrowthFactor, 4.0);
end

function testExplicitUpperFactorIsPreserved(testCase)
settings = ResolveAdaptivePorosityStepSettings(struct('porosityStepUpperFactor', 2.0));

verifyEqual(testCase, settings.porosityStepUpperFactor, 2.0);
end

function testExplicitGrowthFactorIsPreserved(testCase)
settings = ResolveAdaptivePorosityStepSettings(struct( ...
    'targetDissolutionSlices', 100, ...
    'adaptiveGrowthFactor', 1.5));

verifyEqual(testCase, settings.adaptiveGrowthFactor, 1.5);
end
