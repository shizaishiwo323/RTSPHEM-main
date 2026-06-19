function tests = test_ResolveConcentrationLimitedStep
tests = functiontests(localfunctions);
end

function testLegacyModeShrinksOnConfiguredOvershoot(testCase)
[nextStep, limited] = ResolveConcentrationLimitedStep( ...
    0.4, 2.0e-4, 1.0e-4, 1.05, 10.0, 0.5, 1e-5, false);

verifyEqual(testCase, nextStep, 0.2, 'AbsTol', 1e-12);
verifyTrue(testCase, limited);
end

function testTargetSliceModeIgnoresMildOvershoot(testCase)
[nextStep, limited] = ResolveConcentrationLimitedStep( ...
    0.4, 1.164e-4, 1.0e-4, 1.05, 10.0, 0.5, 1e-5, true);

verifyEqual(testCase, nextStep, 0.4, 'AbsTol', 1e-12);
verifyFalse(testCase, limited);
end

function testTargetSliceModeStillShrinksSevereOvershoot(testCase)
[nextStep, limited] = ResolveConcentrationLimitedStep( ...
    0.4, 2.0e-3, 1.0e-4, 1.05, 10.0, 0.5, 1e-5, true);

verifyEqual(testCase, nextStep, 0.2, 'AbsTol', 1e-12);
verifyTrue(testCase, limited);
end
