function tests = test_LimitTargetDissolutionStepSize
tests = functiontests(localfunctions);
end

function testLimitsStepToNextTargetBoundary(testCase)
[dt, wasLimited, boundary] = LimitTargetDissolutionStepSize( ...
    2.0, 1.0, 0.01, 0.55, 0.545, 1, 100, 1e-5, 60, 1.05);

expectedBoundary = 0.545 + (1 - 0.545) * 3 / 100;
expectedDt = 1.05 * (expectedBoundary - 0.55) / 0.01;

verifyTrue(testCase, wasLimited);
verifyEqual(testCase, boundary, expectedBoundary, 'AbsTol', 1e-12);
verifyEqual(testCase, dt, expectedDt, 'RelTol', 1e-12);
end

function testKeepsSmallerProposedStep(testCase)
[dt, wasLimited] = LimitTargetDissolutionStepSize( ...
    0.1, 1.0, 0.01, 0.55, 0.545, 1, 100, 1e-5, 60, 1.05);

verifyFalse(testCase, wasLimited);
verifyEqual(testCase, dt, 0.1);
end

function testNoLimitWithoutPorositySlope(testCase)
[dt, wasLimited] = LimitTargetDissolutionStepSize( ...
    2.0, 1.0, 0.0, 0.55, 0.545, 1, 100, 1e-5, 60, 1.05);

verifyFalse(testCase, wasLimited);
verifyEqual(testCase, dt, 2.0);
end
