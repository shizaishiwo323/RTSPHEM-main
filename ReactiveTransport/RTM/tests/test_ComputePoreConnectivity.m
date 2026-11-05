function tests = test_ComputePoreConnectivity
tests = functiontests(localfunctions);
end

function testConnectedAndIsolatedPores(testCase)
grid = struct('numT', 4, 'T0E', [1 2; 2 3; 1 0; 4 0]);
[numComponents, largestFraction] = ComputePoreConnectivity(grid, true(4, 1));

verifyEqual(testCase, numComponents, 2);
verifyEqual(testCase, largestFraction, 0.75);
end

function testNoFluid(testCase)
grid = struct('numT', 2, 'T0E', [1 2; 1 0; 2 0]);
[numComponents, largestFraction] = ComputePoreConnectivity(grid, false(2, 1));

verifyEqual(testCase, numComponents, 0);
verifyEqual(testCase, largestFraction, 0);
end
