function tests = test_ComputeOutletConcentration
tests = functiontests(localfunctions);
end

function testAveragesOutletTriangles(testCase)
concentration = [1; 2; 3; 4];
outletMask = [false; true; true; false];

value = ComputeOutletConcentration(concentration, outletMask);

verifyEqual(testCase, value, 2.5);
end

function testReturnsNaNWhenNoOutletTriangles(testCase)
value = ComputeOutletConcentration([1; 2], [false; false]);

verifyTrue(testCase, isnan(value));
end
