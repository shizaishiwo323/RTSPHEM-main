function tests = test_ResolveMeshPartitions
tests = functiontests(localfunctions);
end

function testUsesTargetElementSize(testCase)
[nxParts, nyParts, mode] = ResolveMeshPartitions(0.15, 0.266, 128, [], [], 0.01);

verifyEqual(testCase, nxParts, 15);
verifyEqual(testCase, nyParts, 27);
verifyEqual(testCase, mode, "target_element_size");
end

function testUsesExplicitPartitions(testCase)
[nxParts, nyParts, mode] = ResolveMeshPartitions(0.15, 0.266, 128, 80, 120, []);

verifyEqual(testCase, nxParts, 80);
verifyEqual(testCase, nyParts, 120);
verifyEqual(testCase, mode, "explicit_xy");
end

function testCompletesOneExplicitAxisWithPhysicalAspect(testCase)
[nxParts, nyParts, mode] = ResolveMeshPartitions(0.15, 0.30, 128, 60, [], []);

verifyEqual(testCase, nxParts, 60);
verifyEqual(testCase, nyParts, 120);
verifyEqual(testCase, mode, "explicit_x");
end
