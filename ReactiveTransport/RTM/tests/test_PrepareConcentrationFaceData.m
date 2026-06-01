function tests = test_PrepareConcentrationFaceData
tests = functiontests(localfunctions);
end

function testMasksNonPoreTriangles(testCase)
triangles = [1 2 3; 1 3 4; 2 3 4];
levelSetData = [-1; -1; -1; 1];
concentration = [0.1; 0.2; 0.3];

faceData = PrepareConcentrationFaceData(concentration, triangles, levelSetData);

verifyEqual(testCase, faceData(1), 0.1);
verifyTrue(testCase, isnan(faceData(2)));
verifyTrue(testCase, isnan(faceData(3)));
end
