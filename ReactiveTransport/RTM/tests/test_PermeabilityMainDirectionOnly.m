function tests = test_PermeabilityMainDirectionOnly
tests = functiontests(localfunctions);
end

function testKxxMatchesFullTensor(testCase)
mesh = FoldedCartesianGrid(2, [0, 1, 0, 1], [24, 24]);
grid = Grid(mesh.coordinates, mesh.triangles);
level = 0.18 - hypot(grid.coordV(:, 1) - 0.5, grid.coordV(:, 2) - 0.5);

fullK = computePhysicalPermeabilityTensor(grid.coordV, grid.V0T, level);
mainK = computePhysicalPermeabilityTensor( ...
    grid.coordV, grid.V0T, level, 'MainDirectionOnly');

verifyEqual(testCase, mainK(1, 1), fullK(1, 1), 'RelTol', 1e-12);
verifyTrue(testCase, all(isnan(mainK(2, :))));
end

function testNestedDissectionLargeSystemMatchesReference(testCase)
mesh = FoldedCartesianGrid(2, [0, 1, 0, 1], [192, 192]);
grid = Grid(mesh.coordinates, mesh.triangles);
level = 0.18 - hypot(grid.coordV(:, 1) - 0.5, grid.coordV(:, 2) - 0.5);

mainK = computePhysicalPermeabilityTensor( ...
    grid.coordV, grid.V0T, level, 'MainDirectionOnly');

verifyEqual(testCase, mainK(1, 1), 0.040361019031682203, ...
    'RelTol', 1e-10);
end
