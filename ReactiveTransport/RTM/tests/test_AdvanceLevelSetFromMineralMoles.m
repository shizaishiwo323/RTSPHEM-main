function tests = test_AdvanceLevelSetFromMineralMoles
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testUniformCircleDissolutionShrinksSignedDistanceRadius(testCase)
mesh = structuredTriangleMesh(-1, 1, -1, 1, 0.05);
radiusBefore = 0.5;
normalDisplacement = 0.001;
levelSet = radiusBefore - hypot(mesh.vertices_cm(:, 1), mesh.vertices_cm(:, 2));
geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet);
realizedMoles = normalDisplacement .* geometry.interface_area_cm2;
options = struct('molarVolume_cm3_mol', 1, ...
    'maxDisplacementOverH', Inf, 'thickness_cm', 1);

[newLevelSet, moveInfo] = rtm.geometry.AdvanceLevelSetFromMineralMoles( ...
    mesh, levelSet, geometry, realizedMoles, 1, options);

radiusAfter = estimateRightAxisRadius(mesh, newLevelSet);
verifyTrue(testCase, moveInfo.accepted);
verifyEqual(testCase, radiusAfter, radiusBefore - normalDisplacement, ...
    'AbsTol', 5e-3);
verifyLessThan(testCase, sum(moveInfo.new_geometry.solid_volume_cm3), ...
    sum(geometry.solid_volume_cm3));
end

function testRebuiltLevelSetSolidVolumeClosesAgainstRealizedMoles(testCase)
mesh.vertices_cm = [0 0; 1 0; 1 1; 0 1];
mesh.triangles = [1 2 4; 2 3 4];
levelSet = mesh.vertices_cm(:, 1) - 0.5;
geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet);
normalDisplacement = 0.1;
realizedMoles = normalDisplacement .* geometry.interface_area_cm2;
options = struct('molarVolume_cm3_mol', 1, ...
    'maxDisplacementOverH', Inf, 'thickness_cm', 1);

[~, moveInfo] = rtm.geometry.AdvanceLevelSetFromMineralMoles( ...
    mesh, levelSet, geometry, realizedMoles, 1, options);

oldSolidVolume = sum(geometry.solid_volume_cm3);
newSolidVolume = sum(moveInfo.new_geometry.solid_volume_cm3);
expectedDissolvedMoles = sum(realizedMoles);
actualDissolvedMoles = -(newSolidVolume - oldSolidVolume) ./ ...
    options.molarVolume_cm3_mol;
verifyTrue(testCase, moveInfo.accepted);
verifyEqual(testCase, actualDissolvedMoles, expectedDissolvedMoles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-14);
verifyLessThan(testCase, abs(moveInfo.mineral_volume_closure_relative_residual), ...
    1e-12);
end

function testRejectsMineralChangeWithoutInterfaceArea(testCase)
mesh = structuredTriangleMesh(0, 1, 0, 1, 1);
levelSet = -ones(size(mesh.vertices_cm, 1), 1);
geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet);
realizedMoles = ones(size(geometry.water_volume_cm3));

[newLevelSet, moveInfo] = rtm.geometry.AdvanceLevelSetFromMineralMoles( ...
    mesh, levelSet, geometry, realizedMoles, 1);

verifyFalse(testCase, moveInfo.accepted);
verifyEqual(testCase, moveInfo.reject_reason, ...
    "mineral change without interface area");
verifyEqual(testCase, newLevelSet, levelSet, 'AbsTol', 0);
end

function testRejectsDisplacementAboveConfiguredLimit(testCase)
mesh.vertices_cm = [0 0; 1 0; 1 1; 0 1];
mesh.triangles = [1 2 4; 2 3 4];
levelSet = mesh.vertices_cm(:, 1) - 0.5;
geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet);
realizedMoles = 0.2 .* geometry.interface_area_cm2;
options = struct('molarVolume_cm3_mol', 1, ...
    'maxDisplacementOverH', 0.05);

[newLevelSet, moveInfo] = rtm.geometry.AdvanceLevelSetFromMineralMoles( ...
    mesh, levelSet, geometry, realizedMoles, 1, options);

verifyFalse(testCase, moveInfo.accepted);
verifyEqual(testCase, moveInfo.reject_reason, ...
    "geometry displacement exceeds tolerance");
verifyEqual(testCase, newLevelSet, levelSet, 'AbsTol', 0);
end

function mesh = structuredTriangleMesh(xMin, xMax, yMin, yMax, spacing)
x = (xMin:spacing:xMax).';
y = (yMin:spacing:yMax).';
[X, Y] = meshgrid(x, y);
mesh.vertices_cm = [X(:), Y(:)];
nx = numel(x);
ny = numel(y);
triangles = zeros(2 * (nx - 1) * (ny - 1), 3);
k = 0;
for j = 1:(ny - 1)
    for i = 1:(nx - 1)
        v1 = sub2ind([ny, nx], j, i);
        v2 = sub2ind([ny, nx], j, i + 1);
        v3 = sub2ind([ny, nx], j + 1, i + 1);
        v4 = sub2ind([ny, nx], j + 1, i);
        k = k + 1;
        triangles(k, :) = [v1, v2, v4];
        k = k + 1;
        triangles(k, :) = [v2, v3, v4];
    end
end
mesh.triangles = triangles;
end

function radius = estimateRightAxisRadius(mesh, levelSet)
y = mesh.vertices_cm(:, 2);
x = mesh.vertices_cm(:, 1);
axisMask = abs(y) <= 10 * eps & x >= 0;
axisX = x(axisMask);
axisPhi = levelSet(axisMask);
[axisX, order] = sort(axisX);
axisPhi = axisPhi(order);
inside = find(axisPhi >= 0, 1, 'last');
outside = find(axisPhi < 0 & axisX > axisX(inside), 1, 'first');
x1 = axisX(inside);
x2 = axisX(outside);
p1 = axisPhi(inside);
p2 = axisPhi(outside);
radius = x1 - p1 .* (x2 - x1) ./ (p2 - p1);
end
