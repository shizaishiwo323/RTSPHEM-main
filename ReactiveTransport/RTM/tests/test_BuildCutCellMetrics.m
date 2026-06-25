function tests = test_BuildCutCellMetrics
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
addpath(fullfile(fileparts(rtmDir), 'HyPHM', 'classes'));
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testAllFluidTriangleUsesFullCellVolume(testCase)
mesh = singleRightTriangleMesh();
levelSet = [-1; -2; -3];

geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet, struct('thickness_cm', 2));

verifyEqual(testCase, geometry.cell_volume_cm3, 1, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.water_volume_cm3, 1, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.solid_volume_cm3, 0, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.fluid_fraction, 1, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.interface_area_cm2, 0, 'AbsTol', 1e-14);
verifyTrue(testCase, geometry.active_fluid_cell);
verifyFalse(testCase, geometry.cut_cell);
end

function testAllSolidTriangleUsesZeroWaterVolume(testCase)
mesh = singleRightTriangleMesh();
levelSet = [1; 2; 3];

geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet);

verifyEqual(testCase, geometry.cell_volume_cm3, 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.water_volume_cm3, 0, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.solid_volume_cm3, 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.fluid_fraction, 0, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.interface_area_cm2, 0, 'AbsTol', 1e-14);
verifyFalse(testCase, geometry.active_fluid_cell);
verifyFalse(testCase, geometry.cut_cell);
end

function testVerticalCutThroughUnitSquareConservesVolumeAndInterfaceArea(testCase)
mesh.vertices_cm = [0 0; 1 0; 1 1; 0 1];
mesh.triangles = [1 2 4; 2 3 4];
levelSet = mesh.vertices_cm(:, 1) - 0.5;

geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet);

verifyEqual(testCase, sum(geometry.cell_volume_cm3), 1, 'AbsTol', 1e-14);
verifyEqual(testCase, sum(geometry.water_volume_cm3), 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, sum(geometry.solid_volume_cm3), 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, sum(geometry.interface_area_cm2), 1, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.water_volume_cm3 + geometry.solid_volume_cm3, ...
    geometry.cell_volume_cm3, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.fluid_fraction, [0.75; 0.25], 'AbsTol', 1e-14);
verifyTrue(testCase, all(geometry.active_fluid_cell));
verifyTrue(testCase, all(geometry.cut_cell));
end

function testInterfaceLengthScaleUsesLocalWaterVolumeAndInterfaceArea(testCase)
mesh.vertices_cm = [0 0; 1 0; 1 1; 0 1];
mesh.triangles = [1 2 4; 2 3 4];
levelSet = mesh.vertices_cm(:, 1) - 0.5;

geometry = rtm.geometry.BuildCutCellMetrics(mesh, levelSet);

expectedScale = min(sqrt(geometry.water_volume_cm3 ./ ...
    max(geometry.interface_area_cm2, eps)), sqrt(geometry.cell_volume_cm3));
verifyTrue(testCase, isfield(geometry, 'interface_length_scale_cm'));
verifyEqual(testCase, geometry.interface_length_scale_cm, expectedScale, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyGreaterThanOrEqual(testCase, geometry.interface_length_scale_cm, ...
    zeros(size(geometry.interface_length_scale_cm)));
end

function testHyphmGridFieldNamesAreAccepted(testCase)
grid.coordV = [0 0; 1 0; 0 1];
grid.V0T = [1 2 3];
grid.areaT = 0.5;
levelSet = [-1; 1; -1];

geometry = rtm.geometry.BuildCutCellMetrics(grid, levelSet);

verifyEqual(testCase, geometry.cell_volume_cm3, 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.water_volume_cm3 + geometry.solid_volume_cm3, ...
    geometry.cell_volume_cm3, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.interface_area_cm2, 0.5, 'AbsTol', 1e-14);
verifyTrue(testCase, geometry.cut_cell);
end

function testWrappedHyphmRectGridFieldNamesAreAccepted(testCase)
grid.rectGrid.coordV = [0 0; 1 0; 0 1];
grid.rectGrid.V0T = [1 2 3];
grid.rectGrid.areaT = 0.5;
levelSet = [-1; 1; -1];

geometry = rtm.geometry.BuildCutCellMetrics(grid, levelSet);

verifyEqual(testCase, geometry.cell_volume_cm3, 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.water_volume_cm3 + geometry.solid_volume_cm3, ...
    geometry.cell_volume_cm3, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.interface_area_cm2, 0.5, 'AbsTol', 1e-14);
verifyTrue(testCase, geometry.cut_cell);
end

function testHyphmGridObjectPropertiesAreAccepted(testCase)
grid = Grid([0 0; 1 0; 0 1], [1 2 3]);
levelSet = [-1; 1; -1];

geometry = rtm.geometry.BuildCutCellMetrics(grid, levelSet);

verifyEqual(testCase, geometry.cell_volume_cm3, 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, geometry.interface_area_cm2, 0.5, 'AbsTol', 1e-14);
verifyTrue(testCase, geometry.cut_cell);
end

function mesh = singleRightTriangleMesh()
mesh.vertices_cm = [0 0; 1 0; 0 1];
mesh.triangles = [1 2 3];
end
