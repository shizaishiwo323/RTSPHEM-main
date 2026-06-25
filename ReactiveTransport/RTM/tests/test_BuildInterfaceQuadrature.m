function tests = test_BuildInterfaceQuadrature
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

function testOnePointQuadratureUsesInterfaceCentroidAndArea(testCase)
geometry = struct();
geometry.interface_area_cm2 = [0.2; 0; 0.5];
geometry.interface_centroid_cm = [1 2; NaN NaN; 3 4];
geometry.interface_normal = [0 1; NaN NaN; -1 0];

quadrature = rtm.geometry.BuildInterfaceQuadrature(geometry);

verifyEqual(testCase, quadrature.cell_id, [1; 3]);
verifyEqual(testCase, quadrature.point_cm, [1 2; 3 4], 'AbsTol', 1e-14);
verifyEqual(testCase, quadrature.weight_cm2, [0.2; 0.5], 'AbsTol', 1e-14);
verifyEqual(testCase, quadrature.normal, [0 1; -1 0], 'AbsTol', 1e-14);
verifyEqual(testCase, sum(quadrature.weight_cm2), sum(geometry.interface_area_cm2), ...
    'AbsTol', 1e-14);
end

function testExcludesDryInterfaceCellsWhenWaterVolumeIsAvailable(testCase)
geometry = struct();
geometry.interface_area_cm2 = [0.2; 0.4; 0.5];
geometry.interface_centroid_cm = [1 2; 2 3; 3 4];
geometry.interface_normal = [0 1; 1 0; -1 0];
geometry.water_volume_cm3 = [0; 0.1; 0];

quadrature = rtm.geometry.BuildInterfaceQuadrature(geometry);

verifyEqual(testCase, quadrature.cell_id, 2);
verifyEqual(testCase, quadrature.point_cm, [2 3], 'AbsTol', 1e-14);
verifyEqual(testCase, quadrature.weight_cm2, 0.4, 'AbsTol', 1e-14);
verifyEqual(testCase, quadrature.normal, [1 0], 'AbsTol', 1e-14);
end

function testEmptyQuadratureWhenNoInterfaceArea(testCase)
geometry = struct();
geometry.interface_area_cm2 = [0; 0];
geometry.interface_centroid_cm = [NaN NaN; NaN NaN];
geometry.interface_normal = [NaN NaN; NaN NaN];

quadrature = rtm.geometry.BuildInterfaceQuadrature(geometry);

verifyEmpty(testCase, quadrature.cell_id);
verifyEmpty(testCase, quadrature.weight_cm2);
verifyEqual(testCase, size(quadrature.point_cm), [0 2]);
end

function testRejectsMissingCentroids(testCase)
geometry = struct();
geometry.interface_area_cm2 = 1;

verifyError(testCase, @() rtm.geometry.BuildInterfaceQuadrature(geometry), ...
    'RTSPHEM:Geometry:MissingInterfaceCentroid');
end

function testRejectsNegativeInterfaceArea(testCase)
geometry = struct();
geometry.interface_area_cm2 = [0.2; -0.1];
geometry.interface_centroid_cm = [1 2; 3 4];

verifyError(testCase, @() rtm.geometry.BuildInterfaceQuadrature(geometry), ...
    'RTSPHEM:Geometry:NegativeInterfaceArea');
end
