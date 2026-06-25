function tests = test_ApplyMineralVolumeChangeToGeometry
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

function testAcceptedMineralChangeUpdatesWaterAndSolidVolumes(testCase)
geometry = struct();
geometry.water_volume_cm3 = [0.6; 0.2];
geometry.solid_volume_cm3 = [0.4; 0.8];
geometry.cell_volume_cm3 = [1.0; 1.0];
geometry.fluid_fraction = [0.6; 0.2];
geometry.interface_area_cm2 = [2; 3];
geometryInfo = struct();
geometryInfo.cell_solid_volume_change_cm3 = [-0.1; -0.25];

updated = rtm.geometry.ApplyMineralVolumeChange(geometry, geometryInfo);

verifyEqual(testCase, updated.solid_volume_cm3, [0.3; 0.55], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.water_volume_cm3, [0.7; 0.45], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.water_volume_cm3 + updated.solid_volume_cm3, ...
    updated.cell_volume_cm3, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.fluid_fraction, [0.7; 0.45], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.final_surface_area_cm2, 5);
verifyEqual(testCase, updated.final_solid_volume_cm3, 0.85, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testRejectsSolidVolumeOvershootInsteadOfClipping(testCase)
geometry = struct();
geometry.water_volume_cm3 = 0.9;
geometry.solid_volume_cm3 = 0.1;
geometry.cell_volume_cm3 = 1.0;
geometry.fluid_fraction = 0.9;
geometryInfo = struct();
geometryInfo.cell_solid_volume_change_cm3 = -0.2;

verifyError(testCase, ...
    @() rtm.geometry.ApplyMineralVolumeChange(geometry, geometryInfo), ...
    'RTSPHEM:Geometry:NegativeSolidVolume');
end

function testPreservesTinyPositiveSolidVolumeAfterSmallChange(testCase)
geometry = struct();
geometry.water_volume_cm3 = 9e-7;
geometry.solid_volume_cm3 = 1e-7;
geometry.cell_volume_cm3 = 1e-6;
geometry.fluid_fraction = 0.9;
geometryInfo = struct();
geometryInfo.cell_solid_volume_change_cm3 = -9.99995e-8;

updated = rtm.geometry.ApplyMineralVolumeChange(geometry, geometryInfo);

verifyEqual(testCase, updated.solid_volume_cm3, 5e-13, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.water_volume_cm3, geometry.cell_volume_cm3 - 5e-13, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end
