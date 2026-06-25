function tests = test_AdvanceGeometryFromMineralMoles
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

function testComputesSolidVolumeChangeFromRealizedMineralMoles(testCase)
geometry = geometryFixture();
realizedMoles = [1e-6; 0];
options = struct('molarVolume_cm3_mol', 2, 'maxDisplacementOverH', 0.25);

info = rtm.geometry.AdvanceGeometryFromMineralMoles(geometry, realizedMoles, options);

verifyEqual(testCase, info.expected_solid_volume_change_cm3, -2e-6, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, info.cell_solid_volume_change_cm3, [-2e-6; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, info.normal_displacement_cm, [1e-6; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, info.max_displacement_over_h, 0.1, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, info.accepted);
end

function testRejectsDisplacementAboveConfiguredLimit(testCase)
geometry = geometryFixture();
realizedMoles = [1e-4; 0];
options = struct('molarVolume_cm3_mol', 2, 'maxDisplacementOverH', 0.25);

info = rtm.geometry.AdvanceGeometryFromMineralMoles(geometry, realizedMoles, options);

verifyFalse(testCase, info.accepted);
verifyEqual(testCase, info.reject_reason, "geometry displacement exceeds tolerance");
verifyGreaterThan(testCase, info.max_displacement_over_h, 0.25);
end

function testDisplacementCflUsesInterfaceLengthScaleWhenAvailable(testCase)
geometry = geometryFixture();
geometry.interface_h_cm = [1; 1];
geometry.interface_length_scale_cm = [1e-5; 1];
realizedMoles = [1e-6; 0];
options = struct('molarVolume_cm3_mol', 2, 'maxDisplacementOverH', 0.25);

info = rtm.geometry.AdvanceGeometryFromMineralMoles(geometry, realizedMoles, options);

verifyEqual(testCase, info.displacement_over_h(1), 0.1, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, info.accepted);
end

function testRejectsSolidVolumeOvershootBeforeClipping(testCase)
geometry = geometryFixture();
geometry.solid_volume_cm3 = [1e-6; 2e-3];
realizedMoles = [2e-6; 0];
options = struct('molarVolume_cm3_mol', 1, 'maxDisplacementOverH', 1e6);

info = rtm.geometry.AdvanceGeometryFromMineralMoles(geometry, realizedMoles, options);

verifyFalse(testCase, info.accepted);
verifyEqual(testCase, info.reject_reason, "solid volume would become negative");
verifyLessThan(testCase, min(info.cell_solid_volume_after_cm3), 0);
end

function testMissingInterfaceAreaGivesZeroDisplacement(testCase)
geometry = geometryFixture();
geometry.interface_area_cm2 = [0; 0];
realizedMoles = [1e-6; 2e-6];

info = rtm.geometry.AdvanceGeometryFromMineralMoles(geometry, realizedMoles);

verifyEqual(testCase, info.normal_displacement_cm, [0; 0], 'AbsTol', 1e-18);
verifyFalse(testCase, info.accepted);
verifyEqual(testCase, info.reject_reason, "mineral change without interface area");
end

function geometry = geometryFixture()
geometry = struct();
geometry.solid_volume_cm3 = [1e-3; 2e-3];
geometry.interface_area_cm2 = [2; 0];
geometry.interface_h_cm = [1e-5; 1e-5];
end
