function tests = test_ComputeFullMassLedger
tests = functiontests(localfunctions);
end

function setupOnce(~)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rtmDir);
end

function testComputesFullCoupledLedgerFields(testCase)
stepInfo = balancedStepInfo();

ledger = rtm.diagnostics.ComputeFullMassLedger(stepInfo);

verifyEqual(testCase, ledger.component_names, ["Ca", "C"]);
verifyEqual(testCase, ledger.initial_component_moles_total, [10, 20], ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.transport_source_delta_moles_total, [1, 0], ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.transport_internal_flux_delta_moles_total, ...
    [-0.5, 0.5], 'AbsTol', 0);
verifyEqual(testCase, ledger.transport_boundary_delta_moles_total, [2, 3], ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.reaction_delta_moles_total, [4, 4], ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.remap_delta_moles_total, [-0.25, 0.25], ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.transport_roundoff_suppressed_moles_total, ...
    5e-16, 'AbsTol', 0);
verifyEqual(testCase, ledger.final_component_moles_total, ...
    [16.25, 27.75], 'AbsTol', 0);
verifyEqual(testCase, ledger.component_residual_moles, [0, 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.initial_mineral_moles_total, 50, 'AbsTol', 0);
verifyEqual(testCase, ledger.reaction_mineral_delta_moles_total, -2, ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.final_mineral_moles_total, 48, 'AbsTol', 0);
verifyEqual(testCase, ledger.mineral_residual_moles, 0, 'AbsTol', 1e-18);
verifyEqual(testCase, ledger.initial_solid_volume_cm3, 100, 'AbsTol', 0);
verifyEqual(testCase, ledger.expected_solid_volume_change_cm3, -2, ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.actual_solid_volume_change_cm3, -2, ...
    'AbsTol', 0);
verifyEqual(testCase, ledger.solid_volume_residual_cm3, 0, 'AbsTol', 0);
end

function testValidateFullCoupledStepRejectsMassAndSolidResiduals(testCase)
stepInfo = balancedStepInfo();
stepInfo.mass.final_component_moles_total = [16.2501, 27.75];
stepInfo.geometry.actual_solid_volume_change_cm3 = -2.1;

diagnostics = rtm.diagnostics.ValidateFullCoupledStep(stepInfo, ...
    struct('mass_absolute_tolerance_mol', 1e-8, ...
    'mass_relative_tolerance', 1e-8, ...
    'solid_absolute_tolerance_cm3', 1e-8, ...
    'solid_relative_tolerance', 1e-8));

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "component mass residual exceeds tolerance"));
verifyTrue(testCase, any(diagnostics.reasons == ...
    "solid volume residual exceeds tolerance"));
verifyGreaterThan(testCase, diagnostics.ledger.max_abs_component_residual_moles, ...
    1e-8);
verifyGreaterThan(testCase, abs(diagnostics.ledger.solid_volume_residual_cm3), ...
    1e-8);
end

function testValidateFullCoupledStepRejectsMineralResidual(testCase)
stepInfo = balancedStepInfo();
stepInfo.reaction.final_mineral_moles_total = 48.25;

diagnostics = rtm.diagnostics.ValidateFullCoupledStep(stepInfo, ...
    struct('mineral_absolute_tolerance_mol', 1e-8, ...
    'mineral_relative_tolerance', 1e-8));

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "mineral inventory residual exceeds tolerance"));
verifyGreaterThan(testCase, diagnostics.mineral_absolute_residual_moles, 1e-8);
verifyEqual(testCase, diagnostics.ledger.mineral_residual_moles, 0.25, ...
    'AbsTol', 1e-18);
end

function stepInfo = balancedStepInfo()
stepInfo = struct();
stepInfo.mass = struct();
stepInfo.mass.component_names = {'Ca', 'C'};
stepInfo.mass.initial_component_moles_total = [10, 20];
stepInfo.mass.final_component_moles_total = [16.25, 27.75];

stepInfo.transport = struct();
stepInfo.transport.source_delta_moles_total = [1, 0];
stepInfo.transport.internal_flux_delta_moles_total = [-0.5, 0.5];
stepInfo.transport.boundary_delta_moles_total = [2, 3];
stepInfo.transport.roundoff_suppressed_moles_total = 5e-16;

stepInfo.reaction = struct();
stepInfo.reaction.component_delta_moles_total = [4, 4];
stepInfo.reaction.initial_mineral_moles_total = 50;
stepInfo.reaction.mineral_delta_moles_total = -2;
stepInfo.reaction.final_mineral_moles_total = 48;

stepInfo.remap = struct();
stepInfo.remap.component_delta_moles_total = [-0.25, 0.25];

stepInfo.geometry = struct();
stepInfo.geometry.solid_volume_before_cm3 = 100;
stepInfo.geometry.solid_volume_after_cm3 = 98;
stepInfo.geometry.expected_solid_volume_change_cm3 = -2;
stepInfo.geometry.actual_solid_volume_change_cm3 = -2;
stepInfo.geometry.max_displacement_over_h = 0.1;

stepInfo.chemistry = struct();
stepInfo.chemistry.converged = true;
stepInfo.chemistry.failed_cells = [];
stepInfo.chemistry.error_message = "";
end
