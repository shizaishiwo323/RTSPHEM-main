function tests = test_ValidateAcceptedStep
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

function testAcceptedWhenAllResidualsWithinTolerance(testCase)
diagnostics = rtm.diagnostics.ValidateAcceptedStep( ...
    baseStepInfo(), toleranceOptions());

verifyTrue(testCase, diagnostics.accepted);
verifyEqual(testCase, diagnostics.reasons, strings(0, 1));
verifyLessThanOrEqual(testCase, diagnostics.component_relative_residual, 1e-8);
verifyLessThanOrEqual(testCase, diagnostics.solid_volume_relative_residual, 1e-8);
end

function testRejectsLargeComponentResidual(testCase)
info = baseStepInfo();
info.mass.component_residual_moles = [3e-8, 0];

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == "component mass residual exceeds tolerance"));
verifyGreaterThan(testCase, diagnostics.component_relative_residual, 1e-8);
end

function testRejectsFullLedgerResidualEvenWhenLegacyMassResidualIsClean(testCase)
info = baseStepInfo();
info.mass.initial_component_moles_total = [1, 2];
info.mass.final_component_moles_total = [1.1, 2];
info.mass.component_residual_moles = [0, 0];
info.transport = struct();
info.reaction = struct();
info.remap = struct();

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "component mass residual exceeds tolerance"));
verifyTrue(testCase, isfield(diagnostics, 'full_mass_ledger'));
verifyEqual(testCase, ...
    diagnostics.full_mass_ledger.component_residual_moles, [0.1, 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testExposesFullLedgerMineralResidualOnAcceptedStepDiagnostics(testCase)
info = baseStepInfo();
info.mass.initial_component_moles_total = [1, 2];
info.mass.final_component_moles_total = [1, 2];
info.mass.component_residual_moles = [0, 0];
info.transport = struct();
info.reaction = struct( ...
    'initial_mineral_moles_total', 1, ...
    'mineral_delta_moles_total', -0.25, ...
    'final_mineral_moles_total', 0.80);
info.remap = struct();

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "mineral inventory residual exceeds tolerance"));
verifyEqual(testCase, diagnostics.mineral_absolute_residual_moles, 0.05, ...
    'AbsTol', 1e-15);
verifyEqual(testCase, diagnostics.mineral_relative_residual, 0.05, ...
    'AbsTol', 1e-15);
end

function testRejectsSolidVolumeMineralMismatch(testCase)
info = baseStepInfo();
info.geometry.actual_solid_volume_change_cm3 = -2e-6;

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == "solid volume residual exceeds tolerance"));
verifyGreaterThan(testCase, diagnostics.solid_volume_relative_residual, 1e-8);
end

function testRejectsChemistryFailure(testCase)
info = baseStepInfo();
info.chemistry.converged = false;
info.chemistry.failed_cells = [3; 4];
info.chemistry.error_message = "PHREEQC failed";

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == "chemistry backend failed"));
verifyEqual(testCase, diagnostics.failed_cells, [3; 4]);
verifyEqual(testCase, diagnostics.error_message, "PHREEQC failed");
end

function testRejectsFailedCellsEvenWhenChemistryReportsConverged(testCase)
info = baseStepInfo();
info.chemistry.converged = true;
info.chemistry.failed_cells = [2; 7];

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "chemistry has failed cells"));
verifyEqual(testCase, diagnostics.failed_cells, [2; 7]);
end

function testRejectsChargeBalanceResidualAboveTolerance(testCase)
info = baseStepInfo();
info.chemistry.charge_balance_residual_eq = [2e-7; -1e-9];

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "chemistry charge balance residual exceeds tolerance"));
verifyEqual(testCase, diagnostics.charge_absolute_residual_eq, 2e-7, ...
    'AbsTol', 1e-18);
end

function testRejectsExcessGeometryDisplacement(testCase)
info = baseStepInfo();
info.geometry.max_displacement_over_h = 0.5;

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == "geometry displacement exceeds tolerance"));
end

function testRejectsGeometryBackendFailureReason(testCase)
info = baseStepInfo();
info.geometry.accepted = false;
info.geometry.reject_reason = "solid volume would become negative";

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "solid volume would become negative"));
end

function testRejectsFlowDiagnosticsFailure(testCase)
info = baseStepInfo();
info.flow = struct();
info.flow.accepted = false;
info.flow.failure_reasons = "inlet/outlet flow imbalance exceeds tolerance";
info.flow.inlet_outlet_relative_residual = 1e-3;
info.flow.max_abs_cell_divergence_cm3_s = 2e-9;
info.flow.global_residual_cm3_s = 0;

diagnostics = rtm.diagnostics.ValidateAcceptedStep(info, toleranceOptions());

verifyFalse(testCase, diagnostics.accepted);
verifyTrue(testCase, any(diagnostics.reasons == ...
    "flow diagnostics failed"));
verifyTrue(testCase, any(diagnostics.reasons == ...
    "inlet/outlet flow imbalance exceeds tolerance"));
verifyEqual(testCase, diagnostics.flow_inlet_outlet_relative_residual, ...
    1e-3, 'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.flow_max_abs_cell_divergence_cm3_s, ...
    2e-9, 'AbsTol', 1e-18);
end

function info = baseStepInfo()
info = struct();
info.mass = struct();
info.mass.initial_component_moles_total = [1, 2];
info.mass.final_component_moles_total = [1 + 1e-12, 2 - 1e-12];
info.mass.component_residual_moles = [1e-12, -1e-12];

info.geometry = struct();
info.geometry.solid_volume_before_cm3 = 1;
info.geometry.solid_volume_after_cm3 = 1 - 1e-6;
info.geometry.mineral_dissolved_moles = 1e-6;
info.geometry.molar_volume_cm3_mol = 1;
info.geometry.actual_solid_volume_change_cm3 = -1e-6;
info.geometry.max_displacement_over_h = 0.1;

info.chemistry = struct();
info.chemistry.converged = true;
info.chemistry.failed_cells = [];
info.chemistry.error_message = "";
end

function options = toleranceOptions()
options = struct();
options.mass_relative_tolerance = 1e-8;
options.mass_absolute_tolerance_mol = 1e-14;
options.solid_relative_tolerance = 1e-8;
options.solid_absolute_tolerance_cm3 = 1e-14;
options.charge_absolute_tolerance_eq = 1e-8;
options.max_displacement_over_h = 0.25;
end
