function tests = test_WriteStepDiagnosticTables
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

function testWritesMassAndSolidDiagnosticCsvFiles(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupOutput(outputDir));

paths = rtm.diagnostics.WriteStepDiagnosticTables(outputDir, sampleStepResults());

verifyEqual(testCase, paths.mass_balance_components, ...
    string(fullfile(outputDir, 'mass_balance_components.csv')));
verifyEqual(testCase, paths.solid_geometry_balance, ...
    string(fullfile(outputDir, 'solid_geometry_balance.csv')));
verifyEqual(testCase, paths.chemistry_step_status, ...
    string(fullfile(outputDir, 'chemistry_step_status.csv')));
verifyEqual(testCase, paths.step_rejection_log, ...
    string(fullfile(outputDir, 'step_rejection_log.csv')));
verifyEqual(testCase, paths.transport_flux_balance, ...
    string(fullfile(outputDir, 'transport_flux_balance.csv')));
verifyEqual(testCase, paths.reaction_cluster_diagnostics, ...
    string(fullfile(outputDir, 'reaction_cluster_diagnostics.csv')));
verifyTrue(testCase, exist(paths.mass_balance_components, 'file') == 2);
verifyTrue(testCase, exist(paths.solid_geometry_balance, 'file') == 2);
verifyTrue(testCase, exist(paths.chemistry_step_status, 'file') == 2);
verifyTrue(testCase, exist(paths.step_rejection_log, 'file') == 2);
verifyTrue(testCase, exist(paths.transport_flux_balance, 'file') == 2);
verifyTrue(testCase, exist(paths.reaction_cluster_diagnostics, 'file') == 2);

massTable = readtable(paths.mass_balance_components, 'TextType', 'string');
solidTable = readtable(paths.solid_geometry_balance, 'TextType', 'string');
chemistryTable = readtable(paths.chemistry_step_status, 'TextType', 'string');
rejectionTable = readtable(paths.step_rejection_log, 'TextType', 'string');
fluxTable = readtable(paths.transport_flux_balance, 'TextType', 'string');
clusterTable = readtable(paths.reaction_cluster_diagnostics, 'TextType', 'string');

verifyEqual(testCase, height(massTable), 2);
verifyEqual(testCase, massTable.step_index, [1; 2]);
verifyEqual(testCase, logical(massTable.accepted), [true; false]);
verifyEqual(testCase, massTable.component_name, ["H_reactant"; "H_reactant"]);
verifyEqual(testCase, massTable.component_residual_moles, [0; 1e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, massTable.chemistry_delta_moles, [-1e-7; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, massTable.boundary_delta_moles, [5e-8; 0], ...
    'AbsTol', 1e-18);
verifyTrue(testCase, any(string(massTable.Properties.VariableNames) == ...
    "remap_delta_moles"));
verifyEqual(testCase, massTable.remap_delta_moles, [4e-8; -5e-8], ...
    'AbsTol', 1e-18);

verifyEqual(testCase, height(solidTable), 2);
verifyEqual(testCase, logical(solidTable.accepted), [true; false]);
verifyEqual(testCase, solidTable.realized_mineral_moles, [1e-7; 0], ...
    'AbsTol', 1e-18);
verifyTrue(testCase, any(string(solidTable.Properties.VariableNames) == ...
    "initial_mineral_moles"));
verifyEqual(testCase, solidTable.initial_mineral_moles, [1; 1], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, solidTable.reaction_mineral_delta_moles, [-1e-7; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, solidTable.final_mineral_moles, [1 - 1e-7; 1], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, solidTable.mineral_residual_moles, [0; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, solidTable.expected_solid_volume_change_cm3, [-1e-7; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, solidTable.max_displacement_over_h, [0.1; 0.5], ...
    'AbsTol', 1e-18);

verifyEqual(testCase, height(chemistryTable), 2);
verifyEqual(testCase, logical(chemistryTable.converged), [true; false]);
verifyEqual(testCase, chemistryTable.failed_cell_count, [0; 2]);
verifyTrue(testCase, ismissing(chemistryTable.error_message(1)) || ...
    chemistryTable.error_message(1) == "");
verifyEqual(testCase, chemistryTable.error_message(2), "PHREEQC failed");
verifyEqual(testCase, chemistryTable.charge_balance_max_abs_eq, [0; 2e-6], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, chemistryTable.transaction_status, ...
    ["committed"; "rolled_back"]);
verifyTrue(testCase, any(string(chemistryTable.Properties.VariableNames) == ...
    "phreeqc_run_status"));
verifyEqual(testCase, chemistryTable.phreeqc_run_status, [0; 9]);

verifyEqual(testCase, height(rejectionTable), 1);
verifyEqual(testCase, rejectionTable.step_index, 2);
verifyEqual(testCase, rejectionTable.reason, "chemistry failed; mass residual");
verifyEqual(testCase, rejectionTable.dt_s, 0.1, 'AbsTol', 1e-18);
verifyTrue(testCase, any(string(rejectionTable.Properties.VariableNames) == ...
    "component_absolute_residual_moles"));
verifyEqual(testCase, rejectionTable.component_absolute_residual_moles, ...
    1e-9, 'AbsTol', 1e-18);
verifyEqual(testCase, rejectionTable.solid_volume_absolute_residual_cm3, ...
    2e-10, 'AbsTol', 1e-18);
verifyEqual(testCase, rejectionTable.mineral_absolute_residual_moles, ...
    3e-10, 'AbsTol', 1e-18);
verifyEqual(testCase, rejectionTable.charge_absolute_residual_eq, ...
    2e-6, 'AbsTol', 1e-18);

verifyEqual(testCase, height(fluxTable), 2);
verifyEqual(testCase, fluxTable.component_name, ["H_reactant"; "H_reactant"]);
verifyEqual(testCase, fluxTable.source_delta_moles, [2e-7; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, fluxTable.internal_flux_delta_moles, [0; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, fluxTable.boundary_delta_moles, [5e-8; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, fluxTable.boundary_advective_delta_moles, [3e-8; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, fluxTable.boundary_diffusive_delta_moles, [2e-8; 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, fluxTable.roundoff_suppressed_moles, [3e-16; 0], ...
    'AbsTol', 1e-24);
verifyEqual(testCase, fluxTable.inlet_flow_cm3_s, [1; 1], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, fluxTable.outlet_flow_cm3_s, [0.99; 0.99], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, fluxTable.inlet_outlet_relative_residual, [0.01; 0.01], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);

verifyEqual(testCase, height(clusterTable), 2);
verifyEqual(testCase, clusterTable.cluster_count, [1; 0]);
verifyEqual(testCase, clusterTable.realized_interface_moles_total, ...
    [1e-7; 0], 'AbsTol', 1e-18);
verifyEqual(testCase, clusterTable.candidate_interface_moles_total, ...
    [2e-7; 0], 'AbsTol', 1e-18);
end

function testRejectsMissingOutputDirectory(testCase)
missingDir = fullfile(tempname, 'missing');
verifyError(testCase, ...
    @() rtm.diagnostics.WriteStepDiagnosticTables(missingDir, sampleStepResults()), ...
    'RTSPHEM:Diagnostics:MissingOutputDirectory');
end

function testMissingFullMineralLedgerWritesNaNInsteadOfZero(testCase)
outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() cleanupOutput(outputDir));

stepResult = baseStepResult();
stepResult.diagnostics = rmfield(stepResult.diagnostics, 'full_mass_ledger');

paths = rtm.diagnostics.WriteStepDiagnosticTables(outputDir, stepResult);
solidTable = readtable(paths.solid_geometry_balance, 'TextType', 'string');

verifyTrue(testCase, isnan(solidTable.initial_mineral_moles));
verifyTrue(testCase, isnan(solidTable.reaction_mineral_delta_moles));
verifyTrue(testCase, isnan(solidTable.final_mineral_moles));
verifyTrue(testCase, isnan(solidTable.mineral_residual_moles));
end

function stepResults = sampleStepResults()
stepResults = repmat(baseStepResult(), 1, 2);

stepResults(2).diagnostics.accepted = false;
stepResults(2).transport_ledger.final_moles_total = 9.5e-7;
stepResults(2).transport_ledger.source_delta_moles_total = 0;
stepResults(2).transport_ledger.boundary_delta_moles_total = 0;
stepResults(2).transport_ledger.boundary_advective_delta_moles_total = 0;
stepResults(2).transport_ledger.boundary_diffusive_delta_moles_total = 0;
stepResults(2).transport_ledger.roundoff_suppressed_moles_total = 0;
stepResults(2).chemistry_ledger.component_delta_moles_total = 0;
stepResults(2).diagnostics.component_absolute_residual_moles = 1e-9;
stepResults(2).diagnostics.solid_volume_absolute_residual_cm3 = 2e-10;
stepResults(2).diagnostics.mineral_absolute_residual_moles = 3e-10;
stepResults(2).diagnostics.full_mass_ledger.final_component_moles_total = 9.5e-7;
stepResults(2).diagnostics.full_mass_ledger.transport_source_delta_moles_total = 0;
stepResults(2).diagnostics.full_mass_ledger.transport_boundary_delta_moles_total = 0;
stepResults(2).diagnostics.full_mass_ledger.reaction_delta_moles_total = 0;
stepResults(2).diagnostics.full_mass_ledger.remap_delta_moles_total = -5e-8;
stepResults(2).diagnostics.full_mass_ledger.component_residual_moles = 1e-9;
stepResults(2).diagnostics.full_mass_ledger.reaction_mineral_delta_moles_total = 0;
stepResults(2).diagnostics.full_mass_ledger.final_mineral_moles_total = 1;
stepResults(2).diagnostics.full_mass_ledger.mineral_residual_moles = 0;
stepResults(2).diagnostics.charge_absolute_residual_eq = 2e-6;
stepResults(2).geometry_info.realized_mineral_moles = 0;
stepResults(2).geometry_info.expected_solid_volume_change_cm3 = 0;
stepResults(2).geometry_info.solid_volume_after_cm3 = 1;
stepResults(2).geometry_info.max_displacement_over_h = 0.5;
stepResults(2).diagnostics.reasons = ["chemistry failed", "mass residual"];
stepResults(2).reaction_result.converged = false;
stepResults(2).reaction_result.failed_cells = [3; 5];
stepResults(2).reaction_result.error_message = "PHREEQC failed";
stepResults(2).reaction_result.aux.phreeqc_run_status = 9;
stepResults(2).reaction_result.realized_interface_moles = 0;
stepResults(2).reaction_result.candidate_interface_moles = 0;
stepResults(2).transaction_status = 'rolled_back';
end

function result = baseStepResult()
result = struct();
result.diagnostics = struct('accepted', true, ...
    'component_absolute_residual_moles', 0, ...
    'charge_absolute_residual_eq', 0, ...
    'reasons', strings(0, 1));
result.diagnostics.full_mass_ledger = struct( ...
    'component_names', "H_reactant", ...
    'initial_component_moles_total', 1e-6, ...
    'final_component_moles_total', 1.19e-6, ...
    'transport_source_delta_moles_total', 2e-7, ...
    'transport_internal_flux_delta_moles_total', 0, ...
    'transport_boundary_delta_moles_total', 5e-8, ...
    'reaction_delta_moles_total', -1e-7, ...
    'remap_delta_moles_total', 4e-8, ...
    'component_residual_moles', 0, ...
    'initial_mineral_moles_total', 1, ...
    'reaction_mineral_delta_moles_total', -1e-7, ...
    'final_mineral_moles_total', 1 - 1e-7, ...
    'mineral_residual_moles', 0);
result.transaction_status = 'committed';
result.flow_diagnostics = struct( ...
    'inlet_flow_cm3_s', 1, ...
    'outlet_flow_cm3_s', 0.99, ...
    'inlet_outlet_relative_residual', 0.01);
result.transport_ledger = struct();
result.transport_ledger.dt_s = 0.1;
result.transport_ledger.component_names = {'H_reactant'};
result.transport_ledger.initial_moles_total = 1e-6;
result.transport_ledger.final_moles_total = 1.15e-6;
result.transport_ledger.source_delta_moles_total = 2e-7;
result.transport_ledger.internal_flux_delta_moles_total = 0;
result.transport_ledger.boundary_delta_moles_total = 5e-8;
result.transport_ledger.boundary_advective_delta_moles_total = 3e-8;
result.transport_ledger.boundary_diffusive_delta_moles_total = 2e-8;
result.transport_ledger.roundoff_suppressed_moles_total = 3e-16;
result.chemistry_ledger = struct();
result.chemistry_ledger.component_delta_moles_total = -1e-7;
result.geometry_info = struct();
result.geometry_info.solid_volume_before_cm3 = 1;
result.geometry_info.solid_volume_after_cm3 = 1 - 1e-7;
result.geometry_info.expected_solid_volume_change_cm3 = -1e-7;
result.geometry_info.actual_solid_volume_change_cm3 = -1e-7;
result.geometry_info.realized_mineral_moles = 1e-7;
result.geometry_info.max_displacement_over_h = 0.1;
result.reaction_result = struct();
result.reaction_result.converged = true;
result.reaction_result.failed_cells = [];
result.reaction_result.error_message = "";
result.reaction_result.realized_interface_moles = 1e-7;
result.reaction_result.candidate_interface_moles = 2e-7;
result.reaction_result.aux = struct('chemistry_mode', "strict_molins");
result.reaction_result.aux.phreeqc_run_status = 0;
end

function cleanupOutput(outputDir)
massPath = fullfile(outputDir, 'mass_balance_components.csv');
solidPath = fullfile(outputDir, 'solid_geometry_balance.csv');
chemistryPath = fullfile(outputDir, 'chemistry_step_status.csv');
rejectionPath = fullfile(outputDir, 'step_rejection_log.csv');
fluxPath = fullfile(outputDir, 'transport_flux_balance.csv');
clusterPath = fullfile(outputDir, 'reaction_cluster_diagnostics.csv');
if exist(massPath, 'file') == 2
    delete(massPath);
end
if exist(solidPath, 'file') == 2
    delete(solidPath);
end
if exist(chemistryPath, 'file') == 2
    delete(chemistryPath);
end
if exist(rejectionPath, 'file') == 2
    delete(rejectionPath);
end
if exist(fluxPath, 'file') == 2
    delete(fluxPath);
end
if exist(clusterPath, 'file') == 2
    delete(clusterPath);
end
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end
