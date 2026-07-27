function tests = test_WritePhreeqcSpeciesTable
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
helperDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.helperDir = helperDir;
addpath(helperDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testSpeciesTableWritesMassLedgerFields(testCase)
outputDir = tempname;
mkdir(outputDir);
testCase.addTeardown(@() cleanupOutput(outputDir));

grid = struct();
grid.numT = 2;
grid.baryT = [0 0; 1 1];

speciesData = struct();
speciesData.h_mol_cm3 = [1e-7; 2e-7];
speciesData.ca_total_mol_cm3 = [1e-8; 2e-8];
speciesData.c_total_mol_cm3 = [1e-8; 2e-8];
speciesData.calciteDissolvedMoles = [1e-9; 2e-9];
speciesData.water_phase_ca_delta_moles = [1e-9; 2e-9];
speciesData.water_phase_c_delta_moles = [1e-9; 2e-9];
speciesData.calcite_ca_stoich_residual_moles = [0; 0];
speciesData.calcite_c_stoich_residual_moles = [0; 0];
speciesData.interfaceRealizedCalciteMoles = [1e-9; 2e-9];
speciesData.ionicStrength_mol_kgw = [0.01; 0.02];
speciesData.fluidConductivity_S_m = [0.1; 0.2];
speciesData.h_activity_dimensionless = [1e-4; 2e-4];
speciesData.water_volume_cm3 = [1; 2];
speciesData.interface_area_cm2 = [3; 4];

WritePhreeqcSpeciesTable(outputDir, 3, 12.5, grid, speciesData);
T = readtable(fullfile(outputDir, 'phreeqc_species_0003.csv'));

verifyTrue(testCase, any(strcmp(T.Properties.VariableNames, ...
    'water_phase_ca_delta_moles')));
verifyTrue(testCase, any(strcmp(T.Properties.VariableNames, ...
    'water_phase_c_delta_moles')));
verifyTrue(testCase, any(strcmp(T.Properties.VariableNames, ...
    'calcite_ca_stoich_residual_moles')));
verifyTrue(testCase, any(strcmp(T.Properties.VariableNames, ...
    'calcite_c_stoich_residual_moles')));
verifyEqual(testCase, T.water_phase_ca_delta_moles, [1e-9; 2e-9], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, T.fluidConductivity_S_m, [0.1; 0.2], 'AbsTol', 1e-12);
verifyEqual(testCase, T.interface_area_cm2, [3; 4], 'AbsTol', 1e-12);
end

function testSeismoelectricSummaryUsesPhysicalWeightsAndLeavesZetaUncalibrated(testCase)
outputDir = tempname;
mkdir(outputDir);
csvFile = fullfile(outputDir, 'seismoelectric_chemistry.csv');
testCase.addTeardown(@() cleanupSummaryOutput(outputDir, csvFile));

rtState = struct('porosity', 0.2, 'permeability_mD', 10, 'k_k0', 2, ...
    'surface_area_cm2', 5, 'grain_volume_cm3', 6, ...
    'outlet_h_conc_mol_cm3', 1e-6, 'hydraulic_tortuosity', 1.5);
speciesData = struct();
speciesData.water_volume_cm3 = [1; 3];
speciesData.interface_area_cm2 = [2; 0];
speciesData.pH = [2; 4];
speciesData.ionicStrength_mol_kgw = [0.1; 0.3];
speciesData.fluidConductivity_S_m = [1; 3];
speciesData.phreeqcPercentError = [1e-4; 2e-4];
speciesData.chargeBalance = [1e-8; 2e-8];
speciesData.h_mol_cm3 = [1e-5; 1e-7];
speciesData.ca_total_mol_cm3 = [1e-6; 3e-6];
speciesData.c_total_mol_cm3 = [2e-6; 4e-6];
speciesData.na_total_mol_cm3 = [1e-5; 1e-5];
speciesData.cl_total_mol_cm3 = [1.2e-5; 1.4e-5];
speciesData.ca_mol_cm3 = [0.8e-6; 2.4e-6];
speciesData.hco3_mol_cm3 = [0.5e-6; 1.5e-6];
speciesData.co3_mol_cm3 = [1e-9; 3e-9];
speciesData.h_activity_dimensionless = [1e-2; 1e-4];
speciesData.ca_activity_dimensionless = [1e-3; 3e-3];
speciesData.hco3_activity_dimensionless = [5e-4; 1.5e-3];
speciesData.co3_activity_dimensionless = [1e-7; 3e-7];
speciesData.calciteSI = [-2; -1];

AppendSeismoelectricChemistrySummary(csvFile, 3, 12.5, rtState, speciesData);
T = readtable(csvFile);

verifyEqual(testCase, T.pH_pore_volume_weighted, 3.5, 'AbsTol', 1e-12);
verifyEqual(testCase, T.pH_interface_area_weighted, 2, 'AbsTol', 1e-12);
verifyEqual(testCase, T.IonicStrength_mol_kgw, 0.25, 'AbsTol', 1e-12);
verifyEqual(testCase, T.FluidConductivity_S_m, 2.5, 'AbsTol', 1e-12);
verifyEqual(testCase, T.OutletHConc, 1e-6, 'AbsTol', 1e-15);
verifyTrue(testCase, isnan(T.ZetaPotential_V));
verifyTrue(testCase, isnan(T.Tortuosity));
verifyEqual(testCase, string(T.ChemistryStatus), "ready_for_calibrated_zeta_model");
verifyEqual(testCase, string(T.ZetaModelStatus), ...
    "not_computed_requires_calcite_scm_or_measurement");
verifyEqual(testCase, string(T.SEInputStatus), ...
    "blocked_requires_zeta_and_alpha_inf_calibration");
end

function cleanupOutput(outputDir)
deleteIfFile(fullfile(outputDir, 'phreeqc_species_0003.csv'));
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end

function deleteIfFile(pathValue)
if exist(pathValue, 'file') == 2
    delete(pathValue);
end
end

function cleanupSummaryOutput(outputDir, csvFile)
deleteIfFile(csvFile);
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end
