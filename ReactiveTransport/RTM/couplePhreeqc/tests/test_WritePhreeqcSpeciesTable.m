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
