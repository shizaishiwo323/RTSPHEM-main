function tests = test_PackagePhreeqcIO
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

function testPackageBuildCalciteInputDelegatesToRuntimeBuilder(testCase)
state = minimalPhreeqcState();
text = string(rtm.phreeqc.BuildCalcitePhreeqcInput(state, struct( ...
    'rateLaw', 'external_tst_phreeqc', ...
    'timeStepSize', 1)));

verifyTrue(testCase, contains(text, "SELECTED_OUTPUT"));
verifyTrue(testCase, contains(text, "KIN_DELTA"));
verifyTrue(testCase, contains(text, "Alkalinity"));
end

function testPackageParseSelectedOutputDelegatesToRuntimeParser(testCase)
raw = {
    'soln', 'pH', 'Ca(mol/kgw)', 'C(mol/kgw)', 'Na(mol/kgw)', 'Cl(mol/kgw)', 'Alkalinity(mol/kgw)', 'KIN_DELTA_Calcite', 'RATE_Calcite';
    1, 7.1, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, -6e-9, 7e-10
    };

result = rtm.phreeqc.ParsePhreeqcSelectedOutput(raw, 1);

verifyEqual(testCase, result.ca_total_mol_cm3, 1e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, result.c_total_mol_cm3, 2e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, result.alkalinity_mol_cm3, 5e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, result.calciteDeltaMoles, -6e-9, 'AbsTol', 1e-18);
verifyEqual(testCase, result.calciteDissolvedMoles, 6e-9, 'AbsTol', 1e-18);
verifyEqual(testCase, result.calcite_cell_rate_mol_s, 7e-10, 'AbsTol', 1e-18);
end

function state = minimalPhreeqcState()
state = struct();
state.h_mol_cm3 = 1e-10;
state.ca_total_mol_cm3 = 0;
state.c_total_mol_cm3 = 0;
state.na_total_mol_cm3 = 1e-8;
state.cl_total_mol_cm3 = 1e-8;
state.alkalinity_mol_cm3 = 2e-9;
state.interface_area_cm2 = 1;
state.water_volume_cm3 = 1;
state.reaction_water_volume_cm3 = 1;
state.calcite_moles = 1;
state.prescribed_calcite_dissolved_moles = 1e-12;
end
