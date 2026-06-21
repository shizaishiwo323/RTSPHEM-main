function tests = test_CalcitePhreeqcHelpers
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
helperDir = fileparts(testDir);
testCase.TestData.helperDir = helperDir;
addpath(helperDir);
end

function teardownOnce(testCase)
rmpath(testCase.TestData.helperDir);
end

function testBuildInputContainsCalciteHclNaclChemistry(testCase)
state = struct();
state.h_mol_cm3 = [1e-4; 5e-5];
state.ca_mol_cm3 = [0; 1e-8];
state.c_mol_cm3 = [0; 2e-8];
state.na_mol_cm3 = [1e-5; 1e-5];
state.cl_mol_cm3 = [1.1e-4; 6e-5];
state.interface_area_cm2 = [2e-4; 1e-4];
state.water_volume_cm3 = [1e-3; 2e-3];
state.calcite_moles = [3e-6; 4e-6];

options = struct();
options.timeStepSize = 0.25;
options.temperatureC = 25;
options.solutionUnits = 'mol/kgw';
options.mineralName = 'Calcite';
options.kineticsCorrectionFactor = 1;

text = BuildCalcitePhreeqcInput(state, options);

verifyFalse(testCase, contains(string(text), 'RATES'), ...
    'Calcite kinetics must come from phreeqc-m.dat, not an inline RATES block.');
verifyFalse(testCase, contains(string(text), 'Calcite_RTSPHEM'), ...
    'The simplified Calcite_RTSPHEM rate law must not be used.');
verifyTextContains(testCase, text, 'SOLUTION 1');
verifyTextContains(testCase, text, 'pH 1');
verifyFalse(testCase, contains(string(text), 'water 1e-06'), ...
    'Active PHREEQC cells should use the reference solution water mass, not tiny PNM water masses.');
verifyTextContains(testCase, text, 'Na 0.01');
verifyTextContains(testCase, text, 'Cl 0.11');
verifyTextContains(testCase, text, 'KINETICS 1');
verifyTextContains(testCase, text, 'Calcite');
verifyTextContains(testCase, text, '-formula  CaCO3  1');
verifyTextContains(testCase, text, '-m        1');
verifyTextContains(testCase, text, '-m0       1');
verifyTextContains(testCase, text, '-parms    0.02  1');
verifyTextContains(testCase, text, '-bad_step_max 5000');
verifyTextContains(testCase, text, 'RUN_CELLS');
verifyTextContains(testCase, text, '-cells 1-2');
verifyTextContains(testCase, text, '-time_step 0.25');
verifyTextContains(testCase, text, '-molalities H+ Ca+2 HCO3- CO3-2 Cl- Na+');
verifyTextContains(testCase, text, 'KIN_DELTA("Calcite")/KIN_TIME');
end

function testBuildInputCapsSpecificSurfaceAreaForThinWaterCells(testCase)
state = struct();
state.h_mol_cm3 = 1e-12;
state.interface_area_cm2 = 1e-4;
state.water_volume_cm3 = 1e-12;
state.calcite_moles = 1e-8;

options = struct();
options.maxSpecificSurfaceArea = 100;
options.badStepMax = 1234;

text = BuildCalcitePhreeqcInput(state, options);

verifyTextContains(testCase, text, '-parms    100  1');
verifyTextContains(testCase, text, '-bad_step_max 1234');
end

function testBuildInputUsesNeutralWaterFloorForZeroAcid(testCase)
state = struct();
state.h_mol_cm3 = 0;
state.interface_area_cm2 = 0;
state.water_volume_cm3 = 1e-6;
state.calcite_moles = 0;

text = BuildCalcitePhreeqcInput(state, struct());

verifyTextContains(testCase, text, 'pH 7');
verifyFalse(testCase, contains(string(text), 'pH 14'));
end

function testBuildInputDeactivatesKineticsWhenNoInterfaceExists(testCase)
state = struct();
state.h_mol_cm3 = 0;
state.interface_area_cm2 = 0;
state.water_volume_cm3 = 1e-12;
state.calcite_moles = 1e-6;

text = BuildCalcitePhreeqcInput(state, struct());

verifyTextContains(testCase, text, '-m        1e-30');
verifyTextContains(testCase, text, '-m0       1e-30');
verifyTextContains(testCase, text, '-parms    0  1');
verifyFalse(testCase, contains(string(text), '-m        1e-06'));
end

function testRunBatchSkipsComWhenNoActiveCellsExist(testCase)
state = struct();
state.h_mol_cm3 = [0; 0];
state.ca_mol_cm3 = [0; 0];
state.c_mol_cm3 = [0; 0];
state.na_mol_cm3 = [0; 0];
state.cl_mol_cm3 = [0; 0];
state.interface_area_cm2 = [0; 0];
state.water_volume_cm3 = [0; 0];
state.calcite_moles = [1e-6; 1e-6];

options = struct();
options.workDir = tempdir;
options.timeStepIndex = 9001;
options.databasePath = 'missing-database-for-inactive-test.dat';

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.pH, [7; 7], 'AbsTol', 1e-12);
verifyEqual(testCase, result.calciteDissolvedMoles, [0; 0], 'AbsTol', 0);
verifyEqual(testCase, result.solutionNumber, [1; 2]);
end

function testRunBatchSkipsNeutralInterfaceCellsByDefault(testCase)
state = struct();
state.h_mol_cm3 = 0;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 0;
state.interface_area_cm2 = 1e-4;
state.water_volume_cm3 = 1e-6;
state.calcite_moles = 1e-6;

options = struct();
options.workDir = tempdir;
options.timeStepIndex = 9002;
options.databasePath = 'missing-database-for-neutral-interface-test.dat';

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.pH, 7, 'AbsTol', 1e-12);
verifyEqual(testCase, result.calciteDissolvedMoles, 0, 'AbsTol', 0);
end

function testParseSelectedOutputKeepsLastRowsAndNamedSpecies(testCase)
raw = {
    'sim', 'state', 'soln', 'pH', 'charge(eq)', 'Ca(mol/kgw)', 'C(mol/kgw)', 'Na(mol/kgw)', 'Cl(mol/kgw)', 'm_H+(mol/kgw)', 'm_Ca+2(mol/kgw)', 'm_HCO3-(mol/kgw)', 'm_CO3-2(mol/kgw)', 'm_Cl-(mol/kgw)', 'm_Na+(mol/kgw)', 'si_Calcite', 'KIN_DELTA_Calcite', 'RATE_Calcite';
    1, 'initial', 1, 2.0, 0.01, 0, 0, 0.01, 0.11, 0.1, 0, 0, 0, 0.11, 0.01, -5, 0, 1e-5;
    1, 'initial', 2, 3.0, 0.02, 0, 0, 0.01, 0.06, 0.001, 0, 0, 0, 0.06, 0.01, -4, 0, 2e-5;
    1, 'react', 1, 2.1, 0.03, 1e-4, 1e-4, 0.01, 0.11, 0.08, 1e-4, 8e-5, 1e-8, 0.11, 0.01, -3, -3e-6, 3e-5;
    1, 'react', 2, 3.1, 0.04, 2e-4, 2e-4, 0.01, 0.06, 8e-4, 2e-4, 1.8e-4, 2e-8, 0.06, 0.01, -2, -4e-6, 4e-5
};

result = ParsePhreeqcSelectedOutput(raw, 2);

verifyEqual(testCase, result.pH, [2.1; 3.1], 'AbsTol', 1e-12);
verifyEqual(testCase, result.h_mol_cm3, [0.08; 8e-4] / 1000, 'AbsTol', 1e-15);
verifyEqual(testCase, result.ca_mol_cm3, [1e-4; 2e-4] / 1000, 'AbsTol', 1e-15);
verifyEqual(testCase, result.hco3_mol_cm3, [8e-5; 1.8e-4] / 1000, 'AbsTol', 1e-15);
verifyEqual(testCase, result.co3_mol_cm3, [1e-8; 2e-8] / 1000, 'AbsTol', 1e-15);
verifyEqual(testCase, result.cl_mol_cm3, [0.11; 0.06] / 1000, 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteRate_mol_dm2_s, [3e-5; 4e-5], 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteDeltaMoles, [-3e-6; -4e-6], 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteDissolvedMoles, [3e-6; 4e-6], 'AbsTol', 1e-15);
verifyEqual(testCase, result.solutionNumber, [1; 2]);
end

function testInferCalciteDissolutionFromTotalsWhenKinDeltaIsUnavailable(testCase)
state = struct();
state.ca_mol_cm3 = [0; 1e-8];
state.c_mol_cm3 = [0; 2e-8];
state.water_volume_cm3 = [1e-3; 2e-3];

result = struct();
result.ca_total_mol_cm3 = [6e-5; 5e-8];
result.c_total_mol_cm3 = [5e-5; 6e-8];
result.calciteDeltaMoles = [0; 0];
result.calciteDissolvedMoles = [0; 0];
result.calciteRate_mol_s = [NaN; NaN];
result.calciteRate_mol_dm2_s = [NaN; NaN];

updated = InferCalciteDissolutionFromTotals(result, state, 2);

verifyEqual(testCase, updated.calciteDissolvedMoles, [5e-8; 8e-11], 'AbsTol', 1e-15);
verifyEqual(testCase, updated.calciteDeltaMoles, [-5e-8; -8e-11], 'AbsTol', 1e-15);
verifyEqual(testCase, updated.calciteRate_mol_s, [2.5e-8; 4e-11], 'AbsTol', 1e-15);
end

function verifyTextContains(testCase, text, pattern)
verifyTrue(testCase, contains(string(text), pattern), sprintf('Expected text to contain: %s', pattern));
end
