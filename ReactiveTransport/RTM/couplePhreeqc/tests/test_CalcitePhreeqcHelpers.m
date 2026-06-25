function tests = test_CalcitePhreeqcHelpers
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
helperDir = fileparts(testDir);
testCase.TestData.helperDir = helperDir;
addpath(helperDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
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
verifyTextContains(testCase, text, 'pH 1 charge');
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

function testBuildInputUsesPrescribedCalciteReactionForTstMatch(testCase)
state = struct();
state.h_mol_cm3 = 1e-7;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 1e-7;
state.interface_area_cm2 = 2e-4;
state.water_volume_cm3 = 1e-3;
state.calcite_moles = 3e-6;
state.prescribed_calcite_dissolved_moles = 2e-9;

options = struct();
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 2.5e-4;
options.maxSpecificSurfaceArea = 100;

text = string(BuildCalcitePhreeqcInput(state, options));

verifyFalse(testCase, contains(text, 'RATES'), ...
    'Unified-H TST-match must not let PHREEQC run a second H-driven kinetic rate.');
verifyFalse(testCase, contains(text, 'KINETICS'), ...
    'Unified-H TST-match must receive prescribed CaCO3 dissolution instead of solving kinetics.');
verifyTextContains(testCase, text, 'REACTION 1');
verifyTextContains(testCase, text, 'CaCO3 1');
verifyTextContains(testCase, text, '0.002 moles');
verifyTextContains(testCase, text, '-headings KIN_DELTA_Calcite RATE_Calcite');
verifyTextContains(testCase, text, 'PUNCH 0 0');
end

function testBuildInputUsesReactionWaterVolumeForCutCellMixing(testCase)
state = struct();
state.h_mol_cm3 = 1e-7;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 1e-7;
state.interface_area_cm2 = 2e-4;
state.water_volume_cm3 = 1e-6;
state.reaction_water_volume_cm3 = 1e-3;
state.calcite_moles = 3e-6;
state.prescribed_calcite_dissolved_moles = 2e-9;

options = struct();
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 2.5e-4;

text = string(BuildCalcitePhreeqcInput(state, options));

verifyTextContains(testCase, text, 'REACTION 1');
verifyTextContains(testCase, text, '0.002 moles');
verifyFalse(testCase, contains(text, newline + "2 moles"), ...
    'Cut-cell reaction mixing volume should prevent extreme mol/kgw reaction inputs.');
end

function testBuildInputPrefersTotalComponentsOverFreeIonMolality(testCase)
state = struct();
state.h_mol_cm3 = 1e-6;
state.ca_mol_cm3 = 2e-8;
state.ca_total_mol_cm3 = 8e-8;
state.c_mol_cm3 = 3e-8;
state.c_total_mol_cm3 = 9e-8;
state.na_mol_cm3 = 4e-8;
state.na_total_mol_cm3 = 1.1e-7;
state.cl_mol_cm3 = 5e-8;
state.cl_total_mol_cm3 = 1.2e-7;
state.interface_area_cm2 = 1;
state.water_volume_cm3 = 1000;
state.calcite_moles = 1;

options = struct();
options.rateLaw = 'tst_match';

text = string(BuildCalcitePhreeqcInput(state, options));

verifyTextContains(testCase, text, 'Ca 8e-05');
verifyTextContains(testCase, text, 'C 9e-05');
verifyTextContains(testCase, text, 'Na 0.00011');
verifyTextContains(testCase, text, 'Cl 0.00012');
end

function testTstMatchBatchReturnsPrescribedDissolution(testCase)
databasePath = findPhreeqcDatabase(testCase);
assumeTrue(testCase, ~isempty(databasePath), 'No PHREEQC database found for COM integration test.');
assumeTrue(testCase, canCreateIPhreeqcCom(), 'IPhreeqcCOM is not available on this machine.');

state = struct();
state.h_mol_cm3 = 1e-7;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 1e-7;
state.interface_area_cm2 = 1;
state.water_volume_cm3 = 1000;
state.calcite_moles = 1;
state.prescribed_calcite_dissolved_moles = 3e-6;

options = struct();
options.databasePath = databasePath;
options.workDir = tempdir;
options.timeStepIndex = 9101;
options.timeStepSize = 1;
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 1e-4;
options.maxSpecificSurfaceArea = 10;
options.minActiveWaterVolumeFraction = 0;

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.calciteDissolvedMoles, state.prescribed_calcite_dissolved_moles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteRate_mol_s, ...
    state.prescribed_calcite_dissolved_moles ./ options.timeStepSize, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyGreaterThan(testCase, result.ca_total_mol_cm3, state.ca_mol_cm3);
verifyGreaterThan(testCase, result.c_total_mol_cm3, state.c_mol_cm3);
end

function testTstMatchPrescribedReactionActivatesReceiverWithoutCalciteInventory(testCase)
databasePath = findPhreeqcDatabase(testCase);
assumeTrue(testCase, ~isempty(databasePath), 'No PHREEQC database found for COM integration test.');
assumeTrue(testCase, canCreateIPhreeqcCom(), 'IPhreeqcCOM is not available on this machine.');

state = struct();
state.h_mol_cm3 = 0;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 0;
state.interface_area_cm2 = 0;
state.water_volume_cm3 = 1000;
state.calcite_moles = 0;
state.prescribed_calcite_dissolved_moles = 1e-6;

options = struct();
options.databasePath = databasePath;
options.workDir = tempdir;
options.timeStepIndex = 9104;
options.timeStepSize = 1;
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 1e-4;
options.minActiveWaterVolumeFraction = 0;

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.calciteDissolvedMoles, ...
    state.prescribed_calcite_dissolved_moles, 'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyGreaterThan(testCase, result.ca_total_mol_cm3, 0);
verifyGreaterThan(testCase, result.c_total_mol_cm3, 0);
end

function testTstMatchPrescribedReactionScalesThinWaterCells(testCase)
databasePath = findPhreeqcDatabase(testCase);
assumeTrue(testCase, ~isempty(databasePath), 'No PHREEQC database found for COM integration test.');
assumeTrue(testCase, canCreateIPhreeqcCom(), 'IPhreeqcCOM is not available on this machine.');

state = struct();
state.h_mol_cm3 = 2e-8;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 2e-8;
state.interface_area_cm2 = 1;
state.water_volume_cm3 = 1e-6;
state.calcite_moles = 1;
state.prescribed_calcite_dissolved_moles = 5e-13;

options = struct();
options.databasePath = databasePath;
options.workDir = tempdir;
options.timeStepIndex = 9102;
options.timeStepSize = 1;
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 1e-4;
options.maxSpecificSurfaceArea = Inf;
options.kineticsReservoirMoles = 1;
options.minActiveWaterVolumeFraction = 0;

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.calciteDissolvedMoles, state.prescribed_calcite_dissolved_moles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyGreaterThan(testCase, result.ca_total_mol_cm3, state.ca_mol_cm3);
verifyGreaterThan(testCase, result.c_total_mol_cm3, state.c_mol_cm3);
end

function testAgglomeratedDeltaUpdatesAllMembersConservatively(testCase)
base = struct();
base.h_mol_cm3 = [1e-7; 1e-7; 1e-7];
base.ca_total_mol_cm3 = [0; 0; 0];
base.c_total_mol_cm3 = [0; 0; 0];
base.na_total_mol_cm3 = [1e-5; 1e-5; 1e-5];
base.cl_total_mol_cm3 = [1.1e-4; 1.1e-4; 1.1e-4];

reactionInput = base;
reactionInput.ca_total_mol_cm3(1) = 0;
reactionInput.c_total_mol_cm3(1) = 0;

reactionOutput = reactionInput;
reactionOutput.ca_total_mol_cm3(1) = 2e-8;
reactionOutput.c_total_mol_cm3(1) = 2e-8;
reactionOutput.h_mol_cm3(1) = 8e-8;

waterVolumeCm3 = [1; 3; 6];
weightMatrix = sparse([1; 2; 3], [1; 1; 1], waterVolumeCm3 ./ sum(waterVolumeCm3), 3, 3);

updated = ApplyPhreeqcAgglomeratedDelta(base, reactionInput, reactionOutput, ...
    waterVolumeCm3, weightMatrix);

verifyEqual(testCase, updated.ca_total_mol_cm3, [2e-8; 2e-8; 2e-8], ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, updated.c_total_mol_cm3, [2e-8; 2e-8; 2e-8], ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, updated.h_mol_cm3, [8e-8; 8e-8; 8e-8], ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, sum(updated.ca_total_mol_cm3 .* waterVolumeCm3), ...
    2e-8 * sum(waterVolumeCm3), 'RelTol', 1e-12, 'AbsTol', 1e-20);
end

function testAgglomeratedStateUsesWaterWeightedInitialSolution(testCase)
state = struct();
state.h_mol_cm3 = [1e-7; 3e-7; 5e-7];
state.ca_mol_cm3 = [0; 2e-8; 4e-8];
state.c_mol_cm3 = [0; 1e-8; 3e-8];
state.na_mol_cm3 = [1e-5; 2e-5; 4e-5];
state.cl_mol_cm3 = [1.1e-4; 1.2e-4; 1.4e-4];
state.na_total_mol_cm3 = state.na_mol_cm3;
state.cl_total_mol_cm3 = state.cl_mol_cm3;
state.water_volume_cm3 = [1; 3; 6];
state.reaction_water_volume_cm3 = [10; 3; 6];

weightMatrix = sparse([1; 2; 3], [1; 1; 1], ...
    state.water_volume_cm3 ./ sum(state.water_volume_cm3), 3, 3);

mixed = BuildPhreeqcAgglomeratedState(state, state.water_volume_cm3, weightMatrix);

verifyEqual(testCase, mixed.h_mol_cm3(1), 4e-7, ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, mixed.ca_mol_cm3(1), 3e-8, ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, mixed.c_mol_cm3(1), 2.1e-8, ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, mixed.na_mol_cm3(1), 3.1e-5, ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, mixed.cl_mol_cm3(1), 1.31e-4, ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, mixed.na_total_mol_cm3(1), 3.1e-5, ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, mixed.cl_total_mol_cm3(1), 1.31e-4, ...
    'RelTol', 1e-12, 'AbsTol', 1e-20);
verifyEqual(testCase, mixed.water_volume_cm3(1), 10);
verifyEqual(testCase, mixed.h_mol_cm3(2:3), state.h_mol_cm3(2:3));
end

function testTstMatchPrescribedChemistryDoesNotOverdrivePhreeqc(testCase)
databasePath = findPhreeqcDatabase(testCase);
assumeTrue(testCase, ~isempty(databasePath), 'No PHREEQC database found for COM integration test.');
assumeTrue(testCase, canCreateIPhreeqcCom(), 'IPhreeqcCOM is not available on this machine.');

state = struct();
state.h_mol_cm3 = 0.0313951397792009 / 1000;
state.ca_mol_cm3 = 0.430233798081538 / 1000;
state.c_mol_cm3 = 0.43023379825878 / 1000;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 0.031300500717835 / 1000;
state.water_volume_cm3 = 1000;
state.interface_area_cm2 = 2373.2443460919 * 1e4;
state.calcite_moles = 1000;
state.prescribed_calcite_dissolved_moles = 1e-4;

options = struct();
options.databasePath = databasePath;
options.workDir = tempdir;
options.timeStepIndex = 9103;
options.timeStepSize = 1.6;
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 1e-4;
options.maxSpecificSurfaceArea = Inf;
options.kineticsReservoirMoles = 1;
options.minActiveWaterVolumeFraction = 0;

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.calciteDissolvedMoles, state.prescribed_calcite_dissolved_moles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyGreaterThan(testCase, result.pH, 1.5);
end

function testTstMatchReturnsHydrogenActivityForRateState(testCase)
databasePath = findPhreeqcDatabase(testCase);
assumeTrue(testCase, ~isempty(databasePath), 'No PHREEQC database found for COM integration test.');
assumeTrue(testCase, canCreateIPhreeqcCom(), 'IPhreeqcCOM is not available on this machine.');

state = struct();
state.h_mol_cm3 = 1e-5;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 1e-5;
state.interface_area_cm2 = 1;
state.water_volume_cm3 = 1000;
state.reaction_water_volume_cm3 = 1000;
state.calcite_moles = 1;
state.prescribed_calcite_dissolved_moles = 1e-4;

options = struct();
options.databasePath = databasePath;
options.workDir = tempdir;
options.timeStepIndex = 9105;
options.timeStepSize = 1;
options.rateLaw = 'tst_match';
options.minActiveWaterVolumeFraction = 0;

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.h_activity_mol_cm3, 10.^(-result.pH) ./ 1000, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyTrue(testCase, isfinite(result.h_mol_cm3));
verifyLessThan(testCase, result.h_activity_mol_cm3, state.h_mol_cm3, ...
    'TST-match rate state should use H+ activity, not raw molality.');
end

function testConfigureRunGroupKeepsDatabaseCalciteAsDefault(testCase)
state = oneCellAcidCalciteState();
cfg = struct();
cfg.rateCoefficientTST = 2e-4;
cfg.phreeqcKineticsCorrectionFactor = 7;

cfg = ConfigurePhreeqcRunGroup(cfg, 'phreeqc_database_calcite', 'unitstamp');
text = BuildCalcitePhreeqcInput(state, cfg);

verifyEqual(testCase, cfg.reactionModel, 'phreeqc');
verifyEqual(testCase, cfg.phreeqcRunGroup, 'phreeqc_database_calcite');
verifyEqual(testCase, cfg.phreeqcRateLaw, 'database_calcite');
verifyEqual(testCase, cfg.runName, 'phreeqc_database_calcite_unitstamp');
verifyFalse(testCase, contains(string(text), 'RATES'), ...
    'Database calcite group must keep using phreeqc-m.dat RATES.');
verifyTextContains(testCase, text, 'Calcite');
verifyTextContains(testCase, text, '-parms    0.0001  7');
end

function testConfigureRunGroupPassesTstMatchToInputBuilder(testCase)
state = oneCellAcidCalciteState();
cfg = struct();
cfg.rateCoefficientTST = 2e-4;
cfg.phreeqcKineticsCorrectionFactor = 7;
cfg.phreeqcMaxSpecificSurfaceArea = 10;

cfg = ConfigurePhreeqcRunGroup(cfg, 'external_tst_phreeqc', 'unitstamp');
text = BuildCalcitePhreeqcInput(state, cfg);

verifyEqual(testCase, cfg.reactionModel, 'phreeqc');
verifyEqual(testCase, cfg.phreeqcRunGroup, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.phreeqcRateLaw, 'tst_match');
verifyEqual(testCase, cfg.phreeqcTstRateCoefficient, cfg.rateCoefficientTST);
verifyTrue(testCase, isinf(cfg.phreeqcMaxSpecificSurfaceArea), ...
    'TST-match must not cap specific surface area if it is meant to reproduce legacy TST numerically.');
verifyEqual(testCase, cfg.runName, 'external_tst_phreeqc_unitstamp');
verifyTextContains(testCase, text, 'REACTION 1');
verifyFalse(testCase, contains(string(text), 'Calcite_TST_Match'));
end

function testConfigureRunGroupUsesExternalTstPhreeqcAsCanonicalName(testCase)
state = oneCellAcidCalciteState();
cfg = struct();
cfg.rateCoefficientTST = 2e-4;
cfg.phreeqcMaxSpecificSurfaceArea = 10;

cfg = ConfigurePhreeqcRunGroup(cfg, 'external_tst_phreeqc', 'unitstamp');
text = BuildCalcitePhreeqcInput(state, cfg);

verifyEqual(testCase, cfg.reactionModel, 'phreeqc');
verifyEqual(testCase, cfg.phreeqcRunGroup, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.phreeqcRateLaw, 'tst_match');
verifyEqual(testCase, cfg.chemistryMode, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.chemistrySemantics, ...
    'explicit external TST rate + PHREEQC equilibrium closure');
verifyFalse(testCase, isfield(cfg, 'legacyPhreeqcRunGroupAlias'));
verifyEqual(testCase, cfg.runName, 'external_tst_phreeqc_unitstamp');
verifyTrue(testCase, isinf(cfg.phreeqcMaxSpecificSurfaceArea), ...
    'External TST PHREEQC must keep the legacy TST-match no-cap behavior.');
verifyTextContains(testCase, text, 'REACTION 1');
verifyFalse(testCase, contains(string(text), 'KINETICS'), ...
    'External TST PHREEQC must prescribe realized calcite moles instead of running PHREEQC kinetics.');
end

function testConfigureRunGroupKeepsPhreeqcTstMatchAsWarningAlias(testCase)
cfg = struct();
cfg.rateCoefficientTST = 2e-4;

cfg = verifyWarning(testCase, ...
    @() ConfigurePhreeqcRunGroup(cfg, 'phreeqc_tst_match', 'unitstamp'), ...
    'RTSPHEM:Phreeqc:DeprecatedRunGroup');

verifyEqual(testCase, cfg.phreeqcRunGroup, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.legacyPhreeqcRunGroupAlias, 'phreeqc_tst_match');
verifyEqual(testCase, cfg.phreeqcRateLaw, 'tst_match');
verifyEqual(testCase, cfg.chemistryMode, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.runName, 'external_tst_phreeqc_unitstamp');
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

function testExternalTstOnlyRunsPhreeqcForPrescribedReactionCells(testCase)
databasePath = findPhreeqcDatabase(testCase);
assumeTrue(testCase, ~isempty(databasePath), 'No PHREEQC database found for mock integration test.');

raw = {
    'soln', 'pH', 'charge(eq)', 'Ca', 'C', 'Na', 'Cl', ...
        'm_H+', 'm_Ca+2', 'm_HCO3-', 'm_CO3-2', 'm_Cl-', 'm_Na+', ...
        'si_Calcite', 'KIN_DELTA_Calcite', 'RATE_Calcite';
    1, 7.0, 0, 1e-6, 1e-6, 1e-5, 1e-5, ...
        1e-7, 1e-6, 1e-6, 0, 1e-5, 1e-5, ...
        0, 0, 0
    };
mockEngine = MockIPhreeqcEngine(raw);
session = rtm.phreeqc.PhreeqcSession(struct('engine', mockEngine));
cleanupSession = onCleanup(@() session.close());

state = struct();
state.h_mol_cm3 = [1e-7; 1e-7; 1e-7];
state.ca_mol_cm3 = [0; 0; 0];
state.c_mol_cm3 = [0; 0; 0];
state.na_mol_cm3 = [1e-8; 1e-8; 1e-8];
state.cl_mol_cm3 = [1e-8; 1e-8; 1e-8];
state.interface_area_cm2 = [1; 0; 0];
state.water_volume_cm3 = [1; 1; 1];
state.calcite_moles = [1; 1; 1];
state.prescribed_calcite_dissolved_moles = [1e-9; 0; 0];

options = struct();
options.databasePath = databasePath;
options.phreeqcSession = session;
options.writeInputFiles = false;
options.rateLaw = 'external_tst_phreeqc';
options.timeStepSize = 1;

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, mockEngine.RunStringCount, 1);
verifyTextContains(testCase, mockEngine.LastInputText, '-cells 1-1');
verifyEqual(testCase, result.calciteDissolvedMoles, [1e-9; 0; 0], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, result.ca_total_mol_cm3(2:3), state.ca_mol_cm3(2:3), ...
    'AbsTol', 0);
verifyEqual(testCase, result.na_total_mol_cm3(2:3), state.na_mol_cm3(2:3), ...
    'AbsTol', 0);

clear cleanupSession;
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
verifyEqual(testCase, result.calcite_cell_rate_mol_s, [3e-5; 4e-5], 'AbsTol', 1e-15);
verifyFalse(testCase, isfield(result, 'calciteRate_mol_dm2_s'));
verifyEqual(testCase, result.calciteDeltaMoles, [-3e-6; -4e-6], 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteDissolvedMoles, [3e-6; 4e-6], 'AbsTol', 1e-15);
verifyEqual(testCase, result.solutionNumber, [1; 2]);
end

function testDiagnosticInferCalciteDissolutionFromTotalsWhenKinDeltaIsUnavailable(testCase)
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
result.calcite_cell_rate_mol_s = [NaN; NaN];

updated = DiagnosticInferCalciteDissolutionFromTotals(result, state, 2);

verifyEqual(testCase, updated.diagnostic_inferred_calcite_dissolved_moles, ...
    [5e-8; 8e-11], 'AbsTol', 1e-15);
verifyEqual(testCase, updated.diagnostic_inferred_calcite_delta_moles, ...
    [-5e-8; -8e-11], 'AbsTol', 1e-15);
verifyEqual(testCase, updated.diagnostic_inferred_calcite_cell_rate_mol_s, ...
    [2.5e-8; 4e-11], 'AbsTol', 1e-15);
verifyEqual(testCase, updated.calciteDissolvedMoles, [0; 0], 'AbsTol', 0, ...
    'Diagnostic inference must not overwrite production PHREEQC mineral deltas.');
end

function testTstMatchInterfaceRateUsesLegacyPerAreaFormula(testCase)
result = struct();
result.calciteDissolvedMoles = [1e-20; 2e-20; 3e-20];
interfaceAreaCm2 = [1; 2; 0];
timeStepSize = 54;
legacyHMolCm3 = [1e-5; 2e-5; 3e-5];
options = struct();
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 8.9125093813374588e-5;

ratePerArea = ComputePhreeqcInterfaceRatePerArea( ...
    result, interfaceAreaCm2, timeStepSize, options, legacyHMolCm3);

expected = legacyHMolCm3(:) * 1000 * options.rateCoefficientTST;
expected(interfaceAreaCm2(:) <= 0) = 0;
verifyEqual(testCase, ratePerArea, expected, 'RelTol', 1e-12, 'AbsTol', 1e-15);
end

function testTstMatchInterfaceRatePrefersHydrogenActivity(testCase)
result = struct();
result.h_mol_cm3 = [1e-5; 2e-5];
result.h_activity_mol_cm3 = [4e-6; 6e-6];
interfaceAreaCm2 = [1; 1];
timeStepSize = 1;
options = struct();
options.rateLaw = 'tst_match';
options.rateCoefficientTST = 2e-4;

ratePerArea = ComputePhreeqcInterfaceRatePerArea( ...
    result, interfaceAreaCm2, timeStepSize, options);

expected = result.h_activity_mol_cm3(:) * 1000 * options.rateCoefficientTST;
verifyEqual(testCase, ratePerArea, expected, 'RelTol', 1e-12, 'AbsTol', 1e-15);
end

function testLegacyEquivalentTstTriangleRateUsesVertexMinimum(testCase)
gridV0T = [1 2 3; 2 4 3];
vertexTriMatrix = cell(1, 4);
vertexTriMatrix{1} = 1;
vertexTriMatrix{2} = [1 2];
vertexTriMatrix{3} = [1 2];
vertexTriMatrix{4} = 2;
triangleHydrogen = [1e-5; 4e-6];
interfaceArea = [1; 1];
rateCoefficient = 2e-4;

ratePerArea = ComputeLegacyEquivalentTstTriangleRate( ...
    triangleHydrogen, interfaceArea, gridV0T, vertexTriMatrix, 4, rateCoefficient);

vertexRate = [1e-5; 4e-6; 4e-6; 4e-6] * 1000 * rateCoefficient;
expected = [mean(vertexRate([1 2 3])); mean(vertexRate([2 4 3]))];
verifyEqual(testCase, ratePerArea, expected, 'RelTol', 1e-12, 'AbsTol', 1e-15);
end

function testTstMatchStepDissolutionUsesOneDriverForMotionAndMoles(testCase)
gridV0T = [1 2 3; 2 4 3];
vertexTriMatrix = cell(1, 4);
vertexTriMatrix{1} = 1;
vertexTriMatrix{2} = [1 2];
vertexTriMatrix{3} = [1 2];
vertexTriMatrix{4} = 2;
driverMolCm3 = [1e-5; 4e-6];
interfaceAreaCm2 = [2; 3];
calciteMoles = [1e-4; 1e-12];
dt = 5;
rateCoefficient = 2e-4;
molarVolume = 36.9;

step = ComputeTstMatchStepDissolution(driverMolCm3, interfaceAreaCm2, ...
    calciteMoles, dt, gridV0T, vertexTriMatrix, 4, rateCoefficient, molarVolume);

vertexRate = [1e-5; 4e-6; 4e-6; 4e-6] * 1000 * rateCoefficient;
expectedRate = [mean(vertexRate([1 2 3])); mean(vertexRate([2 4 3]))];
expectedMoles = min(expectedRate .* interfaceAreaCm2(:) .* dt, calciteMoles(:));

expectedRealizedRate = expectedMoles ./ interfaceAreaCm2(:) ./ dt;

verifyEqual(testCase, step.candidateRatePerArea_mol_cm2_s, expectedRate, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, step.ratePerArea_mol_cm2_s, expectedRealizedRate, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, step.prescribedMoles, expectedMoles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, step.candidateMoles, expectedRate .* interfaceAreaCm2(:) .* dt, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, step.realizedMoles, expectedMoles, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, step.inventoryLimited, [false; true]);
verifyEqual(testCase, step.normalSpeed_cm_s, expectedInventoryLimitedVertexSpeed( ...
    gridV0T, vertexTriMatrix, expectedMoles, interfaceAreaCm2, dt, molarVolume), ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
end

function testTstMatchStepDissolutionStopsMotionWhenInventoryIsExhausted(testCase)
gridV0T = [1 2 3];
vertexTriMatrix = cell(1, 3);
vertexTriMatrix{1} = 1;
vertexTriMatrix{2} = 1;
vertexTriMatrix{3} = 1;
driverMolCm3 = 1e-5;
interfaceAreaCm2 = 2;
calciteMoles = 0;
dt = 5;
rateCoefficient = 2e-4;
molarVolume = 36.9;

step = ComputeTstMatchStepDissolution(driverMolCm3, interfaceAreaCm2, ...
    calciteMoles, dt, gridV0T, vertexTriMatrix, 3, rateCoefficient, molarVolume);

verifyGreaterThan(testCase, step.candidateMoles, 0);
verifyEqual(testCase, step.realizedMoles, 0, 'AbsTol', 0);
verifyEqual(testCase, step.prescribedMoles, 0, 'AbsTol', 0);
verifyTrue(testCase, all(step.inventoryLimited));
verifyEqual(testCase, step.ratePerArea_mol_cm2_s, 0, 'AbsTol', 0);
verifyEqual(testCase, step.normalSpeed_cm_s, zeros(3, 1), 'AbsTol', 0);
end

function testTstMatchInterfaceRateAcceptsRunGroupAlias(testCase)
result = struct();
result.calciteDissolvedMoles = [9e-8; 8e-8];
interfaceAreaCm2 = [1; 1];
timeStepSize = 10;
legacyHMolCm3 = [1e-5; 2e-5];
options = struct();
options.rateLaw = 'phreeqc_tst_match';
options.rateCoefficientTST = 8.9125093813374588e-5;

ratePerArea = ComputePhreeqcInterfaceRatePerArea( ...
    result, interfaceAreaCm2, timeStepSize, options, legacyHMolCm3);

expected = legacyHMolCm3(:) * 1000 * options.rateCoefficientTST;
verifyEqual(testCase, ratePerArea, expected, 'RelTol', 1e-12, 'AbsTol', 1e-15);
end

function testDatabaseCalciteInterfaceRateUsesPhreeqcDissolution(testCase)
result = struct();
result.calciteDissolvedMoles = [2e-8; 6e-8; 3e-8];
interfaceAreaCm2 = [1; 2; 0];
timeStepSize = 2;
legacyHMolCm3 = [1e-5; 2e-5; 3e-5];
options = struct();
options.rateLaw = 'database_calcite';
options.rateCoefficientTST = 8.9125093813374588e-5;

ratePerArea = ComputePhreeqcInterfaceRatePerArea( ...
    result, interfaceAreaCm2, timeStepSize, options, legacyHMolCm3);

verifyEqual(testCase, ratePerArea, [1e-8; 1.5e-8; 0], 'AbsTol', 1e-15);
end

function vertexSpeed = expectedInventoryLimitedVertexSpeed( ...
    gridV0T, vertexTriMatrix, realizedMoles, interfaceAreaCm2, dt, molarVolume)
realizedRate = zeros(size(realizedMoles(:)));
active = interfaceAreaCm2(:) > 0 & dt > 0;
realizedRate(active) = realizedMoles(active) ./ interfaceAreaCm2(active) ./ dt;
vertexSpeed = zeros(numel(vertexTriMatrix), 1);
for iVertex = 1:numel(vertexTriMatrix)
    triIndices = vertexTriMatrix{1, iVertex};
    values = realizedRate(triIndices);
    values = values(isfinite(values));
    if isempty(values)
        vertexSpeed(iVertex) = 0;
    else
        vertexSpeed(iVertex) = mean(values) * molarVolume;
    end
end
end

function verifyTextContains(testCase, text, pattern)
verifyTrue(testCase, contains(string(text), pattern), sprintf('Expected text to contain: %s', pattern));
end

function state = oneCellAcidCalciteState()
state = struct();
state.h_mol_cm3 = 1e-7;
state.ca_mol_cm3 = 0;
state.c_mol_cm3 = 0;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = 1e-7;
state.interface_area_cm2 = 1;
state.water_volume_cm3 = 1000;
state.calcite_moles = 1;
end

function databasePath = findPhreeqcDatabase(testCase)
helperDir = testCase.TestData.helperDir;
rtmDir = fileparts(helperDir);
candidates = {
    fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc-m.dat')
    fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc.dat')
    fullfile(helperDir, 'phreeqc-m.dat')
    fullfile(rtmDir, 'couplePhreeqc', 'phreeqc-m.dat')
    'C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat'
    };
databasePath = '';
for iCandidate = 1:numel(candidates)
    if exist(candidates{iCandidate}, 'file') == 2
        databasePath = candidates{iCandidate};
        return;
    end
end
end

function ok = canCreateIPhreeqcCom()
ok = false;
try
    iphreeqc = actxserver('IPhreeqcCOM.Object');
    cleanupObj = onCleanup(@() delete(iphreeqc));
    ok = true;
    clear cleanupObj;
catch
    ok = false;
end
end
