function tests = test_ChemistryModeSemantics
tests = functiontests(localfunctions);
end

function setupOnce(~)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));
end

function testStrictMolinsConfigDisablesPhreeqc(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');

verifyEqual(testCase, cfg.chemistry.mode, 'strict_molins');
verifyEqual(testCase, cfg.chemistry.semantics, ...
    "strict Molins benchmark without PHREEQC");
verifyEqual(testCase, cfg.phreeqc.engine, 'none');
verifyEqual(testCase, cfg.phreeqc.databasePolicy, 'not_used');
end

function testExternalTstPhreeqcUsesReactionClosureNotKinetics(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('integration_phreeqc');
state = phreeqcStateWithPrescribedReaction();
options = struct('rateLaw', 'external_tst_phreeqc', ...
    'timeStepSize', 1, 'rateCoefficientTST', 0.1);

text = string(BuildCalcitePhreeqcInput(state, options));

verifyEqual(testCase, cfg.chemistry.mode, 'external_tst_phreeqc');
verifyEqual(testCase, cfg.chemistry.semantics, ...
    "explicit external TST rate + PHREEQC equilibrium closure");
verifyTrue(testCase, contains(text, 'REACTION 1'));
verifyFalse(testCase, contains(text, 'KINETICS'));
end

function testPhreeqcKineticsUsesKineticsNotPrescribedReaction(testCase)
cfg = rtm.config.ValidateReactiveTransportConfig(struct( ...
    'chemistry', struct('mode', 'phreeqc_kinetics'), ...
    'phreeqc', struct('engine', 'iphreeqc_com', ...
        'databasePolicy', 'exact_local')));
state = phreeqcStateWithPrescribedReaction();
state = rmfield(state, 'prescribed_calcite_dissolved_moles');
options = struct('rateLaw', 'database_calcite', 'timeStepSize', 1);

text = string(BuildCalcitePhreeqcInput(state, options));

verifyEqual(testCase, cfg.chemistry.mode, 'phreeqc_kinetics');
verifyEqual(testCase, cfg.chemistry.semantics, ...
    "PHREEQC native KINETICS/RATES");
verifyTrue(testCase, contains(text, 'KINETICS 1'));
verifyFalse(testCase, contains(text, 'REACTION 1'));
end

function state = phreeqcStateWithPrescribedReaction()
state = struct();
state.h_mol_cm3 = 1e-7;
state.ca_total_mol_cm3 = 0;
state.c_total_mol_cm3 = 0;
state.na_total_mol_cm3 = 1e-8;
state.cl_total_mol_cm3 = 1e-8;
state.water_volume_cm3 = 1;
state.reaction_water_volume_cm3 = 1;
state.interface_area_cm2 = 1;
state.calcite_moles = 1;
state.prescribed_calcite_dissolved_moles = 1e-10;
end
