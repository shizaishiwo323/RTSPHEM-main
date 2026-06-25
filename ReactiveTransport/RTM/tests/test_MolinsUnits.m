function tests = test_MolinsUnits
tests = functiontests(localfunctions);
end

function setupOnce(~)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rtmDir);
end

function testMolinsStrictRateUsesCentimeterSecondUnits(testCase)
params = rtm.units.MolinsStrictParameters();

rate = params.rate_constant_cm_s .* params.h_inlet_mol_cm3;

verifyEqual(testCase, params.units.concentration, "mol/cm^3");
verifyEqual(testCase, params.units.rate_per_area, "mol/cm^2/s");
verifyEqual(testCase, params.units.velocity, "cm/s");
verifyEqual(testCase, rate, 1.255e-7, 'RelTol', 1e-12);
end

function testCalciteMolarVolumeIsPhysicalDefault(testCase)
verifyEqual(testCase, rtm.units.CalciteMolarVolumeCm3Mol(), 36.9, ...
    'RelTol', 1e-12);
end

function testMolinsBenchmarkConfigRecordsUnitConventions(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');

verifyEqual(testCase, cfg.units.concentration, "mol/cm^3");
verifyEqual(testCase, cfg.units.rate_per_area, "mol/cm^2/s");
verifyEqual(testCase, cfg.units.velocity, "cm/s");
verifyEqual(testCase, cfg.units.rate_law, "r_mol_cm2_s = k_cm_s * C_mol_cm3");
end

function testRuntimeManifestUsesConfigUnitConventions(testCase)
cfg = rtm.config.CreateMolinsBenchmarkConfig('partI_strict');
cfg.units.rate_constant = "cm/s";

manifest = rtm.diagnostics.CreateRuntimeManifest(cfg);

verifyEqual(testCase, manifest.units.concentration, "mol/cm^3");
verifyEqual(testCase, manifest.units.rate_per_area, "mol/cm^2/s");
verifyEqual(testCase, manifest.units.rate_constant, "cm/s");
verifyEqual(testCase, manifest.units.rate_law, ...
    "r_mol_cm2_s = k_cm_s * C_mol_cm3");
end

function testLegacyDissolutionPathDoesNotExposeMisleadingDm2RateFields(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
sourceText = string(fileread(fullfile(rtmDir, 'PNM_beauty3.m')));

verifyFalse(testCase, contains(sourceText, 'calciteRate_mol_dm2_s'), ...
    'Legacy PHREEQC fallback data must use calciteRatePerArea_mol_cm2_s for per-area rates.');
verifyFalse(testCase, contains(sourceText, 'rateCoefficientTST_mol_dm2_s'), ...
    'TST rate constants must not be serialized with a mol/dm2/s suffix.');
verifyFalse(testCase, contains(sourceText, 'phreeqcTstRateCoefficient_mol_dm2_s'), ...
    'PHREEQC external TST rate constants must not be serialized with a mol/dm2/s suffix.');
end

function testExternalTstRateConstantCmSDoesNotApplyLiterConversion(testCase)
state = carbonateState();
geometry = singleCellGeometry();
options = struct();
options.h_mol_cm3 = 1.255e-6;
options.h_activity_mol_cm3 = 1.255e-6;
options.rate_constant_cm_s = 0.1;
options.runBatchFunction = @mockRunBatch;

result = rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options);

verifyEqual(testCase, result.candidate_interface_rate_mol_cm2_s, ...
    1.255e-7, 'RelTol', 1e-12);
verifyEqual(testCase, result.candidate_interface_moles, ...
    1.255e-7, 'RelTol', 1e-12);

    function batchResult = mockRunBatch(batchState, ~)
        dissolved = batchState.prescribed_calcite_dissolved_moles(:);
        batchResult = struct();
        batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
            dissolved ./ batchState.water_volume_cm3(:);
        batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
        batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
        batchResult.calciteDissolvedMoles = dissolved;
    end
end

function state = carbonateState()
state = struct();
state.component_names = {'Ca', 'C', 'Na', 'Cl'};
state.component_moles = [0, 0, 1e-8, 1e-8];
state.mineral_names = {'Calcite'};
state.mineral_moles = 1;
state.temperature_C = 25;
state.pressure_atm = 1;
state.time_s = 0;
end

function geometry = singleCellGeometry()
geometry = struct();
geometry.water_volume_cm3 = 1;
geometry.interface_area_cm2 = 1;
geometry.solid_volume_cm3 = 1;
geometry.interface_h_cm = 1;
end
