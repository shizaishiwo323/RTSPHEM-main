function tests = test_RunPhreeqcSeBatchConfig
tests = functiontests(localfunctions);
end

function testThreeRegimesAndSliceBudget(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rtmDir);
design = run_phreeqc_se_batch('design');

verifyEqual(testCase, string(design.Mode), ...
    ["face"; "uniform"; "ramified_wormhole"]);
verifyEqual(testCase, design.Pe, ...
    design.InletVelocity_cm_s .* design.CharacteristicLength_cm ./ ...
    design.DiffusionCoefficient_cm2_s, 'AbsTol', 1e-12);
computedDaII = design.AcidConcentration_mol_cm3 .* ...
    design.MolarVolume_cm3_mol .* design.RateCoefficient_mol_dm2_s .* ...
    1000 .* design.CharacteristicLength_cm ./ ...
    design.DiffusionCoefficient_cm2_s;
verifyEqual(testCase, design.DaII, computedDaII, 'AbsTol', 1e-12);
verifyEqual(testCase, design.DaI, design.DaII ./ design.Pe, ...
    'AbsTol', 1e-12);

verifyTrue(testCase, design.Pe(1) <= 0.01 && design.DaI(1) > 1);
verifyTrue(testCase, design.DaI(2) <= 1);
verifyTrue(testCase, design.Pe(3) > 10 && design.DaI(3) > 1);
verifyEqual(testCase, design.SmokeTargetSlices, repmat(10, 3, 1));
verifyGreaterThanOrEqual(testCase, sum(design.ProductionTargetSlices), 100);
verifyEqual(testCase, numel(unique(design.InletVelocity_cm_s)), 3);
verifyEqual(testCase, numel(unique(design.AcidConcentration_mol_L)), 3);
verifyEqual(testCase, numel(unique(design.RateCoefficient_mol_dm2_s)), 3);
initialGridCourant = design.InitialTimeStep_s .* ...
    design.InletVelocity_cm_s ./ design.CharacteristicLength_cm;
verifyLessThanOrEqual(testCase, initialGridCourant, repmat(0.05, 3, 1));
end

function testBatchRequiresEarlyChemistryRows(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
text = string(fileread(fullfile(rtmDir, 'run_phreeqc_se_batch.m')));

verifyTrue(testCase, contains(text, "cfg.phreeqcEarlyExportSteps = earlyExportSteps;"));
verifyTrue(testCase, contains(text, ...
    "hasEarlyChemistryRows(chemistryCsv, earlyExportSteps)"));
verifyTrue(testCase, contains(text, ...
    "all(ismember((1:requiredSteps)', chemistry.TimeStep(:)))"));
end

function testBatchHasExistingDirectoryResumeMode(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
text = string(fileread(fullfile(rtmDir, 'run_phreeqc_se_batch.m')));

verifyTrue(testCase, contains(text, 'stage == "resume"'));
verifyTrue(testCase, contains(text, "Existing batch result directory: "));
verifyTrue(testCase, contains(text, "strip(string(options.batchDir), '""')"));
verifyTrue(testCase, contains(text, "isfolder(options.batchDir)"));
verifyTrue(testCase, contains(text, "SkippedComplete"));
end
