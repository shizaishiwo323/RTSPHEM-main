function tests = test_RunBenchmarkMolinsPartIIPhreeqcTstConfig
tests = functiontests(localfunctions);
end

function testRunnerUsesExternalTstPhreeqcWithoutChangingLegacyTst(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
runnerPath = fullfile(rtmDir, 'run_benchmark_molins_partII_phreeqc_tst.m');
text = string(fileread(runnerPath));

verifyTrue(testCase, contains(text, ...
    "ConfigurePhreeqcRunGroup(baseCfg, 'external_tst_phreeqc')"), ...
    'The PHREEQC branch should use the explicit external TST + PHREEQC closure mode.');
verifyTrue(testCase, contains(text, "tstCfg.reactionModel = 'tst';"), ...
    'The legacy TST branch must stay on the original TST reaction model.');
verifyFalse(testCase, contains(text, ...
    "ConfigurePhreeqcRunGroup(baseCfg, 'phreeqc_tst_match')"), ...
    'The Part II runner should not use the deprecated PHREEQC TST-match group directly.');
end
