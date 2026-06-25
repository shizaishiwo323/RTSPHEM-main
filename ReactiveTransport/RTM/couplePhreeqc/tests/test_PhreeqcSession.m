function tests = test_PhreeqcSession
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
rtmDir = fileparts(fileparts(testDir));
testCase.TestData.testDir = testDir;
testCase.TestData.rtmDir = rtmDir;
addpath(testDir);
addpath(rtmDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testSessionLoadsDatabaseOnceAndRunsString(testCase)
databasePath = createTempDatabase(testCase);
rawOutput = minimalSelectedOutput();
mockEngine = MockIPhreeqcEngine(rawOutput);

session = rtm.phreeqc.PhreeqcSession(struct( ...
    'engineType', 'mock', ...
    'engineFactory', @() mockEngine));
cleanupObj = onCleanup(@() session.close());

session.loadDatabaseExact(databasePath);
session.loadDatabaseExact(databasePath);
raw = session.runString("SOLUTION 1" + newline + "END");
status = session.getLastStatus();
manifest = session.getDatabaseManifest();

verifyEqual(testCase, raw, rawOutput);
verifyEqual(testCase, mockEngine.LoadDatabaseCount, 1);
verifyEqual(testCase, mockEngine.RunStringCount, 1);
verifyEqual(testCase, status.loadStatus, 0);
verifyEqual(testCase, status.runStatus, 0);
verifyEqual(testCase, manifest.databasePath, string(databasePath));
verifyEqual(testCase, strlength(manifest.databaseSha256), 64);
clear cleanupObj;
end

function testRunBatchUsesInjectedSessionAndRunString(testCase)
databasePath = createTempDatabase(testCase);
rawOutput = minimalSelectedOutput();
mockEngine = MockIPhreeqcEngine(rawOutput);
session = rtm.phreeqc.PhreeqcSession(struct( ...
    'engineType', 'mock', ...
    'engineFactory', @() mockEngine));
cleanupObj = onCleanup(@() session.close());

workDir = tempname;
mkdir(workDir);
state = oneCellAcidCalciteState();
options = struct();
options.databasePath = databasePath;
options.workDir = workDir;
options.timeStepIndex = 77;
options.timeStepSize = 1;
options.rateLaw = 'external_tst_phreeqc';
options.phreeqcSession = session;

result = RunPhreeqcCalciteBatch(state, options);
inputPath = fullfile(workDir, 'phreeqc_calcite_step_0077.phr');

verifyEqual(testCase, mockEngine.LoadDatabaseCount, 1);
verifyEqual(testCase, mockEngine.RunStringCount, 1);
verifyFalse(testCase, exist(inputPath, 'file') == 2, ...
    'Injected PHREEQC sessions should use RunString without writing per-step .phr files.');
verifyEqual(testCase, result.phreeqcRunMethod, "RunString");
verifyTrue(testCase, result.phreeqcSessionReused);
verifyEqual(testCase, result.databasePath, string(databasePath));
verifyEqual(testCase, result.pH, 2.1, 'AbsTol', 1e-12);
verifyGreaterThan(testCase, strlength(result.databaseSha256), 0);
verifyEqual(testCase, result.water_phase_ca_delta_moles, 1e-4, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, result.water_phase_c_delta_moles, 1e-4, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, result.water_phase_ca_delta_moles_total, 1e-4, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, result.water_phase_c_delta_moles_total, 1e-4, ...
    'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, result.calcite_ca_stoich_residual_moles, 0, ...
    'AbsTol', 1e-15);
verifyEqual(testCase, result.calcite_c_stoich_residual_moles, 0, ...
    'AbsTol', 1e-15);
verifyEqual(testCase, result.calcite_stoich_max_abs_residual_moles, 0, ...
    'AbsTol', 1e-15);
clear cleanupObj;
end

function testRunBatchDoesNotInferMineralDeltaFromTotals(testCase)
databasePath = createTempDatabase(testCase);
rawOutput = selectedOutputWithTotalsButNoKinDelta();
mockEngine = MockIPhreeqcEngine(rawOutput);
session = rtm.phreeqc.PhreeqcSession(struct( ...
    'engineType', 'mock', ...
    'engineFactory', @() mockEngine));
cleanupObj = onCleanup(@() session.close());

workDir = tempname;
mkdir(workDir);
testCase.addTeardown(@() cleanupWorkDir(workDir));

state = oneCellAcidCalciteState();
state = rmfield(state, 'prescribed_calcite_dissolved_moles');
options = struct();
options.databasePath = databasePath;
options.workDir = workDir;
options.timeStepIndex = 78;
options.timeStepSize = 1;
options.rateLaw = 'database_calcite';
options.phreeqcSession = session;

result = RunPhreeqcCalciteBatch(state, options);

verifyEqual(testCase, result.calciteDissolvedMoles, 0, 'AbsTol', 0);
verifyEqual(testCase, result.calciteDeltaMoles, 0, 'AbsTol', 0);
verifyEqual(testCase, result.diagnostic_inferred_calcite_dissolved_moles, ...
    1e-4, 'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, result.calcite_ca_stoich_residual_moles, ...
    1e-4, 'RelTol', 1e-12, 'AbsTol', 1e-15);
verifyEqual(testCase, result.calcite_c_stoich_residual_moles, ...
    1e-4, 'RelTol', 1e-12, 'AbsTol', 1e-15);
clear cleanupObj;
end

function testRunBatchChecksRunFileStatusWhenNoSession(testCase)
databasePath = createTempDatabase(testCase);
mockEngine = MockIPhreeqcEngine(minimalSelectedOutput());
mockEngine.RunStatus = 7;
mockEngine.ErrorString = "mock runfile failure";

workDir = tempname;
mkdir(workDir);
testCase.addTeardown(@() cleanupWorkDir(workDir));

state = oneCellAcidCalciteState();
options = struct();
options.databasePath = databasePath;
options.workDir = workDir;
options.timeStepIndex = 78;
options.timeStepSize = 1;
options.rateLaw = 'external_tst_phreeqc';
options.engineFactory = @() mockEngine;
options.writeInputFiles = true;

verifyError(testCase, @() RunPhreeqcCalciteBatch(state, options), ...
    'RTSPHEM:Phreeqc:RunFileFailed');
verifyEqual(testCase, mockEngine.LoadDatabaseCount, 1);
verifyEqual(testCase, mockEngine.RunFileCount, 1);
verifyTrue(testCase, contains(mockEngine.LastInputPath, ...
    'phreeqc_calcite_step_0078.phr'));
end

function databasePath = createTempDatabase(testCase)
databasePath = [tempname, '.dat'];
fid = fopen(databasePath, 'w');
if fid == -1
    error('test_PhreeqcSession:TempFileOpenFailed', 'Cannot create temp database.');
end
fprintf(fid, 'TITLE mock database\n');
fclose(fid);
testCase.addTeardown(@() deleteExistingFile(databasePath));
end

function deleteExistingFile(pathValue)
if exist(pathValue, 'file') == 2
    delete(pathValue);
end
end

function cleanupWorkDir(workDir)
deleteExistingFile(fullfile(workDir, 'phreeqc_calcite_step_0078.phr'));
if exist(workDir, 'dir') == 7
    rmdir(workDir);
end
end

function rawOutput = minimalSelectedOutput()
rawOutput = {
    'sim', 'state', 'soln', 'pH', 'charge(eq)', 'Ca(mol/kgw)', 'C(mol/kgw)', 'Na(mol/kgw)', 'Cl(mol/kgw)', 'm_H+(mol/kgw)', 'm_Ca+2(mol/kgw)', 'm_HCO3-(mol/kgw)', 'm_CO3-2(mol/kgw)', 'm_Cl-(mol/kgw)', 'm_Na+(mol/kgw)', 'si_Calcite', 'KIN_DELTA_Calcite', 'RATE_Calcite';
    1, 'react', 1, 2.1, 0.03, 1e-4, 1e-4, 0.01, 0.11, 0.08, 1e-4, 8e-5, 1e-8, 0.11, 0.01, -3, -3e-6, 3e-5
    };
end

function rawOutput = selectedOutputWithTotalsButNoKinDelta()
rawOutput = {
    'sim', 'state', 'soln', 'pH', 'charge(eq)', 'Ca(mol/kgw)', 'C(mol/kgw)', 'Na(mol/kgw)', 'Cl(mol/kgw)', 'm_H+(mol/kgw)', 'm_Ca+2(mol/kgw)', 'm_HCO3-(mol/kgw)', 'm_CO3-2(mol/kgw)', 'm_Cl-(mol/kgw)', 'm_Na+(mol/kgw)', 'si_Calcite', 'KIN_DELTA_Calcite', 'RATE_Calcite';
    1, 'react', 1, 2.1, 0.03, 1e-4, 1e-4, 0.01, 0.11, 0.08, 1e-4, 8e-5, 1e-8, 0.11, 0.01, -3, 0, 0
    };
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
state.prescribed_calcite_dissolved_moles = 1e-4;
end
