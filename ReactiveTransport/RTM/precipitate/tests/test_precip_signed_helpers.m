function tests = test_precip_signed_helpers
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
moduleDir = fileparts(testDir);
rtmDir = fileparts(moduleDir);
reactiveRoot = fileparts(rtmDir);
levelSetDir = fullfile(reactiveRoot, 'src', 'LevelSetSolver2ndOrder');
testCase.TestData.moduleDir = moduleDir;
testCase.TestData.levelSetDir = levelSetDir;
addpath(moduleDir, '-begin');
addpath(levelSetDir);
end

function teardownOnce(testCase)
rmpath(testCase.TestData.moduleDir);
rmpath(testCase.TestData.levelSetDir);
end

function testConfigureBenchmarkKeepsPnmCompatibilityFields(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4, 'thickness', 0.003));

verifyEqual(testCase, cfg.thicknessCm, 0.003, 'AbsTol', 0);
verifyEqual(testCase, cfg.splitInletY, 0.2, 'AbsTol', 0);
verifyEqual(testCase, cfg.targetLengthYAxis, cfg.lengthYAxis, 'AbsTol', 0);
verifyEqual(testCase, cfg.targetAspectRatio, cfg.lengthXAxis / cfg.lengthYAxis, 'AbsTol', 0);
verifyEqual(testCase, cfg.inletCalciumConcentration, 12.5e-6, 'AbsTol', 1e-15);
verifyEqual(testCase, cfg.inletCarbonConcentration, 12.5e-6, 'AbsTol', 1e-15);
verifyEqual(testCase, cfg.mineralEvolutionMode, 'signed_calcite_surface');
end

function testConfigureBenchmarkPreservesExplicitSplitInletRanges(testCase)
overrides = struct();
overrides.lengthYAxis = 0.4;
overrides.inletA = struct('yRange', [0.05, 0.18]);
overrides.inletB = struct('yRange', [0.18, 0.35]);

cfg = precip_ConfigureZhangYoonBenchmark(overrides);

verifyEqual(testCase, cfg.inletA.yRange, [0.05, 0.18], 'AbsTol', 0);
verifyEqual(testCase, cfg.inletB.yRange, [0.18, 0.35], 'AbsTol', 0);
verifyEqual(testCase, cfg.inletA.Ca_total, 25e-6, 'AbsTol', 1e-15);
verifyEqual(testCase, cfg.inletB.C_total, 25e-6, 'AbsTol', 1e-15);
end

function testComputePhreeqcTransportUpperBoundsUsesSplitInletMaxima(testCase)
cfg = precip_ConfigureZhangYoonBenchmark();

[upperBounds, speciesNames] = precip_ComputePhreeqcTransportUpperBounds(cfg, 2);

verifyEqual(testCase, speciesNames, {'H', 'Ca_total', 'C_total', 'Na_total', 'Cl_total'});
verifyEqual(testCase, upperBounds(1), max([cfg.initialHydrogenConcentration, ...
    cfg.inletA.H_total, cfg.inletB.H_total]) * 2, 'AbsTol', 1e-18);
verifyEqual(testCase, upperBounds(2), 25e-6 * 2, 'AbsTol', 1e-15);
verifyEqual(testCase, upperBounds(3), 25e-6 * 2, 'AbsTol', 1e-15);
verifyEqual(testCase, upperBounds(4), 50e-6 * 2, 'AbsTol', 1e-15);
verifyEqual(testCase, upperBounds(5), 50e-6 * 2, 'AbsTol', 1e-15);
end

function testComputePhreeqcTransportUpperBoundsAllowsInfiniteLimiter(testCase)
cfg = precip_ConfigureZhangYoonBenchmark();

upperBounds = precip_ComputePhreeqcTransportUpperBounds(cfg, Inf);

verifyEqual(testCase, upperBounds, Inf(5, 1));
end

function testComputeBoundaryInletMassFluxUsesUniformFallback(testCase)
cfg = struct();

massFlux = precip_ComputeBoundaryInletMassFlux('uniform_left_inlet', cfg, ...
    'H_total', 3e-6, 0.25, 0.01, 0.4);

verifyEqual(testCase, massFlux, 3e-6 * 0.25 * 0.4 * 0.01, 'AbsTol', 1e-18);
end

function testComputeBoundaryInletMassFluxIntegratesSplitInletRanges(testCase)
cfg = struct();
cfg.inletA = struct('H_total', 2e-6, 'yRange', [0.00, 0.15]);
cfg.inletB = struct('H_total', 8e-6, 'yRange', [0.15, 0.40]);

massFlux = precip_ComputeBoundaryInletMassFlux('split_left_inlet', cfg, ...
    'H_total', 99e-6, 0.25, 0.01, 0.4);

expected = (2e-6 * 0.15 + 8e-6 * 0.25) * 0.25 * 0.01;
verifyEqual(testCase, massFlux, expected, 'AbsTol', 1e-18);
end

function testComputeBoundaryInletMassFluxRejectsInvalidSplitConfig(testCase)
cfg = struct();
cfg.inletA = struct('H_total', 2e-6, 'yRange', [0.00, 0.20]);
cfg.inletB = struct('H_total', 8e-6, 'yRange', [0.25, 0.40]);

verifyError(testCase, @() precip_ComputeBoundaryInletMassFlux('split_left_inlet', ...
    cfg, 'H_total', 1e-6, 0.25, 0.01, 0.4), ...
    'RTSPHEM:Precipitate:InvalidSplitInletRange');
end

function testComputeBoundaryInletConcentrationScaleUsesUniformFallback(testCase)
cfg = struct();

scale = precip_ComputeBoundaryInletConcentrationScale('uniform_left_inlet', ...
    cfg, 'H_total', 3e-6);

verifyEqual(testCase, scale, 3e-6, 'AbsTol', 0);
end

function testComputeBoundaryInletConcentrationScaleUsesSplitMaximum(testCase)
cfg = struct();
cfg.inletA = struct('H_total', 2e-6, 'yRange', [0.00, 0.15]);
cfg.inletB = struct('H_total', 8e-6, 'yRange', [0.15, 0.40]);

scale = precip_ComputeBoundaryInletConcentrationScale('split_left_inlet', ...
    cfg, 'H_total', 99e-6);

verifyEqual(testCase, scale, 8e-6, 'AbsTol', 0);
end

function testPrecipPnmConcentrationLimiterUsesInletConcentrationScale(testCase)
solverPath = fullfile(testCase.TestData.moduleDir, 'precip_PNM_beauty3.m');
solverText = string(fileread(solverPath));
callMatch = regexp(solverText, ...
    'ResolveConcentrationLimitedStep\([^\)]*\)', 'match', 'once');

verifyNotEmpty(testCase, callMatch);
verifyTrue(testCase, contains(callMatch, 'diagnosticInletConcentrationScale'));
verifyFalse(testCase, contains(callMatch, 'initialHydrogenConcentration'));
end

function testResolveAdaptiveTimeGridExtensionHonorsExplicitOverride(testCase)
cfg = struct('allowAdaptiveTimeGridExtension', true);

tf = precip_ResolveAdaptiveTimeGridExtension(cfg, false);

verifyTrue(testCase, tf);
end

function testResolveAdaptiveTimeGridExtensionUsesDefaultWhenUnset(testCase)
cfg = struct();

verifyFalse(testCase, precip_ResolveAdaptiveTimeGridExtension(cfg, false));
verifyTrue(testCase, precip_ResolveAdaptiveTimeGridExtension(cfg, true));
end

function testResolveAdaptiveTimeGridExtensionParsesFalseString(testCase)
cfg = struct('allowAdaptiveTimeGridExtension', "false");

tf = precip_ResolveAdaptiveTimeGridExtension(cfg, true);

verifyFalse(testCase, tf);
end

function testResolveAdaptiveTimeGridExtensionRejectsNonScalar(testCase)
cfg = struct('allowAdaptiveTimeGridExtension', [true, false]);

verifyError(testCase, @() precip_ResolveAdaptiveTimeGridExtension(cfg, false), ...
    'RTSPHEM:Precipitate:InvalidAdaptiveTimeGridExtension');
end

function testResolveAdaptiveTimeGridExtensionRejectsAmbiguousNumericAndText(testCase)
cfgNumeric = struct('allowAdaptiveTimeGridExtension', 2);
cfgText = struct('allowAdaptiveTimeGridExtension', "yes");

verifyError(testCase, @() precip_ResolveAdaptiveTimeGridExtension(cfgNumeric, false), ...
    'RTSPHEM:Precipitate:InvalidAdaptiveTimeGridExtension');
verifyError(testCase, @() precip_ResolveAdaptiveTimeGridExtension(cfgText, false), ...
    'RTSPHEM:Precipitate:InvalidAdaptiveTimeGridExtension');
end

function testSignedBatchUsesPrecipLocalPhreeqcInputBuilder(testCase)
runnerPath = fullfile(testCase.TestData.moduleDir, ...
    'precip_RunPhreeqcCalciteBatchSigned.m');
runnerText = string(fileread(runnerPath));

localBuilderPath = which('precip_BuildCalcitePhreeqcInputSigned');
verifyEqual(testCase, string(localBuilderPath), ...
    string(fullfile(testCase.TestData.moduleDir, 'precip_BuildCalcitePhreeqcInputSigned.m')));
verifyNotEmpty(testCase, regexp(runnerText, ...
    '\<precip_BuildCalcitePhreeqcInputSigned\s*\(', 'once'));
verifyEmpty(testCase, regexp(runnerText, ...
    '(?<!precip_)\<BuildCalcitePhreeqcInput\s*(\(|;|$)', 'once'));
end

function testPrecipBuildCalcitePhreeqcInputSignedWritesSignedKinetics(testCase)
state = makeMinimalPhreeqcBuilderState();
options = struct('timeStepSize', 2, 'rateLaw', 'database_calcite');

inputText = string(precip_BuildCalcitePhreeqcInputSigned(state, options));

verifyTrue(testCase, contains(inputText, 'SOLUTION 1'));
verifyTrue(testCase, contains(inputText, 'KINETICS 1'));
verifyTrue(testCase, contains(inputText, 'USER_PUNCH'));
verifyTrue(testCase, contains(inputText, 'KIN_DELTA("Calcite")'));
verifyTrue(testCase, contains(inputText, '-time_step 2'));
end

function testPrecipBuildCalcitePhreeqcInputSignedUsesPrecipitateErrorIds(testCase)
state = makeMinimalPhreeqcBuilderState();
state = rmfield(state, 'h_mol_cm3');

verifyError(testCase, @() precip_BuildCalcitePhreeqcInputSigned(state, struct()), ...
    'RTSPHEM:Precipitate:MissingStateField');
end

function testPrecipPnmUsesLocalPhreeqcOutputHelpers(testCase)
solverPath = fullfile(testCase.TestData.moduleDir, 'precip_PNM_beauty3.m');
solverText = string(fileread(solverPath));

verifyNotEmpty(testCase, regexp(solverText, ...
    '\<precip_WritePhreeqcSpeciesTable\s*\(', 'once'));
verifyNotEmpty(testCase, regexp(solverText, ...
    '\<precip_ExportPhreeqcSpeciesPlots\s*\(', 'once'));
verifyEmpty(testCase, regexp(solverText, ...
    '(?<!precip_)\<WritePhreeqcSpeciesTable\s*\(', 'once'));
verifyEmpty(testCase, regexp(solverText, ...
    '(?<!precip_)\<ExportPhreeqcSpeciesPlots\s*\(', 'once'));
end

function testPrecipWritePhreeqcSpeciesTableKeepsSignedCalciteColumns(testCase)
workDir = tempname;
mkdir(workDir);
grid = struct();
grid.numT = 2;
grid.baryT = [0.1, 0.2; 0.3, 0.4];
speciesData = struct();
speciesData.h_mol_cm3 = [1e-10; 2e-10];
speciesData.calciteDeltaMoles = [2e-9; -3e-9];
speciesData.calcitePrecipitatedMoles = [2e-9; 0];
speciesData.calciteDissolvedMoles = [0; 3e-9];
speciesData.calciteSignedRate_mol_s = [1e-9; -1.5e-9];
speciesData.calciteRatePerArea_mol_cm2_s = [4e-12; -5e-12];

precip_WritePhreeqcSpeciesTable(workDir, 7, 12.5, grid, speciesData);

tableData = readtable(fullfile(workDir, 'phreeqc_species_0007.csv'));
verifyTrue(testCase, ismember('calciteDeltaMoles', tableData.Properties.VariableNames));
verifyTrue(testCase, ismember('calcitePrecipitatedMoles', tableData.Properties.VariableNames));
verifyTrue(testCase, ismember('calciteSignedRate_mol_s', tableData.Properties.VariableNames));
verifyEqual(testCase, tableData.calciteDeltaMoles, speciesData.calciteDeltaMoles, 'AbsTol', 0);
verifyEqual(testCase, tableData.calciteSignedRate_mol_s, ...
    speciesData.calciteSignedRate_mol_s, 'AbsTol', 0);
end

function testPrecipPrepareConcentrationFaceDataMasksMixedAndSolidTriangles(testCase)
concentration = [1; 2; 3];
triangles = [1, 2, 3; 3, 4, 5; 5, 6, 7];
levelSetData = [-1; -1; -1; -1; 1; 1; 1];

faceData = precip_PrepareConcentrationFaceData(concentration, triangles, levelSetData);

verifyEqual(testCase, faceData(1), 1, 'AbsTol', 0);
verifyTrue(testCase, isnan(faceData(2)));
verifyTrue(testCase, isnan(faceData(3)));
verifyEqual(testCase, precip_PrepareConcentrationFaceData(concentration, triangles, []), ...
    concentration, 'AbsTol', 0);
end

function testExtendTransportProblemsToStepperVisitsAllSpeciesAndVariables(testCase)
fieldNames = {'Q', 'U', 'A', 'B', 'C', 'D', 'E', 'F', ...
    'uD', 'gF', 'uR', 'rob', 'L'};
speciesNames = {'H', 'Ca', 'C', 'Na', 'Cl'};
transportProblems = cell(numel(speciesNames), 1);
for iSpecies = 1:numel(speciesNames)
    transportProblems{iSpecies} = makeTransportFixture(speciesNames{iSpecies}, fieldNames);
end
visited = cell(numel(speciesNames) * numel(fieldNames), 1);
visitCount = 0;
extensionFcn = @(value, mode) recordExtension(value, mode);

summary = precip_ExtendTransportProblemsToStepper(transportProblems, extensionFcn);

verifyEqual(testCase, summary.numProblems, 5);
verifyEqual(testCase, summary.numVariablesExtended, 65);
verifyEqual(testCase, visitCount, 65);
verifyEqual(testCase, visited{1}, 'H_Q');
verifyEqual(testCase, visited{13}, 'H_L');
verifyEqual(testCase, visited{14}, 'Ca_Q');
verifyEqual(testCase, visited{65}, 'Cl_L');

    function recordExtension(value, mode)
        verifyEqual(testCase, mode, 'copy_previous');
        visitCount = visitCount + 1;
        visited{visitCount} = value;
    end
end

function fixture = makeTransportFixture(speciesName, fieldNames)
fixture = struct();
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    fixture.(fieldName) = sprintf('%s_%s', speciesName, fieldName);
end
end

function state = makeMinimalPhreeqcBuilderState()
state = struct();
state.h_mol_cm3 = 1e-10;
state.ca_total_mol_cm3 = 25e-6;
state.c_total_mol_cm3 = 25e-6;
state.na_total_mol_cm3 = 50e-6;
state.cl_total_mol_cm3 = 50e-6;
state.interface_area_cm2 = 1e-4;
state.water_volume_cm3 = 1e-6;
state.calcite_moles = 1e-9;
end

function testSignedParserKeepsPrecipitationPositiveDelta(testCase)
raw = {
    'soln', 'pH', 'Ca(mol/kgw)', 'C(mol/kgw)', 'm_H+(mol/kgw)', 'm_Ca+2(mol/kgw)', 'si_Calcite', 'KIN_DELTA_Calcite', 'RATE_Calcite';
    1, 8.4, 0.010, 0.010, 1e-8, 0.009, 1.2, 2e-6, -2e-6;
    2, 6.5, 0.020, 0.020, 1e-6, 0.018, -0.5, -3e-6, 3e-6
    };

result = precip_ParsePhreeqcSelectedOutputSigned(raw, 2, 10);

verifyEqual(testCase, result.ca_total_mol_cm3, [1e-5; 2e-5], 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteDeltaMoles, [2e-6; -3e-6], 'AbsTol', 1e-15);
verifyEqual(testCase, result.calcitePrecipitatedMoles, [2e-6; 0], 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteDissolvedMoles, [0; 3e-6], 'AbsTol', 1e-15);
verifyEqual(testCase, result.calciteSignedRate_mol_s, [2e-7; -3e-7], 'AbsTol', 1e-15);
end

function testSignedParserRequiresCalciteDeltaHeadings(testCase)
raw = {
    'soln', 'pH', 'Ca(mol/kgw)', 'C(mol/kgw)';
    1, 8.4, 0.010, 0.010
    };

verifyError(testCase, @() precip_ParsePhreeqcSelectedOutputSigned(raw, 1, 10), ...
    'RTSPHEM:Precipitate:MissingSelectedOutputHeading');
end

function testSignedParserRequiresCalciteRateHeading(testCase)
raw = {
    'soln', 'pH', 'Ca(mol/kgw)', 'C(mol/kgw)', 'KIN_DELTA_Calcite';
    1, 8.4, 0.010, 0.010, 2e-6
    };

verifyError(testCase, @() precip_ParsePhreeqcSelectedOutputSigned(raw, 1, 10), ...
    'RTSPHEM:Precipitate:MissingSelectedOutputHeading');
end

function testScaleSignedDeltaDoesNotCapPrecipitationByDefault(testCase)
result = struct();
result.calciteDeltaMoles = [2; -3];
result.ca_total_mol_cm3 = [1e-5; 1e-5];
result.c_total_mol_cm3 = [1e-5; 1e-5];

state = struct();
state.water_volume_cm3 = [1e-3; 1e-3];
state.calcite_moles = [Inf; 1e-9];
state.ca_total_mol_cm3 = [1e-6; 1e-6];
state.c_total_mol_cm3 = [1e-6; 1e-6];

updated = precip_ScaleSignedCalciteDeltaToCellInventory(result, state, struct('timeStepSize', 2));

verifyEqual(testCase, updated.calciteDeltaMoles(1), 2e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.calciteDeltaMoles(2), -1e-9, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.calciteSignedRate_mol_s(1), 1e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.calciteRate_mol_s(2), 0.5e-9, 'AbsTol', 1e-18);
end

function testOptionalPrecipitationCapBlendsAqueousFields(testCase)
result = struct();
result.calciteDeltaMoles = 2;
result.ca_total_mol_cm3 = 1e-5;
result.c_total_mol_cm3 = 1e-5;

state = struct();
state.water_volume_cm3 = 1e-3;
state.calcite_moles = Inf;
state.ca_total_mol_cm3 = 1e-6;
state.c_total_mol_cm3 = 1e-6;

updated = precip_ScaleSignedCalciteDeltaToCellInventory(result, state, ...
    struct('timeStepSize', 1, 'capPrecipitationByInitialAqueousInventory', true));

verifyEqual(testCase, updated.calciteDeltaMoles, 1e-9, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.ca_total_mol_cm3, 1.0045e-6, 'AbsTol', 1e-12);
end

function testLimitedDeltaBlendsTotalsFromLegacyStateFields(testCase)
result = struct();
result.calciteDeltaMoles = -2;
result.ca_total_mol_cm3 = 1.0e-5;
result.c_total_mol_cm3 = 9.0e-6;
result.na_total_mol_cm3 = 5.0e-6;
result.cl_total_mol_cm3 = 6.0e-6;

state = struct();
state.water_volume_cm3 = 1e-3;
state.calcite_moles = 1e-9;
state.ca_mol_cm3 = 2.0e-6;
state.c_mol_cm3 = 3.0e-6;
state.na_mol_cm3 = 4.0e-6;
state.cl_mol_cm3 = 5.0e-6;

updated = precip_ScaleSignedCalciteDeltaToCellInventory(result, state, struct('timeStepSize', 1));

verifyEqual(testCase, updated.calciteDeltaMoles, -1e-9, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.ca_total_mol_cm3, 2.004e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.c_total_mol_cm3, 3.003e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.na_total_mol_cm3, 4.0005e-6, 'AbsTol', 1e-18);
verifyEqual(testCase, updated.cl_total_mol_cm3, 5.0005e-6, 'AbsTol', 1e-18);
end

function testPrescribedDissolutionIsRestoredForSignedTstMatch(testCase)
result = struct();
result.calciteDeltaMoles = [0; 0];

state = struct();
state.calcite_moles = [1e-8; 2e-8];
state.prescribed_calcite_dissolved_moles = [3e-9; 4e-8];

updated = precip_ApplyPrescribedCalciteDissolutionSigned(result, state, ...
    struct('timeStepSize', 2));

verifyEqual(testCase, updated.calciteDeltaMoles, [-3e-9; -2e-8], 'AbsTol', 1e-18);
verifyEqual(testCase, updated.calciteDissolvedMoles, [3e-9; 2e-8], 'AbsTol', 1e-18);
verifyEqual(testCase, updated.calciteSignedRate_mol_s, [-1.5e-9; -1e-8], 'AbsTol', 1e-18);
verifyEqual(testCase, updated.calciteRate_mol_s, [1.5e-9; 1e-8], 'AbsTol', 1e-18);
end

function testSignedInterfaceRateMakesPrecipitationNegative(testCase)
result = struct();
result.calciteDeltaMoles = [2e-9; -3e-9; 1e-9];
rate = precip_ComputeSignedCalciteInterfaceRatePerArea(result, [1e-3; 2e-3; 0], 10);

verifyEqual(testCase, rate, [-2e-7; 1.5e-7; 0], 'AbsTol', 1e-18);
end

function testSignedLevelSetSmokeMakesNegativeSpeedGrowSolid(testCase)
smoke = precip_RunLevelSetSignedVelocitySmoke();
expectedInitialArea = pi * smoke.initialRadius_cm^2;
expectedPrecipitationArea = pi * smoke.precipitationExpectedRadius_cm^2;
expectedDissolutionArea = pi * smoke.dissolutionExpectedRadius_cm^2;

verifyGreaterThan(testCase, smoke.precipitationSolidAreaAfter_cm2, ...
    smoke.initialSolidArea_cm2);
verifyLessThan(testCase, smoke.dissolutionSolidAreaAfter_cm2, ...
    smoke.initialSolidArea_cm2);
verifyLessThan(testCase, smoke.precipitationNormalSpeed_cm_s, 0);
verifyGreaterThan(testCase, smoke.dissolutionNormalSpeed_cm_s, 0);
verifyEqual(testCase, smoke.levelSetStepper, 'levelSetEquationTimeStep');
verifyEqual(testCase, smoke.initialSolidArea_cm2, expectedInitialArea, 'AbsTol', 5e-3);
verifyEqual(testCase, smoke.precipitationSolidAreaAfter_cm2, expectedPrecipitationArea, 'AbsTol', 5e-3);
verifyEqual(testCase, smoke.dissolutionSolidAreaAfter_cm2, expectedDissolutionArea, 'AbsTol', 5e-3);
end

function testSignedLevelSetSmokeRejectsInvalidOptions(testCase)
verifyError(testCase, @() precip_RunLevelSetSignedVelocitySmoke( ...
    struct('numSamples', 20)), 'RTSPHEM:Precipitate:InvalidSmokeGrid');
verifyError(testCase, @() precip_RunLevelSetSignedVelocitySmoke( ...
    struct('initialRadiusCm', 0.5)), 'RTSPHEM:Precipitate:InvalidSmokeRadius');
verifyError(testCase, @() precip_RunLevelSetSignedVelocitySmoke( ...
    struct('timeStepSize', 0)), 'RTSPHEM:Precipitate:InvalidSmokeStep');
verifyError(testCase, @() precip_RunLevelSetSignedVelocitySmoke( ...
    struct('speedMagnitudeCmS', 0)), 'RTSPHEM:Precipitate:InvalidSmokeStep');
end

function testSplitInletFluxSeparatesCalciumAndCarbonate(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4));
inletVelocity = 0.02;
epsValue = 1e-9;
points = [0, 0, 0.01; 0.01, 0.20, 0.20];

caFlux = precip_CreateTransportMultiInlet('split_left_inlet', cfg, ...
    'Ca_total', 0, inletVelocity, epsValue);
cFlux = precip_CreateTransportMultiInlet('split_left_inlet', cfg, ...
    'C_total', 0, inletVelocity, epsValue);

verifyEqual(testCase, caFlux(0, points), ...
    [-inletVelocity * cfg.inletA.Ca_total, 0, 0], 'AbsTol', 1e-18);
verifyEqual(testCase, cFlux(0, points), ...
    [0, -inletVelocity * cfg.inletB.C_total, 0], 'AbsTol', 1e-18);
end

function testSplitInletFluxCoversHydrogenSodiumAndChloride(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4));
inletVelocity = 0.02;
epsValue = 1e-9;
points = [0, 0; 0.01, 0.20];

hFlux = precip_CreateTransportMultiInlet('split_left_inlet', cfg, ...
    'H_total', 0, inletVelocity, epsValue);
naFlux = precip_CreateTransportMultiInlet('split_left_inlet', cfg, ...
    'Na_total', 0, inletVelocity, epsValue);
clFlux = precip_CreateTransportMultiInlet('split_left_inlet', cfg, ...
    'Cl_total', 0, inletVelocity, epsValue);

verifyEqual(testCase, hFlux(0, points), ...
    -inletVelocity * [cfg.inletA.H_total, cfg.inletB.H_total], 'AbsTol', 1e-18);
verifyEqual(testCase, naFlux(0, points), ...
    [0, -inletVelocity * cfg.inletB.Na_total], 'AbsTol', 1e-18);
verifyEqual(testCase, clFlux(0, points), ...
    [-inletVelocity * cfg.inletA.Cl_total, 0], 'AbsTol', 1e-18);
end

function testUniformInletFluxKeepsLegacyScalarBehaviorAndNx2Coordinates(testCase)
flux = precip_CreateTransportMultiInlet('uniform_left_inlet', struct(), ...
    '', 3e-6, 0.02, 1e-9);
points = [0, 0.01; 0, 0.20; 0.01, 0.20];

verifyEqual(testCase, flux(0, points), [-6e-8, -6e-8, 0], 'AbsTol', 1e-18);
end

function testSplitInletRejectsMissingConfigOrSpecies(testCase)
cfg = precip_ConfigureZhangYoonBenchmark();

verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    struct(), 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:MissingSplitInletConfig');

cfg.inletA = rmfield(cfg.inletA, 'Ca_total');
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfg, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:MissingSplitInletSpecies');
end

function testSplitInletRejectsInvalidConcentrations(testCase)
cfg = precip_ConfigureZhangYoonBenchmark();

cfgEmpty = cfg;
cfgEmpty.inletA.Ca_total = [];
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfgEmpty, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:InvalidSplitInletConcentration');

cfgNan = cfg;
cfgNan.inletA.Ca_total = NaN;
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfgNan, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:InvalidSplitInletConcentration');

cfgVector = cfg;
cfgVector.inletA.Ca_total = [1e-6, 2e-6];
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfgVector, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:InvalidSplitInletConcentration');
end

function testSplitInletRejectsInvalidRanges(testCase)
cfg = precip_ConfigureZhangYoonBenchmark();

cfgEmpty = cfg;
cfgEmpty.inletA.yRange = [];
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfgEmpty, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:InvalidSplitInletRange');

cfgNan = cfg;
cfgNan.inletA.yRange = [0, NaN];
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfgNan, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:InvalidSplitInletRange');

cfgZeroWidth = cfg;
cfgZeroWidth.inletA.yRange = [0.1, 0.1];
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfgZeroWidth, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:InvalidSplitInletRange');

cfgOverlap = cfg;
cfgOverlap.inletA.yRange = [0, 0.12];
cfgOverlap.inletB.yRange = [0.10, cfg.lengthYAxis];
verifyError(testCase, @() precip_CreateTransportMultiInlet('split_left_inlet', ...
    cfgOverlap, 'Ca_total', 0, 0.02, 1e-9), ...
    'RTSPHEM:Precipitate:InvalidSplitInletRange');
end

function testSplitInletMeshAlignmentRejectsMissingLeftBoundaryEndpoint(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4));
leftBoundaryNodes = [ ...
    0, 0.0; ...
    0, 0.1; ...
    0, 0.3; ...
    0, 0.4];

verifyError(testCase, @() precip_ValidateSplitInletMeshAlignment( ...
    leftBoundaryNodes, cfg, 1e-12), ...
    'RTSPHEM:Precipitate:SplitInletNotMeshAligned');
end

function testSplitInletMeshAlignmentAcceptsRequiredEndpoints(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4));
leftBoundaryNodes = [ ...
    0, 0.0; ...
    0, 0.2; ...
    0, 0.4];

verifyWarningFree(testCase, @() precip_ValidateSplitInletMeshAlignment( ...
    leftBoundaryNodes, cfg, 1e-12));
end

function testSplitInletMeshAlignmentAcceptsCoordVGridStruct(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4));
grid = struct('coordV', [ ...
    0, 0.0; ...
    0, 0.2; ...
    0, 0.4; ...
    0.1, 0.2]);

verifyWarningFree(testCase, @() precip_ValidateSplitInletMeshAlignment( ...
    grid, cfg, 1e-12));
end

function testSplitInletMeshAlignmentRejectsCoordVGridStructMissingEndpoint(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4));
grid = struct('coordV', [ ...
    0, 0.0; ...
    0, 0.1; ...
    0, 0.4]);

verifyError(testCase, @() precip_ValidateSplitInletMeshAlignment( ...
    grid, cfg, 1e-12), ...
    'RTSPHEM:Precipitate:SplitInletNotMeshAligned');
end

function testSplitInletMeshAlignmentAcceptsTwoByNCoordinates(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4));
leftBoundaryNodes = [ ...
    0, 0, 0; ...
    0.0, 0.2, 0.4];

verifyWarningFree(testCase, @() precip_ValidateSplitInletMeshAlignment( ...
    leftBoundaryNodes, cfg, 1e-12));
end

function testPrecipitationAreaMetricsComputesTotalAndWindowNetArea(testCase)
grid = localUnitSquareGrid();
initialLevels = [-1; -1; -1; -1];
currentLevels = [1; 1; 1; 1];
config = struct();
config.benchmarkFirstPoreXMaxCm = 0.5;
config.benchmarkFirstThreePoresXMaxCm = 1.0;

metrics = precip_ComputePrecipitationAreaMetrics(grid, currentLevels, initialLevels, config);

verifyEqual(testCase, metrics.totalSolidArea_cm2, 1.0, 'AbsTol', 1e-12);
verifyEqual(testCase, metrics.totalNetSolidArea_cm2, 1.0, 'AbsTol', 1e-12);
verifyEqual(testCase, metrics.firstPoreSolidArea_cm2, 0.5, 'AbsTol', 1e-12);
verifyEqual(testCase, metrics.firstPoreNetSolidArea_cm2, 0.5, 'AbsTol', 1e-12);
verifyEqual(testCase, metrics.firstThreePoresSolidArea_cm2, 1.0, 'AbsTol', 1e-12);
verifyEqual(testCase, metrics.firstThreePoresNetSolidArea_cm2, 1.0, 'AbsTol', 1e-12);
end

function testPrecipitationAreaMetricsUsesZhangWindowDefaults(testCase)
grid = localUnitSquareGrid();
levels = [1; -1; 1; -1];

metrics = precip_ComputePrecipitationAreaMetrics(grid, levels, levels, ...
    struct('layoutType', 'zhang2010_micromodel_local'));

verifyEqual(testCase, metrics.totalNetSolidArea_cm2, 0, 'AbsTol', 1e-12);
verifyTrue(testCase, isfinite(metrics.firstPoreSolidArea_cm2));
verifyTrue(testCase, isfinite(metrics.firstThreePoresSolidArea_cm2));
verifyGreaterThanOrEqual(testCase, metrics.firstThreePoresSolidArea_cm2, ...
    metrics.firstPoreSolidArea_cm2);
end

function testPrecipitationAreaMetricsComputesClippedTriangleArea(testCase)
grid = localSingleTriangleGrid();
initialLevels = [-1; -1; -1];
currentLevels = [-1; 1; 1];
config = struct('benchmarkFirstPoreXMaxCm', 1, 'benchmarkFirstThreePoresXMaxCm', 1);

metrics = precip_ComputePrecipitationAreaMetrics(grid, currentLevels, initialLevels, config);

verifyEqual(testCase, metrics.totalSolidArea_cm2, 0.375, 'AbsTol', 1e-12);
verifyEqual(testCase, metrics.totalNetSolidArea_cm2, 0.375, 'AbsTol', 1e-12);
end

function testPrecipitationAreaMetricsRejectsMalformedGrid(testCase)
grid = localSingleTriangleGrid();
levels = [-1; -1; -1];

badCoord = grid;
badCoord.coordV = [0; 1; 2];
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    badCoord, levels, levels, struct()), 'RTSPHEM:Precipitate:InvalidAreaGrid');

badV0T = grid;
badV0T.V0T = [1, 2];
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    badV0T, levels, levels, struct()), 'RTSPHEM:Precipitate:InvalidAreaGrid');

badIndex = grid;
badIndex.V0T = [1, 2, 4];
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    badIndex, levels, levels, struct()), 'RTSPHEM:Precipitate:InvalidAreaGrid');

badArea = grid;
badArea.areaT = -0.5;
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    badArea, levels, levels, struct()), 'RTSPHEM:Precipitate:InvalidAreaGrid');

badBary = grid;
badBary.baryT = 0.25;
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    badBary, levels, levels, struct()), 'RTSPHEM:Precipitate:InvalidAreaGrid');
end

function testPrecipitationAreaMetricsRejectsInvalidLevelSets(testCase)
grid = localSingleTriangleGrid();
levels = [-1; -1; -1];

nanLevels = levels;
nanLevels(2) = NaN;
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    grid, nanLevels, levels, struct()), 'RTSPHEM:Precipitate:InvalidAreaLevels');

infLevels = levels;
infLevels(2) = Inf;
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    grid, levels, infLevels, struct()), 'RTSPHEM:Precipitate:InvalidAreaLevels');

complexLevels = levels;
complexLevels(2) = 1i;
verifyError(testCase, @() precip_ComputePrecipitationAreaMetrics( ...
    grid, complexLevels, levels, struct()), 'RTSPHEM:Precipitate:InvalidAreaLevels');
end

function testPrecipitationAreaTimeseriesWritesRequiredColumns(testCase)
workDir = tempname;
mkdir(workDir);
csvFile = fullfile(workDir, 'precipitation_area_timeseries.csv');
metrics = localAreaMetricsStruct();

precip_InitializePrecipitationAreaTimeseries(csvFile);
precip_AppendPrecipitationAreaTimeseries(csvFile, 3, 12.5, metrics);

tableData = readtable(csvFile);
expectedNames = {'timestep', 'time_s', ...
    'total_net_solid_area_cm2', 'first_pore_net_solid_area_cm2', ...
    'first_three_pores_net_solid_area_cm2', 'total_solid_area_cm2', ...
    'first_pore_solid_area_cm2', 'first_three_pores_solid_area_cm2'};

verifyEqual(testCase, tableData.Properties.VariableNames, expectedNames);
verifyEqual(testCase, tableData.timestep, 3);
verifyEqual(testCase, tableData.time_s, 12.5, 'AbsTol', 1e-12);
verifyEqual(testCase, tableData.total_net_solid_area_cm2, ...
    metrics.totalNetSolidArea_cm2, 'AbsTol', 1e-12);
verifyEqual(testCase, tableData.first_three_pores_solid_area_cm2, ...
    metrics.firstThreePoresSolidArea_cm2, 'AbsTol', 1e-12);
end

function testPrecipitationAreaTimeseriesPlotCreatesPng(testCase)
workDir = tempname;
mkdir(workDir);
csvFile = fullfile(workDir, 'precipitation_area_timeseries.csv');
pngFile = fullfile(workDir, 'precipitation_area_timeseries.png');

precip_InitializePrecipitationAreaTimeseries(csvFile);
metrics = localAreaMetricsStruct();
precip_AppendPrecipitationAreaTimeseries(csvFile, 1, 60, metrics);
metrics.totalNetSolidArea_cm2 = metrics.totalNetSolidArea_cm2 * 2;
metrics.firstPoreNetSolidArea_cm2 = metrics.firstPoreNetSolidArea_cm2 * 2;
metrics.firstThreePoresNetSolidArea_cm2 = metrics.firstThreePoresNetSolidArea_cm2 * 2;
precip_AppendPrecipitationAreaTimeseries(csvFile, 2, 120, metrics);

precip_PlotPrecipitationAreaTimeseries(csvFile, pngFile, [60, 120]);

verifyTrue(testCase, isfile(pngFile));
verifyGreaterThan(testCase, dir(pngFile).bytes, 0);
end

function testPrecipitationAreaTimeseriesRejectsNonIntegerTimestep(testCase)
workDir = tempname;
mkdir(workDir);
csvFile = fullfile(workDir, 'precipitation_area_timeseries.csv');
precip_InitializePrecipitationAreaTimeseries(csvFile);

verifyError(testCase, @() precip_AppendPrecipitationAreaTimeseries( ...
    csvFile, 1.5, 60, localAreaMetricsStruct()), ...
    'RTSPHEM:Precipitate:InvalidAreaTimeseriesStep');
end

function testBenchmarkSnapshotsSkipUntilDue(testCase)
workDir = tempname;
mkdir(workDir);
grid = localUnitSquareGrid();
initialLevels = localVerticalInterfaceLevels(0.50);
currentLevels = localVerticalInterfaceLevels(0.55);
exported = false(2, 1);

exported = precip_ExportBenchmarkSnapshots(workDir, [60; 120], exported, ...
    1, 30, grid, currentLevels, initialLevels, localSnapshotConfig());

verifyEqual(testCase, exported, false(2, 1));
verifyEmpty(testCase, dir(fullfile(workDir, 'benchmark_snapshot_*.png')));
end

function testBenchmarkSnapshotsExportMinuteNamesWhenDue(testCase)
workDir = tempname;
mkdir(workDir);
grid = localUnitSquareGrid();
initialLevels = localVerticalInterfaceLevels(0.50);
currentLevels = localVerticalInterfaceLevels(0.55);
exported = false(3, 1);

exported = precip_ExportBenchmarkSnapshots(workDir, [780; 1080; 7080], exported, ...
    5, 7080, grid, currentLevels, initialLevels, localSnapshotConfig());

verifyEqual(testCase, exported, true(3, 1));
expectedFiles = { ...
    fullfile(workDir, 'benchmark_snapshot_013min.png')
    fullfile(workDir, 'benchmark_snapshot_018min.png')
    fullfile(workDir, 'benchmark_snapshot_118min.png')};
for iFile = 1:numel(expectedFiles)
    verifyTrue(testCase, isfile(expectedFiles{iFile}));
    verifyGreaterThan(testCase, dir(expectedFiles{iFile}).bytes, 0);
end
end

function testBenchmarkSnapshotsUseSecondNamesForArtificialSmokeTimes(testCase)
workDir = tempname;
mkdir(workDir);
grid = localUnitSquareGrid();
initialLevels = localVerticalInterfaceLevels(0.50);
currentLevels = localVerticalInterfaceLevels(0.60);

exported = precip_ExportBenchmarkSnapshots(workDir, [1; 1.5; 2], false(3, 1), ...
    3, 2, grid, currentLevels, initialLevels, localSnapshotConfig());

verifyEqual(testCase, exported, true(3, 1));
verifyTrue(testCase, isfile(fullfile(workDir, 'benchmark_snapshot_001s.png')));
verifyTrue(testCase, isfile(fullfile(workDir, 'benchmark_snapshot_1p5s.png')));
verifyTrue(testCase, isfile(fullfile(workDir, 'benchmark_snapshot_002s.png')));
end

function testBenchmarkSnapshotsRejectInvalidGrid(testCase)
workDir = tempname;
mkdir(workDir);
grid = localUnitSquareGrid();
grid.V0T = [1, 2, 5];

verifyError(testCase, @() precip_ExportBenchmarkSnapshots(workDir, 60, false, ...
    1, 60, grid, [-1; 1; 1; -1], [-1; 1; 1; -1], localSnapshotConfig()), ...
    'RTSPHEM:Precipitate:InvalidSnapshotGrid');
end

function testBenchmarkSnapshotFilenameRejectsInvalidTarget(testCase)
verifyError(testCase, @() precip_BenchmarkSnapshotFilename(-1), ...
    'RTSPHEM:Precipitate:InvalidSnapshotTime');
verifyError(testCase, @() precip_BenchmarkSnapshotFilename([1, 2]), ...
    'RTSPHEM:Precipitate:InvalidSnapshotTime');
end

function testCompareZhangYoonBenchmarkWritesOutputs(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'reference_curves.csv');
writeLocalReferenceCurveCsv(referenceCsv);

report = precip_CompareZhangYoonBenchmark(workDir, referenceCsv);

verifyTrue(testCase, isfile(fullfile(workDir, 'zhang_yoon_area_comparison.csv')));
verifyTrue(testCase, isfile(fullfile(workDir, 'zhang_yoon_area_comparison.png')));
verifyTrue(testCase, isfile(fullfile(workDir, 'benchmark_comparison_report.md')));
verifyEqual(testCase, string(report.runDir), string(workDir));
verifyEqual(testCase, report.numSimulationRows, 6);
verifyEqual(testCase, report.numReferenceRows, 6);
verifyTrue(testCase, all(ismember( ...
    ["entire_domain", "first_pore", "first_three_pores", "zhang_upgradient", "zhang_middle"], ...
    report.regions)));

combined = readtable(fullfile(workDir, 'zhang_yoon_area_comparison.csv'), ...
    'Delimiter', ',', 'VariableNamingRule', 'preserve');
verifyTrue(testCase, all(ismember( ...
    {'source', 'case', 'region', 'time_min', 'precipitated_area_cm2', ...
    'precipitated_area_norm', 'data_role', 'note'}, ...
    combined.Properties.VariableNames)));
verifyTrue(testCase, any(strcmp(combined.source, 'RTSPHEM')));
verifyTrue(testCase, any(strcmp(combined.source, 'Zhang2010')));
verifyTrue(testCase, any(strcmp(combined.case, 'case_5')));
end

function testCompareZhangYoonBenchmarkReportsLimiterAndStabilityFlags(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
writeLocalRunDiagnostics(workDir, 2, "overshoot_c;advective_cfl_gt_1;mass_balance_drift");
referenceCsv = fullfile(workDir, 'reference_curves.csv');
writeLocalReferenceCurveCsv(referenceCsv);

precip_CompareZhangYoonBenchmark(workDir, referenceCsv);

reportText = string(fileread(fullfile(workDir, 'benchmark_comparison_report.md')));
verifyTrue(testCase, contains(reportText, "phreeqcTransportMaxFactor = 2"));
verifyTrue(testCase, contains(reportText, "finite PHREEQC transport limiter"));
verifyTrue(testCase, contains(reportText, "advective_cfl_gt_1"));
verifyTrue(testCase, contains(reportText, "mass_balance_drift"));
end

function testCompareZhangYoonBenchmarkSupportsEntireDomainYoonCase1(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'entire_domain_reference_curves.csv');
referenceTable = table( ...
    "Yoon2012", "case_1", "entire_domain", 2, 0.95, 3.6e4, ...
    "synthetic unit-test value for Yoon case 1 entire-domain curve", ...
    'VariableNames', {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note'});
writetable(referenceTable, referenceCsv);

report = precip_CompareZhangYoonBenchmark(workDir, referenceCsv);

verifyTrue(testCase, any(report.regions == "entire_domain"));
verifyEqual(testCase, report.numSimulationRows, 6);
combined = readtable(fullfile(workDir, 'zhang_yoon_area_comparison.csv'), ...
    'Delimiter', ',', 'VariableNamingRule', 'preserve');
simulationEntireDomain = strcmp(combined.source, 'RTSPHEM') & ...
    strcmp(combined.region, 'entire_domain');
referenceEntireDomain = strcmp(combined.source, 'Yoon2012') & ...
    strcmp(combined.case, 'case_1') & strcmp(combined.region, 'entire_domain');
verifyEqual(testCase, sum(simulationEntireDomain), 2);
verifyEqual(testCase, sum(referenceEntireDomain), 1);
verifyEqual(testCase, combined.precipitated_area_cm2(simulationEntireDomain), ...
    [0.30; 0.30], 'AbsTol', 1e-12);
end

function testCompareZhangYoonBenchmarkRejectsZhangEntireDomainReference(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'zhang_entire_domain_reference_curves.csv');
referenceTable = table( ...
    "Zhang2010", "experiment_25mM", "entire_domain", 2, 0.95, 3.6e4, ...
    "synthetic unit-test value for invalid Zhang entire-domain curve", ...
    'VariableNames', {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note'});
writetable(referenceTable, referenceCsv);

verifyError(testCase, @() precip_CompareZhangYoonBenchmark(workDir, referenceCsv), ...
    'RTSPHEM:Precipitate:InvalidReferenceCurves');
end

function testCompareZhangYoonBenchmarkSupportsZhangSelectedPoreRegions(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'zhang_selected_pore_reference_curves.csv');
referenceTable = table( ...
    ["Zhang2010"; "Zhang2010"; "Zhang2010"], ...
    ["experiment_25mM"; "experiment_25mM"; "experiment_25mM"], ...
    ["zhang_upgradient"; "zhang_middle"; "zhang_downgradient"], ...
    [0.25; 0.25; 0.25], ...
    [0.40; 0.30; 0.20], ...
    [4000; 3000; 2000], ...
    ["synthetic unit-test value for Zhang SI Fig S3b region"; ...
     "synthetic unit-test value for Zhang SI Fig S3b region"; ...
     "synthetic unit-test value for Zhang SI Fig S3b region"], ...
    'VariableNames', {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note'});
writetable(referenceTable, referenceCsv);

report = precip_CompareZhangYoonBenchmark(workDir, referenceCsv);

verifyTrue(testCase, all(ismember( ...
    ["zhang_upgradient", "zhang_middle", "zhang_downgradient"], ...
    report.regions)));
verifyEqual(testCase, report.numReferenceRows, 3);
reportText = fileread(fullfile(workDir, 'benchmark_comparison_report.md'));
verifyTrue(testCase, contains(reportText, ...
    'Zhang selected-pore regions are reference-only'));
verifyTrue(testCase, contains(reportText, ...
    'no matching RTSPHEM simulation rows'));
end

function testCompareZhangYoonBenchmarkRejectsYoonZhangSelectedPoreRegion(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'yoon_zhang_region_reference_curves.csv');
referenceTable = table( ...
    "Yoon2012", "case_1", "zhang_upgradient", 2, 0.95, 3.6e4, ...
    "synthetic unit-test value for invalid Yoon Zhang-region curve", ...
    'VariableNames', {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note'});
writetable(referenceTable, referenceCsv);

verifyError(testCase, @() precip_CompareZhangYoonBenchmark(workDir, referenceCsv), ...
    'RTSPHEM:Precipitate:InvalidReferenceCurves');
end

function testCompareZhangYoonBenchmarkSupportsLegacyAreaCsvWithoutTotal(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalLegacySimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'legacy_reference_curves.csv');
writeLocalReferenceCurveCsv(referenceCsv);

report = precip_CompareZhangYoonBenchmark(workDir, referenceCsv);

verifyEqual(testCase, report.numSimulationRows, 4);
verifyTrue(testCase, all(ismember( ...
    ["first_pore", "first_three_pores", "zhang_upgradient", "zhang_middle"], ...
    report.regions)));
reportText = fileread(fullfile(workDir, 'benchmark_comparison_report.md'));
verifyTrue(testCase, contains(reportText, 'Legacy simulation area CSV lacks total_net_solid_area_cm2'));
combined = readtable(fullfile(workDir, 'zhang_yoon_area_comparison.csv'), ...
    'Delimiter', ',', 'VariableNamingRule', 'preserve');
verifyFalse(testCase, any(strcmp(combined.source, 'RTSPHEM') & ...
    strcmp(combined.region, 'entire_domain')));
end

function testCompareZhangYoonBenchmarkRejectsBadReferenceCsv(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
badReferenceCsv = fullfile(workDir, 'bad_reference_curves.csv');
writetable(table("Zhang2010", "case_1", 'VariableNames', {'source', 'case'}), ...
    badReferenceCsv);

verifyError(testCase, @() precip_CompareZhangYoonBenchmark(workDir, badReferenceCsv), ...
    'RTSPHEM:Precipitate:InvalidReferenceCurves');
end

function testCompareZhangYoonBenchmarkFillsReferenceNormFromCm2(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'cm2_only_reference_curves.csv');
referenceTable = localReferenceCurveTable();
referenceTable.precipitated_area_norm(:) = NaN;
writetable(referenceTable, referenceCsv);

precip_CompareZhangYoonBenchmark(workDir, referenceCsv);

combined = readtable(fullfile(workDir, 'zhang_yoon_area_comparison.csv'), ...
    'Delimiter', ',', 'VariableNamingRule', 'preserve');
dataRoleColumn = strcmp(combined.Properties.VariableNames, 'data_role');
verifyTrue(testCase, any(dataRoleColumn));
referenceRows = strcmp(combined{:, dataRoleColumn}, 'reference');
verifyTrue(testCase, all(isfinite(combined.precipitated_area_norm(referenceRows))));
reportText = fileread(fullfile(workDir, 'benchmark_comparison_report.md'));
verifyTrue(testCase, contains(reportText, 'simulation-internal normalization'));
verifyTrue(testCase, contains(reportText, 'missing reference normalized values are filled from cm2'));
end

function testCompareZhangYoonBenchmarkRejectsInvalidSourceCasePair(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'invalid_pair_reference_curves.csv');
referenceTable = localReferenceCurveTable();
referenceTable.("case")(1) = "case_5";
writetable(referenceTable, referenceCsv);

verifyError(testCase, @() precip_CompareZhangYoonBenchmark(workDir, referenceCsv), ...
    'RTSPHEM:Precipitate:InvalidReferenceCurves');
end

function testCompareZhangYoonBenchmarkAcceptsRepositoryReferenceCurves(testCase)
moduleDir = testCase.TestData.moduleDir;
referenceCsv = fullfile(moduleDir, 'reference_data', 'zhang_yoon_reference_curves.csv');
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);

report = precip_CompareZhangYoonBenchmark(workDir, referenceCsv);

verifyGreaterThanOrEqual(testCase, report.numReferenceRows, 18);
verifyTrue(testCase, isfile(fullfile(workDir, 'zhang_yoon_area_comparison.csv')));
verifyTrue(testCase, isfile(fullfile(workDir, 'benchmark_comparison_report.md')));
combined = readtable(fullfile(workDir, 'zhang_yoon_area_comparison.csv'), ...
    'Delimiter', ',', 'VariableNamingRule', 'preserve');
referenceRows = strcmp(combined.data_role, 'reference');
yoonRows = referenceRows & strcmp(combined.source, 'Yoon2012');
zhangRows = referenceRows & strcmp(combined.source, 'Zhang2010');
verifyEqual(testCase, sum(referenceRows), 18);
verifyTrue(testCase, all(isfinite(combined.precipitated_area_norm(yoonRows))));
verifyTrue(testCase, any(abs(combined.precipitated_area_cm2(yoonRows) - 3.62e-4) < 1e-12));
verifyTrue(testCase, any(strcmp(combined.region(zhangRows), 'zhang_upgradient')));
reportText = fileread(fullfile(workDir, 'benchmark_comparison_report.md'));
verifyTrue(testCase, contains(reportText, ...
    'Zhang selected-pore regions are reference-only'));
verifyTrue(testCase, contains(reportText, ...
    'missing reference normalized values are filled from cm2'));
end

function testCompareZhangYoonBenchmarkRejectsHeaderOnlyReferenceCsv(testCase)
workDir = tempname;
mkdir(workDir);
writeLocalSimulationAreaCsv(workDir);
writeLocalBenchmarkComparisonCsv(workDir);
referenceCsv = fullfile(workDir, 'header_only_reference_curves.csv');
fid = fopen(referenceCsv, 'w');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, 'source,case,region,time_min,precipitated_area_norm,precipitated_area_cm2,note\n');
clear cleanupObj;

verifyError(testCase, @() precip_CompareZhangYoonBenchmark(workDir, referenceCsv), ...
    'RTSPHEM:Precipitate:MissingReferenceCurves');
end

function testSignedPhreeqcSingleCellLocalDatabaseSmoke(testCase)
databasePath = findInstalledPhreeqcDatabase(testCase);
assumeTrue(testCase, ~isempty(databasePath), 'No PHREEQC database found for COM integration test.');
assumeTrue(testCase, canCreateIPhreeqcCom(), 'IPhreeqcCOM is not available on this machine.');

options = struct();
options.databasePath = databasePath;
options.timeStepIndex = 1;
options.timeStepSize = 1000;
options.rateLaw = 'database_calcite';
options.maxSpecificSurfaceArea = 1e12;
options.kineticsReservoirMoles = 10;

premixedState = localDirectionalPhreeqcState(10^(-8.3) * 1e-3, 12.5e-6, 12.5e-6);
options.workDir = tempname;
premixedResult = precip_RunPhreeqcCalciteBatchSigned(premixedState, options);

verifyEqual(testCase, string(premixedResult.databasePath), string(databasePath));
verifyTrue(testCase, isfile(char(premixedResult.inputPath)));
[inputDir, ~, ~] = fileparts(char(premixedResult.inputPath));
verifyEqual(testCase, string(inputDir), string(options.workDir));
inputText = string(fileread(char(premixedResult.inputPath)));
verifyTrue(testCase, contains(inputText, "USER_PUNCH"));
verifyTrue(testCase, contains(inputText, "-headings KIN_DELTA_Calcite RATE_Calcite"));
verifyTrue(testCase, all(isfinite(premixedResult.calciteDeltaMoles)));
verifyTrue(testCase, all(isfinite(premixedResult.pH)));
verifyTrue(testCase, all(isfinite(premixedResult.ca_total_mol_cm3)));
verifyTrue(testCase, all(isfinite(premixedResult.c_total_mol_cm3)));
verifyTrue(testCase, all(isfinite(premixedResult.calciteSI)));
verifyEqual(testCase, premixedResult.calciteSignedRate_mol_s, ...
    premixedResult.calciteDeltaMoles ./ options.timeStepSize, 'AbsTol', 1e-18);

supersaturatedState = localDirectionalPhreeqcState(1e-10, 1e-3, 1e-3);
options.workDir = tempname;
options.timeStepIndex = 2;
supersaturatedResult = precip_RunPhreeqcCalciteBatchSigned(supersaturatedState, options);

undersaturatedState = localDirectionalPhreeqcState(1e-6, 0, 0);
options.workDir = tempname;
options.timeStepIndex = 3;
undersaturatedResult = precip_RunPhreeqcCalciteBatchSigned(undersaturatedState, options);

verifyGreaterThan(testCase, supersaturatedResult.calciteDeltaMoles, 0);
verifyLessThan(testCase, undersaturatedResult.calciteDeltaMoles, 0);
verifyEqual(testCase, supersaturatedResult.calcitePrecipitatedMoles, ...
    supersaturatedResult.calciteDeltaMoles, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, undersaturatedResult.calciteDissolvedMoles, ...
    -undersaturatedResult.calciteDeltaMoles, 'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function databasePath = findInstalledPhreeqcDatabase(testCase)
moduleDir = testCase.TestData.moduleDir;
rtmDir = fileparts(moduleDir);
envDatabase = getenv('PHREEQC_DATABASE_PATH');
candidates = {
    envDatabase
    fullfile(rtmDir, 'couplePhreeqc', 'phreeqc-m.dat')
    'C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat'
    };
databasePath = '';
for iCandidate = 1:numel(candidates)
    if exist(candidates{iCandidate}, 'file')
        databasePath = candidates{iCandidate};
        return;
    end
end
end

function ok = canCreateIPhreeqcCom()
ok = false;
try
    iphreeqc = actxserver('IPhreeqcCOM.Object');
    cleanupObj = onCleanup(@() safeDeleteComObject(iphreeqc));
    ok = true;
    clear cleanupObj;
catch
end
end

function safeDeleteComObject(obj)
try
    delete(obj);
catch
end
end

function state = localDirectionalPhreeqcState(hMolCm3, caMolCm3, cMolCm3)
state = struct();
state.h_mol_cm3 = hMolCm3;
state.ca_mol_cm3 = caMolCm3;
state.c_mol_cm3 = cMolCm3;
state.na_mol_cm3 = 0;
state.cl_mol_cm3 = max(2 * caMolCm3, 1e-8);
state.interface_area_cm2 = 1e4;
state.water_volume_cm3 = 1000;
state.calcite_moles = 10;
end

function grid = localUnitSquareGrid()
grid = struct();
grid.coordV = [ ...
    0, 0; ...
    1, 0; ...
    1, 1; ...
    0, 1];
grid.V0T = [1, 2, 3; 1, 3, 4];
grid.areaT = [0.5; 0.5];
grid.baryT = [2/3, 1/3; 1/3, 2/3];
end

function grid = localSingleTriangleGrid()
grid = struct();
grid.coordV = [ ...
    0, 0; ...
    1, 0; ...
    0, 1];
grid.V0T = [1, 2, 3];
grid.areaT = 0.5;
grid.baryT = [1/3, 1/3];
end

function metrics = localAreaMetricsStruct()
metrics = struct();
metrics.totalNetSolidArea_cm2 = 0.30;
metrics.firstPoreNetSolidArea_cm2 = 0.10;
metrics.firstThreePoresNetSolidArea_cm2 = 0.20;
metrics.totalSolidArea_cm2 = 0.80;
metrics.firstPoreSolidArea_cm2 = 0.25;
metrics.firstThreePoresSolidArea_cm2 = 0.55;
end

function levels = localVerticalInterfaceLevels(xZero)
grid = localUnitSquareGrid();
levels = grid.coordV(:, 1) - xZero;
end

function config = localSnapshotConfig()
config = struct();
config.inletA = struct('name', 'CaCl2');
config.inletB = struct('name', 'Na2CO3');
end

function writeLocalSimulationAreaCsv(workDir)
metrics = localAreaMetricsStruct();
csvFile = fullfile(workDir, 'precipitation_area_timeseries.csv');
precip_InitializePrecipitationAreaTimeseries(csvFile);
precip_AppendPrecipitationAreaTimeseries(csvFile, 1, 60, metrics);
metrics.firstPoreNetSolidArea_cm2 = 0.20;
metrics.firstThreePoresNetSolidArea_cm2 = 0.35;
precip_AppendPrecipitationAreaTimeseries(csvFile, 2, 120, metrics);
end

function writeLocalLegacySimulationAreaCsv(workDir)
csvFile = fullfile(workDir, 'precipitation_area_timeseries.csv');
fid = fopen(csvFile, 'w');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, ['timestep,time_s,first_pore_net_solid_area_cm2,', ...
    'first_three_pores_net_solid_area_cm2,first_pore_solid_area_cm2,', ...
    'first_three_pores_solid_area_cm2\n']);
fprintf(fid, '1,60,0.10,0.20,0.25,0.55\n');
fprintf(fid, '2,120,0.20,0.35,0.35,0.65\n');
clear cleanupObj;
end

function writeLocalBenchmarkComparisonCsv(workDir)
csvFile = fullfile(workDir, 'benchmark_comparison_times_log.csv');
fid = fopen(csvFile, 'w');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, ['target_time_s,target_time_min,captured_timestep,captured_time_s,', ...
    'porosity,permeability_mD,k_k0,rate_mol_cm2_s,surface_area_cm2,grain_volume_cm3,', ...
    'first_pore_net_solid_area_cm2,first_three_pores_net_solid_area_cm2,', ...
    'first_pore_solid_area_cm2,first_three_pores_solid_area_cm2\n']);
fprintf(fid, '60,1,1,60,0.4,10,1,0,1,1,0.1,0.2,0.3,0.4\n');
fprintf(fid, '120,2,2,120,0.39,9,0.9,0,1,1,0.2,0.35,0.4,0.55\n');
clear cleanupObj;
end

function writeLocalRunDiagnostics(workDir, limiterFactor, diagnosticFlag)
metadataFile = fullfile(workDir, 'run_metadata.json');
fidMeta = fopen(metadataFile, 'w');
cleanupMeta = onCleanup(@() fclose(fidMeta));
fprintf(fidMeta, '{"parameters":{"phreeqcTransportMaxFactor":%.15g}}', limiterFactor);
clear cleanupMeta;

stabilityFile = fullfile(workDir, 'stability_diagnostics_log.csv');
fidStab = fopen(stabilityFile, 'w');
cleanupStab = onCleanup(@() fclose(fidStab));
fprintf(fidStab, ['timestep,time_s,dt_s,c_min,c_max,total_mass_mol,', ...
    'mass_balance_residual_mol,mass_balance_relative,local_u_max_cm_s,', ...
    'local_pe_max,advective_cfl,reaction_cfl,normal_speed_max_cm_s,', ...
    'pore_components,largest_pore_component_fraction,porosity,permeability_mD,', ...
    'k_k0,transport_upwind,diagnostic_flag\n']);
fprintf(fidStab, '1,60,60,0,1e-9,1e-12,1e-13,1,0.1,1,320,0,0,1,1,0.5,10,1,exp,"%s"\n', diagnosticFlag);
fprintf(fidStab, '2,120,60,0,2e-9,1e-12,1e-13,2,0.1,1,120,0,0,1,1,0.5,10,1,exp,"overshoot_c"\n');
clear cleanupStab;
end

function writeLocalReferenceCurveCsv(referenceCsv)
referenceTable = localReferenceCurveTable();
writetable(referenceTable, referenceCsv);
end

function referenceTable = localReferenceCurveTable()
referenceTable = table( ...
    ["Zhang2010"; "Zhang2010"; "Yoon2012"; "Yoon2012"; "Yoon2012"; "Yoon2012"], ...
    ["experiment_25mM"; "experiment_25mM"; "case_1"; "case_1"; "case_5"; "case_5"], ...
    ["zhang_upgradient"; "zhang_middle"; "first_pore"; "first_three_pores"; "first_pore"; "first_three_pores"], ...
    [1; 1; 2; 2; 2; 2], ...
    [0.10; 0.20; 0.15; 0.22; 0.18; 0.28], ...
    [0.010; 0.020; 0.015; 0.022; 0.018; 0.028], ...
    ["synthetic unit-test value"; "synthetic unit-test value"; ...
     "synthetic unit-test value"; "synthetic unit-test value"; ...
     "synthetic unit-test value"; "synthetic unit-test value"], ...
    'VariableNames', {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note'});
end
