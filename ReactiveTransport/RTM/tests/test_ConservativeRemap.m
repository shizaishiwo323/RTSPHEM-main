function tests = test_ConservativeRemap
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

function testIdentityRemapPreservesComponentAndMineralMoles(testCase)
state = validState([1; 2], [0.5; 0.25]);
oldGeometry = geometryWithWater([1; 2]);
newGeometry = geometryWithWater([1; 2]);

[newState, ledger] = rtm.geometry.ConservativeRemap(state, oldGeometry, newGeometry);

verifyEqual(testCase, newState.component_moles, state.component_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, newState.mineral_moles, state.mineral_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.unassigned_component_moles, [0 0], ...
    'AbsTol', 1e-18);
end

function testDefaultSameGridRemapUsesSparseOverlapLedger(testCase)
numCells = 2000;
waterVolume = ones(numCells, 1);
state = validState(waterVolume, ones(numCells, 1));
oldGeometry = geometryWithWater(waterVolume);
newGeometry = geometryWithWater(waterVolume);

[newState, ledger] = rtm.geometry.ConservativeRemap( ...
    state, oldGeometry, newGeometry);

verifyEqual(testCase, newState.component_moles, state.component_moles, ...
    'AbsTol', 1e-18);
verifyTrue(testCase, issparse(ledger.overlap_volume_cm3));
verifyEqual(testCase, nnz(ledger.overlap_volume_cm3), numCells);
end

function testSplitCellRemapConservesComponentTotals(testCase)
state = validState(1, 0.5);
oldGeometry = geometryWithWater(1);
newGeometry = geometryWithWater([0.25; 0.75]);
options.overlap_volume_cm3 = [0.25, 0.75];
options.new_mineral_moles = [0.1; 0.4];

[newState, ledger] = rtm.geometry.ConservativeRemap(state, oldGeometry, newGeometry, options);

verifyEqual(testCase, newState.component_moles, ...
    [2.5e-7, 5e-7; 7.5e-7, 1.5e-6], ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, sum(newState.component_moles, 1), ...
    sum(state.component_moles, 1), 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, newState.mineral_moles, [0.1; 0.4], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], ...
    'AbsTol', 1e-18);
end

function testNewWaterCellWithoutOverlapReceivesZeroComponentMoles(testCase)
state = validState(1, 0.5);
oldGeometry = geometryWithWater(1);
newGeometry = geometryWithWater([1; 0.2]);
options.overlap_volume_cm3 = [1, 0];
options.new_mineral_moles = [0.5; 0.2];

[newState, ledger] = rtm.geometry.ConservativeRemap(state, oldGeometry, newGeometry, options);

verifyEqual(testCase, newState.component_moles(1, :), state.component_moles, ...
    'AbsTol', 1e-18);
verifyEqual(testCase, newState.component_moles(2, :), [0 0], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.new_water_without_overlap_cells, 2);
verifyEqual(testCase, ledger.component_residual_moles, [0 0], ...
    'AbsTol', 1e-18);
end

function testUnmappedOldMolesAreReportedAsResidual(testCase)
state = validState([1; 1], [0.5; 0.5]);
oldGeometry = geometryWithWater([1; 1]);
newGeometry = geometryWithWater(1);
options.overlap_volume_cm3 = [1; 0];
options.new_mineral_moles = 0.5;

[newState, ledger] = rtm.geometry.ConservativeRemap(state, oldGeometry, newGeometry, options);

verifyEqual(testCase, newState.component_moles, [1e-6, 2e-6], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.unassigned_component_moles, [2e-6, 4e-6], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, ledger.component_residual_moles, -[2e-6, 4e-6], ...
    'AbsTol', 1e-18);
end

function testRejectsNegativeNewMineralMolesInsteadOfClipping(testCase)
state = validState(1, 0.5);
oldGeometry = geometryWithWater(1);
newGeometry = geometryWithWater(1);
options.new_mineral_moles = -0.1;

verifyError(testCase, ...
    @() rtm.geometry.ConservativeRemap(state, oldGeometry, newGeometry, options), ...
    'RTSPHEM:Geometry:NegativeMineralMoles');
end

function testRejectsNegativeWaterVolumeInsteadOfClipping(testCase)
state = validState(1, 0.5);
oldGeometry = geometryWithWater(-1);
newGeometry = geometryWithWater(1);

verifyError(testCase, ...
    @() rtm.geometry.ConservativeRemap(state, oldGeometry, newGeometry), ...
    'RTSPHEM:Geometry:NegativeWaterVolume');
end

function state = validState(waterVolume, mineralMoles)
numCells = numel(waterVolume);
state = struct();
state.component_names = {'Ca_total', 'C_total'};
state.component_moles = [1e-6 * (1:numCells)', 2e-6 * (1:numCells)'];
state.mineral_names = {'Calcite'};
state.mineral_moles = mineralMoles(:);
state.temperature_C = 25 * ones(numCells, 1);
state.pressure_atm = ones(numCells, 1);
state.time_s = 12;
end

function geometry = geometryWithWater(waterVolume)
waterVolume = waterVolume(:);
geometry = struct();
geometry.water_volume_cm3 = waterVolume;
geometry.cell_volume_cm3 = waterVolume;
geometry.solid_volume_cm3 = zeros(size(waterVolume));
geometry.fluid_fraction = double(waterVolume > 0);
end
