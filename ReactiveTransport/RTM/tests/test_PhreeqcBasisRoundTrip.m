function tests = test_PhreeqcBasisRoundTrip
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testExternalTstConfigDeclaresConservedBasisAndDerivedSpecies(testCase)
cfg = struct();
cfg.chemistry.mode = 'external_tst_phreeqc';
cfg.phreeqc.engine = 'mock';

validated = rtm.config.ValidateReactiveTransportConfig(cfg);

verifyEqual(testCase, validated.chemistry.basis, ...
    {'Ca', 'C', 'Na', 'Cl', 'Alkalinity'});
verifyEqual(testCase, validated.chemistry.derived, ...
    {'pH', 'H+', 'HCO3-', 'CO3-2', 'SI_Calcite'});
end

function testPhreeqcBasisRejectsFreeHydrogenAsTransportComponent(testCase)
state = carbonateState({'Ca', 'C', 'H+'});

verifyError(testCase, ...
    @() rtm.chemistry.ValidatePhreeqcTransportBasis(state), ...
    'RTSPHEM:Chemistry:ForbiddenTransportComponent');
end

function testExternalTstBackendRejectsDerivedSpeciesInTransportState(testCase)
state = carbonateState({'Ca', 'C', 'pH'});
geometry = singleInterfaceGeometry();
options = struct('h_mol_cm3', 1e-7, ...
    'h_activity_mol_cm3', 1e-7, ...
    'runBatchFunction', @mockRunBatch);

verifyError(testCase, ...
    @() rtm.chemistry.ExternalTstPhreeqcBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:ForbiddenTransportComponent');
end

function testPhreeqcKineticsBackendRejectsCarbonateSpeciesTransportState(testCase)
state = carbonateState({'Ca', 'C', 'HCO3-'});
geometry = singleInterfaceGeometry();
options = struct('h_mol_cm3', 1e-7, 'runBatchFunction', @mockRunBatch);

verifyError(testCase, ...
    @() rtm.chemistry.PhreeqcKineticsBackend(state, geometry, 1, options), ...
    'RTSPHEM:Chemistry:ForbiddenTransportComponent');
end

function testBuildPhreeqcBasisStateConvertsConservedMolesToTotals(testCase)
state = carbonateState({'Ca', 'C', 'Na', 'Cl', 'Alkalinity'});
state.component_moles = [
    2e-9, 3e-9, 4e-9, 5e-9, 6e-9
    8e-9, 9e-9, 1e-8, 1.1e-8, 1.2e-8
    ];
state.mineral_moles = [1e-3; 2e-3];
state.temperature_C = [25; 30];
state.pressure_atm = [1; 1];
geometry = singleInterfaceGeometry();
geometry.water_volume_cm3 = [1; 2];
geometry.interface_area_cm2 = [0.1; 0.2];

basisState = rtm.chemistry.BuildPhreeqcBasisState(state, geometry, ...
    struct('h_mol_cm3', [1e-7; 2e-7]));

verifyEqual(testCase, basisState.ca_total_mol_cm3, [2e-9; 4e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, basisState.c_total_mol_cm3, [3e-9; 4.5e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, basisState.na_total_mol_cm3, [4e-9; 5e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, basisState.cl_total_mol_cm3, [5e-9; 5.5e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, basisState.alkalinity_mol_cm3, [6e-9; 6e-9], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, basisState.h_mol_cm3, [1e-7; 2e-7], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, basisState.water_volume_cm3, [1; 2], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, basisState.calcite_moles, [1e-3; 2e-3], ...
    'AbsTol', 1e-18);
end

function testBuildPhreeqcBasisStateCanScaleComponentTotalsWithReactionWater(testCase)
state = carbonateState({'Ca', 'C'});
state.component_moles = [2e-9, 4e-9; 8e-9, 1.2e-8];
state.mineral_moles = [1e-3; 2e-3];
state.temperature_C = [25; 25];
state.pressure_atm = [1; 1];
geometry = singleInterfaceGeometry();
geometry.water_volume_cm3 = [1e-6; 2e-6];
reactionWaterVolume = [1000; 1000];

basisState = rtm.chemistry.BuildPhreeqcBasisState(state, geometry, ...
    struct('reactionWaterVolumeCm3', reactionWaterVolume, ...
    'componentWaterVolumeCm3', reactionWaterVolume));

verifyEqual(testCase, basisState.water_volume_cm3, geometry.water_volume_cm3, ...
    'AbsTol', 0);
verifyEqual(testCase, basisState.reaction_water_volume_cm3, ...
    reactionWaterVolume, 'AbsTol', 0);
verifyEqual(testCase, basisState.ca_total_mol_cm3, ...
    state.component_moles(:, 1) ./ reactionWaterVolume, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyEqual(testCase, basisState.c_total_mol_cm3, ...
    state.component_moles(:, 2) ./ reactionWaterVolume, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
end

function testBuildPhreeqcBasisStateRejectsDerivedTransportComponents(testCase)
state = carbonateState({'Ca', 'C', 'pH'});
geometry = singleInterfaceGeometry();

verifyError(testCase, ...
    @() rtm.chemistry.BuildPhreeqcBasisState(state, geometry, struct()), ...
    'RTSPHEM:Chemistry:ForbiddenTransportComponent');
end

function state = carbonateState(componentNames)
state = struct();
state.component_names = componentNames;
state.component_moles = repmat(1e-8, 1, numel(componentNames));
state.mineral_names = {'Calcite'};
state.mineral_moles = 1e-3;
state.temperature_C = 25;
state.pressure_atm = 1;
state.time_s = 0;
end

function geometry = singleInterfaceGeometry()
geometry = struct();
geometry.water_volume_cm3 = 1;
geometry.solid_volume_cm3 = 1e-3;
geometry.cell_volume_cm3 = 1.001;
geometry.fluid_fraction = 1;
geometry.interface_area_cm2 = 1;
geometry.interface_h_cm = 1;
geometry.interface_centroid_cm = [0 0];
geometry.cell_centroid_cm = [0 0];
end

function batchResult = mockRunBatch(batchState, ~)
dissolved = zeros(size(batchState.water_volume_cm3(:)));
if isfield(batchState, 'prescribed_calcite_dissolved_moles')
    dissolved = batchState.prescribed_calcite_dissolved_moles(:);
end
batchResult = struct();
batchResult.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + ...
    dissolved ./ max(batchState.water_volume_cm3(:), eps);
batchResult.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + ...
    dissolved ./ max(batchState.water_volume_cm3(:), eps);
batchResult.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
batchResult.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
batchResult.calciteDissolvedMoles = dissolved;
end
