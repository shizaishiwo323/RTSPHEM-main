function tests = test_FluidSideReconstructionOneSided
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

function testOneSidedStencilExcludesNeighborsBehindFluidNormal(testCase)
geometry = geometryFixture();
concentration = [1; 1; 1; 1; 100; 100] .* 1e-6;
state = stateFromConcentration(concentration, geometry);
connectivity.cell_neighbors = {[2 3 4 5 6], 1, 1, 1, 1, 1};
options = struct('useLimiter', false);

interfaceState = rtm.chemistry.ReconstructFluidSideState( ...
    state, geometry, connectivity, options);

verifyEqual(testCase, interfaceState.component_concentration_mol_cm3(1), ...
    1e-6, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyFalse(testCase, interfaceState.fallback_used(1));
verifyEqual(testCase, interfaceState.quality_flag(1), "linear_limited");
end

function testOneSidedStencilFallsBackWhenFluidSideIsUnderdetermined(testCase)
geometry = geometryFixture();
concentration = [2; 4; 4; 4; 100; 100] .* 1e-6;
state = stateFromConcentration(concentration, geometry);
connectivity.cell_neighbors = {[2], 1, 1, 1, 1, 1};

interfaceState = rtm.chemistry.ReconstructFluidSideState( ...
    state, geometry, connectivity);

verifyEqual(testCase, interfaceState.component_concentration_mol_cm3(1), ...
    2e-6, 'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, interfaceState.fallback_used(1));
verifyEqual(testCase, interfaceState.quality_flag(1), "cell_average");
end

function geometry = geometryFixture()
geometry = struct();
geometry.water_volume_cm3 = ones(6, 1);
geometry.cell_centroid_cm = [
    0 0
    1 0
    1 1
    1 -1
    -1 0
    -1 1];
geometry.interface_centroid_cm = [
    0.5 0
    NaN NaN
    NaN NaN
    NaN NaN
    NaN NaN
    NaN NaN];
geometry.interface_area_cm2 = [1; 0; 0; 0; 0; 0];
geometry.interface_normal = [
    1 0
    NaN NaN
    NaN NaN
    NaN NaN
    NaN NaN
    NaN NaN];
geometry.active_fluid_cell = true(6, 1);
end

function state = stateFromConcentration(concentration, geometry)
state = struct();
state.component_names = {'H_reactant'};
state.component_moles = concentration(:) .* geometry.water_volume_cm3(:);
state.mineral_names = {'Calcite'};
state.mineral_moles = ones(numel(concentration), 1);
state.temperature_C = 25 * ones(numel(concentration), 1);
state.pressure_atm = ones(numel(concentration), 1);
state.time_s = 0;
end
