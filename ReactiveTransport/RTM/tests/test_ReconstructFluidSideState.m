function tests = test_ReconstructFluidSideState
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

function testFallbackUsesCellAverageWhenNeighborsAreMissing(testCase)
geometry = geometryFixture([1 1], [0 0; 1 0], [0.25 0; NaN NaN], [1; 0]);
state = stateFromConcentration([2e-6; 5e-6], geometry);
connectivity.cell_neighbors = {[], []};

interfaceState = rtm.chemistry.ReconstructFluidSideState(state, geometry, connectivity);

verifyEqual(testCase, interfaceState.component_concentration_mol_cm3(1), 2e-6, ...
    'AbsTol', 1e-18);
verifyTrue(testCase, interfaceState.fallback_used(1));
verifyEqual(testCase, interfaceState.quality_flag(1), "cell_average");
verifyTrue(testCase, isnan(interfaceState.component_concentration_mol_cm3(2)));
end

function testLinearReconstructionEvaluatesAtInterfaceCentroid(testCase)
cellCentroids = [0 0; 1 0; 0 1; -1 0; 0 -1];
interfaceCentroids = [0.5 0.5; NaN NaN; NaN NaN; NaN NaN; NaN NaN];
geometry = geometryFixture(ones(5, 1), cellCentroids, interfaceCentroids, [1; 0; 0; 0; 0]);
concentration = 1e-6 + 2e-7 .* cellCentroids(:, 1) + 3e-7 .* cellCentroids(:, 2);
state = stateFromConcentration(concentration, geometry);
connectivity.cell_neighbors = {[2 3 4 5], 1, 1, 1, 1};

interfaceState = rtm.chemistry.ReconstructFluidSideState(state, geometry, connectivity);

verifyEqual(testCase, interfaceState.component_concentration_mol_cm3(1), 1.25e-6, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyFalse(testCase, interfaceState.fallback_used(1));
verifyEqual(testCase, interfaceState.quality_flag(1), "linear_limited");
end

function testLimiterBoundsExtrapolatedStateByLocalStencil(testCase)
geometry = geometryFixture(ones(4, 1), ...
    [0 0; 0.1 0; 0 0.1; 0.1 0.1], ...
    [10 10; NaN NaN; NaN NaN; NaN NaN], [1; 0; 0; 0]);
state = stateFromConcentration([1e-6; 2e-6; 2e-6; 2e-6], geometry);
connectivity.cell_neighbors = {[2 3 4], 1, 1, 1};

interfaceState = rtm.chemistry.ReconstructFluidSideState(state, geometry, connectivity);

verifyLessThanOrEqual(testCase, interfaceState.component_concentration_mol_cm3(1), 2e-6);
verifyGreaterThanOrEqual(testCase, interfaceState.component_concentration_mol_cm3(1), 1e-6);
verifyFalse(testCase, interfaceState.fallback_used(1));
end

function testRejectsNegativeGeometryMeasure(testCase)
geometry = geometryFixture([1; -0.1], [0 0; 1 0], [0.25 0; NaN NaN], [1; 0]);
state = stateFromConcentration([2e-6; 5e-6], geometryFixture([1; 0.1], ...
    [0 0; 1 0], [0.25 0; NaN NaN], [1; 0]));
connectivity.cell_neighbors = {[], []};

verifyError(testCase, ...
    @() rtm.chemistry.ReconstructFluidSideState(state, geometry, connectivity), ...
    'RTSPHEM:Chemistry:NegativeGeometryMeasure');
end

function geometry = geometryFixture(waterVolume, cellCentroids, interfaceCentroids, interfaceArea)
geometry = struct();
geometry.water_volume_cm3 = waterVolume(:);
geometry.cell_centroid_cm = cellCentroids;
geometry.interface_centroid_cm = interfaceCentroids;
geometry.interface_area_cm2 = interfaceArea(:);
geometry.active_fluid_cell = waterVolume(:) > 0;
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
