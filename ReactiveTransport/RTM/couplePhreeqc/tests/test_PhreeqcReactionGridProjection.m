function tests = test_PhreeqcReactionGridProjection
tests = functiontests(localfunctions);
end

function testAggregateFineStateToReactionGridUsesWaterWeightedConcentrations(testCase)
fineState = struct();
fineState.h_mol_cm3 = [1; 3; 10] * 1e-6;
fineState.ca_total_mol_cm3 = [2; 4; 20] * 1e-9;
fineState.c_total_mol_cm3 = [5; 7; 30] * 1e-9;
fineState.na_total_mol_cm3 = [0; 0; 0];
fineState.cl_total_mol_cm3 = [1; 3; 10] * 1e-6;
fineState.ca_mol_cm3 = fineState.ca_total_mol_cm3;
fineState.c_mol_cm3 = fineState.c_total_mol_cm3;
fineState.na_mol_cm3 = fineState.na_total_mol_cm3;
fineState.cl_mol_cm3 = fineState.cl_total_mol_cm3;
fineState.interface_area_cm2 = [2; 3; 5] * 1e-4;
fineState.water_volume_cm3 = [2; 6; 4] * 1e-6;
fineState.calcite_moles = [11; 13; 17] * 1e-9;

projection = struct();
projection.fineToReactionCell = [1; 1; 2];
projection.numReactionCells = 2;

coarse = AggregatePhreeqcStateToReactionGrid(fineState, projection);

verifyEqual(testCase, coarse.h_mol_cm3, [2.5e-6; 10e-6], 'RelTol', 1e-12);
verifyEqual(testCase, coarse.ca_total_mol_cm3, [3.5e-9; 20e-9], 'RelTol', 1e-12);
verifyEqual(testCase, coarse.c_total_mol_cm3, [6.5e-9; 30e-9], 'RelTol', 1e-12);
verifyEqual(testCase, coarse.interface_area_cm2, [5; 5] * 1e-4, 'RelTol', 1e-12);
verifyEqual(testCase, coarse.water_volume_cm3, [8; 4] * 1e-6, 'RelTol', 1e-12);
verifyEqual(testCase, coarse.calcite_moles, [24; 17] * 1e-9, 'RelTol', 1e-12);
end

function testExpandReactionGridResultDistributesKineticRatesByFineInterface(testCase)
fineTemplate = struct();
fineTemplate.h_mol_cm3 = zeros(3, 1);
fineTemplate.ca_total_mol_cm3 = zeros(3, 1);
fineTemplate.c_total_mol_cm3 = zeros(3, 1);
fineTemplate.na_total_mol_cm3 = zeros(3, 1);
fineTemplate.cl_total_mol_cm3 = zeros(3, 1);
fineTemplate.interface_area_cm2 = [1; 0; 2] * 1e-4;

coarseResult = struct();
coarseResult.h_mol_cm3 = [4; 8] * 1e-6;
coarseResult.ca_total_mol_cm3 = [1; 2] * 1e-9;
coarseResult.c_total_mol_cm3 = [3; 4] * 1e-9;
coarseResult.na_total_mol_cm3 = [0; 0];
coarseResult.cl_total_mol_cm3 = [4; 8] * 1e-6;
coarseResult.pH = [5; 6];
coarseResult.calciteKinDeltaRate_mol_s = [9; 7] * 1e-8;
coarseResult.calciteRate_mol_s = [1; 2] * 1e-10;
coarseResult.calciteDissolvedMoles = [5; 6] * 1e-12;

projection = struct();
projection.fineToReactionCell = [1; 1; 2];
projection.numReactionCells = 2;

fineResult = ExpandPhreeqcReactionGridResult(coarseResult, fineTemplate, projection);

verifyEqual(testCase, fineResult.h_mol_cm3, [4; 4; 8] * 1e-6, 'RelTol', 1e-12);
verifyEqual(testCase, fineResult.ca_total_mol_cm3, [1; 1; 2] * 1e-9, 'RelTol', 1e-12);
verifyEqual(testCase, fineResult.calciteKinDeltaRate_mol_s, [9; 0; 7] * 1e-8, 'RelTol', 1e-12);
verifyEqual(testCase, fineResult.calciteRate_mol_s, [1; 0; 2] * 1e-10, 'RelTol', 1e-12);
verifyEqual(testCase, fineResult.calciteDissolvedMoles, [5; 0; 6] * 1e-12, 'RelTol', 1e-12);
verifyEqual(testCase, fineResult.calciteKinDeltaRate_mol_s(2), 0, ...
    'AbsTol', 1e-18, 'Fine cells without interface area must not receive kinetic rate.');
end
