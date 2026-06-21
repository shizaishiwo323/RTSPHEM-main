function text = BuildCalcitePhreeqcInput(state, options)
% BuildCalcitePhreeqcInput - Build PHREEQC input for single-mineral calcite coupling.
%
% Inputs are stored in mol/cm^3 to match the RTM code and converted to
% PHREEQC mol/kgw-style numbers with the common dilute-solution approximation
% 1 mol/L = 1e-3 mol/cm^3. Following the reference Live Script coupling,
% the PHREEQC batch solution uses a reference water mass by default; the PNM
% caller scales KIN_DELTA back to each cell's real water volume after the run.

arguments
    state struct
    options struct
end

h = requireColumn(state, 'h_mol_cm3');
numCells = numel(h);
ca = optionalColumn(state, 'ca_mol_cm3', numCells, 0);
c = optionalColumn(state, 'c_mol_cm3', numCells, 0);
na = optionalColumn(state, 'na_mol_cm3', numCells, 0);
cl = optionalColumn(state, 'cl_mol_cm3', numCells, 0);
interfaceAreaCm2 = optionalColumn(state, 'interface_area_cm2', numCells, 1);
waterVolumeCm3 = optionalColumn(state, 'water_volume_cm3', numCells, 1000);
calciteMoles = optionalColumn(state, 'calcite_moles', numCells, 1);

dt = getOption(options, 'timeStepSize', 1);
temperatureC = getOption(options, 'temperatureC', 25);
solutionUnits = char(getOption(options, 'solutionUnits', 'mol/kgw'));
minWaterKg = getOption(options, 'minWaterKg', 1e-12);
solutionWaterKg = getOption(options, 'solutionWaterKg', 1);
writeSolutionWaterLine = logical(getOption(options, 'writeSolutionWaterLine', false));
minCalciteMoles = getOption(options, 'minCalciteMoles', 0);
minSpecificSurfaceArea = getOption(options, 'minSpecificSurfaceArea', 0);
minKineticsMoles = getOption(options, 'minKineticsMoles', 1e-30);
minKineticsM0 = getOption(options, 'minKineticsM0', minKineticsMoles);
kineticsReservoirMoles = getOption(options, 'kineticsReservoirMoles', 1);
minKineticsSurfaceArea = getOption(options, 'minKineticsSurfaceArea', 0);
maxSpecificSurfaceArea = getOption(options, 'maxSpecificSurfaceArea', 10);
minSolutionWaterKg = getOption(options, 'minSolutionWaterKg', minWaterKg);
minHForPH = getOption(options, 'minHForPHMolL', 1e-7);
kineticsCorrectionFactor = getOption(options, 'kineticsCorrectionFactor', 1);
kineticsTolerance = getOption(options, 'kineticsTolerance', 1e-8);
badStepMax = getOption(options, 'badStepMax', 5000);
mineralName = char(getOption(options, 'mineralName', 'Calcite'));
mineralFormula = char(getOption(options, 'mineralFormula', 'CaCO3'));

interfaceAreaCm2 = max(interfaceAreaCm2(:), 0);
waterVolumeCm3 = max(waterVolumeCm3(:), 0);
waterKg = max(waterVolumeCm3 * 1e-3, minWaterKg);
solutionWaterKg = max(solutionWaterKg, minSolutionWaterKg);
specificSurfaceArea = max((interfaceAreaCm2 * 1e-4) ./ waterKg, minSpecificSurfaceArea);
specificSurfaceArea = min(specificSurfaceArea, maxSpecificSurfaceArea);
calciteMoles = max(calciteMoles(:), minCalciteMoles);
activeKinetics = interfaceAreaCm2 > 0 & waterVolumeCm3 > 0 & calciteMoles > 0;
kineticsMoles = repmat(max(kineticsReservoirMoles, minKineticsMoles), numCells, 1);
kineticsM0 = repmat(max(kineticsReservoirMoles, minKineticsM0), numCells, 1);
kineticsSurfaceArea = specificSurfaceArea;
kineticsMoles(~activeKinetics) = minKineticsMoles;
kineticsM0(~activeKinetics) = minKineticsM0;
kineticsSurfaceArea(~activeKinetics) = minKineticsSurfaceArea;

lines = strings(0, 1);
lines(end + 1) = "TITLE RTSPHEM single-calcite PHREEQC coupling";

for iCell = 1:numCells
    hMolL = max(h(iCell) * 1000, minHForPH);
    lines(end + 1) = sprintf('SOLUTION %d', iCell);
    lines(end + 1) = sprintf('temp %.15g', temperatureC);
    lines(end + 1) = sprintf('pH %.15g', -log10(hMolL));
    lines(end + 1) = ['units ', solutionUnits];
    if writeSolutionWaterLine || abs(solutionWaterKg - 1) > eps
        lines(end + 1) = sprintf('water %.15g', solutionWaterKg);
    end
    lines(end + 1) = sprintf('Ca %.15g', max(ca(iCell), 0) * 1000);
    lines(end + 1) = sprintf('C %.15g', max(c(iCell), 0) * 1000);
    lines(end + 1) = sprintf('Na %.15g', max(na(iCell), 0) * 1000);
    lines(end + 1) = sprintf('Cl %.15g', max(cl(iCell), 0) * 1000);
    lines(end + 1) = sprintf('KINETICS %d', iCell);
    lines(end + 1) = mineralName;
    lines(end + 1) = sprintf('-formula  %s  1', mineralFormula);
    lines(end + 1) = sprintf('-m        %.15g', kineticsMoles(iCell));
    lines(end + 1) = sprintf('-m0       %.15g', kineticsM0(iCell));
    lines(end + 1) = sprintf('-parms    %.15g  %.15g', kineticsSurfaceArea(iCell), kineticsCorrectionFactor);
    lines(end + 1) = sprintf('-tol %.15g', kineticsTolerance);
    lines(end + 1) = sprintf('-bad_step_max %.15g', badStepMax);
end

lines(end + 1) = "SELECTED_OUTPUT";
lines(end + 1) = "-reset false";
lines(end + 1) = "-simulation true";
lines(end + 1) = "-state true";
lines(end + 1) = "-solution true";
lines(end + 1) = "-pH true";
lines(end + 1) = "-charge_balance true";
lines(end + 1) = "-totals Ca C Na Cl";
lines(end + 1) = "-molalities H+ Ca+2 HCO3- CO3-2 Cl- Na+";
lines(end + 1) = "-saturation_indices Calcite";
lines(end + 1) = "USER_PUNCH";
lines(end + 1) = "-headings KIN_DELTA_Calcite RATE_Calcite";
lines(end + 1) = "-start";
lines(end + 1) = sprintf('10 PUNCH KIN_DELTA("%s") -KIN_DELTA("%s")/KIN_TIME', mineralName, mineralName);
lines(end + 1) = "-end";
lines(end + 1) = "END";
lines(end + 1) = "RUN_CELLS";
lines(end + 1) = sprintf('-cells 1-%d', numCells);
lines(end + 1) = sprintf('-time_step %.15g', dt);
lines(end + 1) = "END";

text = char(strjoin(lines, newline));
end

function values = requireColumn(state, fieldName)
if ~isfield(state, fieldName)
    error('RTSPHEM:Phreeqc:MissingStateField', 'Missing PHREEQC state field: %s', fieldName);
end
values = state.(fieldName)(:);
end

function values = optionalColumn(state, fieldName, numCells, defaultValue)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    values = state.(fieldName)(:);
else
    values = repmat(defaultValue, numCells, 1);
end
if numel(values) ~= numCells
    error('RTSPHEM:Phreeqc:StateSizeMismatch', ...
        'State field %s has %d values, expected %d.', fieldName, numel(values), numCells);
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
