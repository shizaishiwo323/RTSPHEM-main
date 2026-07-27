function row = AppendSeismoelectricChemistrySummary( ...
    csvFile, stepIndex, timeSeconds, rtState, speciesData)
% AppendSeismoelectricChemistrySummary - Export chemistry needed by SE models.
%
% Fluid quantities are pore-water-volume weighted. Interface quantities are
% interface-area weighted over cells that contain water. Zeta potential,
% high-frequency tortuosity, and viscous characteristic length are left NaN
% because PHREEQC aqueous speciation and hydraulic tortuosity do not determine
% those quantities.

waterVolume = requireVector(speciesData, 'water_volume_cm3');
numCells = numel(waterVolume);
interfaceArea = optionalVector(speciesData, 'interface_area_cm2', numCells, 0);
fluidMask = isfinite(waterVolume) & waterVolume > 0;
poreWeights = zeros(numCells, 1);
poreWeights(fluidMask) = waterVolume(fluidMask);
interfaceWeights = zeros(numCells, 1);
interfaceWeights(fluidMask) = max(interfaceArea(fluidMask), 0);

record = struct();
record.TimeStep = stepIndex;
record.Time_s = timeSeconds;
record.Porosity = scalarField(rtState, 'porosity');
record.Permeability_mD = scalarField(rtState, 'permeability_mD');
record.k_k0 = scalarField(rtState, 'k_k0');
record.SurfaceArea_cm2 = scalarField(rtState, 'surface_area_cm2');
record.GrainVolume_cm3 = scalarField(rtState, 'grain_volume_cm3');
record.OutletHConc = scalarField(rtState, 'outlet_h_conc_mol_cm3');
record.OutletHConc_unit = "mol_cm3";
record.HydraulicTortuosity = scalarField(rtState, 'hydraulic_tortuosity');
record.PoreWaterVolume_cm3 = sum(poreWeights);
record.InterfaceArea_cm2 = sum(interfaceWeights);
record.PoreCellCount = nnz(poreWeights > 0);
record.InterfaceCellCount = nnz(interfaceWeights > 0);

record.pH_pore_volume_weighted = weightedField(speciesData, 'pH', poreWeights);
record.pH_pore_min = finiteExtreme(speciesData, 'pH', poreWeights, 'min');
record.pH_pore_max = finiteExtreme(speciesData, 'pH', poreWeights, 'max');
record.pH_interface_area_weighted = weightedField(speciesData, 'pH', interfaceWeights);
record.IonicStrength_mol_kgw = weightedField( ...
    speciesData, 'ionicStrength_mol_kgw', poreWeights);
record.IonicStrength_interface_mol_kgw = weightedField( ...
    speciesData, 'ionicStrength_mol_kgw', interfaceWeights);
record.FluidConductivity_S_m = weightedField( ...
    speciesData, 'fluidConductivity_S_m', poreWeights);
record.FluidConductivity_interface_S_m = weightedField( ...
    speciesData, 'fluidConductivity_S_m', interfaceWeights);
record.PHREEQCPercentError_abs_max = maxFiniteAbs( ...
    speciesData, 'phreeqcPercentError', poreWeights);
record.ChargeBalance_abs_max_eq = maxFiniteAbs( ...
    speciesData, 'chargeBalance', poreWeights);

record.H_free_mol_L = 1000 * weightedField(speciesData, 'h_mol_cm3', poreWeights);
record.CaTotal_mol_L = 1000 * weightedField(speciesData, 'ca_total_mol_cm3', poreWeights);
record.CTotal_mol_L = 1000 * weightedField(speciesData, 'c_total_mol_cm3', poreWeights);
record.NaTotal_mol_L = 1000 * weightedField(speciesData, 'na_total_mol_cm3', poreWeights);
record.ClTotal_mol_L = 1000 * weightedField(speciesData, 'cl_total_mol_cm3', poreWeights);
record.Alkalinity_mol_L = 1000 * weightedField( ...
    speciesData, 'alkalinity_mol_cm3', poreWeights);
record.Ca2_mol_L = 1000 * weightedField(speciesData, 'ca_mol_cm3', poreWeights);
record.HCO3_mol_L = 1000 * weightedField(speciesData, 'hco3_mol_cm3', poreWeights);
record.CO3_mol_L = 1000 * weightedField(speciesData, 'co3_mol_cm3', poreWeights);
record.H_activity = weightedField(speciesData, 'h_activity_dimensionless', poreWeights);
record.Ca_activity = weightedField(speciesData, 'ca_activity_dimensionless', poreWeights);
record.HCO3_activity = weightedField(speciesData, 'hco3_activity_dimensionless', poreWeights);
record.CO3_activity = weightedField(speciesData, 'co3_activity_dimensionless', poreWeights);
record.Na_activity = weightedField(speciesData, 'na_activity_dimensionless', poreWeights);
record.Cl_activity = weightedField(speciesData, 'cl_activity_dimensionless', poreWeights);
record.CalciteSI_interface_area_weighted = weightedField( ...
    speciesData, 'calciteSI', interfaceWeights);

record.ZetaPotential_V = NaN;
record.Tortuosity = NaN;
record.HighFrequencyTortuosity = NaN;
record.ViscousCharacteristicLength_m = NaN;
record.UpperFluidConductivity_S_m = NaN;
record.ChemistryStatus = chemistryStatus(record);
record.ZetaModelStatus = "not_computed_requires_calcite_scm_or_measurement";
record.DynamicGeometryStatus = ...
    "not_computed_hydraulic_tortuosity_is_not_alpha_inf";
record.UpperFluidStatus = "not_computed_requires_separate_upper_fluid_chemistry";
record.SEInputStatus = "blocked_requires_zeta_and_alpha_inf_calibration";

row = struct2table(record, 'AsArray', true);
outputDir = fileparts(csvFile);
if ~isempty(outputDir) && exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end
% ponytail: one durable CSV row per exported step; batch only if I/O becomes measurable.
if exist(csvFile, 'file') == 2
    writetable(row, csvFile, 'WriteMode', 'append', 'WriteVariableNames', false);
else
    writetable(row, csvFile);
end
end

function values = requireVector(data, fieldName)
if ~isfield(data, fieldName) || isempty(data.(fieldName))
    error('RTSPHEM:SeismoelectricExport:MissingField', ...
        'Missing required chemistry field: %s.', fieldName);
end
values = data.(fieldName)(:);
end

function values = optionalVector(data, fieldName, numCells, defaultValue)
if isfield(data, fieldName) && ~isempty(data.(fieldName))
    values = data.(fieldName)(:);
else
    values = repmat(defaultValue, numCells, 1);
end
if numel(values) ~= numCells
    error('RTSPHEM:SeismoelectricExport:SizeMismatch', ...
        'Field %s has %d values; expected %d.', fieldName, numel(values), numCells);
end
end

function value = scalarField(data, fieldName)
if isfield(data, fieldName) && isscalar(data.(fieldName))
    value = double(data.(fieldName));
else
    value = NaN;
end
end

function value = weightedField(data, fieldName, weights)
values = optionalVector(data, fieldName, numel(weights), NaN);
valid = isfinite(values) & isfinite(weights) & weights > 0;
if ~any(valid)
    value = NaN;
else
    value = sum(values(valid) .* weights(valid)) / sum(weights(valid));
end
end

function value = finiteExtreme(data, fieldName, weights, mode)
values = optionalVector(data, fieldName, numel(weights), NaN);
values = values(isfinite(values) & weights > 0);
if isempty(values)
    value = NaN;
elseif strcmp(mode, 'min')
    value = min(values);
else
    value = max(values);
end
end

function value = maxFiniteAbs(data, fieldName, weights)
values = optionalVector(data, fieldName, numel(weights), NaN);
values = abs(values(isfinite(values) & weights > 0));
if isempty(values)
    value = NaN;
else
    value = max(values);
end
end

function status = chemistryStatus(record)
required = [record.pH_pore_volume_weighted, record.IonicStrength_mol_kgw, ...
    record.FluidConductivity_S_m, record.PHREEQCPercentError_abs_max];
if record.PoreWaterVolume_cm3 <= 0 || any(~isfinite(required))
    status = "incomplete";
elseif record.FluidConductivity_S_m <= 0 || record.IonicStrength_mol_kgw < 0
    status = "invalid";
elseif record.PHREEQCPercentError_abs_max > 1e-3
    status = "charge_balance_warning";
else
    status = "ready_for_calibrated_zeta_model";
end
end
