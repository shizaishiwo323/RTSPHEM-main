function diagnostics = ValidateAcceptedStep(stepInfo, options)
%VALIDATEACCEPTEDSTEP Evaluate mass, geometry, and chemistry acceptance gates.

if nargin < 2 || isempty(options)
    options = struct();
end

diagnostics = struct();
diagnostics.accepted = true;
diagnostics.reasons = strings(0, 1);

[componentAbsResidual, componentRelResidual] = componentResidual(stepInfo);
diagnostics.component_absolute_residual_moles = componentAbsResidual;
diagnostics.component_relative_residual = componentRelResidual;
try
    diagnostics.full_mass_ledger = ...
        rtm.diagnostics.ComputeFullMassLedger(stepInfo);
catch
    diagnostics.full_mass_ledger = struct();
end

[solidAbsResidual, solidRelResidual] = solidVolumeResidual(stepInfo);
diagnostics.solid_volume_absolute_residual_cm3 = solidAbsResidual;
diagnostics.solid_volume_relative_residual = solidRelResidual;
diagnostics.mineral_absolute_residual_moles = 0;
diagnostics.mineral_relative_residual = 0;
chargeAbsResidual = chargeResidual(stepInfo);
diagnostics.charge_absolute_residual_eq = chargeAbsResidual;

massAbsTol = getOption(options, 'mass_absolute_tolerance_mol', 1e-14);
massRelTol = getOption(options, 'mass_relative_tolerance', 1e-8);
solidAbsTol = getOption(options, 'solid_absolute_tolerance_cm3', 1e-14);
solidRelTol = getOption(options, 'solid_relative_tolerance', 1e-8);
chargeAbsTol = getOption(options, 'charge_absolute_tolerance_eq', Inf);
maxDisplacement = getOption(options, 'max_displacement_over_h', 0.25);

diagnostics = evaluateFullCoupledLedger(diagnostics, options);

if componentAbsResidual > massAbsTol && componentRelResidual > massRelTol
    diagnostics = reject(diagnostics, "component mass residual exceeds tolerance");
end
if solidAbsResidual > solidAbsTol && solidRelResidual > solidRelTol
    diagnostics = reject(diagnostics, "solid volume residual exceeds tolerance");
end
if chargeAbsResidual > chargeAbsTol
    diagnostics = reject(diagnostics, ...
        "chemistry charge balance residual exceeds tolerance");
end

[chemistryConverged, failedCells, errorMessage] = chemistryStatus(stepInfo);
diagnostics.failed_cells = failedCells;
diagnostics.error_message = errorMessage;
if ~chemistryConverged
    diagnostics = reject(diagnostics, "chemistry backend failed");
end
if ~isempty(failedCells)
    diagnostics = reject(diagnostics, "chemistry has failed cells");
end

geometryDisplacement = getNestedField(stepInfo, {'geometry', 'max_displacement_over_h'}, 0);
diagnostics.max_displacement_over_h = geometryDisplacement;
if geometryDisplacement > maxDisplacement
    diagnostics = reject(diagnostics, "geometry displacement exceeds tolerance");
end

geometryAccepted = getNestedField(stepInfo, {'geometry', 'accepted'}, true);
if ~logical(geometryAccepted)
    reason = string(getNestedField(stepInfo, {'geometry', 'reject_reason'}, ...
        "geometry backend rejected step"));
    if strlength(reason) == 0
        reason = "geometry backend rejected step";
    end
    diagnostics = reject(diagnostics, reason);
end

diagnostics = evaluateFlowDiagnostics(diagnostics, stepInfo);
end

function diagnostics = evaluateFullCoupledLedger(diagnostics, options)
if ~isfield(diagnostics, 'full_mass_ledger') || ...
        ~isstruct(diagnostics.full_mass_ledger) || ...
        isempty(fieldnames(diagnostics.full_mass_ledger))
    return;
end
fullDiagnostics = rtm.diagnostics.ValidateFullCoupledStep( ...
    diagnostics.full_mass_ledger, options);
diagnostics.full_mass_diagnostics = fullDiagnostics;
diagnostics.component_absolute_residual_moles = max( ...
    diagnostics.component_absolute_residual_moles, ...
    fullDiagnostics.component_absolute_residual_moles);
diagnostics.component_relative_residual = max( ...
    diagnostics.component_relative_residual, ...
    fullDiagnostics.component_relative_residual);
diagnostics.solid_volume_absolute_residual_cm3 = max( ...
    diagnostics.solid_volume_absolute_residual_cm3, ...
    fullDiagnostics.solid_volume_absolute_residual_cm3);
diagnostics.solid_volume_relative_residual = max( ...
    diagnostics.solid_volume_relative_residual, ...
    fullDiagnostics.solid_volume_relative_residual);
diagnostics.mineral_absolute_residual_moles = max( ...
    diagnostics.mineral_absolute_residual_moles, ...
    fullDiagnostics.mineral_absolute_residual_moles);
diagnostics.mineral_relative_residual = max( ...
    diagnostics.mineral_relative_residual, ...
    fullDiagnostics.mineral_relative_residual);
if fullDiagnostics.accepted
    return;
end
for iReason = 1:numel(fullDiagnostics.reasons)
    diagnostics = reject(diagnostics, fullDiagnostics.reasons(iReason));
end
end

function [absoluteResidual, relativeResidual] = componentResidual(stepInfo)
mass = requireField(stepInfo, 'mass');
residual = getRequiredField(mass, 'component_residual_moles');
initialTotal = getRequiredField(mass, 'initial_component_moles_total');
finalTotal = getRequiredField(mass, 'final_component_moles_total');

absoluteResidual = max(abs(residual(:)));
scale = max([abs(initialTotal(:)); abs(finalTotal(:)); 1e-300]);
relativeResidual = absoluteResidual ./ max(scale, eps);
end

function [absoluteResidual, relativeResidual] = solidVolumeResidual(stepInfo)
geometry = requireField(stepInfo, 'geometry');
if isfield(geometry, 'expected_solid_volume_change_cm3')
    expectedChange = geometry.expected_solid_volume_change_cm3;
else
    dissolvedMoles = getRequiredField(geometry, 'mineral_dissolved_moles');
    molarVolume = getRequiredField(geometry, 'molar_volume_cm3_mol');
    expectedChange = -dissolvedMoles .* molarVolume;
end
actualChange = getRequiredField(geometry, 'actual_solid_volume_change_cm3');
absoluteResidual = abs(actualChange - expectedChange);
solidBefore = getRequiredField(geometry, 'solid_volume_before_cm3');
solidAfter = getRequiredField(geometry, 'solid_volume_after_cm3');
scale = max([abs(solidBefore), abs(solidAfter), abs(expectedChange), 1e-300]);
relativeResidual = absoluteResidual ./ max(scale, eps);
end

function absoluteResidual = chargeResidual(stepInfo)
if ~isfield(stepInfo, 'chemistry') || isempty(stepInfo.chemistry) || ...
        ~isfield(stepInfo.chemistry, 'charge_balance_residual_eq') || ...
        isempty(stepInfo.chemistry.charge_balance_residual_eq)
    absoluteResidual = 0;
    return;
end
values = stepInfo.chemistry.charge_balance_residual_eq(:);
values = values(isfinite(values));
if isempty(values)
    absoluteResidual = 0;
else
    absoluteResidual = max(abs(values));
end
end

function [converged, failedCells, errorMessage] = chemistryStatus(stepInfo)
if ~isfield(stepInfo, 'chemistry') || isempty(stepInfo.chemistry)
    converged = true;
    failedCells = [];
    errorMessage = "";
    return;
end
chemistry = stepInfo.chemistry;
converged = getFieldOrDefault(chemistry, 'converged', true);
failedCells = getFieldOrDefault(chemistry, 'failed_cells', []);
errorMessage = string(getFieldOrDefault(chemistry, 'error_message', ""));
end

function diagnostics = evaluateFlowDiagnostics(diagnostics, stepInfo)
flow = getNestedField(stepInfo, {'flow'}, struct());
if ~isstruct(flow) || isempty(fieldnames(flow))
    diagnostics.flow_inlet_outlet_relative_residual = NaN;
    diagnostics.flow_max_abs_cell_divergence_cm3_s = NaN;
    diagnostics.flow_global_residual_cm3_s = NaN;
    return;
end

diagnostics.flow_inlet_outlet_relative_residual = getFieldOrDefault( ...
    flow, 'inlet_outlet_relative_residual', NaN);
diagnostics.flow_max_abs_cell_divergence_cm3_s = getFieldOrDefault( ...
    flow, 'max_abs_cell_divergence_cm3_s', NaN);
diagnostics.flow_global_residual_cm3_s = getFieldOrDefault( ...
    flow, 'global_residual_cm3_s', NaN);

flowAccepted = logical(getFieldOrDefault(flow, 'accepted', true));
flowReasons = string(getFieldOrDefault(flow, 'failure_reasons', strings(0, 1)));
flowReasons = flowReasons(:);
flowReasons = flowReasons(strlength(strtrim(flowReasons)) > 0);
if ~flowAccepted
    diagnostics = reject(diagnostics, "flow diagnostics failed");
end
for iReason = 1:numel(flowReasons)
    diagnostics = reject(diagnostics, flowReasons(iReason));
end
end

function diagnostics = reject(diagnostics, reason)
diagnostics.accepted = false;
diagnostics.reasons(end + 1, 1) = reason;
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function value = requireField(structValue, fieldName)
if ~isfield(structValue, fieldName)
    error('RTSPHEM:Diagnostics:MissingField', ...
        'Missing required field: %s.', fieldName);
end
value = structValue.(fieldName);
end

function value = getRequiredField(structValue, fieldName)
if ~isfield(structValue, fieldName)
    error('RTSPHEM:Diagnostics:MissingField', ...
        'Missing required field: %s.', fieldName);
end
value = structValue.(fieldName);
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function value = getNestedField(structValue, path, defaultValue)
value = structValue;
for iPath = 1:numel(path)
    if ~isstruct(value) || ~isfield(value, path{iPath})
        value = defaultValue;
        return;
    end
    value = value.(path{iPath});
end
end
