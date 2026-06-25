function diagnostics = ValidateFullCoupledStep(stepInfoOrLedger, options)
%VALIDATEFULLCOUPLEDSTEP Validate full coupled mass and geometry ledgers.

if nargin < 2 || isempty(options)
    options = struct();
end
if isstruct(stepInfoOrLedger) && isfield(stepInfoOrLedger, ...
        'component_residual_moles') && isfield(stepInfoOrLedger, ...
        'solid_volume_residual_cm3')
    ledger = stepInfoOrLedger;
else
    ledger = rtm.diagnostics.ComputeFullMassLedger(stepInfoOrLedger);
end

diagnostics = struct();
diagnostics.accepted = true;
diagnostics.reasons = strings(0, 1);
diagnostics.ledger = ledger;
diagnostics.component_absolute_residual_moles = ...
    ledger.max_abs_component_residual_moles;
diagnostics.component_relative_residual = ledger.component_relative_residual;
diagnostics.solid_volume_absolute_residual_cm3 = ...
    abs(ledger.solid_volume_residual_cm3);
diagnostics.solid_volume_relative_residual = ...
    ledger.solid_volume_relative_residual;
diagnostics.mineral_absolute_residual_moles = ...
    getLedgerScalar(ledger, 'max_abs_mineral_residual_moles', 0);
diagnostics.mineral_relative_residual = ...
    getLedgerScalar(ledger, 'mineral_relative_residual', 0);

massAbsTol = getOption(options, 'mass_absolute_tolerance_mol', 1e-14);
massRelTol = getOption(options, 'mass_relative_tolerance', 1e-8);
solidAbsTol = getOption(options, 'solid_absolute_tolerance_cm3', 1e-14);
solidRelTol = getOption(options, 'solid_relative_tolerance', 1e-8);
mineralAbsTol = getOption(options, 'mineral_absolute_tolerance_mol', massAbsTol);
mineralRelTol = getOption(options, 'mineral_relative_tolerance', massRelTol);

if diagnostics.component_absolute_residual_moles > massAbsTol && ...
        diagnostics.component_relative_residual > massRelTol
    diagnostics = reject(diagnostics, ...
        "component mass residual exceeds tolerance");
end
if diagnostics.solid_volume_absolute_residual_cm3 > solidAbsTol && ...
        diagnostics.solid_volume_relative_residual > solidRelTol
    diagnostics = reject(diagnostics, ...
        "solid volume residual exceeds tolerance");
end
if diagnostics.mineral_absolute_residual_moles > mineralAbsTol && ...
        diagnostics.mineral_relative_residual > mineralRelTol
    diagnostics = reject(diagnostics, ...
        "mineral inventory residual exceeds tolerance");
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

function value = getLedgerScalar(ledger, fieldName, defaultValue)
if isstruct(ledger) && isfield(ledger, fieldName) && ~isempty(ledger.(fieldName))
    values = ledger.(fieldName);
    value = values(1);
else
    value = defaultValue;
end
if ~isfinite(value)
    value = defaultValue;
end
end
