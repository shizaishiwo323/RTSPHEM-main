function result = precip_RunTransportSubcycle(initialComponents, spec, dt_s, candidateFcn, options)
% precip_RunTransportSubcycle - Reject/shrink/retry conservative transport.
%
% Inputs:
%   initialComponents - component concentration fields.
%   spec              - benchmark spec with componentNames.
%   dt_s              - requested transport substep.
%   candidateFcn      - function handle @(components, dt_s) trialComponents.
%   options           - shrinkFactor, maxRejectedSteps, massTolerance.
%
% Output:
%   result            - accepted components and rejection diagnostics.

if nargin < 5
    options = struct();
end
if dt_s <= 0
    error('RTSPHEM:Precipitate:InvalidTransportSubstep', ...
        'Transport substep must be positive.');
end
if ~isa(candidateFcn, 'function_handle')
    error('RTSPHEM:Precipitate:InvalidTransportCandidate', ...
        'candidateFcn must be a function handle.');
end

shrinkFactor = getFieldOrDefault(options, 'shrinkFactor', 0.5);
maxRejectedSteps = getFieldOrDefault(options, 'maxRejectedSteps', 8);
massTolerance = getFieldOrDefault(options, 'massTolerance', 1e-4);
if shrinkFactor <= 0 || shrinkFactor >= 1
    error('RTSPHEM:Precipitate:InvalidTransportShrinkFactor', ...
        'shrinkFactor must be between 0 and 1.');
end

trialDt = dt_s;
rejectedStepCount = 0;
diagnostics = repmat(makeDiagnostic(NaN, false, "not_attempted", NaN), ...
    maxRejectedSteps + 1, 1);

for iAttempt = 1:(maxRejectedSteps + 1)
    trialComponents = candidateFcn(initialComponents, trialDt);
    audit = auditTransportTrial(initialComponents, trialComponents, spec, massTolerance);
    diagnostics(iAttempt) = makeDiagnostic(trialDt, audit.accepted, ...
        audit.reason, audit.maxRelativeMassError);
    if audit.accepted
        result = struct();
        result.components = trialComponents;
        result.acceptedDt_s = trialDt;
        result.rejectedStepCount = rejectedStepCount;
        result.diagnostics = diagnostics(1:iAttempt);
        return;
    end
    rejectedStepCount = rejectedStepCount + 1;
    trialDt = trialDt * shrinkFactor;
end

error('RTSPHEM:Precipitate:TransportSubcycleFailed', ...
    'Transport subcycle failed after %d rejected attempts.', maxRejectedSteps);
end

function audit = auditTransportTrial(initialComponents, trialComponents, spec, massTolerance)
audit = struct();
audit.accepted = true;
audit.reason = "accepted";
audit.maxRelativeMassError = 0;

nonnegativeComponents = {'Ca_total', 'C_total', 'Na_total', 'Cl_total'};
for iComponent = 1:numel(nonnegativeComponents)
    fieldName = nonnegativeComponents{iComponent};
    if min(trialComponents.(fieldName)(:)) < -1e-15
        audit.accepted = false;
        audit.reason = "negative_component";
        return;
    end
end

ledger = precip_ComputeComponentMassLedger(initialComponents, trialComponents, spec);
relativeErrors = zeros(numel(spec.componentNames), 1);
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    relativeErrors(iComponent) = abs(ledger.relative.(fieldName));
end
audit.maxRelativeMassError = max(relativeErrors);
if audit.maxRelativeMassError > massTolerance
    audit.accepted = false;
    audit.reason = "mass_balance_drift";
end
end

function diagnostic = makeDiagnostic(dt_s, accepted, reason, maxRelativeMassError)
diagnostic = struct();
diagnostic.dt_s = dt_s;
diagnostic.accepted = accepted;
diagnostic.reason = reason;
diagnostic.maxRelativeMassError = maxRelativeMassError;
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
