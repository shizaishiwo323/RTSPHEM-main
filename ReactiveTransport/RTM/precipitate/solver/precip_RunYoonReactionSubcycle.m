function [state, diagnostics] = precip_RunYoonReactionSubcycle(state, spec, totalDtS, options)
% precip_RunYoonReactionSubcycle - Advance Yoon reaction with Vm-limited substeps.

if nargin < 4
    options = struct();
end

remainingDt = totalDtS;
advancedDt = 0;
acceptedSubstepCount = 0;
lastAcceptedDt = 0;
maxAcceptedVmChange = 0;
stableDtMin = Inf;
limiter = "none";

while remainingDt > max(10 * eps(totalDtS), 1e-15)
    samples = stateComponentsToSamples(state, spec);
    [chem, chemistryMetadata] = precip_SpeciateYoonComponents(samples, spec);
    rate = reshape(precip_YoonVateriteRate(chem, spec), size(state.Vm));
    area = precip_ComputeYoonReactiveArea(state, spec);
    [stableDt, stableDiagnostics] = precip_ComputeYoonReactionStableDt( ...
        state, rate, area, spec);
    attemptedDt = min(remainingDt, stableDt);
    if ~isfinite(attemptedDt) || attemptedDt <= 0
        attemptedDt = remainingDt;
    end

    beforeVm = state.Vm;
    state = precip_AdvanceReaction(state, spec, attemptedDt);
    acceptedSubstepCount = acceptedSubstepCount + 1;
    lastAcceptedDt = attemptedDt;
    advancedDt = advancedDt + attemptedDt;
    remainingDt = max(totalDtS - advancedDt, 0);
    maxAcceptedVmChange = max(maxAcceptedVmChange, ...
        max(abs(state.Vm(:) - beforeVm(:))));
    stableDtMin = min(stableDtMin, stableDt);
    if stableDiagnostics.limiter ~= "none"
        limiter = stableDiagnostics.limiter;
    end
    state.chemistryBackend = chemistryMetadata.backend;
    state.chemistryBackendMetadata = chemistryMetadata;
end

diagnostics = struct();
diagnostics.requestedDt_s = totalDtS;
diagnostics.acceptedDt_s = lastAcceptedDt;
diagnostics.totalAdvancedDt_s = advancedDt;
diagnostics.acceptedSubstepCount = acceptedSubstepCount;
diagnostics.maxAcceptedVmChange = maxAcceptedVmChange;
diagnostics.stableDt_s = stableDtMin;
diagnostics.stabilityLimiter = limiter;
end

function samples = stateComponentsToSamples(state, spec)
samples = struct();
samples.fixedPH = nan(numel(state.Vm), 1);
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    samples.(fieldName) = state.components.(fieldName)(:);
end
end
