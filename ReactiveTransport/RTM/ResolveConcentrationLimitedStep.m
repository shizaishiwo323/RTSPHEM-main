function [nextStepSize, wasLimited, overshootRatio, activeThreshold] = ...
    ResolveConcentrationLimitedStep(nextStepSize, cMax, inletConcentration, ...
    overshootThreshold, hardOvershootThreshold, shrinkFactor, minStepSize, targetSliceMode)
% ResolveConcentrationLimitedStep applies concentration overshoot dt limiting.
% In target-slice mode, mild bounded overshoot should not freeze the adaptive
% porosity stepper; only severe overshoot should force shrinkage.

overshootRatio = cMax / max(inletConcentration, eps);
if targetSliceMode
    activeThreshold = max(overshootThreshold, hardOvershootThreshold);
else
    activeThreshold = overshootThreshold;
end

wasLimited = overshootRatio > activeThreshold;
if wasLimited
    nextStepSize = max(minStepSize, nextStepSize * shrinkFactor);
end
end
