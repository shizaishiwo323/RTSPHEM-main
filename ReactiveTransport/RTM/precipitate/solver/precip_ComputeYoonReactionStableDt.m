function [stableDt, diagnostics] = precip_ComputeYoonReactionStableDt(state, rate, areaCm2, spec)
% precip_ComputeYoonReactionStableDt - Limit reaction step by Vm change.

if isfield(spec, 'maxVmChangePerStep') && isfinite(spec.maxVmChangePerStep)
    maxVmChange = max(spec.maxVmChangePerStep, 0);
else
    maxVmChange = Inf;
end

vmRate = spec.vateriteMolarVolume_cm3_mol .* rate .* areaCm2 ./ ...
    spec.cellVolume_cm3;
maxAbsVmRate = max(abs(vmRate(:)));
if maxAbsVmRate <= 0 || ~isfinite(maxAbsVmRate) || isinf(maxVmChange)
    stableDt = Inf;
    limiter = "none";
else
    stableDt = maxVmChange ./ maxAbsVmRate;
    limiter = "max_vm_change";
end

diagnostics = struct();
diagnostics.maxVmChange = maxVmChange;
diagnostics.maxAbsVmRate_s = maxAbsVmRate;
diagnostics.limiter = limiter;
end
