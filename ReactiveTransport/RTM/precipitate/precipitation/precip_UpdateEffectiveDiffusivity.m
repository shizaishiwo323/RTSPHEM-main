function deff = precip_UpdateEffectiveDiffusivity(state, spec)
% precip_UpdateEffectiveDiffusivity - Compute Yoon D_eff = D(1 - Vm)^n.
%
% Inputs:
%   state - Yoon state with Vm.
%   spec  - benchmark spec with diffusion coefficient and exponent.
%
% Output:
%   deff  - effective diffusivity field, cm2/s.

vm = min(max(state.Vm, 0), 1);
deff = spec.diffusionCoefficient_cm2_s .* (1 - vm) .^ spec.diffusionExponent;
if isfield(state, 'substrateMask')
    deff(state.substrateMask) = 0;
end
end
