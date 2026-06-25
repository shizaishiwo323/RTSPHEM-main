function [stableDt, diagnostics] = precip_ComputeMicrocontinuumStableDt(state, spec, options)
% precip_ComputeMicrocontinuumStableDt - State-aware transport CFL wrapper.

if nargin < 3
    options = struct();
end
if isfield(state, 'substrateMask') && ~isempty(state.substrateMask)
    options.substrateMask = state.substrateMask;
end
if isfield(state, 'blockedMask') && ~isempty(state.blockedMask)
    options.blockedMask = state.blockedMask;
end
if isfield(state, 'effectiveDiffusivity_cm2_s') && ...
        ~isempty(state.effectiveDiffusivity_cm2_s)
    options.effectiveDiffusivity_cm2_s = state.effectiveDiffusivity_cm2_s;
end
[stableDt, diagnostics] = precip_ComputeTransportStableDt(spec, options);
end
