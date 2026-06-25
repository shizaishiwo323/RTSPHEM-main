function [updated, diagnostics] = precip_UpdateFlowMask(state, spec)
% precip_UpdateFlowMask - Update Yoon blocked-cell mask from Vm threshold.
%
% Inputs:
%   state - Yoon state with Vm, blockedMask, and optional substrateMask.
%   spec  - benchmark spec with blockedVmThreshold.
%
% Output:
%   updated     - state with refreshed blockedMask.
%   diagnostics - counts and topology-change flag.

if ~isfield(state, 'Vm')
    error('RTSPHEM:Precipitate:MissingVm', 'state.Vm is required.');
end
previousMask = false(size(state.Vm));
if isfield(state, 'blockedMask') && ~isempty(state.blockedMask)
    previousMask = logical(state.blockedMask);
end

newMask = state.Vm >= spec.blockedVmThreshold;
if isfield(state, 'substrateMask')
    newMask = newMask & ~state.substrateMask;
end

updated = state;
updated.blockedMask = newMask;
diagnostics = struct();
diagnostics.topologyChanged = any(newMask(:) ~= previousMask(:));
diagnostics.numBlockedCells = nnz(newMask);
diagnostics.numNewBlockedCells = nnz(newMask & ~previousMask);
diagnostics.blockedVmThreshold = spec.blockedVmThreshold;
end
