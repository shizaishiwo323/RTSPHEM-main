function areaCm2 = precip_ComputeYoonReactiveArea(state, spec)
% precip_ComputeYoonReactiveArea - Yoon top/bottom reactive area per cell.
%
% Inputs:
%   state - Yoon state with substrateMask.
%   spec  - benchmark spec with dx/dy.
%
% Output:
%   areaCm2 - per-cell reactive area, cm2.

if ~isfield(state, 'substrateMask')
    error('RTSPHEM:Precipitate:MissingSubstrateMask', ...
        'state.substrateMask is required.');
end

areaCm2 = ones(size(state.substrateMask)) .* (2 * spec.dx_cm * spec.dy_cm);
areaCm2(state.substrateMask) = 0;
end
