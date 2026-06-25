function state = precip_UpdateCase5BoostState(state, spec, currentTime_s)
% precip_UpdateCase5BoostState - Activate Case 5 boost after center clogging.
%
% The Case 5 dissolution multiplier is delayed until a blocked cell appears in
% the center mixing band. Once activated, it remains active for later steps.

if nargin < 3
    currentTime_s = NaN;
end
if ~isfield(state, 'case5BoostActive') || isempty(state.case5BoostActive)
    state.case5BoostActive = false;
end
if ~isfield(state, 'case5ActivationTime_s') || ...
        isempty(state.case5ActivationTime_s)
    state.case5ActivationTime_s = NaN;
end
if state.case5BoostActive || getFieldOrDefault(spec, 'dissolutionFactor', 1) <= 1
    return;
end
if ~isfield(state, 'blockedMask') || isempty(state.blockedMask)
    return;
end

centerBandMask = resolveCenterBandMask(state, spec);
if any(state.blockedMask(centerBandMask))
    state.case5BoostActive = true;
    state.case5ActivationTime_s = currentTime_s;
end
end

function centerBandMask = resolveCenterBandMask(state, spec)
if isfield(spec, 'centerBandMask') && ~isempty(spec.centerBandMask)
    centerBandMask = logical(spec.centerBandMask);
    if ~isequal(size(centerBandMask), size(state.blockedMask))
        error('RTSPHEM:Precipitate:InvalidCenterBandMask', ...
            'spec.centerBandMask must match state.blockedMask size.');
    end
    return;
end

yCenters = ((1:size(state.blockedMask, 1))' - 0.5) .* spec.dy_cm;
halfWidth = getFieldOrDefault(spec, 'centerBandHalfWidth_cm', spec.dy_cm);
rowMask = abs(yCenters - spec.splitInletY_cm) <= halfWidth;
centerBandMask = repmat(rowMask, 1, size(state.blockedMask, 2));
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
