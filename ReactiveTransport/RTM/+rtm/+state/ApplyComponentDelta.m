function updated = ApplyComponentDelta(state, deltaComponentNames, deltaMoles)
% ApplyComponentDelta - Add named component mole increments to conserved state.

rtm.state.ValidateState(state);
deltaMoles = deltaMoles(:,:);
numCells = size(state.component_moles, 1);
if size(deltaMoles, 1) ~= numCells
    error('RTSPHEM:State:StateSizeMismatch', ...
        'deltaMoles rows (%d) must match state cells (%d).', ...
        size(deltaMoles, 1), numCells);
end
if numel(deltaComponentNames) ~= size(deltaMoles, 2)
    error('RTSPHEM:State:ComponentNameMismatch', ...
        'deltaComponentNames must match deltaMoles columns.');
end

updated = state;
for iDelta = 1:numel(deltaComponentNames)
    componentName = char(deltaComponentNames{iDelta});
    idx = find(strcmp(updated.component_names, componentName), 1);
    if isempty(idx)
        error('RTSPHEM:State:UnknownComponent', ...
            'Unknown component in delta: %s.', componentName);
    end
    updated.component_moles(:, idx) = updated.component_moles(:, idx) + ...
        deltaMoles(:, iDelta);
end

rtm.state.ValidateState(updated);
end
