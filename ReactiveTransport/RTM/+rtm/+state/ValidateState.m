function ValidateState(state)
% ValidateState - Check conserved RTM state structure invariants.

requiredFields = {'component_names', 'component_moles', 'mineral_names', ...
    'mineral_moles', 'temperature_C', 'pressure_atm', 'time_s'};
for iField = 1:numel(requiredFields)
    if ~isfield(state, requiredFields{iField})
        error('RTSPHEM:State:MissingField', ...
            'Missing state field: %s.', requiredFields{iField});
    end
end

componentMoles = state.component_moles;
if ~isnumeric(componentMoles) || ~ismatrix(componentMoles)
    error('RTSPHEM:State:InvalidComponentMoles', ...
        'component_moles must be a numeric matrix.');
end
if numel(state.component_names) ~= size(componentMoles, 2)
    error('RTSPHEM:State:ComponentNameMismatch', ...
        'component_names must match component_moles columns.');
end
if any(componentMoles(:) < 0)
    error('RTSPHEM:State:NegativeComponentMoles', ...
        'component_moles must be nonnegative.');
end
if any(~isfinite(componentMoles(:)))
    error('RTSPHEM:State:NonfiniteComponentMoles', ...
        'component_moles must be finite.');
end

mineralMoles = state.mineral_moles;
if ~isnumeric(mineralMoles) || ~ismatrix(mineralMoles)
    error('RTSPHEM:State:InvalidMineralMoles', ...
        'mineral_moles must be a numeric matrix.');
end
if size(mineralMoles, 1) ~= size(componentMoles, 1)
    error('RTSPHEM:State:StateSizeMismatch', ...
        'mineral_moles rows must match component_moles rows.');
end
if numel(state.mineral_names) ~= size(mineralMoles, 2)
    error('RTSPHEM:State:MineralNameMismatch', ...
        'mineral_names must match mineral_moles columns.');
end
if any(mineralMoles(:) < 0)
    error('RTSPHEM:State:NegativeMineralMoles', ...
        'mineral_moles must be nonnegative.');
end

numCells = size(componentMoles, 1);
if numel(state.temperature_C) ~= numCells || numel(state.pressure_atm) ~= numCells
    error('RTSPHEM:State:StateSizeMismatch', ...
        'temperature_C and pressure_atm must match state cell count.');
end
end
