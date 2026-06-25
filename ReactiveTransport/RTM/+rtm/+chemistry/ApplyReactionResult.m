function [updated, ledger] = ApplyReactionResult(state, reactionResult)
%APPLYREACTIONRESULT Apply chemistry deltas to conserved state with ledger.

rtm.state.ValidateState(state);
componentDelta = requireDelta(reactionResult, 'component_delta_moles', ...
    size(state.component_moles));
mineralDelta = requireDelta(reactionResult, 'mineral_delta_moles', ...
    size(state.mineral_moles));

updated = state;
updated.component_moles = state.component_moles + componentDelta;
updated.mineral_moles = state.mineral_moles + mineralDelta;
rtm.state.ValidateState(updated);

ledger = struct();
ledger.component_names = state.component_names;
ledger.mineral_names = state.mineral_names;
ledger.initial_component_moles_total = sum(state.component_moles, 1);
ledger.final_component_moles_total = sum(updated.component_moles, 1);
ledger.component_delta_moles_total = sum(componentDelta, 1);
ledger.initial_mineral_moles_total = sum(state.mineral_moles, 1);
ledger.final_mineral_moles_total = sum(updated.mineral_moles, 1);
ledger.mineral_delta_moles_total = sum(mineralDelta, 1);
ledger.realized_interface_moles_total = realizedInterfaceTotal(reactionResult);
ledger.converged = logical(getFieldOrDefault(reactionResult, 'converged', true));
ledger.failed_cells = getFieldOrDefault(reactionResult, 'failed_cells', []);
ledger.error_message = string(getFieldOrDefault(reactionResult, 'error_message', ""));
end

function values = requireDelta(reactionResult, fieldName, expectedSize)
if ~isstruct(reactionResult) || ~isfield(reactionResult, fieldName)
    error('RTSPHEM:Chemistry:MissingReactionResultField', ...
        'reactionResult.%s is required.', fieldName);
end
values = reactionResult.(fieldName);
values = values(:,:);
if ~isequal(size(values), expectedSize)
    error('RTSPHEM:Chemistry:ReactionResultSizeMismatch', ...
        'reactionResult.%s has size [%s], expected [%s].', ...
        fieldName, num2str(size(values)), num2str(expectedSize));
end
if any(~isfinite(values(:)))
    error('RTSPHEM:Chemistry:NonfiniteReactionDelta', ...
        'reactionResult.%s must contain finite values.', fieldName);
end
end

function total = realizedInterfaceTotal(reactionResult)
if isfield(reactionResult, 'realized_interface_moles') && ...
        ~isempty(reactionResult.realized_interface_moles)
    total = sum(reactionResult.realized_interface_moles(:), 'omitnan');
else
    total = 0;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
