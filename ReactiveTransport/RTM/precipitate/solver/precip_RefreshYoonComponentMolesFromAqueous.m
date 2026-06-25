function state = precip_RefreshYoonComponentMolesFromAqueous(state, spec)
% precip_RefreshYoonComponentMolesFromAqueous - Store aqueous fields as moles.
%
% The Yoon micro-continuum state keeps componentMoles as the conservative
% inventory. components and aqueousConcentration remain water-phase
% concentrations for chemistry, output, and older callers.

state = ensureFluidVolumeFraction(state);
waterVolumeCm3 = max(state.fluidVolumeFraction, 0) .* spec.cellVolume_cm3;
state.componentMoles = struct();
state.componentMoles.names = spec.componentNames;
state.aqueousConcentration = struct();
state.aqueousConcentration.names = spec.componentNames;

for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    concentration = state.components.(fieldName);
    state.aqueousConcentration.(fieldName) = concentration;
    state.componentMoles.(fieldName) = concentration .* waterVolumeCm3;
end
end

function state = ensureFluidVolumeFraction(state)
if ~isfield(state, 'fluidVolumeFraction') || isempty(state.fluidVolumeFraction)
    state.fluidVolumeFraction = double(~state.substrateMask) .* (1 - state.Vm);
end
end
