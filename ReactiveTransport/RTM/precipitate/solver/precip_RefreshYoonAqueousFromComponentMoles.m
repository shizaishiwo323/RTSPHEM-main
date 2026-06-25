function state = precip_RefreshYoonAqueousFromComponentMoles(state, spec)
% precip_RefreshYoonAqueousFromComponentMoles - Recover aqueous concentration.

state = ensureFluidVolumeFraction(state);
waterVolumeCm3 = max(state.fluidVolumeFraction, 0) .* spec.cellVolume_cm3;
minWaterVolumeCm3 = getFieldOrDefault(spec, 'minWaterVolume_cm3', 1e-30);
aqueousMask = waterVolumeCm3 > minWaterVolumeCm3;

state.components = struct();
state.components.names = spec.componentNames;
state.aqueousConcentration = struct();
state.aqueousConcentration.names = spec.componentNames;

for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    concentration = zeros(size(waterVolumeCm3));
    concentration(aqueousMask) = state.componentMoles.(fieldName)(aqueousMask) ./ ...
        waterVolumeCm3(aqueousMask);
    state.aqueousConcentration.(fieldName) = concentration;
    state.components.(fieldName) = concentration;
end
end

function state = ensureFluidVolumeFraction(state)
if ~isfield(state, 'fluidVolumeFraction') || isempty(state.fluidVolumeFraction)
    state.fluidVolumeFraction = double(~state.substrateMask) .* (1 - state.Vm);
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
