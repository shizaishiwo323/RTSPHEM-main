function units = NormalizeRtmUnits(units)
%NORMALIZERTMUNITS Fill canonical RTM cgs + mol unit conventions.

if nargin < 1 || isempty(units)
    units = struct();
end
units.length = string(getFieldOrDefault(units, 'length', "cm"));
units.area = string(getFieldOrDefault(units, 'area', "cm^2"));
units.volume = string(getFieldOrDefault(units, 'volume', "cm^3"));
units.time = string(getFieldOrDefault(units, 'time', "s"));
units.component_state = string(getFieldOrDefault(units, ...
    'component_state', "mol"));
units.concentration = string(getFieldOrDefault(units, ...
    'concentration', "mol/cm^3"));
units.rate_per_area = string(getFieldOrDefault(units, ...
    'rate_per_area', "mol/cm^2/s"));
units.interface_rate = string(getFieldOrDefault(units, ...
    'interface_rate', units.rate_per_area));
units.cell_rate = string(getFieldOrDefault(units, 'cell_rate', "mol/s"));
units.velocity = string(getFieldOrDefault(units, 'velocity', "cm/s"));
units.rate_constant = string(getFieldOrDefault(units, ...
    'rate_constant', "cm/s"));
units.rate_law = string(getFieldOrDefault(units, ...
    'rate_law', "r_mol_cm2_s = k_cm_s * C_mol_cm3"));
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
