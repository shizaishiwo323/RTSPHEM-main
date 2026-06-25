function state = CreateConservedState(componentConcentrationMolCm3, geometry, componentNames, options)
% CreateConservedState - Build a conserved molar RTM state from concentrations.
%
% componentConcentrationMolCm3 is nCells-by-nComponents in mol/cm^3.
% The stored primary state is component_moles = C * water_volume_cm3.

if nargin < 4 || isempty(options)
    options = struct();
end

waterVolume = requireWaterVolume(geometry);
componentConcentrationMolCm3 = componentConcentrationMolCm3(:,:);
numCells = numel(waterVolume);
if size(componentConcentrationMolCm3, 1) ~= numCells
    error('RTSPHEM:State:StateSizeMismatch', ...
        'Concentration rows (%d) must match water volumes (%d).', ...
        size(componentConcentrationMolCm3, 1), numCells);
end
if numel(componentNames) ~= size(componentConcentrationMolCm3, 2)
    error('RTSPHEM:State:ComponentNameMismatch', ...
        'componentNames must match concentration columns.');
end

waterVolume = waterVolume(:);
validateConcentration(componentConcentrationMolCm3);
state = struct();
state.component_names = reshape(componentNames, 1, []);
state.component_moles = componentConcentrationMolCm3 .* waterVolume;
state.component_moles(waterVolume == 0, :) = 0;
state.mineral_names = getOption(options, 'mineralNames', {});
state.mineral_moles = getOption(options, 'mineralMoles', zeros(numCells, 0));
state.temperature_C = getOption(options, 'temperature_C', 25 * ones(numCells, 1));
state.pressure_atm = getOption(options, 'pressure_atm', ones(numCells, 1));
state.chemistry_aux = getOption(options, 'chemistry_aux', struct());
state.time_s = getOption(options, 'time_s', 0);

rtm.state.ValidateState(state);
end

function waterVolume = requireWaterVolume(geometry)
if ~isstruct(geometry) || ~isfield(geometry, 'water_volume_cm3')
    error('RTSPHEM:State:MissingWaterVolume', ...
        'geometry.water_volume_cm3 is required.');
end
waterVolume = geometry.water_volume_cm3(:);
if any(~isfinite(waterVolume))
    error('RTSPHEM:State:InvalidWaterVolume', ...
        'geometry.water_volume_cm3 must contain finite values.');
end
if any(waterVolume < 0)
    error('RTSPHEM:State:NegativeWaterVolume', ...
        'geometry.water_volume_cm3 must be nonnegative.');
end
end

function validateConcentration(concentration)
if any(~isfinite(concentration(:)))
    error('RTSPHEM:State:NonfiniteConcentration', ...
        'component concentrations must be finite.');
end
if any(concentration(:) < 0)
    error('RTSPHEM:State:NegativeConcentration', ...
        'component concentrations must be nonnegative.');
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
