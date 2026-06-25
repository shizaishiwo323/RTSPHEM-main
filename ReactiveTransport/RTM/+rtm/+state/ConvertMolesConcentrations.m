function values = ConvertMolesConcentrations(stateOrConcentration, geometry, mode)
% ConvertMolesConcentrations - Convert between conserved moles and concentration.

if nargin < 3 || isempty(mode)
    mode = 'moles_to_concentration';
end
mode = lower(strrep(strtrim(char(mode)), '-', '_'));
waterVolume = requireWaterVolume(geometry);

switch mode
    case {'moles_to_concentration', 'moles_to_concentrations'}
        state = stateOrConcentration;
        rtm.state.ValidateState(state);
        componentMoles = state.component_moles;
        values = zeros(size(componentMoles));
        active = waterVolume(:) > 0;
        values(active, :) = componentMoles(active, :) ./ waterVolume(active);
    case {'concentration_to_moles', 'concentrations_to_moles'}
        concentration = stateOrConcentration(:,:);
        if size(concentration, 1) ~= numel(waterVolume)
            error('RTSPHEM:State:StateSizeMismatch', ...
                'Concentration rows (%d) must match water volumes (%d).', ...
                size(concentration, 1), numel(waterVolume));
        end
        validateConcentration(concentration);
        values = concentration .* waterVolume(:);
        values(waterVolume(:) <= 0, :) = 0;
    otherwise
        error('RTSPHEM:State:UnknownConversionMode', ...
            'Unknown conversion mode: %s.', char(mode));
end
if any(~isfinite(values(:)))
    error('RTSPHEM:State:NonfiniteConversionResult', ...
        'Converted state values must be finite.');
end
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
