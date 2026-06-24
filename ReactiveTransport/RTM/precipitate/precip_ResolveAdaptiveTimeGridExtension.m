function tf = precip_ResolveAdaptiveTimeGridExtension(config, defaultValue)
% precip_ResolveAdaptiveTimeGridExtension - Resolve adaptive time-grid extension.
%
% The copied PNM solver can shrink future time steps after CFL or concentration
% diagnostics. When no preplanned step remains, this option allows the internal
% stepper to append additional steps so the run can still reach endTime.

if nargin < 2
    defaultValue = false;
end

tf = defaultValue;
if isstruct(config) && isfield(config, 'allowAdaptiveTimeGridExtension') && ...
        ~isempty(config.allowAdaptiveTimeGridExtension)
    tf = parseScalarBoolean(config.allowAdaptiveTimeGridExtension);
end
end

function tf = parseScalarBoolean(value)
if islogical(value) || isnumeric(value)
    if ~isscalar(value)
        error('RTSPHEM:Precipitate:InvalidAdaptiveTimeGridExtension', ...
            'allowAdaptiveTimeGridExtension must be a scalar boolean value.');
    end
    if isnumeric(value) && ~(value == 0 || value == 1)
        error('RTSPHEM:Precipitate:InvalidAdaptiveTimeGridExtension', ...
            'Numeric allowAdaptiveTimeGridExtension must be 0 or 1.');
    end
    tf = logical(value);
    return;
end

if ischar(value) || (isstring(value) && isscalar(value))
    valueText = lower(strtrim(char(value)));
    switch valueText
        case {'true', '1'}
            tf = true;
        case {'false', '0'}
            tf = false;
        otherwise
            error('RTSPHEM:Precipitate:InvalidAdaptiveTimeGridExtension', ...
                'Invalid allowAdaptiveTimeGridExtension value: %s.', char(value));
    end
    return;
end

error('RTSPHEM:Precipitate:InvalidAdaptiveTimeGridExtension', ...
    'allowAdaptiveTimeGridExtension must be scalar logical, numeric, or true/false text.');
end
