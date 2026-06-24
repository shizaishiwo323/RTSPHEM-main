function massFlux = precip_ComputeBoundaryInletMassFlux(mode, config, fieldName, ...
    fallbackConcentration, inletVelocity, thickness, crossSectionLength)
% precip_ComputeBoundaryInletMassFlux - Diagnostic inlet solute flux [mol/s].

if nargin < 7 || isempty(crossSectionLength)
    crossSectionLength = 0;
end
if nargin < 6 || isempty(thickness)
    thickness = 1;
end
if nargin < 5 || isempty(inletVelocity)
    inletVelocity = 0;
end
if nargin < 4 || isempty(fallbackConcentration)
    fallbackConcentration = 0;
end

normalizedMode = normalizeMode(mode);
if ~strcmp(normalizedMode, 'split_left_inlet')
    massFlux = fallbackConcentration * abs(inletVelocity) * crossSectionLength * thickness;
    return;
end

validateSplitConfig(config, fieldName);
aRange = splitRange(config.inletA, 'yRange', [], 'inletA.yRange');
bRange = splitRange(config.inletB, 'yRange', [], 'inletB.yRange');
validateSplitRangePair(aRange, bRange);
aLength = max(aRange) - min(aRange);
bLength = max(bRange) - min(bRange);
aConcentration = validateScalarFinite(config.inletA.(fieldName), ['inletA.', fieldName]);
bConcentration = validateScalarFinite(config.inletB.(fieldName), ['inletB.', fieldName]);

massFlux = (aConcentration * aLength + bConcentration * bLength) * ...
    abs(inletVelocity) * thickness;
end

function mode = normalizeMode(mode)
if isstring(mode)
    mode = char(mode);
end
if ~ischar(mode)
    mode = 'uniform_left_inlet';
    return;
end
mode = lower(strrep(strtrim(mode), '-', '_'));
end

function validateSplitConfig(config, fieldName)
if ~isstruct(config) || ~isfield(config, 'inletA') || ~isfield(config, 'inletB')
    error('RTSPHEM:Precipitate:MissingSplitInletConfig', ...
        'split_left_inlet requires config.inletA and config.inletB.');
end
if ~isfield(config.inletA, fieldName) || ~isfield(config.inletB, fieldName)
    error('RTSPHEM:Precipitate:MissingSplitInletSpecies', ...
        'split_left_inlet requires %s in both inletA and inletB.', fieldName);
end
end

function value = validateScalarFinite(value, label)
if isempty(value) || ~isscalar(value) || ~isnumeric(value) || ~isfinite(value)
    error('RTSPHEM:Precipitate:InvalidSplitInletConcentration', ...
        '%s must be a finite numeric scalar.', label);
end
end

function range = splitRange(config, fieldName, defaultValue, label)
if isstruct(config) && isfield(config, fieldName)
    range = config.(fieldName);
else
    range = defaultValue;
end
if ~isnumeric(range) || numel(range) ~= 2 || any(~isfinite(range))
    error('RTSPHEM:Precipitate:InvalidSplitInletRange', ...
        '%s must be a finite numeric two-element vector.', label);
end
range = reshape(range, 1, 2);
if range(2) <= range(1)
    error('RTSPHEM:Precipitate:InvalidSplitInletRange', ...
        '%s must have positive width.', label);
end
end

function validateSplitRangePair(aRange, bRange)
if abs(aRange(2) - bRange(1)) > max(eps(max(abs([aRange, bRange]))), 1e-12)
    error('RTSPHEM:Precipitate:InvalidSplitInletRange', ...
        'inletA.yRange and inletB.yRange must be contiguous and non-overlapping.');
end
end
