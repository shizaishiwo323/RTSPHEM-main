function concentrationScale = precip_ComputeBoundaryInletConcentrationScale(mode, config, ...
    fieldName, fallbackConcentration)
% precip_ComputeBoundaryInletConcentrationScale - Diagnostic inlet scale [mol/cm3].

if nargin < 4 || isempty(fallbackConcentration)
    fallbackConcentration = 0;
end

normalizedMode = normalizeMode(mode);
if ~strcmp(normalizedMode, 'split_left_inlet')
    concentrationScale = fallbackConcentration;
    return;
end

validateSplitConfig(config, fieldName);
aConcentration = validateScalarFinite(config.inletA.(fieldName), ['inletA.', fieldName]);
bConcentration = validateScalarFinite(config.inletB.(fieldName), ['inletB.', fieldName]);
concentrationScale = max(abs([aConcentration, bConcentration]));
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
