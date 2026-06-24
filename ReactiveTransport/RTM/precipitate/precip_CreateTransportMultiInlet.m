function fluxFunction = precip_CreateTransportMultiInlet(mode, config, fieldName, fallbackConcentration, inletVelocity, epsValue)
% precip_CreateTransportMultiInlet - Build inlet flux function for split-flow cases.
%
% Inputs:
%   mode                  Boundary mode. 'split_left_inlet' uses config.inletA/B.
%   config                Benchmark configuration struct.
%   fieldName             Species field in inletA/inletB, e.g. 'Ca_total'.
%   fallbackConcentration Scalar inlet concentration for legacy uniform mode.
%   inletVelocity         Left boundary Darcy velocity [cm/s].
%   epsValue              Coordinate tolerance for the left boundary [cm].
%
% Output:
%   fluxFunction          Function handle compatible with HyPHM gF callbacks.

if nargin < 6 || isempty(epsValue)
    epsValue = eps;
end
if nargin < 5 || isempty(inletVelocity)
    inletVelocity = 0;
end
if nargin < 4 || isempty(fallbackConcentration)
    fallbackConcentration = 0;
end

normalizedMode = normalizeMode(mode);
isSplitLeft = strcmpi(normalizedMode, 'split_left_inlet');
if ~isSplitLeft
    fluxFunction = @(~, x) uniformLeftFlux(x, fallbackConcentration, inletVelocity, epsValue);
    return;
end

validateSplitConfig(config, fieldName);
inletA = config.inletA;
inletB = config.inletB;
splitInletY = cfgget(config, 'splitInletY', 0.5 * cfgget(config, 'lengthYAxis', 1));
aRange = splitRange(inletA, 'yRange', [0, splitInletY], 'inletA.yRange');
bRange = splitRange(inletB, 'yRange', [splitInletY, cfgget(config, 'lengthYAxis', Inf)], 'inletB.yRange');
validateSplitRangePair(aRange, bRange);
aConcentration = inletA.(fieldName);
bConcentration = inletB.(fieldName);
validateScalarFinite(aConcentration, ['inletA.', fieldName]);
validateScalarFinite(bConcentration, ['inletB.', fieldName]);

fluxFunction = @(~, x) splitLeftFlux(x, aRange, bRange, ...
    aConcentration, bConcentration, inletVelocity, epsValue);
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

function flux = uniformLeftFlux(x, concentration, inletVelocity, epsValue)
[xCoord, ~] = splitCoordinates(x);
flux = -inletVelocity .* (xCoord < epsValue) .* concentration;
end

function flux = splitLeftFlux(x, aRange, bRange, aConcentration, bConcentration, inletVelocity, epsValue)
[xCoord, yCoord] = splitCoordinates(x);
onLeft = xCoord < epsValue;
lowerInlet = yCoord >= min(aRange) & yCoord < max(aRange);
upperInlet = yCoord >= min(bRange) & yCoord <= max(bRange);
concentration = aConcentration .* lowerInlet + bConcentration .* upperInlet;
flux = -inletVelocity .* onLeft .* concentration;
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

function validateScalarFinite(value, label)
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

function [xCoord, yCoord] = splitCoordinates(x)
if size(x, 2) == 2 && size(x, 1) ~= 2
    xCoord = x(:, 1).';
    yCoord = x(:, 2).';
elseif size(x, 1) >= 2
    xCoord = x(1, :);
    yCoord = x(2, :);
elseif size(x, 2) >= 2
    xCoord = x(:, 1).';
    yCoord = x(:, 2).';
else
    error('RTSPHEM:Precipitate:InvalidFluxCoordinates', ...
        'Flux coordinates must contain x and y components.');
end
end

function value = cfgget(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end
