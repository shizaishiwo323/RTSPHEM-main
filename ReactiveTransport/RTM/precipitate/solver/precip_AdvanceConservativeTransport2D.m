function [updated, ledger] = precip_AdvanceConservativeTransport2D(components, spec, dt, options)
% precip_AdvanceConservativeTransport2D - Explicit conservative transport step.
%
% Inputs:
%   components - struct of conservative component concentration fields.
%   spec       - Yoon benchmark spec with grid, velocity, D, and inlets.
%   dt         - time step in seconds.
%   options    - boundaryMode: 'closed' or 'split_inlet'.
%
% Output:
%   updated - component fields after one no-clipping transport step.

if nargin < 4
    options = struct();
end
if dt < 0
    error('RTSPHEM:Precipitate:InvalidTransportStep', ...
        'Transport time step must be nonnegative.');
end

boundaryMode = getFieldOrDefault(options, 'boundaryMode', 'closed');
flowField = getFieldOrDefault(options, 'flowField', []);
transport = buildTransportContext(spec, options);
updated = components;
ledger = initializeLedger(spec);
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    c = components.(fieldName);
    validateComponentField(c, spec, fieldName);
    [updated.(fieldName), scalarLedger] = advanceScalar(c, fieldName, ...
        spec, dt, boundaryMode, flowField, transport);
    ledger.massInitial.(fieldName) = scalarLedger.massInitial;
    ledger.massFinal.(fieldName) = scalarLedger.massFinal;
    ledger.massChange.(fieldName) = scalarLedger.massChange;
    ledger.boundaryNetFlux.(fieldName) = scalarLedger.boundaryNetFlux;
    ledger.boundaryClosureError.(fieldName) = scalarLedger.boundaryClosureError;
end
ledger.maxBoundaryClosureError = max(abs(struct2array( ...
    ledger.boundaryClosureError)));
end

function [cNew, ledger] = advanceScalar(c, fieldName, spec, dt, boundaryMode, flowField, transport)
if dt == 0
    cNew = c;
    ledger = scalarLedger(c, cNew, 0, spec, dt);
    return;
end

u = getFieldOrDefault(spec, 'darcyVelocity_cm_s', 0);
d = max(getFieldOrDefault(spec, 'diffusionCoefficient_cm2_s', 0), 0);
dx = spec.dx_cm;
dy = spec.dy_cm;

if ~isempty(flowField)
    [advInc, boundaryNetFlux] = velocityFieldAdvectionIncrement(c, ...
        fieldName, spec, flowField, dx, dy, dt, boundaryMode, transport);
else
    [advInc, boundaryNetFlux] = advectionIncrement(c, fieldName, spec, u, ...
        dx, dt, boundaryMode, transport);
end
cNew = c + diffusionIncrement(c, d, dx, dy, dt, transport) + advInc;
ledger = scalarLedger(c, cNew, boundaryNetFlux, spec, dt);
end

function inc = diffusionIncrement(c, d, dx, dy, dt, transport)
dCell = transport.effectiveDiffusivity_cm2_s;
if isempty(dCell)
    dCell = ones(size(c)) .* d;
end
dCell(~transport.diffusionMask) = 0;
if all(dCell(:) == 0)
    inc = zeros(size(c));
    return;
end

inc = zeros(size(c));
if size(c, 2) > 1
    dFaceX = harmonicMeanNonnegative(dCell(:, 1:end-1), dCell(:, 2:end));
    deltaX = dFaceX .* dt .* (c(:, 2:end) - c(:, 1:end-1)) ./ dx.^2;
    inc(:, 1:end-1) = inc(:, 1:end-1) + deltaX;
    inc(:, 2:end) = inc(:, 2:end) - deltaX;
end
if size(c, 1) > 1
    dFaceY = harmonicMeanNonnegative(dCell(1:end-1, :), dCell(2:end, :));
    deltaY = dFaceY .* dt .* (c(2:end, :) - c(1:end-1, :)) ./ dy.^2;
    inc(1:end-1, :) = inc(1:end-1, :) + deltaY;
    inc(2:end, :) = inc(2:end, :) - deltaY;
end
end

function [inc, boundaryNetFlux] = advectionIncrement(c, fieldName, spec, u, dx, dt, boundaryMode, transport)
if u == 0
    inc = zeros(size(c));
    boundaryNetFlux = 0;
    return;
end
if u < 0
    error('RTSPHEM:Precipitate:NegativeDarcyVelocityUnsupported', ...
        'This smoke transport operator currently supports u >= 0 only.');
end

advectiveMask = transport.advectiveMask;
inlet = inletColumn(fieldName, spec, c, boundaryMode);
fluxX = zeros(size(c, 1), size(c, 2) + 1);
if ~strcmp(boundaryMode, 'closed')
    fluxX(advectiveMask(:, 1), 1) = u .* inlet(advectiveMask(:, 1));
end
if size(c, 2) > 1
    faceMask = advectiveMask(:, 1:end-1) & advectiveMask(:, 2:end);
    fluxX(:, 2:end-1) = u .* c(:, 1:end-1) .* faceMask;
end
if ~strcmp(boundaryMode, 'closed')
    fluxX(advectiveMask(:, end), end) = u .* c(advectiveMask(:, end), end);
end
inc = -dt ./ dx .* (fluxX(:, 2:end) - fluxX(:, 1:end-1));
leftFlux = sum(fluxX(:, 1)) .* spec.dy_cm .* spec.thickness_cm;
rightFlux = sum(fluxX(:, end)) .* spec.dy_cm .* spec.thickness_cm;
if strcmp(boundaryMode, 'closed')
    leftFlux = 0;
    rightFlux = 0;
end
boundaryNetFlux = leftFlux - rightFlux;
end

function [inc, boundaryNetFlux] = velocityFieldAdvectionIncrement(c, ...
    fieldName, spec, flowField, dx, dy, dt, boundaryMode, transport)
uCell = requireVelocityField(flowField, 'velocityX_cm_s', spec);
if isfield(flowField, 'velocityY_cm_s') && ~isempty(flowField.velocityY_cm_s)
    vCell = requireVelocityField(flowField, 'velocityY_cm_s', spec);
else
    vCell = zeros(size(uCell));
end

[numY, numX] = size(c);
advectiveMask = transport.advectiveMask;
fluxX = zeros(numY, numX + 1);
inlet = inletColumn(fieldName, spec, c, boundaryMode);
for ixFace = 1:(numX + 1)
    if ixFace == 1
        uFace = uCell(:, 1);
        cUpwind = c(:, 1);
        if ~strcmp(boundaryMode, 'closed')
            cUpwind(uFace >= 0) = inlet(uFace >= 0);
        end
        faceMask = advectiveMask(:, 1);
    elseif ixFace == numX + 1
        uFace = uCell(:, end);
        cUpwind = c(:, end);
        faceMask = advectiveMask(:, end);
    else
        uFace = 0.5 .* (uCell(:, ixFace - 1) + uCell(:, ixFace));
        cUpwind = c(:, ixFace - 1);
        reverse = uFace < 0;
        cUpwind(reverse) = c(reverse, ixFace);
        faceMask = advectiveMask(:, ixFace - 1) & advectiveMask(:, ixFace);
    end
    if strcmp(boundaryMode, 'closed') && (ixFace == 1 || ixFace == numX + 1)
        faceMask(:) = false;
    end
    fluxX(:, ixFace) = uFace .* cUpwind .* faceMask;
end

fluxY = zeros(numY + 1, numX);
for iyFace = 2:numY
    vFace = 0.5 .* (vCell(iyFace - 1, :) + vCell(iyFace, :));
    cUpwind = c(iyFace - 1, :);
    reverse = vFace < 0;
    cUpwind(reverse) = c(iyFace, reverse);
    faceMask = advectiveMask(iyFace - 1, :) & advectiveMask(iyFace, :);
    fluxY(iyFace, :) = vFace .* cUpwind .* faceMask;
end

inc = -dt ./ dx .* (fluxX(:, 2:end) - fluxX(:, 1:end-1)) - ...
    dt ./ dy .* (fluxY(2:end, :) - fluxY(1:end-1, :));
if strcmp(boundaryMode, 'closed')
    boundaryNetFlux = 0;
else
    boundaryNetFlux = (sum(fluxX(:, 1)) - sum(fluxX(:, end))) .* ...
        spec.dy_cm .* spec.thickness_cm;
end
end

function transport = buildTransportContext(spec, options)
substrateMask = getLogicalMaskOption(options, 'substrateMask', spec, false);
blockedMask = getLogicalMaskOption(options, 'blockedMask', spec, false);
transport.diffusionMask = ~substrateMask;
transport.advectiveMask = ~(substrateMask | blockedMask);
transport.effectiveDiffusivity_cm2_s = getDiffusivityFieldOption( ...
    options, spec);
end

function mask = getLogicalMaskOption(options, fieldName, spec, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    mask = logical(options.(fieldName));
    if ~isequal(size(mask), [spec.numY, spec.numX])
        error('RTSPHEM:Precipitate:InvalidTransportMaskSize', ...
            '%s must have size [numY, numX].', fieldName);
    end
else
    mask = false(spec.numY, spec.numX);
    mask(:) = defaultValue;
end
end

function dCell = getDiffusivityFieldOption(options, spec)
if isfield(options, 'effectiveDiffusivity_cm2_s') && ...
        ~isempty(options.effectiveDiffusivity_cm2_s)
    dCell = options.effectiveDiffusivity_cm2_s;
elseif isfield(options, 'effectiveDiffusivity') && ...
        ~isempty(options.effectiveDiffusivity)
    dCell = options.effectiveDiffusivity;
else
    dCell = [];
    return;
end
if ~isequal(size(dCell), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidDiffusivityFieldSize', ...
        'effectiveDiffusivity_cm2_s must have size [numY, numX].');
end
if any(~isfinite(dCell(:))) || any(dCell(:) < 0)
    error('RTSPHEM:Precipitate:InvalidDiffusivityFieldValue', ...
        'effectiveDiffusivity_cm2_s must contain finite nonnegative values.');
end
end

function h = harmonicMeanNonnegative(a, b)
denom = a + b;
h = zeros(size(denom));
mask = denom > 0 & a > 0 & b > 0;
h(mask) = 2 .* a(mask) .* b(mask) ./ denom(mask);
end

function velocity = requireVelocityField(flowField, fieldName, spec)
if ~isfield(flowField, fieldName) || isempty(flowField.(fieldName))
    error('RTSPHEM:Precipitate:MissingVelocityField', ...
        'flowField.%s is required for velocity-field transport.', fieldName);
end
velocity = flowField.(fieldName);
if ~isequal(size(velocity), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidVelocityFieldSize', ...
        'flowField.%s must have size [numY, numX].', fieldName);
end
if any(~isfinite(velocity(:)))
    error('RTSPHEM:Precipitate:InvalidVelocityFieldValue', ...
        'flowField.%s contains non-finite values.', fieldName);
end
end

function inlet = inletColumn(fieldName, spec, c, boundaryMode)
switch boundaryMode
    case 'closed'
        inlet = c(:, 1);
    case 'split_inlet'
        yCenters = ((1:spec.numY)' - 0.5) .* spec.dy_cm;
        lowerMask = yCenters < spec.splitInletY_cm;
        inlet = zeros(spec.numY, 1);
        inlet(lowerMask) = spec.inletA.(fieldName);
        inlet(~lowerMask) = spec.inletB.(fieldName);
    otherwise
        error('RTSPHEM:Precipitate:InvalidTransportBoundaryMode', ...
            'Unsupported boundaryMode: %s.', boundaryMode);
end
end

function validateComponentField(c, spec, fieldName)
if ~isequal(size(c), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidComponentFieldSize', ...
        '%s must have size [numY, numX].', fieldName);
end
if any(~isfinite(c(:)))
    error('RTSPHEM:Precipitate:InvalidComponentFieldValue', ...
        '%s contains non-finite values.', fieldName);
end
end

function ledger = scalarLedger(cInitial, cFinal, boundaryNetFlux, spec, dt)
ledger = struct();
ledger.massInitial = sum(cInitial(:)) * spec.cellVolume_cm3;
ledger.massFinal = sum(cFinal(:)) * spec.cellVolume_cm3;
ledger.massChange = ledger.massFinal - ledger.massInitial;
ledger.boundaryNetFlux = boundaryNetFlux;
ledger.boundaryClosureError = ledger.massChange - boundaryNetFlux * dt;
end

function ledger = initializeLedger(spec)
ledger = struct();
ledger.massInitial = struct();
ledger.massFinal = struct();
ledger.massChange = struct();
ledger.boundaryNetFlux = struct();
ledger.boundaryClosureError = struct();
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    ledger.massInitial.(fieldName) = NaN;
    ledger.massFinal.(fieldName) = NaN;
    ledger.massChange.(fieldName) = NaN;
    ledger.boundaryNetFlux.(fieldName) = NaN;
    ledger.boundaryClosureError.(fieldName) = NaN;
end
ledger.maxBoundaryClosureError = NaN;
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
