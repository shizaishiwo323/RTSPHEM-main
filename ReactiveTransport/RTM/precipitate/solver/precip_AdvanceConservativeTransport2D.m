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
updated = components;
ledger = initializeLedger(spec);
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    c = components.(fieldName);
    validateComponentField(c, spec, fieldName);
    [updated.(fieldName), scalarLedger] = advanceScalar(c, fieldName, ...
        spec, dt, boundaryMode, flowField);
    ledger.massInitial.(fieldName) = scalarLedger.massInitial;
    ledger.massFinal.(fieldName) = scalarLedger.massFinal;
    ledger.massChange.(fieldName) = scalarLedger.massChange;
    ledger.boundaryNetFlux.(fieldName) = scalarLedger.boundaryNetFlux;
    ledger.boundaryClosureError.(fieldName) = scalarLedger.boundaryClosureError;
end
ledger.maxBoundaryClosureError = max(abs(struct2array( ...
    ledger.boundaryClosureError)));
end

function [cNew, ledger] = advanceScalar(c, fieldName, spec, dt, boundaryMode, flowField)
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
        fieldName, spec, flowField, dx, dy, dt, boundaryMode);
else
    [advInc, boundaryNetFlux] = advectionIncrement(c, fieldName, spec, u, ...
        dx, dt, boundaryMode);
end
cNew = c + diffusionIncrement(c, d, dx, dy, dt, boundaryMode) + advInc;
ledger = scalarLedger(c, cNew, boundaryNetFlux, spec, dt);
end

function inc = diffusionIncrement(c, d, dx, dy, dt, boundaryMode)
if d == 0
    inc = zeros(size(c));
    return;
end

left = [c(:, 1), c(:, 1:end-1)];
right = [c(:, 2:end), c(:, end)];
down = [c(1, :); c(1:end-1, :)];
up = [c(2:end, :); c(end, :)];

if strcmp(boundaryMode, 'split_inlet')
    left = c;
end

laplacian = (left - 2 .* c + right) ./ dx.^2 + ...
    (down - 2 .* c + up) ./ dy.^2;
inc = d .* dt .* laplacian;
end

function [inc, boundaryNetFlux] = advectionIncrement(c, fieldName, spec, u, dx, dt, boundaryMode)
if u == 0
    inc = zeros(size(c));
    boundaryNetFlux = 0;
    return;
end
if u < 0
    error('RTSPHEM:Precipitate:NegativeDarcyVelocityUnsupported', ...
        'This smoke transport operator currently supports u >= 0 only.');
end

inlet = inletColumn(fieldName, spec, c, boundaryMode);
upwindLeft = [inlet, c(:, 1:end-1)];
inc = -u .* dt ./ dx .* (c - upwindLeft);
leftFlux = u .* sum(inlet) .* spec.dy_cm .* spec.thickness_cm;
rightFlux = u .* sum(c(:, end)) .* spec.dy_cm .* spec.thickness_cm;
if strcmp(boundaryMode, 'closed')
    leftFlux = 0;
    rightFlux = 0;
end
boundaryNetFlux = leftFlux - rightFlux;
end

function [inc, boundaryNetFlux] = velocityFieldAdvectionIncrement(c, ...
    fieldName, spec, flowField, dx, dy, dt, boundaryMode)
uCell = requireVelocityField(flowField, 'velocityX_cm_s', spec);
if isfield(flowField, 'velocityY_cm_s') && ~isempty(flowField.velocityY_cm_s)
    vCell = requireVelocityField(flowField, 'velocityY_cm_s', spec);
else
    vCell = zeros(size(uCell));
end

[numY, numX] = size(c);
fluxX = zeros(numY, numX + 1);
inlet = inletColumn(fieldName, spec, c, boundaryMode);
for ixFace = 1:(numX + 1)
    if ixFace == 1
        uFace = uCell(:, 1);
        cUpwind = c(:, 1);
        if ~strcmp(boundaryMode, 'closed')
            cUpwind(uFace >= 0) = inlet(uFace >= 0);
        end
    elseif ixFace == numX + 1
        uFace = uCell(:, end);
        cUpwind = c(:, end);
    else
        uFace = 0.5 .* (uCell(:, ixFace - 1) + uCell(:, ixFace));
        cUpwind = c(:, ixFace - 1);
        reverse = uFace < 0;
        cUpwind(reverse) = c(reverse, ixFace);
    end
    fluxX(:, ixFace) = uFace .* cUpwind;
end

fluxY = zeros(numY + 1, numX);
for iyFace = 2:numY
    vFace = 0.5 .* (vCell(iyFace - 1, :) + vCell(iyFace, :));
    cUpwind = c(iyFace - 1, :);
    reverse = vFace < 0;
    cUpwind(reverse) = c(iyFace, reverse);
    fluxY(iyFace, :) = vFace .* cUpwind;
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
