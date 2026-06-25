function [fluxes, ledger] = precip_AssembleComponentFaceFluxes(state, fieldName, spec, options)
% precip_AssembleComponentFaceFluxes - Assemble mol/s component face fluxes.

if nargin < 4
    options = struct();
end
boundaryMode = getFieldOrDefault(options, 'boundaryMode', 'closed');
flowField = getFieldOrDefault(options, 'flowField', []);
transport = buildTransportContext(state, spec, options);
c = state.components.(fieldName);

[numY, numX] = size(c);
fluxX = zeros(numY, numX + 1);
fluxY = zeros(numY + 1, numX);

[advX, advY] = assembleAdvection(c, fieldName, spec, boundaryMode, ...
    flowField, transport);
[diffX, diffY] = assembleDiffusion(c, spec, transport);
fluxX = fluxX + advX + diffX;
fluxY = fluxY + advY + diffY;

leftFlux = sum(fluxX(:, 1));
rightFlux = sum(fluxX(:, end));
if strcmp(boundaryMode, 'closed')
    leftFlux = 0;
    rightFlux = 0;
end

fluxes = struct();
fluxes.fluxX_mol_s = fluxX;
fluxes.fluxY_mol_s = fluxY;
fluxes.boundaryNetFlux_mol_s = leftFlux - rightFlux;

ledger = struct();
ledger.massInitial = sum(state.componentMoles.(fieldName)(:));
ledger.massFinal = NaN;
ledger.massChange = NaN;
ledger.boundaryNetFlux = fluxes.boundaryNetFlux_mol_s;
ledger.boundaryClosureError = NaN;
end

function [fluxX, fluxY] = assembleAdvection(c, fieldName, spec, ...
    boundaryMode, flowField, transport)
if ~isempty(flowField)
    [fluxX, fluxY] = velocityFieldAdvection(c, fieldName, spec, ...
        boundaryMode, flowField, transport);
    return;
end

u = getFieldOrDefault(spec, 'darcyVelocity_cm_s', 0);
if u < 0
    error('RTSPHEM:Precipitate:NegativeDarcyVelocityUnsupported', ...
        'This transport operator currently supports u >= 0 only.');
end
[numY, numX] = size(c);
fluxX = zeros(numY, numX + 1);
fluxY = zeros(numY + 1, numX);
if u == 0
    return;
end

areaX = spec.dy_cm * spec.thickness_cm;
advectiveMask = transport.advectiveMask;
inlet = inletColumn(fieldName, spec, c, boundaryMode);
if ~strcmp(boundaryMode, 'closed')
    fluxX(advectiveMask(:, 1), 1) = u .* inlet(advectiveMask(:, 1)) .* areaX;
end
if numX > 1
    faceMask = advectiveMask(:, 1:end-1) & advectiveMask(:, 2:end);
    fluxX(:, 2:end-1) = u .* c(:, 1:end-1) .* areaX .* faceMask;
end
if ~strcmp(boundaryMode, 'closed')
    fluxX(advectiveMask(:, end), end) = u .* c(advectiveMask(:, end), end) .* areaX;
end
end

function [fluxX, fluxY] = velocityFieldAdvection(c, fieldName, spec, ...
    boundaryMode, flowField, transport)
uCell = requireVelocityField(flowField, 'velocityX_cm_s', spec);
if isfield(flowField, 'velocityY_cm_s') && ~isempty(flowField.velocityY_cm_s)
    vCell = requireVelocityField(flowField, 'velocityY_cm_s', spec);
else
    vCell = zeros(size(uCell));
end

[numY, numX] = size(c);
areaX = spec.dy_cm * spec.thickness_cm;
areaY = spec.dx_cm * spec.thickness_cm;
advectiveMask = transport.advectiveMask;
inlet = inletColumn(fieldName, spec, c, boundaryMode);
fluxX = zeros(numY, numX + 1);
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
    fluxX(:, ixFace) = uFace .* cUpwind .* areaX .* faceMask;
end

fluxY = zeros(numY + 1, numX);
for iyFace = 2:numY
    vFace = 0.5 .* (vCell(iyFace - 1, :) + vCell(iyFace, :));
    cUpwind = c(iyFace - 1, :);
    reverse = vFace < 0;
    cUpwind(reverse) = c(iyFace, reverse);
    faceMask = advectiveMask(iyFace - 1, :) & advectiveMask(iyFace, :);
    fluxY(iyFace, :) = vFace .* cUpwind .* areaY .* faceMask;
end
end

function [fluxX, fluxY] = assembleDiffusion(c, spec, transport)
[numY, numX] = size(c);
areaX = spec.dy_cm * spec.thickness_cm;
areaY = spec.dx_cm * spec.thickness_cm;
fluxX = zeros(numY, numX + 1);
fluxY = zeros(numY + 1, numX);
dCell = transport.effectiveDiffusivity_cm2_s;
dCell(~transport.diffusionMask) = 0;
if numX > 1
    dFaceX = harmonicMeanNonnegative(dCell(:, 1:end-1), dCell(:, 2:end));
    fluxX(:, 2:end-1) = dFaceX .* (c(:, 1:end-1) - c(:, 2:end)) ./ ...
        spec.dx_cm .* areaX;
end
if numY > 1
    dFaceY = harmonicMeanNonnegative(dCell(1:end-1, :), dCell(2:end, :));
    fluxY(2:end-1, :) = dFaceY .* (c(1:end-1, :) - c(2:end, :)) ./ ...
        spec.dy_cm .* areaY;
end
end

function transport = buildTransportContext(state, spec, options)
substrateMask = getLogicalMask(state, options, 'substrateMask', spec, false);
blockedMask = getLogicalMask(state, options, 'blockedMask', spec, false);
transport.diffusionMask = ~substrateMask;
transport.advectiveMask = ~(substrateMask | blockedMask);
transport.effectiveDiffusivity_cm2_s = getDiffusivityField(state, options, spec);
end

function mask = getLogicalMask(state, options, fieldName, spec, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    mask = logical(options.(fieldName));
elseif isfield(state, fieldName) && ~isempty(state.(fieldName))
    mask = logical(state.(fieldName));
else
    mask = false(spec.numY, spec.numX);
    mask(:) = defaultValue;
end
if ~isequal(size(mask), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidTransportMaskSize', ...
        '%s must have size [numY, numX].', fieldName);
end
end

function dCell = getDiffusivityField(state, options, spec)
if isfield(options, 'effectiveDiffusivity_cm2_s') && ...
        ~isempty(options.effectiveDiffusivity_cm2_s)
    dCell = options.effectiveDiffusivity_cm2_s;
elseif isfield(state, 'effectiveDiffusivity_cm2_s') && ...
        ~isempty(state.effectiveDiffusivity_cm2_s)
    dCell = state.effectiveDiffusivity_cm2_s;
else
    dCell = ones(spec.numY, spec.numX) .* ...
        getFieldOrDefault(spec, 'diffusionCoefficient_cm2_s', 0);
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

function h = harmonicMeanNonnegative(a, b)
denom = a + b;
h = zeros(size(denom));
mask = denom > 0 & a > 0 & b > 0;
h(mask) = 2 .* a(mask) .* b(mask) ./ denom(mask);
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
