function mapped = precip_MapHyPHMFluxToCartesianFaces(hyphmFlow, state, spec, options)
% precip_MapHyPHMFluxToCartesianFaces - Map HyPHM velocity samples to Yoon cells.
%
% This bridge intentionally uses explicit velocity samples from the HyPHM-side
% solver. It keeps the Yoon finite-volume transport API independent of HyPHM
% classes while preserving flow diagnostics needed by readiness gates.

if nargin < 4
    options = struct();
end

[samplePoints, velocitySamples] = readVelocitySamples(hyphmFlow);
[xCenters, yCenters] = cellCenters(state, spec);
[xGrid, yGrid] = meshgrid(xCenters, yCenters);

velocityX = zeros(size(xGrid));
velocityY = zeros(size(yGrid));
for iCell = 1:numel(xGrid)
    iSample = nearestSampleIndex([xGrid(iCell), yGrid(iCell)], samplePoints);
    velocityX(iCell) = velocitySamples(iSample, 1);
    velocityY(iCell) = velocitySamples(iSample, 2);
end

substrateMask = logicalMask(state, 'substrateMask');
blockedMask = logicalMask(state, 'blockedMask');
mobileMask = ~(substrateMask | blockedMask);
velocityX(~mobileMask) = 0;
velocityY(~mobileMask) = 0;

dy = getFieldOrDefault(spec, 'dy_cm', 1);
dx = getFieldOrDefault(spec, 'dx_cm', 1);
depth = getFieldOrDefault(spec, 'thickness_cm', ...
    getFieldOrDefault(spec, 'depth_cm', 1));
inletFlux = sum(max(velocityX(mobileMask(:, 1), 1), 0)) * dy * depth;
outletFlux = sum(max(velocityX(mobileMask(:, end), end), 0)) * dy * depth;
closure = abs(inletFlux - outletFlux) / ...
    max([abs(inletFlux), abs(outletFlux), eps]);

mapped = struct();
mapped.velocityX_cm_s = velocityX;
mapped.velocityY_cm_s = velocityY;
mapped.inletFlux_cm3_s = inletFlux;
mapped.outletFlux_cm3_s = outletFlux;
mapped.boundaryFluxClosureRelativeError = closure;
mapped.maxDivergenceResidual_s_inv = maxDivergenceResidual( ...
    velocityX, velocityY, mobileMask, dx, dy);
mapped.numMobileCells = nnz(mobileMask);
mapped.numBlockedCells = nnz(blockedMask & ~substrateMask);
mapped.numFluidCells = nnz(~substrateMask);
mapped.source = getFieldOrDefault(hyphmFlow, 'source', ...
    'hyphm_velocity_samples');
mapped.mappingMethod = getFieldOrDefault(options, 'mappingMethod', ...
    'nearest_velocity_sample');
end

function [samplePoints, velocitySamples] = readVelocitySamples(hyphmFlow)
if ~isstruct(hyphmFlow)
    error('RTSPHEM:Precipitate:InvalidHyphmFlow', ...
        'HyPHM flow output must be a struct.');
end
if isfield(hyphmFlow, 'velocitySamplePoints_cm')
    samplePoints = hyphmFlow.velocitySamplePoints_cm;
else
    error('RTSPHEM:Precipitate:MissingHyphmVelocitySamples', ...
        'hyphmFlow.velocitySamplePoints_cm is required.');
end
if isfield(hyphmFlow, 'velocitySamples_cm_s')
    velocitySamples = hyphmFlow.velocitySamples_cm_s;
elseif isfield(hyphmFlow, 'velocity_cm_s')
    velocitySamples = hyphmFlow.velocity_cm_s;
else
    error('RTSPHEM:Precipitate:MissingHyphmVelocitySamples', ...
        'hyphmFlow.velocitySamples_cm_s is required.');
end
if size(samplePoints, 2) ~= 2 || size(velocitySamples, 2) ~= 2 || ...
        size(samplePoints, 1) ~= size(velocitySamples, 1)
    error('RTSPHEM:Precipitate:InvalidHyphmVelocitySamples', ...
        'HyPHM velocity sample points and velocities must be N-by-2 arrays.');
end
if any(~isfinite(samplePoints(:))) || any(~isfinite(velocitySamples(:)))
    error('RTSPHEM:Precipitate:InvalidHyphmVelocitySamples', ...
        'HyPHM velocity samples must be finite.');
end
end

function iSample = nearestSampleIndex(point, samplePoints)
dist2 = (samplePoints(:, 1) - point(1)).^2 + ...
    (samplePoints(:, 2) - point(2)).^2;
[~, iSample] = min(dist2);
end

function residual = maxDivergenceResidual(velocityX, velocityY, mobileMask, dx, dy)
[numY, numX] = size(mobileMask);
values = zeros(nnz(mobileMask), 1);
index = 0;
for iy = 1:numY
    for ix = 1:numX
        if ~mobileMask(iy, ix)
            continue;
        end
        div = localDerivative(velocityX, mobileMask, iy, ix, 2, dx) + ...
            localDerivative(velocityY, mobileMask, iy, ix, 1, dy);
        index = index + 1;
        values(index) = div;
    end
end
if index == 0
    residual = 0;
else
    residual = max(abs(values(1:index)));
end
end

function derivative = localDerivative(values, mobileMask, iy, ix, axisIndex, h)
[numY, numX] = size(mobileMask);
if axisIndex == 2
    hasMinus = ix > 1 && mobileMask(iy, ix - 1);
    hasPlus = ix < numX && mobileMask(iy, ix + 1);
    centerValue = values(iy, ix);
    minusValue = values(iy, max(ix - 1, 1));
    plusValue = values(iy, min(ix + 1, numX));
else
    hasMinus = iy > 1 && mobileMask(iy - 1, ix);
    hasPlus = iy < numY && mobileMask(iy + 1, ix);
    centerValue = values(iy, ix);
    minusValue = values(max(iy - 1, 1), ix);
    plusValue = values(min(iy + 1, numY), ix);
end
if hasMinus && hasPlus
    derivative = (plusValue - minusValue) / (2 * h);
elseif hasPlus
    derivative = (plusValue - centerValue) / h;
elseif hasMinus
    derivative = (centerValue - minusValue) / h;
else
    derivative = 0;
end
end

function [xCenters, yCenters] = cellCenters(state, spec)
if isfield(state, 'grid') && isfield(state.grid, 'xCenters_cm') && ...
        isfield(state.grid, 'yCenters_cm')
    xCenters = state.grid.xCenters_cm(:)';
    yCenters = state.grid.yCenters_cm(:);
else
    xCenters = ((1:spec.numX) - 0.5) .* getFieldOrDefault(spec, 'dx_cm', 1);
    yCenters = ((1:spec.numY)' - 0.5) .* getFieldOrDefault(spec, 'dy_cm', 1);
end
end

function mask = logicalMask(state, fieldName)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    mask = logical(state.(fieldName));
elseif isfield(state, 'Vm')
    mask = false(size(state.Vm));
else
    error('RTSPHEM:Precipitate:MissingVm', ...
        'state.Vm is required when %s is absent.', fieldName);
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
