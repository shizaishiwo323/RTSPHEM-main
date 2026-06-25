function flow = precip_SolveYoonDarcyFlow2D(state, spec, options)
% precip_SolveYoonDarcyFlow2D - Topology-aware pressure-flow recomputation.
%
% This is a finite-volume Darcy pressure solve for Yoon blockage feedback.
% It is not a Stokes solver; it is a spatially resolved replacement for the
% previous blocked-fraction proxy while the production Stokes path is pending.

if nargin < 3 || isempty(options)
    options = struct();
end

mobileMask = ~(logicalMask(state, 'substrateMask') | ...
    logicalMask(state, 'blockedMask'));
permeability = computeCellMobility(state, mobileMask);
leftPressure = getOption(options, 'leftPressure', 1);
rightPressure = getOption(options, 'rightPressure', 0);

[pressure, flux, solveInfo] = solvePressureNetwork(mobileMask, permeability, ...
    leftPressure, rightPressure);

baselineMobileMask = ~logicalMask(state, 'substrateMask');
baselinePermeability = computeCellMobility(state, baselineMobileMask);
[~, baselineFlux] = solvePressureNetwork(baselineMobileMask, ...
    baselinePermeability, leftPressure, rightPressure);

outletFlux = max(flux.outlet, 0);
baselineOutletFlux = max(baselineFlux.outlet, eps);
relativePermeability = min(max(outletFlux / baselineOutletFlux, 0), 1);

flow = struct();
flow.solver = 'finite_volume_darcy_pressure';
flow.isProxy = false;
flow.isStokes = false;
flow.pressure = pressure;
flow.velocityX_cm_s = flux.velocityX_cm_s;
flow.velocityY_cm_s = flux.velocityY_cm_s;
flow.inletFlux_cm3_s = flux.inlet;
flow.outletFlux_cm3_s = flux.outlet;
flow.baselineOutletFlux_cm3_s = baselineOutletFlux;
flow.relativePermeability = relativePermeability;
flow.pressureDropRelative = 1 / max(relativePermeability, eps);
flow.flowRateRelative = relativePermeability;
flow.numBlockedCells = nnz(logicalMask(state, 'blockedMask') & ...
    ~logicalMask(state, 'substrateMask'));
flow.numFluidCells = nnz(~logicalMask(state, 'substrateMask'));
flow.blockedFraction = flow.numBlockedCells / max(flow.numFluidCells, 1);
flow.numMobileCells = nnz(mobileMask);
flow.leftPressure = leftPressure;
flow.rightPressure = rightPressure;
flow.linearSystemSize = solveInfo.linearSystemSize;
flow.note = ['Finite-volume Darcy pressure solve; use a Stokes solver ', ...
    'for final Zhang/Yoon benchmark claims.'];
end

function [pressure, flux, solveInfo] = solvePressureNetwork(mobileMask, ...
    permeability, leftPressure, rightPressure)
[numY, numX] = size(mobileMask);
cellIds = zeros(numY, numX);
cellIds(mobileMask) = 1:nnz(mobileMask);
numUnknowns = nnz(mobileMask);

pressure = nan(numY, numX);
velocityX = zeros(numY, numX);
velocityY = zeros(numY, numX);

if numUnknowns == 0
    flux = makeFlux(velocityX, velocityY, 0, 0);
    solveInfo = struct('linearSystemSize', 0);
    return;
end

rows = zeros(0, 1);
cols = zeros(0, 1);
vals = zeros(0, 1);
rhs = zeros(numUnknowns, 1);

for iy = 1:numY
    for ix = 1:numX
        if ~mobileMask(iy, ix)
            continue;
        end
        row = cellIds(iy, ix);
        diagonal = 0;
        [rows, cols, vals, rhs, diagonal] = addHorizontalFace(rows, cols, ...
            vals, rhs, diagonal, row, iy, ix, -1, cellIds, mobileMask, ...
            permeability, leftPressure, rightPressure);
        [rows, cols, vals, rhs, diagonal] = addHorizontalFace(rows, cols, ...
            vals, rhs, diagonal, row, iy, ix, 1, cellIds, mobileMask, ...
            permeability, leftPressure, rightPressure);
        [rows, cols, vals, diagonal] = addVerticalFace(rows, cols, vals, ...
            diagonal, row, iy, ix, -1, cellIds, mobileMask, permeability);
        [rows, cols, vals, diagonal] = addVerticalFace(rows, cols, vals, ...
            diagonal, row, iy, ix, 1, cellIds, mobileMask, permeability);
        if diagonal <= 0
            diagonal = 1;
        end
        rows(end + 1, 1) = row; %#ok<AGROW>
        cols(end + 1, 1) = row; %#ok<AGROW>
        vals(end + 1, 1) = diagonal; %#ok<AGROW>
    end
end

matrix = sparse(rows, cols, vals, numUnknowns, numUnknowns);
solution = matrix \ rhs;
pressure(mobileMask) = solution;

[velocityX, velocityY, inletFlux, outletFlux] = computeFluxes(pressure, ...
    mobileMask, permeability, leftPressure, rightPressure);
flux = makeFlux(velocityX, velocityY, inletFlux, outletFlux);
solveInfo = struct('linearSystemSize', numUnknowns);
end

function [rows, cols, vals, rhs, diagonal] = addHorizontalFace(rows, cols, ...
    vals, rhs, diagonal, row, iy, ix, direction, cellIds, mobileMask, ...
    permeability, leftPressure, rightPressure)
[~, numX] = size(mobileMask);
neighborX = ix + direction;
if neighborX < 1
    conductance = permeability(iy, ix);
    diagonal = diagonal + conductance;
    rhs(row) = rhs(row) + conductance * leftPressure;
elseif neighborX > numX
    conductance = permeability(iy, ix);
    diagonal = diagonal + conductance;
    rhs(row) = rhs(row) + conductance * rightPressure;
elseif mobileMask(iy, neighborX)
    conductance = faceConductance(permeability(iy, ix), ...
        permeability(iy, neighborX));
    diagonal = diagonal + conductance;
    rows(end + 1, 1) = row; %#ok<AGROW>
    cols(end + 1, 1) = cellIds(iy, neighborX); %#ok<AGROW>
    vals(end + 1, 1) = -conductance; %#ok<AGROW>
end
end

function [rows, cols, vals, diagonal] = addVerticalFace(rows, cols, vals, ...
    diagonal, row, iy, ix, direction, cellIds, mobileMask, permeability)
[numY, ~] = size(mobileMask);
neighborY = iy + direction;
if neighborY < 1 || neighborY > numY
    return;
end
if mobileMask(neighborY, ix)
    conductance = faceConductance(permeability(iy, ix), ...
        permeability(neighborY, ix));
    diagonal = diagonal + conductance;
    rows(end + 1, 1) = row; %#ok<AGROW>
    cols(end + 1, 1) = cellIds(neighborY, ix); %#ok<AGROW>
    vals(end + 1, 1) = -conductance; %#ok<AGROW>
end
end

function [velocityX, velocityY, inletFlux, outletFlux] = computeFluxes( ...
    pressure, mobileMask, permeability, leftPressure, rightPressure)
[numY, numX] = size(mobileMask);
velocityX = zeros(numY, numX);
velocityY = zeros(numY, numX);
inletFlux = 0;
outletFlux = 0;

for iy = 1:numY
    for ix = 1:numX
        if ~mobileMask(iy, ix)
            continue;
        end
        fluxLeft = 0;
        fluxRight = 0;
        fluxUp = 0;
        fluxDown = 0;
        if ix == 1
            fluxLeft = permeability(iy, ix) * (leftPressure - pressure(iy, ix));
            inletFlux = inletFlux + max(fluxLeft, 0);
        elseif mobileMask(iy, ix - 1)
            fluxLeft = faceConductance(permeability(iy, ix), ...
                permeability(iy, ix - 1)) * (pressure(iy, ix - 1) - ...
                pressure(iy, ix));
        end
        if ix == numX
            fluxRight = permeability(iy, ix) * (pressure(iy, ix) - rightPressure);
            outletFlux = outletFlux + max(fluxRight, 0);
        elseif mobileMask(iy, ix + 1)
            fluxRight = faceConductance(permeability(iy, ix), ...
                permeability(iy, ix + 1)) * (pressure(iy, ix) - ...
                pressure(iy, ix + 1));
        end
        if iy > 1 && mobileMask(iy - 1, ix)
            fluxUp = faceConductance(permeability(iy, ix), ...
                permeability(iy - 1, ix)) * (pressure(iy - 1, ix) - ...
                pressure(iy, ix));
        end
        if iy < numY && mobileMask(iy + 1, ix)
            fluxDown = faceConductance(permeability(iy, ix), ...
                permeability(iy + 1, ix)) * (pressure(iy, ix) - ...
                pressure(iy + 1, ix));
        end
        velocityX(iy, ix) = 0.5 * (fluxLeft + fluxRight);
        velocityY(iy, ix) = 0.5 * (fluxUp + fluxDown);
    end
end
end

function flux = makeFlux(velocityX, velocityY, inletFlux, outletFlux)
flux = struct();
flux.velocityX_cm_s = velocityX;
flux.velocityY_cm_s = velocityY;
flux.inlet = inletFlux;
flux.outlet = outletFlux;
end

function mobility = computeCellMobility(state, mobileMask)
mobility = zeros(size(mobileMask));
if isfield(state, 'fluidVolumeFraction') && ~isempty(state.fluidVolumeFraction)
    mobility(mobileMask) = max(state.fluidVolumeFraction(mobileMask), 0) .^ 3;
else
    mobility(mobileMask) = 1;
end
mobility(~isfinite(mobility)) = 0;
end

function conductance = faceConductance(leftValue, rightValue)
if leftValue <= 0 || rightValue <= 0
    conductance = 0;
else
    conductance = 2 * leftValue * rightValue / (leftValue + rightValue);
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

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
