function flow = precip_SolveYoonStokesBrinkman2D(state, spec, options)
% precip_SolveYoonStokesBrinkman2D - Cell-centered Stokes-Brinkman flow solve.
%
% This lightweight solver gives the Yoon micro-continuum path a local
% Stokes-like recomputation hook without depending on the legacy level-set
% HyPHM driver. It solves a finite-difference Stokes-Brinkman saddle-point
% system on the Yoon Cartesian cells and marks blocked/substrate cells as
% no-flow. Production benchmark claims still require grid convergence and a
% validated Stokes backend.

if nargin < 3 || isempty(options)
    options = struct();
end

substrateMask = logicalMask(state, 'substrateMask');
blockedMask = logicalMask(state, 'blockedMask');
mobileMask = ~(substrateMask | blockedMask);
baselineMobileMask = ~substrateMask;

[pressure, velocityX, velocityY, solveInfo] = solveStokesSystem( ...
    mobileMask, spec, options);
[~, baselineVelocityX, ~, baselineInfo] = solveStokesSystem( ...
    baselineMobileMask, spec, options);

cellDepth = getSpecValue(spec, 'depth_cm', ...
    getSpecValue(spec, 'thickness_cm', 1));
dy = getSpecValue(spec, 'dy_cm', getSpecValue(spec, 'cellSize_cm', 1));
outletFlux = computeBoundaryFlux(velocityX, mobileMask(:, end), 'right', dy, ...
    cellDepth);
baselineOutletFlux = computeBoundaryFlux(baselineVelocityX, ...
    baselineMobileMask(:, end), 'right', dy, cellDepth);
relativePermeability = min(max(outletFlux / max(baselineOutletFlux, eps), ...
    0), 1);

flow = struct();
flow.solver = 'finite_difference_stokes_brinkman';
flow.isProxy = false;
flow.isStokes = true;
flow.pressure = pressure;
flow.velocityX_cm_s = velocityX;
flow.velocityY_cm_s = velocityY;
flow.inletFlux_cm3_s = computeBoundaryFlux(velocityX, mobileMask(:, 1), ...
    'left', dy, cellDepth);
flow.outletFlux_cm3_s = outletFlux;
flow.baselineOutletFlux_cm3_s = baselineOutletFlux;
flow.relativePermeability = relativePermeability;
flow.pressureDropRelative = max(solveInfo.pressureSpan, eps) / ...
    max(baselineInfo.pressureSpan, eps);
flow.flowRateRelative = relativePermeability;
flow.numBlockedCells = nnz(blockedMask & ~substrateMask);
flow.numFluidCells = nnz(~substrateMask);
flow.blockedFraction = flow.numBlockedCells / max(flow.numFluidCells, 1);
flow.numMobileCells = nnz(mobileMask);
flow.momentumSystemSize = solveInfo.linearSystemSize;
flow.linearSystemSize = solveInfo.linearSystemSize;
flow.linearResidualRelative = solveInfo.linearResidualRelative;
flow.brinkmanViscosity_cm2_s = solveInfo.viscosity;
flow.brinkmanDrag_s_inv = solveInfo.drag;
flow.note = ['Finite-difference Stokes-Brinkman diagnostic; validate ', ...
    'against the production Stokes backend before final benchmark claims.'];
end

function [pressure, velocityX, velocityY, info] = solveStokesSystem( ...
    mobileMask, spec, options)
[numY, numX] = size(mobileMask);
pressure = nan(numY, numX);
velocityX = zeros(numY, numX);
velocityY = zeros(numY, numX);

numMobile = nnz(mobileMask);
viscosity = getOption(options, 'brinkmanViscosity_cm2_s', ...
    getSpecValue(spec, 'kinematicViscosity_cm2_s', 0.01));
drag = getOption(options, 'brinkmanDrag_s_inv', 1);
inletVelocity = getSpecValue(spec, 'darcyVelocity_cm_s', 0);
dx = getSpecValue(spec, 'dx_cm', getSpecValue(spec, 'cellSize_cm', 1));
dy = getSpecValue(spec, 'dy_cm', getSpecValue(spec, 'cellSize_cm', dx));

if numMobile == 0
    info = makeInfo(0, 0, 0, viscosity, drag, pressure);
    return;
end

cellIds = zeros(numY, numX);
cellIds(mobileMask) = 1:numMobile;
uIds = cellIds;
vIds = cellIds + numMobile;
pIds = cellIds + 2 * numMobile;
referenceCell = choosePressureReferenceCell(mobileMask);

rows = zeros(0, 1);
cols = zeros(0, 1);
vals = zeros(0, 1);
rhs = zeros(3 * numMobile, 1);
eq = 0;

for iy = 1:numY
    for ix = 1:numX
        if ~mobileMask(iy, ix)
            continue;
        end

        if ix == 1
            eq = eq + 1;
            [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
                uIds(iy, ix), 1);
            rhs(eq) = inletVelocity;
        else
            eq = eq + 1;
            [rows, cols, vals, rhs] = addMomentumEquation(rows, cols, ...
                vals, rhs, eq, iy, ix, 1, uIds, pIds, mobileMask, ...
                viscosity, drag, dx, dy);
        end

        if ix == 1
            eq = eq + 1;
            [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
                vIds(iy, ix), 1);
        else
            eq = eq + 1;
            [rows, cols, vals, rhs] = addMomentumEquation(rows, cols, ...
                vals, rhs, eq, iy, ix, 2, vIds, pIds, mobileMask, ...
                viscosity, drag, dx, dy);
        end

        eq = eq + 1;
        if iy == referenceCell(1) && ix == referenceCell(2)
            [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
                pIds(iy, ix), 1);
        else
            [rows, cols, vals, rhs] = addContinuityEquation(rows, cols, ...
                vals, rhs, eq, iy, ix, uIds, vIds, mobileMask, ...
                inletVelocity, dx, dy);
        end
    end
end

matrix = sparse(rows, cols, vals, 3 * numMobile, 3 * numMobile);
solution = matrix \ rhs;
velocityX(mobileMask) = solution(uIds(mobileMask));
velocityY(mobileMask) = solution(vIds(mobileMask));
pressure(mobileMask) = solution(pIds(mobileMask));
residual = norm(matrix * solution - rhs);
relativeResidual = residual / max(norm(rhs), eps);
info = makeInfo(3 * numMobile, residual, relativeResidual, viscosity, ...
    drag, pressure);
end

function [rows, cols, vals, rhs] = addMomentumEquation(rows, cols, vals, ...
    rhs, eq, iy, ix, component, velocityIds, pressureIds, mobileMask, ...
    viscosity, drag, dx, dy)
[rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
    velocityIds(iy, ix), drag);
[rows, cols, vals] = addLaplacian(rows, cols, vals, eq, iy, ix, ...
    velocityIds, mobileMask, viscosity, dx, dy);
if component == 1
    [rows, cols, vals] = addPressureGradientX(rows, cols, vals, eq, ...
        iy, ix, pressureIds, mobileMask, dx);
else
    [rows, cols, vals] = addPressureGradientY(rows, cols, vals, eq, ...
        iy, ix, pressureIds, mobileMask, dy);
end
end

function [rows, cols, vals] = addLaplacian(rows, cols, vals, eq, iy, ix, ...
    velocityIds, mobileMask, viscosity, dx, dy)
[numY, numX] = size(mobileMask);
neighbors = [0, -1, dx; 0, 1, dx; -1, 0, dy; 1, 0, dy];
for iNeighbor = 1:size(neighbors, 1)
    ny = iy + neighbors(iNeighbor, 1);
    nx = ix + neighbors(iNeighbor, 2);
    h = neighbors(iNeighbor, 3);
    coeff = viscosity / (h * h);
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        velocityIds(iy, ix), coeff);
    if ny >= 1 && ny <= numY && nx >= 1 && nx <= numX && ...
            mobileMask(ny, nx)
        [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
            velocityIds(ny, nx), -coeff);
    end
end
end

function [rows, cols, vals] = addPressureGradientX(rows, cols, vals, eq, ...
    iy, ix, pressureIds, mobileMask, dx)
[~, numX] = size(mobileMask);
if ix > 1 && mobileMask(iy, ix - 1)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        pressureIds(iy, ix - 1), -1 / (2 * dx));
end
if ix < numX && mobileMask(iy, ix + 1)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        pressureIds(iy, ix + 1), 1 / (2 * dx));
end
end

function [rows, cols, vals] = addPressureGradientY(rows, cols, vals, eq, ...
    iy, ix, pressureIds, mobileMask, dy)
[numY, ~] = size(mobileMask);
if iy > 1 && mobileMask(iy - 1, ix)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        pressureIds(iy - 1, ix), -1 / (2 * dy));
end
if iy < numY && mobileMask(iy + 1, ix)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        pressureIds(iy + 1, ix), 1 / (2 * dy));
end
end

function [rows, cols, vals, rhs] = addContinuityEquation(rows, cols, vals, ...
    rhs, eq, iy, ix, uIds, vIds, mobileMask, inletVelocity, dx, dy)
[numY, numX] = size(mobileMask);
if ix > 1 && mobileMask(iy, ix - 1)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        uIds(iy, ix - 1), -1 / (2 * dx));
elseif ix == 1
    rhs(eq) = rhs(eq) + inletVelocity / (2 * dx);
end
if ix < numX && mobileMask(iy, ix + 1)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        uIds(iy, ix + 1), 1 / (2 * dx));
end
if iy > 1 && mobileMask(iy - 1, ix)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        vIds(iy - 1, ix), -1 / (2 * dy));
end
if iy < numY && mobileMask(iy + 1, ix)
    [rows, cols, vals] = addCoeff(rows, cols, vals, eq, ...
        vIds(iy + 1, ix), 1 / (2 * dy));
end
end

function referenceCell = choosePressureReferenceCell(mobileMask)
[rows, cols] = find(mobileMask);
rightmost = max(cols);
candidate = find(cols == rightmost, 1, 'first');
referenceCell = [rows(candidate), cols(candidate)];
end

function [rows, cols, vals] = addCoeff(rows, cols, vals, row, col, value)
if col <= 0 || value == 0
    return;
end
rows(end + 1, 1) = row; %#ok<AGROW>
cols(end + 1, 1) = col; %#ok<AGROW>
vals(end + 1, 1) = value; %#ok<AGROW>
end

function flux = computeBoundaryFlux(velocityX, boundaryMask, side, dy, depth)
if isempty(velocityX)
    flux = 0;
    return;
end
switch side
    case 'left'
        columnVelocity = velocityX(:, 1);
    case 'right'
        columnVelocity = velocityX(:, end);
    otherwise
        error('RTSPHEM:Precipitate:InvalidFlowBoundarySide', ...
            'Unsupported flow boundary side: %s.', side);
end
flux = sum(max(columnVelocity(boundaryMask), 0)) * dy * depth;
end

function info = makeInfo(linearSystemSize, residual, relativeResidual, ...
    viscosity, drag, pressure)
finitePressure = pressure(isfinite(pressure));
if isempty(finitePressure)
    pressureSpan = 0;
else
    pressureSpan = max(finitePressure) - min(finitePressure);
end
info = struct();
info.linearSystemSize = linearSystemSize;
info.linearResidual = residual;
info.linearResidualRelative = relativeResidual;
info.viscosity = viscosity;
info.drag = drag;
info.pressureSpan = pressureSpan;
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

function value = getSpecValue(spec, fieldName, defaultValue)
if isfield(spec, fieldName) && ~isempty(spec.(fieldName))
    value = spec.(fieldName);
else
    value = defaultValue;
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
