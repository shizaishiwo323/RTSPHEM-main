function diagnostics = ValidateFaceFluxDivergence(flow, numCells, options)
%VALIDATEFACEFLUXDIVERGENCE Check finite-volume face flux closure.
%
% Positive internal velocity moves from internal_face_cells(:,1) to (:,2).
% Positive Dirichlet boundary velocity enters the cell; positive outflow
% boundary velocity leaves the cell.

if nargin < 3 || isempty(options)
    options = struct();
end
if ~(isscalar(numCells) && isfinite(numCells) && numCells == round(numCells) && ...
        numCells >= 1)
    error('RTSPHEM:Flow:InvalidCellCount', ...
        'numCells must be a positive integer scalar.');
end
tolerance = getScalarOption(options, 'absoluteTolerance_cm3_s', 1e-12);
relativeTolerance = getScalarOption(options, 'relativeTolerance', Inf, true);
activeCells = getActiveCells(options, numCells);

cellDivergence = zeros(numCells, 1);
[internalCells, internalFlux] = internalFaceFluxes(flow, numCells);
for iFace = 1:size(internalCells, 1)
    leftCell = internalCells(iFace, 1);
    rightCell = internalCells(iFace, 2);
    flux = internalFlux(iFace);
    cellDivergence(leftCell) = cellDivergence(leftCell) - flux;
    cellDivergence(rightCell) = cellDivergence(rightCell) + flux;
end

[boundaryCells, boundaryFlux, boundaryType] = boundaryFaceFluxes(flow, numCells);
inletFlow = 0;
outletFlow = 0;
for iFace = 1:numel(boundaryCells)
    cellId = boundaryCells(iFace);
    flux = boundaryFlux(iFace);
    if boundaryType(iFace) == "outflow"
        cellDivergence(cellId) = cellDivergence(cellId) - flux;
        outletFlow = outletFlow + flux;
    else
        cellDivergence(cellId) = cellDivergence(cellId) + flux;
        inletFlow = inletFlow + flux;
    end
end

activeDivergence = cellDivergence(activeCells);
maxAbsDivergence = max(abs(activeDivergence), [], 'omitnan');
if isempty(maxAbsDivergence)
    maxAbsDivergence = 0;
end
globalResidual = sum(cellDivergence(activeCells), 'omitnan');
inletOutletRelativeResidual = abs(inletFlow - outletFlow) ./ ...
    max([abs(inletFlow), abs(outletFlow), eps]);
failureReasons = strings(0, 1);
if maxAbsDivergence > tolerance
    failureReasons(end + 1, 1) = "cell divergence exceeds tolerance";
end
if abs(globalResidual) > tolerance
    failureReasons(end + 1, 1) = "global flow residual exceeds tolerance";
end
if inletOutletRelativeResidual > relativeTolerance
    failureReasons(end + 1, 1) = ...
        "inlet/outlet flow imbalance exceeds tolerance";
end

diagnostics = struct();
diagnostics.accepted = isempty(failureReasons);
diagnostics.cell_divergence_cm3_s = cellDivergence;
diagnostics.max_abs_cell_divergence_cm3_s = maxAbsDivergence;
diagnostics.global_residual_cm3_s = globalResidual;
diagnostics.inlet_flow_cm3_s = inletFlow;
diagnostics.outlet_flow_cm3_s = outletFlow;
diagnostics.inlet_outlet_relative_residual = inletOutletRelativeResidual;
diagnostics.absolute_tolerance_cm3_s = tolerance;
diagnostics.relative_tolerance = relativeTolerance;
diagnostics.active_fluid_cell = activeCells;
diagnostics.failure_reasons = failureReasons;
end

function [faceCells, fluxCm3S] = internalFaceFluxes(flow, numCells)
faceCells = zeros(0, 2);
fluxCm3S = zeros(0, 1);
if ~isstruct(flow) || ~isfield(flow, 'internal_face_cells') || ...
        isempty(flow.internal_face_cells)
    return;
end
faceCells = flow.internal_face_cells;
if size(faceCells, 2) ~= 2 || any(faceCells(:) < 1) || ...
        any(faceCells(:) > numCells) || any(faceCells(:) ~= round(faceCells(:)))
    error('RTSPHEM:Flow:InvalidInternalFaces', ...
        'flow.internal_face_cells must reference valid cells.');
end
numFaces = size(faceCells, 1);
if isfield(flow, 'internal_face_flux_cm3_s') && ...
        ~isempty(flow.internal_face_flux_cm3_s)
    fluxCm3S = column(flow.internal_face_flux_cm3_s, numFaces, ...
        'internal_face_flux_cm3_s');
    return;
end
areaCm2 = column(requiredField(flow, 'internal_face_area_cm2'), numFaces, ...
    'internal_face_area_cm2');
velocityCmS = column(requiredField(flow, 'internal_face_velocity_cm_s'), ...
    numFaces, 'internal_face_velocity_cm_s');
fluxCm3S = velocityCmS .* areaCm2;
end

function [faceCells, fluxCm3S, boundaryType] = boundaryFaceFluxes(flow, numCells)
faceCells = zeros(0, 1);
fluxCm3S = zeros(0, 1);
boundaryType = strings(0, 1);
if ~isstruct(flow) || ~isfield(flow, 'boundary_face_cells') || ...
        isempty(flow.boundary_face_cells)
    return;
end
faceCells = column(flow.boundary_face_cells, numel(flow.boundary_face_cells), ...
    'boundary_face_cells');
if any(faceCells < 1) || any(faceCells > numCells) || ...
        any(faceCells ~= round(faceCells))
    error('RTSPHEM:Flow:InvalidBoundaryFaces', ...
        'flow.boundary_face_cells must reference valid cells.');
end
numFaces = numel(faceCells);
if isfield(flow, 'boundary_face_flux_cm3_s') && ...
        ~isempty(flow.boundary_face_flux_cm3_s)
    fluxCm3S = column(flow.boundary_face_flux_cm3_s, numFaces, ...
        'boundary_face_flux_cm3_s');
else
    areaCm2 = column(requiredField(flow, 'boundary_face_area_cm2'), numFaces, ...
        'boundary_face_area_cm2');
    velocityCmS = column(requiredField(flow, 'boundary_face_velocity_cm_s'), ...
        numFaces, 'boundary_face_velocity_cm_s');
    fluxCm3S = velocityCmS .* areaCm2;
end
boundaryType = boundaryTypes(flow, numFaces);
if any(fluxCm3S < 0)
    error('RTSPHEM:Flow:InvalidBoundaryFlux', ...
        'Boundary flux magnitudes must be nonnegative.');
end
end

function values = boundaryTypes(flow, numFaces)
if isfield(flow, 'boundary_type') && ~isempty(flow.boundary_type)
    values = string(flow.boundary_type(:));
else
    values = repmat("dirichlet", numFaces, 1);
end
if isscalar(values) && numFaces > 1
    values = repmat(values, numFaces, 1);
end
if numel(values) ~= numFaces || any(~(values == "dirichlet" | values == "outflow"))
    error('RTSPHEM:Flow:InvalidBoundaryType', ...
        'boundary_type must be dirichlet or outflow.');
end
end

function values = getActiveCells(options, numCells)
if isfield(options, 'active_fluid_cell') && ~isempty(options.active_fluid_cell)
    values = logical(options.active_fluid_cell(:));
else
    values = true(numCells, 1);
end
if numel(values) ~= numCells
    error('RTSPHEM:Flow:InvalidActiveCells', ...
        'options.active_fluid_cell must match numCells.');
end
end

function value = getScalarOption(options, fieldName, defaultValue, allowInf)
if nargin < 4
    allowInf = false;
end
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if allowInf
    valid = isscalar(value) && (isfinite(value) || isinf(value)) && value >= 0;
else
    valid = isscalar(value) && isfinite(value) && value >= 0;
end
if ~valid
    error('RTSPHEM:Flow:InvalidOption', ...
        'options.%s must be a nonnegative scalar.', fieldName);
end
end

function value = requiredField(structValue, fieldName)
if ~isstruct(structValue) || ~isfield(structValue, fieldName) || ...
        isempty(structValue.(fieldName))
    error('RTSPHEM:Flow:MissingFaceField', ...
        'flow.%s is required.', fieldName);
end
value = structValue.(fieldName);
end

function values = column(values, count, fieldName)
values = values(:);
if numel(values) ~= count || any(~isfinite(values))
    error('RTSPHEM:Flow:InvalidFaceField', ...
        '%s must contain %d finite values.', fieldName, count);
end
end
