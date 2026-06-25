function [updated, ledger] = CutCellTransportFV(state, geometry, dtSeconds, options)
%CUTCELLTRANSPORTFV Conservative finite-volume update on component moles.
%
% This first backend layer advances the conserved state in moles. Fluxes are
% extensive mol/s values; positive internal face flux moves from cell(:,1) to
% cell(:,2). Concentrations remain derived quantities outside this function.

if nargin < 4 || isempty(options)
    options = struct();
end
validateInputs(state, geometry, dtSeconds);

componentMoles0 = state.component_moles;
[numCells, numComponents] = size(componentMoles0);

sourceMolS = getMatrixOption(options, 'component_source_mol_s', ...
    zeros(numCells, numComponents), numCells, numComponents);
sourceDeltaMoles = sourceMolS .* dtSeconds;

if strcmp(transportTimeIntegration(options), 'implicit_euler')
    [componentMoles1, internalFluxDeltaMoles, boundaryDeltaMoles, ...
        boundaryAdvectiveDeltaMoles, boundaryDiffusiveDeltaMoles] = ...
        implicitEulerStep(state, geometry, options, sourceDeltaMoles, ...
        numCells, numComponents, dtSeconds);
else
    internalFluxDeltaMoles = zeros(numCells, numComponents);
    if hasAnyInternalFluxOption(options)
        internalFluxDeltaMoles = computeInternalFluxDelta( ...
            options, state, geometry, numCells, numComponents, dtSeconds);
    end

    [boundaryDeltaMoles, boundaryAdvectiveDeltaMoles, boundaryDiffusiveDeltaMoles] = ...
        computeBoundaryDelta(options, state, geometry, numCells, numComponents, dtSeconds);

    componentDeltaMoles = sourceDeltaMoles + internalFluxDeltaMoles + boundaryDeltaMoles;
    componentMoles1 = componentMoles0 + componentDeltaMoles;
end
roundoffSuppressedMoles = zeros(size(componentMoles1));
negativeMask = componentMoles1 < 0;
if any(componentMoles1(:) < -1e-14)
    error('RTSPHEM:Transport:NegativeComponentMoles', ...
        'Cut-cell transport update produced negative component moles.');
end
roundoffSuppressedMoles(negativeMask) = -componentMoles1(negativeMask);
componentMoles1(negativeMask) = 0;

updated = state;
updated.component_moles = componentMoles1;
updated.time_s = state.time_s + dtSeconds;
rtm.state.ValidateState(updated);

ledger = buildLedger(componentMoles0, componentMoles1, sourceDeltaMoles, ...
    internalFluxDeltaMoles, boundaryDeltaMoles, boundaryAdvectiveDeltaMoles, ...
    boundaryDiffusiveDeltaMoles, roundoffSuppressedMoles, dtSeconds, ...
    state.component_names);
end

function value = transportTimeIntegration(options)
if isfield(options, 'time_integration') && ~isempty(options.time_integration)
    value = lower(strrep(strtrim(char(options.time_integration)), '-', '_'));
else
    value = 'explicit_euler';
end
if ~any(strcmp(value, {'explicit_euler', 'implicit_euler'}))
    error('RTSPHEM:Transport:UnknownTimeIntegration', ...
        'Unsupported transport time integration: %s.', value);
end
end

function validateInputs(state, geometry, dtSeconds)
rtm.state.ValidateState(state);
if ~(isscalar(dtSeconds) && isfinite(dtSeconds) && dtSeconds >= 0)
    error('RTSPHEM:Transport:InvalidTimeStep', ...
        'dtSeconds must be a nonnegative finite scalar.');
end
if ~isstruct(geometry) || ~isfield(geometry, 'water_volume_cm3')
    error('RTSPHEM:Transport:MissingWaterVolume', ...
        'geometry.water_volume_cm3 is required.');
end
if numel(geometry.water_volume_cm3) ~= size(state.component_moles, 1)
    error('RTSPHEM:Transport:StateSizeMismatch', ...
        'geometry.water_volume_cm3 must match the state cell count.');
end
end

function values = getMatrixOption(options, fieldName, defaultValue, numRows, numCols)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    values = options.(fieldName);
else
    values = defaultValue;
end
values = values(:,:);
if ~isequal(size(values), [numRows, numCols])
    error('RTSPHEM:Transport:OptionSizeMismatch', ...
        '%s must have size [%d %d].', fieldName, numRows, numCols);
end
if any(~isfinite(values(:)))
    error('RTSPHEM:Transport:NonfiniteOption', ...
        '%s must contain only finite values.', fieldName);
end
end

function [componentMoles1, internalFluxDeltaMoles, boundaryDeltaMoles, ...
        boundaryAdvectiveDeltaMoles, boundaryDiffusiveDeltaMoles] = ...
    implicitEulerStep(state, geometry, options, sourceDeltaMoles, ...
        numCells, numComponents, dtSeconds)
componentMoles1 = zeros(size(state.component_moles));
waterVolume = geometry.water_volume_cm3(:);
if any(waterVolume <= 0 & any(state.component_moles > 0, 2))
    error('RTSPHEM:Transport:InvalidWaterVolume', ...
        'Implicit transport requires positive water volume for cells with component moles.');
end

internal = implicitInternalFaceData(options, numCells, numComponents);
boundary = implicitBoundaryFaceData(options, numCells, numComponents);

for iComponent = 1:numComponents
    [matrixRows, matrixCols, matrixVals, rhs] = implicitSystemBase( ...
        state.component_moles(:, iComponent), sourceDeltaMoles(:, iComponent), ...
        numCells);
    tripletCount = numCells;

    for iFace = 1:size(internal.faceCells, 1)
        leftCell = internal.faceCells(iFace, 1);
        rightCell = internal.faceCells(iFace, 2);
        velocityArea = internal.velocityCmS(iFace) .* internal.areaCm2(iFace);
        diffusiveConductance = internal.diffusionCoefficient(iFace, iComponent) .* ...
            internal.areaCm2(iFace) ./ internal.distanceCm(iFace);
        if velocityArea >= 0
            alphaLeft = (velocityArea + diffusiveConductance) ./ waterVolume(leftCell);
            alphaRight = -diffusiveConductance ./ waterVolume(rightCell);
        else
            alphaLeft = diffusiveConductance ./ waterVolume(leftCell);
            alphaRight = (velocityArea - diffusiveConductance) ./ waterVolume(rightCell);
        end

        [matrixRows, matrixCols, matrixVals, tripletCount] = addTriplet( ...
            matrixRows, matrixCols, matrixVals, tripletCount, ...
            leftCell, leftCell, dtSeconds .* alphaLeft);
        [matrixRows, matrixCols, matrixVals, tripletCount] = addTriplet( ...
            matrixRows, matrixCols, matrixVals, tripletCount, ...
            leftCell, rightCell, dtSeconds .* alphaRight);
        [matrixRows, matrixCols, matrixVals, tripletCount] = addTriplet( ...
            matrixRows, matrixCols, matrixVals, tripletCount, ...
            rightCell, leftCell, -dtSeconds .* alphaLeft);
        [matrixRows, matrixCols, matrixVals, tripletCount] = addTriplet( ...
            matrixRows, matrixCols, matrixVals, tripletCount, ...
            rightCell, rightCell, -dtSeconds .* alphaRight);
    end

    for iFace = 1:numel(boundary.faceCells)
        cellId = boundary.faceCells(iFace);
        velocityArea = boundary.velocityCmS(iFace) .* boundary.areaCm2(iFace);
        diffusiveConductance = boundary.diffusionCoefficient(iFace, iComponent) .* ...
            boundary.areaCm2(iFace) ./ boundary.distanceCm(iFace);
        boundaryConcentration = boundary.concentration(iFace, iComponent);
        if velocityArea >= 0
            constantFlux = (velocityArea + diffusiveConductance) .* boundaryConcentration;
            coefficient = -diffusiveConductance ./ waterVolume(cellId);
        else
            constantFlux = diffusiveConductance .* boundaryConcentration;
            coefficient = (velocityArea - diffusiveConductance) ./ waterVolume(cellId);
        end
        rhs(cellId) = rhs(cellId) + dtSeconds .* constantFlux;
        [matrixRows, matrixCols, matrixVals, tripletCount] = addTriplet( ...
            matrixRows, matrixCols, matrixVals, tripletCount, ...
            cellId, cellId, -dtSeconds .* coefficient);
    end

    A = sparse(matrixRows(1:tripletCount), matrixCols(1:tripletCount), ...
        matrixVals(1:tripletCount), numCells, numCells);
    componentMoles1(:, iComponent) = A \ rhs;
end

finalState = state;
finalState.component_moles = componentMoles1;
internalFluxDeltaMoles = zeros(numCells, numComponents);
if hasAnyInternalFluxOption(options)
    internalFluxDeltaMoles = computeInternalFluxDelta( ...
        options, finalState, geometry, numCells, numComponents, dtSeconds);
end
[boundaryDeltaMoles, boundaryAdvectiveDeltaMoles, boundaryDiffusiveDeltaMoles] = ...
    computeBoundaryDelta(options, finalState, geometry, numCells, numComponents, dtSeconds);
end

function [rows, cols, vals, rhs] = implicitSystemBase(componentMoles, ...
    sourceDeltaMoles, numCells)
maxExtraTriplets = max(1, 8 .* numCells);
rows = zeros(numCells + maxExtraTriplets, 1);
cols = zeros(numCells + maxExtraTriplets, 1);
vals = zeros(numCells + maxExtraTriplets, 1);
rows(1:numCells) = (1:numCells).';
cols(1:numCells) = (1:numCells).';
vals(1:numCells) = 1;
rhs = componentMoles(:) + sourceDeltaMoles(:);
end

function [rows, cols, vals, count] = addTriplet(rows, cols, vals, count, row, col, value)
count = count + 1;
if count > numel(rows)
    growBy = max(numel(rows), 1024);
    rows(end + growBy, 1) = 0;
    cols(end + growBy, 1) = 0;
    vals(end + growBy, 1) = 0;
end
rows(count) = row;
cols(count) = col;
vals(count) = value;
end

function data = implicitInternalFaceData(options, numCells, numComponents)
data = struct();
data.faceCells = zeros(0, 2);
data.areaCm2 = zeros(0, 1);
data.distanceCm = zeros(0, 1);
data.velocityCmS = zeros(0, 1);
data.diffusionCoefficient = zeros(0, numComponents);
if ~hasAnyInternalFluxOption(options) || ...
        (isfield(options, 'internal_face_flux_mol_s') && ...
        ~isempty(options.internal_face_flux_mol_s))
    return;
end
if ~isfield(options, 'internal_face_cells') || isempty(options.internal_face_cells)
    error('RTSPHEM:Transport:MissingInternalFluxOption', ...
        'internal_face_cells is required for implicit internal transport.');
end
faceCells = options.internal_face_cells;
if size(faceCells, 2) ~= 2 || any(faceCells(:) < 1) || ...
        any(faceCells(:) > numCells) || any(faceCells(:) ~= round(faceCells(:)))
    error('RTSPHEM:Transport:InvalidInternalFaces', ...
        'internal_face_cells must reference valid cell indices.');
end
numFaces = size(faceCells, 1);
data.faceCells = faceCells;
data.areaCm2 = columnOption(options.internal_face_area_cm2, numFaces, ...
    'internal_face_area_cm2');
data.distanceCm = columnOption(options.internal_face_distance_cm, numFaces, ...
    'internal_face_distance_cm');
data.velocityCmS = columnOption(options.internal_face_velocity_cm_s, numFaces, ...
    'internal_face_velocity_cm_s');
data.diffusionCoefficient = faceComponentOption(options.diffusion_coefficient_cm2_s, ...
    numFaces, numComponents, 'diffusion_coefficient_cm2_s');
end

function data = implicitBoundaryFaceData(options, numCells, numComponents)
data = struct();
data.faceCells = zeros(0, 1);
data.areaCm2 = zeros(0, 1);
data.distanceCm = zeros(0, 1);
data.velocityCmS = zeros(0, 1);
data.concentration = zeros(0, numComponents);
data.diffusionCoefficient = zeros(0, numComponents);
boundarySpecificFields = {'boundary_face_cells', 'boundary_face_area_cm2', ...
    'boundary_face_distance_cm', 'boundary_face_velocity_cm_s', ...
    'boundary_concentration_mol_cm3'};
hasAnySpecific = false(size(boundarySpecificFields));
for iField = 1:numel(boundarySpecificFields)
    hasAnySpecific(iField) = isfield(options, boundarySpecificFields{iField}) && ...
        ~isempty(options.(boundarySpecificFields{iField}));
end
if ~any(hasAnySpecific)
    return;
end
faceCells = options.boundary_face_cells(:);
numFaces = numel(faceCells);
if any(faceCells < 1) || any(faceCells > numCells) || any(faceCells ~= round(faceCells))
    error('RTSPHEM:Transport:InvalidBoundaryFaces', ...
        'boundary_face_cells must reference valid cell indices.');
end
data.faceCells = faceCells;
data.areaCm2 = columnOption(options.boundary_face_area_cm2, numFaces, ...
    'boundary_face_area_cm2');
data.distanceCm = columnOption(options.boundary_face_distance_cm, numFaces, ...
    'boundary_face_distance_cm');
data.velocityCmS = columnOption(options.boundary_face_velocity_cm_s, numFaces, ...
    'boundary_face_velocity_cm_s');
data.concentration = faceComponentOption(options.boundary_concentration_mol_cm3, ...
    numFaces, numComponents, 'boundary_concentration_mol_cm3');
data.diffusionCoefficient = faceComponentOption(options.diffusion_coefficient_cm2_s, ...
    numFaces, numComponents, 'diffusion_coefficient_cm2_s');
end

function tf = hasAnyInternalFluxOption(options)
internalFields = {'internal_face_cells', 'internal_face_flux_mol_s', ...
    'internal_face_area_cm2', 'internal_face_distance_cm', ...
    'internal_face_velocity_cm_s'};
tf = false;
for iField = 1:numel(internalFields)
    tf = tf || (isfield(options, internalFields{iField}) && ...
        ~isempty(options.(internalFields{iField})));
end
end

function internalFluxDeltaMoles = computeInternalFluxDelta(options, state, geometry, ...
    numCells, numComponents, dtSeconds)
if ~isfield(options, 'internal_face_cells') || isempty(options.internal_face_cells)
    error('RTSPHEM:Transport:MissingInternalFluxOption', ...
        'internal_face_cells is required for internal transport fluxes.');
end

faceCells = options.internal_face_cells;
if size(faceCells, 2) ~= 2
    error('RTSPHEM:Transport:InvalidInternalFaces', ...
        'internal_face_cells must have two columns.');
end
if any(faceCells(:) < 1) || any(faceCells(:) > numCells) || any(faceCells(:) ~= round(faceCells(:)))
    error('RTSPHEM:Transport:InvalidInternalFaces', ...
        'internal_face_cells must reference valid cell indices.');
end
faceFluxMolS = internalFaceFluxMolS(options, state, geometry, faceCells, numComponents);

internalFluxDeltaMoles = zeros(numCells, numComponents);
for iFace = 1:size(faceCells, 1)
    leftCell = faceCells(iFace, 1);
    rightCell = faceCells(iFace, 2);
    fluxDelta = faceFluxMolS(iFace, :) .* dtSeconds;
    internalFluxDeltaMoles(leftCell, :) = internalFluxDeltaMoles(leftCell, :) - fluxDelta;
    internalFluxDeltaMoles(rightCell, :) = internalFluxDeltaMoles(rightCell, :) + fluxDelta;
end
end

function faceFluxMolS = internalFaceFluxMolS(options, state, geometry, faceCells, numComponents)
if isfield(options, 'internal_face_flux_mol_s') && ...
        ~isempty(options.internal_face_flux_mol_s)
    faceFluxMolS = options.internal_face_flux_mol_s;
    if size(faceFluxMolS, 1) ~= size(faceCells, 1) || ...
            size(faceFluxMolS, 2) ~= numComponents
        error('RTSPHEM:Transport:OptionSizeMismatch', ...
            'internal_face_flux_mol_s must have one row per face and one column per component.');
    end
    if any(~isfinite(faceFluxMolS(:)))
        error('RTSPHEM:Transport:NonfiniteOption', ...
            'internal_face_flux_mol_s must contain only finite values.');
    end
    return;
end

requiredFields = {'internal_face_area_cm2', 'internal_face_distance_cm', ...
    'internal_face_velocity_cm_s', 'diffusion_coefficient_cm2_s'};
for iField = 1:numel(requiredFields)
    if ~isfield(options, requiredFields{iField}) || isempty(options.(requiredFields{iField}))
        error('RTSPHEM:Transport:MissingInternalFluxOption', ...
            'Computed internal flux requires %s.', requiredFields{iField});
    end
end

numFaces = size(faceCells, 1);
areaCm2 = columnOption(options.internal_face_area_cm2, numFaces, ...
    'internal_face_area_cm2');
distanceCm = columnOption(options.internal_face_distance_cm, numFaces, ...
    'internal_face_distance_cm');
velocityCmS = columnOption(options.internal_face_velocity_cm_s, numFaces, ...
    'internal_face_velocity_cm_s');
if any(areaCm2 < 0) || any(distanceCm <= 0)
    error('RTSPHEM:Transport:InvalidInternalFaceMeasure', ...
        'Internal face areas must be nonnegative and distances must be positive.');
end
diffusionCoefficient = faceComponentOption(options.diffusion_coefficient_cm2_s, ...
    numFaces, numComponents, 'diffusion_coefficient_cm2_s');
if any(diffusionCoefficient(:) < 0)
    error('RTSPHEM:Transport:InvalidInternalFaceMeasure', ...
        'Internal diffusion coefficients must be nonnegative.');
end

waterVolume = geometry.water_volume_cm3(:);
concentration = zeros(size(state.component_moles));
activeWater = waterVolume > 0;
concentration(activeWater, :) = state.component_moles(activeWater, :) ./ ...
    waterVolume(activeWater);

faceFluxMolS = zeros(numFaces, numComponents);
for iFace = 1:numFaces
    leftCell = faceCells(iFace, 1);
    rightCell = faceCells(iFace, 2);
    if velocityCmS(iFace) >= 0
        upwindConcentration = concentration(leftCell, :);
    else
        upwindConcentration = concentration(rightCell, :);
    end
    advectiveFlux = velocityCmS(iFace) .* areaCm2(iFace) .* upwindConcentration;
    diffusiveFlux = diffusionCoefficient(iFace, :) .* areaCm2(iFace) ./ ...
        distanceCm(iFace) .* ...
        (concentration(leftCell, :) - concentration(rightCell, :));
    faceFluxMolS(iFace, :) = advectiveFlux + diffusiveFlux;
end
end

function [boundaryDeltaMoles, advectiveDeltaMoles, diffusiveDeltaMoles] = ...
    computeBoundaryDelta(options, state, geometry, numCells, numComponents, dtSeconds)
boundaryDeltaMoles = zeros(numCells, numComponents);
advectiveDeltaMoles = zeros(numCells, numComponents);
diffusiveDeltaMoles = zeros(numCells, numComponents);

boundarySpecificFields = {'boundary_face_cells', 'boundary_face_area_cm2', ...
    'boundary_face_distance_cm', 'boundary_face_velocity_cm_s', ...
    'boundary_concentration_mol_cm3'};
hasAnySpecific = false(size(boundarySpecificFields));
for iField = 1:numel(boundarySpecificFields)
    hasAnySpecific(iField) = isfield(options, boundarySpecificFields{iField}) && ...
        ~isempty(options.(boundarySpecificFields{iField}));
end
if ~any(hasAnySpecific)
    return;
end
boundaryFields = [boundarySpecificFields, {'diffusion_coefficient_cm2_s'}];
hasAny = false(size(boundaryFields));
for iField = 1:numel(boundaryFields)
    hasAny(iField) = isfield(options, boundaryFields{iField}) && ...
        ~isempty(options.(boundaryFields{iField}));
end
if ~all(hasAny)
    missing = strjoin(boundaryFields(~hasAny), ', ');
    error('RTSPHEM:Transport:MissingBoundaryFluxOption', ...
        'Boundary Dirichlet flux requires: %s.', missing);
end

faceCells = options.boundary_face_cells(:);
numFaces = numel(faceCells);
if any(faceCells < 1) || any(faceCells > numCells) || any(faceCells ~= round(faceCells))
    error('RTSPHEM:Transport:InvalidBoundaryFaces', ...
        'boundary_face_cells must reference valid cell indices.');
end
areaCm2 = columnOption(options.boundary_face_area_cm2, numFaces, ...
    'boundary_face_area_cm2');
distanceCm = columnOption(options.boundary_face_distance_cm, numFaces, ...
    'boundary_face_distance_cm');
velocityCmS = columnOption(options.boundary_face_velocity_cm_s, numFaces, ...
    'boundary_face_velocity_cm_s');
if any(areaCm2 < 0) || any(distanceCm <= 0)
    error('RTSPHEM:Transport:InvalidBoundaryMeasure', ...
        'Boundary areas must be nonnegative and distances must be positive.');
end

boundaryConcentration = faceComponentOption(options.boundary_concentration_mol_cm3, ...
    numFaces, numComponents, 'boundary_concentration_mol_cm3');
diffusionCoefficient = faceComponentOption(options.diffusion_coefficient_cm2_s, ...
    numFaces, numComponents, 'diffusion_coefficient_cm2_s');
if any(boundaryConcentration(:) < 0) || any(diffusionCoefficient(:) < 0)
    error('RTSPHEM:Transport:InvalidBoundaryMeasure', ...
        'Boundary concentrations and diffusion coefficients must be nonnegative.');
end

waterVolume = geometry.water_volume_cm3(:);
cellConcentration = zeros(numCells, numComponents);
activeWater = waterVolume > 0;
cellConcentration(activeWater, :) = state.component_moles(activeWater, :) ./ ...
    waterVolume(activeWater);

for iFace = 1:numFaces
    cellId = faceCells(iFace);
    advectiveFlux = velocityCmS(iFace) .* areaCm2(iFace) .* ...
        boundaryConcentration(iFace, :);
    diffusiveFlux = diffusionCoefficient(iFace, :) .* areaCm2(iFace) ./ ...
        distanceCm(iFace) .* ...
        (boundaryConcentration(iFace, :) - cellConcentration(cellId, :));
    advectiveDelta = advectiveFlux .* dtSeconds;
    diffusiveDelta = diffusiveFlux .* dtSeconds;
    advectiveDeltaMoles(cellId, :) = advectiveDeltaMoles(cellId, :) + advectiveDelta;
    diffusiveDeltaMoles(cellId, :) = diffusiveDeltaMoles(cellId, :) + diffusiveDelta;
    boundaryDeltaMoles(cellId, :) = boundaryDeltaMoles(cellId, :) + ...
        advectiveDelta + diffusiveDelta;
end
end

function values = columnOption(values, numRows, fieldName)
values = values(:);
if numel(values) == 1
    values = repmat(values, numRows, 1);
end
if numel(values) ~= numRows || any(~isfinite(values))
    error('RTSPHEM:Transport:OptionSizeMismatch', ...
        '%s must be scalar or have one value per boundary face.', fieldName);
end
end

function values = faceComponentOption(values, numRows, numCols, fieldName)
values = values(:,:);
if isscalar(values)
    values = repmat(values, numRows, numCols);
elseif isrow(values) && numel(values) == numCols
    values = repmat(values, numRows, 1);
end
if ~isequal(size(values), [numRows numCols]) || any(~isfinite(values(:)))
    error('RTSPHEM:Transport:OptionSizeMismatch', ...
        '%s must be scalar, one row per component, or size [%d %d].', ...
        fieldName, numRows, numCols);
end
end

function ledger = buildLedger(componentMoles0, componentMoles1, ...
    sourceDeltaMoles, internalFluxDeltaMoles, boundaryDeltaMoles, ...
    boundaryAdvectiveDeltaMoles, boundaryDiffusiveDeltaMoles, ...
    roundoffSuppressedMoles, dtSeconds, componentNames)
sourceDeltaTotal = sum(sourceDeltaMoles, 1);
internalFluxDeltaTotal = sum(internalFluxDeltaMoles, 1);
boundaryDeltaTotal = sum(boundaryDeltaMoles, 1);
boundaryAdvectiveDeltaTotal = sum(boundaryAdvectiveDeltaMoles, 1);
boundaryDiffusiveDeltaTotal = sum(boundaryDiffusiveDeltaMoles, 1);
stateDeltaTotal = sum(componentMoles1 - componentMoles0, 1);
componentResidual = stateDeltaTotal - sourceDeltaTotal - internalFluxDeltaTotal - ...
    boundaryDeltaTotal;

ledger = struct();
ledger.component_names = componentNames;
ledger.dt_s = dtSeconds;
ledger.initial_moles_total = sum(componentMoles0, 1);
ledger.final_moles_total = sum(componentMoles1, 1);
ledger.source_delta_moles_total = sourceDeltaTotal;
ledger.internal_flux_delta_moles_total = internalFluxDeltaTotal;
ledger.boundary_delta_moles_total = boundaryDeltaTotal;
ledger.boundary_advective_delta_moles_total = boundaryAdvectiveDeltaTotal;
ledger.boundary_diffusive_delta_moles_total = boundaryDiffusiveDeltaTotal;
ledger.roundoff_suppressed_moles = roundoffSuppressedMoles;
ledger.roundoff_suppressed_moles_total = sum(roundoffSuppressedMoles(:), 'omitnan');
ledger.roundoff_suppressed_entries = nnz(roundoffSuppressedMoles);
ledger.component_residual_moles = componentResidual;
ledger.max_abs_component_residual_moles = max(abs(componentResidual));
end
