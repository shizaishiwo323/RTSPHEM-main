function flow = MapRt0FluxToCutCellFaces(rt0)
%MAPRT0FLUXTOCUTCELLFACES Convert RT0/HyPHM face fluxes to RTM flow fields.
%
% Input fluxes are volumetric cm^3/s. Positive internal flux moves from
% internal_face_cells(:,1) to internal_face_cells(:,2). Positive boundary
% flux is a magnitude; boundary_type determines whether it enters or leaves.

if nargin < 1 || ~isstruct(rt0)
    error('RTSPHEM:Flow:InvalidRt0FluxInput', ...
        'rt0 must be a struct.');
end

flow = struct();
flow.internal_face_cells = optionalMatrix(rt0, 'internal_face_cells', 2);
numInternalFaces = size(flow.internal_face_cells, 1);
flow.internal_face_area_cm2 = positiveColumn(rt0, ...
    'internal_face_area_cm2', numInternalFaces, 1);
flow.internal_face_distance_cm = positiveColumn(rt0, ...
    'internal_face_distance_cm', numInternalFaces, 1);
internalFlux = column(getFirstField(rt0, ...
    {'internal_face_flux_cm3_s', 'internal_flux_cm3_s', 'rt0_internal_flux_cm3_s'}, ...
    zeros(numInternalFaces, 1)), numInternalFaces, 'internal_face_flux_cm3_s');
flow.internal_face_flux_cm3_s = internalFlux;
flow.internal_face_velocity_cm_s = safeDivide(internalFlux, ...
    flow.internal_face_area_cm2);

flow.boundary_face_cells = optionalColumn(rt0, 'boundary_face_cells');
numBoundaryFaces = numel(flow.boundary_face_cells);
flow.boundary_face_area_cm2 = positiveColumn(rt0, ...
    'boundary_face_area_cm2', numBoundaryFaces, 1);
boundaryFlux = column(getFirstField(rt0, ...
    {'boundary_face_flux_cm3_s', 'boundary_flux_cm3_s', 'rt0_boundary_flux_cm3_s'}, ...
    zeros(numBoundaryFaces, 1)), numBoundaryFaces, 'boundary_face_flux_cm3_s');
flow.boundary_face_flux_cm3_s = boundaryFlux;
flow.boundary_face_velocity_cm_s = safeDivide(boundaryFlux, ...
    flow.boundary_face_area_cm2);
flow.boundary_type = boundaryTypes(rt0, numBoundaryFaces);
if isfield(rt0, 'boundary_face_distance_cm') && ~isempty(rt0.boundary_face_distance_cm)
    flow.boundary_face_distance_cm = positiveColumn(rt0, ...
        'boundary_face_distance_cm', numBoundaryFaces, 1);
end
end

function values = optionalMatrix(structValue, fieldName, numCols)
if isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    values = structValue.(fieldName);
else
    values = zeros(0, numCols);
end
if size(values, 2) ~= numCols || any(~isfinite(values(:)))
    error('RTSPHEM:Flow:InvalidRt0FluxInput', ...
        '%s must have %d columns and finite values.', fieldName, numCols);
end
end

function values = optionalColumn(structValue, fieldName)
if isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    values = structValue.(fieldName)(:);
else
    values = zeros(0, 1);
end
if any(~isfinite(values))
    error('RTSPHEM:Flow:InvalidRt0FluxInput', ...
        '%s must contain finite values.', fieldName);
end
end

function values = positiveColumn(structValue, fieldName, count, defaultValue)
if isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    values = structValue.(fieldName)(:);
else
    values = repmat(defaultValue, count, 1);
end
if numel(values) ~= count || any(~isfinite(values)) || any(values < 0)
    error('RTSPHEM:Flow:InvalidRt0FluxInput', ...
        '%s must contain %d nonnegative finite values.', fieldName, count);
end
end

function value = getFirstField(structValue, fieldNames, defaultValue)
value = defaultValue;
for iField = 1:numel(fieldNames)
    if isfield(structValue, fieldNames{iField}) && ...
            ~isempty(structValue.(fieldNames{iField}))
        value = structValue.(fieldNames{iField});
        return;
    end
end
end

function values = column(values, count, fieldName)
values = values(:);
if numel(values) ~= count || any(~isfinite(values))
    error('RTSPHEM:Flow:InvalidRt0FluxInput', ...
        '%s must contain %d finite values.', fieldName, count);
end
end

function values = safeDivide(flux, area)
values = zeros(size(flux));
active = area > 0;
values(active) = flux(active) ./ area(active);
end

function values = boundaryTypes(structValue, count)
if isfield(structValue, 'boundary_type') && ~isempty(structValue.boundary_type)
    values = string(structValue.boundary_type(:));
else
    values = repmat("dirichlet", count, 1);
end
if isscalar(values) && count > 1
    values = repmat(values, count, 1);
end
if numel(values) ~= count || any(~(values == "dirichlet" | values == "outflow"))
    error('RTSPHEM:Flow:InvalidRt0FluxInput', ...
        'boundary_type must be dirichlet or outflow.');
end
end
