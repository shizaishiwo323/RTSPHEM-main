function options = ApplyFaceFluxesToTransportOptions(options, flow)
%APPLYFACEFLUXESTOTRANSPORTOPTIONS Copy unified flow fields into transport options.

if nargin < 1 || isempty(options)
    options = struct();
end
if nargin < 2 || isempty(flow)
    return;
end
internalFields = {'internal_face_cells', 'internal_face_area_cm2', ...
    'internal_face_distance_cm', 'internal_face_velocity_cm_s'};
for iField = 1:numel(internalFields)
    options = copyField(options, flow, internalFields{iField});
end
boundaryFields = {'boundary_face_cells', 'boundary_face_area_cm2', ...
    'boundary_face_distance_cm', 'boundary_face_velocity_cm_s', ...
    'boundary_type'};
for iField = 1:numel(boundaryFields)
    options = copyField(options, flow, boundaryFields{iField});
end
end

function options = copyField(options, source, fieldName)
if isstruct(source) && isfield(source, fieldName) && ~isempty(source.(fieldName))
    options.(fieldName) = source.(fieldName);
end
end
