function geometry = ApplyMineralVolumeChange(geometry, geometryInfo)
%APPLYMINERALVOLUMECHANGE Update water/solid volumes from accepted mineral change.

if nargin < 2 || isempty(geometryInfo) || ...
        ~isfield(geometryInfo, 'cell_solid_volume_change_cm3')
    error('RTSPHEM:Geometry:MissingGeometryChange', ...
        'geometryInfo.cell_solid_volume_change_cm3 is required.');
end
solidChange = geometryInfo.cell_solid_volume_change_cm3(:);
solidVolume = requireVectorField(geometry, 'solid_volume_cm3', numel(solidChange));
waterVolume = requireVectorField(geometry, 'water_volume_cm3', numel(solidChange));
cellVolume = requireVectorField(geometry, 'cell_volume_cm3', numel(solidChange));

newSolidVolume = solidVolume + solidChange;
solidTolerance = max(1e-14, 1e-12 .* max(solidVolume, 0));
if any(newSolidVolume < -solidTolerance)
    error('RTSPHEM:Geometry:NegativeSolidVolume', ...
        'Mineral volume change would make solid_volume_cm3 negative.');
end
newSolidVolume(abs(newSolidVolume) <= solidTolerance) = 0;
newWaterVolume = cellVolume - newSolidVolume;
waterTolerance = max(1e-14, 1e-12 .* max(cellVolume, 0));
if any(newWaterVolume < -waterTolerance)
    error('RTSPHEM:Geometry:NegativeWaterVolume', ...
        'Mineral volume change would make water_volume_cm3 negative.');
end
newWaterVolume(abs(newWaterVolume) <= waterTolerance) = 0;
geometry.solid_volume_cm3 = reshape(newSolidVolume, size(geometry.solid_volume_cm3));
geometry.water_volume_cm3 = reshape(newWaterVolume, size(geometry.water_volume_cm3));
geometry.fluid_fraction = reshape(newWaterVolume ./ max(cellVolume, eps), ...
    size(geometry.water_volume_cm3));
geometry.final_surface_area_cm2 = sum(getVectorFieldOrDefault(geometry, ...
    'interface_area_cm2', zeros(numel(solidChange), 1), numel(solidChange)));
geometry.final_solid_volume_cm3 = sum(newSolidVolume);
end

function values = requireVectorField(structValue, fieldName, expectedLength)
if ~isstruct(structValue) || ~isfield(structValue, fieldName) || ...
        isempty(structValue.(fieldName))
    error('RTSPHEM:Geometry:MissingGeometryField', ...
        'geometry.%s is required.', fieldName);
end
values = structValue.(fieldName)(:);
if numel(values) ~= expectedLength
    error('RTSPHEM:Geometry:GeometrySizeMismatch', ...
        'geometry.%s must match geometryInfo cell count.', fieldName);
end
if any(~isfinite(values))
    error('RTSPHEM:Geometry:InvalidGeometryField', ...
        'geometry.%s must be finite.', fieldName);
end
end

function values = getVectorFieldOrDefault(structValue, fieldName, defaultValue, expectedLength)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    values = structValue.(fieldName)(:);
else
    values = defaultValue(:);
end
if numel(values) ~= expectedLength
    error('RTSPHEM:Geometry:GeometrySizeMismatch', ...
        'geometry.%s must match geometryInfo cell count.', fieldName);
end
values(~isfinite(values)) = 0;
values = max(values, 0);
end
