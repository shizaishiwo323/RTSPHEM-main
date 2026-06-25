function [newGeometry, rebuildInfo] = RebuildGeometryAfterLevelSetMove( ...
        mesh, oldGeometry, newLevelSet, geometryInfo, options)
%REBUILDGEOMETRYAFTERLEVELSETMOVE Recompute cut-cell metrics after moving level set.

if nargin < 5 || isempty(options)
    options = struct();
end

newGeometry = rtm.geometry.BuildCutCellMetrics(mesh, newLevelSet, options);
oldSolidVolume = sum(requireGeometryVector(oldGeometry, 'solid_volume_cm3'));
newSolidVolume = sum(requireGeometryVector(newGeometry, 'solid_volume_cm3'));
actualVolumeChange = newSolidVolume - oldSolidVolume;
expectedVolumeChange = getFieldOrDefault(geometryInfo, ...
    'expected_solid_volume_change_cm3', 0);
absoluteResidual = actualVolumeChange - expectedVolumeChange;
scale = max(abs(expectedVolumeChange), ...
    max(abs(actualVolumeChange), max(abs(oldSolidVolume), 1) * eps));

rebuildInfo = struct();
rebuildInfo.old_solid_volume_cm3 = oldSolidVolume;
rebuildInfo.new_solid_volume_cm3 = newSolidVolume;
rebuildInfo.actual_solid_volume_change_cm3 = actualVolumeChange;
rebuildInfo.expected_solid_volume_change_cm3 = expectedVolumeChange;
rebuildInfo.mineral_volume_closure_residual_cm3 = absoluteResidual;
rebuildInfo.mineral_volume_closure_relative_residual = ...
    abs(absoluteResidual) ./ scale;
end

function values = requireGeometryVector(geometry, fieldName)
if ~isstruct(geometry) || ~isfield(geometry, fieldName)
    error('RTSPHEM:Geometry:MissingGeometryField', ...
        'geometry.%s is required.', fieldName);
end
values = geometry.(fieldName)(:);
if any(~isfinite(values))
    error('RTSPHEM:Geometry:InvalidGeometryField', ...
        'geometry.%s must be finite.', fieldName);
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
