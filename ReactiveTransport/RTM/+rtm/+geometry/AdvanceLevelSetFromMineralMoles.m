function [newLevelSet, moveInfo] = AdvanceLevelSetFromMineralMoles( ...
        mesh, oldLevelSet, geometry, realizedMoles, dtSeconds, options)
%ADVANCELEVELSETFROMMINERALMOLES Move a signed-distance level set from actual dissolution.
%
% Negative level-set values are fluid and positive values are solid, matching
% BuildCutCellMetrics. Calcite dissolution therefore subtracts the realized
% normal displacement from the signed-distance field.

if nargin < 6 || isempty(options)
    options = struct();
end
if ~(isscalar(dtSeconds) && isfinite(dtSeconds) && dtSeconds > 0)
    error('RTSPHEM:Geometry:InvalidTimeStep', ...
        'dtSeconds must be a positive finite scalar.');
end

oldLevelSet = oldLevelSet(:);
geometryInfo = rtm.geometry.AdvanceGeometryFromMineralMoles( ...
    geometry, realizedMoles, options);
moveInfo = geometryInfo;
moveInfo.dt_s = dtSeconds;
moveInfo.level_set_moved = false;
moveInfo.vertex_normal_speed_cm_s = zeros(size(oldLevelSet));
newLevelSet = oldLevelSet;

if ~geometryInfo.accepted
    return;
end

interfaceArea = requireGeometryVector(geometry, 'interface_area_cm2', numel(realizedMoles));
realizedMoles = max(realizedMoles(:), 0);
ratePerArea = zeros(size(realizedMoles));
active = interfaceArea > 0 & dtSeconds > 0;
ratePerArea(active) = realizedMoles(active) ./ ...
    (interfaceArea(active) .* dtSeconds);
normalSpeed = geometryInfo.molar_volume_cm3_mol .* ratePerArea;
vertexSpeed = rtm.geometry.MapCellNormalSpeedToVertices( ...
    mesh, geometry, normalSpeed, options);
if numel(vertexSpeed) ~= numel(oldLevelSet)
    error('RTSPHEM:Geometry:LevelSetSpeedSizeMismatch', ...
        'Mapped vertex speed must match level-set size.');
end

newLevelSet = oldLevelSet - vertexSpeed(:) .* dtSeconds;
[newGeometry, rebuildInfo] = rtm.geometry.RebuildGeometryAfterLevelSetMove( ...
    mesh, geometry, newLevelSet, geometryInfo, options);

moveInfo.rate_per_area_mol_cm2_s = ratePerArea;
moveInfo.normal_speed_cm_s = normalSpeed;
moveInfo.vertex_normal_speed_cm_s = vertexSpeed(:);
moveInfo.level_set_before = oldLevelSet;
moveInfo.level_set_after = newLevelSet;
moveInfo.level_set_moved = any(abs(newLevelSet - oldLevelSet) > 0);
moveInfo.new_geometry = newGeometry;
moveInfo.cell_actual_solid_volume_change_cm3 = ...
    newGeometry.solid_volume_cm3(:) - geometry.solid_volume_cm3(:);
moveInfo.actual_solid_volume_change_cm3 = ...
    rebuildInfo.actual_solid_volume_change_cm3;
moveInfo.actual_solid_volume_after_cm3 = rebuildInfo.new_solid_volume_cm3;
moveInfo.solid_volume_after_cm3 = rebuildInfo.new_solid_volume_cm3;
moveInfo.final_surface_area_cm2 = sum(newGeometry.interface_area_cm2);
moveInfo.mineral_volume_closure_residual_cm3 = ...
    rebuildInfo.mineral_volume_closure_residual_cm3;
moveInfo.mineral_volume_closure_relative_residual = ...
    rebuildInfo.mineral_volume_closure_relative_residual;
end

function values = requireGeometryVector(geometry, fieldName, expectedLength)
if ~isstruct(geometry) || ~isfield(geometry, fieldName)
    error('RTSPHEM:Geometry:MissingGeometryField', ...
        'geometry.%s is required.', fieldName);
end
values = geometry.(fieldName)(:);
if numel(values) ~= expectedLength
    error('RTSPHEM:Geometry:GeometrySizeMismatch', ...
        'geometry.%s must match realizedMoles length.', fieldName);
end
values(~isfinite(values)) = 0;
values = max(values, 0);
end
