function info = AdvanceGeometryFromMineralMoles(geometry, realizedMoles, options)
%ADVANCEGEOMETRYFROMMINERALMOLES Convert realized mineral dissolution to geometry diagnostics.
%
% This helper does not mutate a mesh. It reports the solid-volume change and
% interface-normal displacement implied by realized mineral moles.

if nargin < 3 || isempty(options)
    options = struct();
end
realizedMoles = realizedMoles(:);
interfaceArea = requireVectorField(geometry, 'interface_area_cm2', numel(realizedMoles));
solidVolume = requireVectorField(geometry, 'solid_volume_cm3', numel(realizedMoles));
interfaceH = getVectorFieldOrDefault(geometry, 'interface_length_scale_cm', ...
    getVectorFieldOrDefault(geometry, 'interface_h_cm', ...
    ones(numel(realizedMoles), 1), numel(realizedMoles)), numel(realizedMoles));

molarVolume = getFieldOrDefault(options, 'molarVolume_cm3_mol', 1);
maxDisplacementOverH = getFieldOrDefault(options, 'maxDisplacementOverH', 0.25);
if ~(isscalar(molarVolume) && isfinite(molarVolume) && molarVolume > 0)
    error('RTSPHEM:Geometry:InvalidMolarVolume', ...
        'molarVolume_cm3_mol must be a positive finite scalar.');
end

realizedMoles(~isfinite(realizedMoles)) = 0;
realizedMoles = max(realizedMoles, 0);
cellSolidVolumeChange = -realizedMoles .* molarVolume;
cellSolidVolumeAfter = solidVolume + cellSolidVolumeChange;
normalDisplacement = zeros(size(realizedMoles));
activeInterface = interfaceArea > 0;
normalDisplacement(activeInterface) = -cellSolidVolumeChange(activeInterface) ./ ...
    interfaceArea(activeInterface);
normalDisplacement(~isfinite(normalDisplacement)) = 0;
displacementOverH = zeros(size(realizedMoles));
activeH = interfaceH > 0;
displacementOverH(activeH) = normalDisplacement(activeH) ./ interfaceH(activeH);

info = struct();
info.realized_mineral_moles = realizedMoles;
info.molar_volume_cm3_mol = molarVolume;
info.solid_volume_before_cm3 = sum(solidVolume);
info.cell_solid_volume_change_cm3 = cellSolidVolumeChange;
info.expected_solid_volume_change_cm3 = sum(cellSolidVolumeChange);
info.cell_solid_volume_after_cm3 = cellSolidVolumeAfter;
info.solid_volume_after_cm3 = sum(max(cellSolidVolumeAfter, 0));
info.normal_displacement_cm = normalDisplacement;
info.displacement_over_h = displacementOverH;
info.max_displacement_over_h = max(displacementOverH);
info.accepted = true;
info.reject_reason = "";

if any(realizedMoles > 0 & interfaceArea <= 0)
    info.accepted = false;
    info.reject_reason = "mineral change without interface area";
elseif any(cellSolidVolumeAfter < -max(1e-14, 1e-12 .* max(solidVolume, 1)))
    info.accepted = false;
    info.reject_reason = "solid volume would become negative";
elseif info.max_displacement_over_h > maxDisplacementOverH
    info.accepted = false;
    info.reject_reason = "geometry displacement exceeds tolerance";
end
end

function values = requireVectorField(structValue, fieldName, expectedLength)
if ~isstruct(structValue) || ~isfield(structValue, fieldName)
    error('RTSPHEM:Geometry:MissingGeometryField', ...
        'geometry.%s is required.', fieldName);
end
values = structValue.(fieldName)(:);
if numel(values) ~= expectedLength
    error('RTSPHEM:Geometry:GeometrySizeMismatch', ...
        'geometry.%s must match realizedMoles length.', fieldName);
end
values(~isfinite(values)) = 0;
values = max(values, 0);
end

function values = getVectorFieldOrDefault(structValue, fieldName, defaultValue, expectedLength)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    values = structValue.(fieldName)(:);
else
    values = defaultValue(:);
end
if numel(values) ~= expectedLength
    error('RTSPHEM:Geometry:GeometrySizeMismatch', ...
        'geometry.%s must match realizedMoles length.', fieldName);
end
values(~isfinite(values)) = 0;
values = max(values, 0);
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
