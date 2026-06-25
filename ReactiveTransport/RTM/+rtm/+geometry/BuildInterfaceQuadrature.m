function quadrature = BuildInterfaceQuadrature(geometry, options)
%BUILDINTERFACEQUADRATURE Build interface quadrature from cut-cell metrics.
%
% The initial implementation uses one quadrature point per cut cell at the
% interface centroid, with weight equal to interface_area_cm2.

if nargin < 2
    options = struct(); %#ok<NASGU>
end
if ~isstruct(geometry) || ~isfield(geometry, 'interface_area_cm2')
    error('RTSPHEM:Geometry:MissingInterfaceArea', ...
        'geometry.interface_area_cm2 is required.');
end
if ~isfield(geometry, 'interface_centroid_cm')
    error('RTSPHEM:Geometry:MissingInterfaceCentroid', ...
        'geometry.interface_centroid_cm is required.');
end

interfaceArea = geometry.interface_area_cm2(:);
if any(~isfinite(interfaceArea))
    error('RTSPHEM:Geometry:InvalidInterfaceArea', ...
        'geometry.interface_area_cm2 must contain finite values.');
end
if any(interfaceArea < 0)
    error('RTSPHEM:Geometry:NegativeInterfaceArea', ...
        'geometry.interface_area_cm2 must be nonnegative.');
end
interfaceCentroid = geometry.interface_centroid_cm;
numCells = numel(interfaceArea);
if size(interfaceCentroid, 1) ~= numCells || size(interfaceCentroid, 2) ~= 2
    error('RTSPHEM:Geometry:InterfaceCentroidSizeMismatch', ...
        'interface_centroid_cm must have size nCells-by-2.');
end

active = interfaceArea > 0 & all(isfinite(interfaceCentroid), 2);
if isfield(geometry, 'water_volume_cm3') && ~isempty(geometry.water_volume_cm3)
    waterVolume = geometry.water_volume_cm3(:);
    if numel(waterVolume) ~= numCells
        error('RTSPHEM:Geometry:WaterVolumeSizeMismatch', ...
            'water_volume_cm3 must match interface_area_cm2 length.');
    end
    if any(~isfinite(waterVolume))
        error('RTSPHEM:Geometry:InvalidWaterVolume', ...
            'geometry.water_volume_cm3 must contain finite values.');
    end
    if any(waterVolume < 0)
        error('RTSPHEM:Geometry:NegativeWaterVolume', ...
            'geometry.water_volume_cm3 must be nonnegative.');
    end
    active = active & waterVolume > 0;
end
cellId = find(active);
quadrature = struct();
quadrature.cell_id = cellId(:);
quadrature.point_cm = interfaceCentroid(active, :);
quadrature.weight_cm2 = interfaceArea(active);

if isfield(geometry, 'interface_normal') && ...
        isequal(size(geometry.interface_normal), [numCells, 2])
    quadrature.normal = geometry.interface_normal(active, :);
else
    quadrature.normal = nan(numel(cellId), 2);
end
end
