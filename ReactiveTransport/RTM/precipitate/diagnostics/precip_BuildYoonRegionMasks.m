function masks = precip_BuildYoonRegionMasks(spec)
% precip_BuildYoonRegionMasks - Build Yoon first-pore region masks.
%
% Inputs:
%   spec - benchmark spec, optionally with regionMasks struct.
%
% Output:
%   masks - firstPoreMask, firstThreePoresMask, and source metadata.

if isfield(spec, 'regionMasks') && ~isempty(spec.regionMasks)
    masks = explicitRegionMasks(spec.regionMasks, spec);
    return;
end
if isfield(spec, 'regionMaskFile') && ~isempty(spec.regionMaskFile)
    masks = precip_LoadYoonRegionMasks(spec.regionMaskFile, spec);
    return;
end

xCenters = ((1:spec.numX) - 0.5) .* spec.dx_cm;
firstPoreXMax = getFieldOrDefault(spec, 'firstPoreXMax_cm', 0.078);
firstThreePoresXMax = getFieldOrDefault(spec, 'firstThreePoresXMax_cm', 0.174);
masks = struct();
masks.firstPoreMask = repmat(xCenters <= firstPoreXMax, spec.numY, 1);
masks.firstThreePoresMask = repmat(xCenters <= firstThreePoresXMax, ...
    spec.numY, 1);
masks.source = 'x_window_fallback';
end

function masks = explicitRegionMasks(inputMasks, spec)
masks = struct();
masks.firstPoreMask = validateMask(inputMasks.firstPoreMask, spec, ...
    'firstPoreMask');
masks.firstThreePoresMask = validateMask(inputMasks.firstThreePoresMask, ...
    spec, 'firstThreePoresMask');
masks.source = getFieldOrDefault(inputMasks, 'source', 'explicit_region_masks');
end

function mask = validateMask(mask, spec, fieldName)
if ~isequal(size(mask), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidYoonRegionMask', ...
        '%s must have size [numY, numX].', fieldName);
end
mask = logical(mask);
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
