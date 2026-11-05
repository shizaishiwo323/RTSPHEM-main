function K_cm2 = computePhysicalPermeabilityTensor(baseCoord_cm, triangles, levelSet, varargin)
% Apply the upstream unit-cell solver, then restore the physical length scale.

domainMin = min(baseCoord_cm, [], 1);
domainSpan = max(baseCoord_cm, [], 1) - domainMin;
physicalArea_cm2 = prod(domainSpan);
assert(isfinite(physicalArea_cm2) && physicalArea_cm2 > 0, ...
    'Permeability domain area must be positive.');

lengthScale_cm = sqrt(physicalArea_cm2);
referenceCoord = (baseCoord_cm - domainMin) / lengthScale_cm;
sourceGrid = Grid(referenceCoord, triangles);
foldedGrid = FoldedGrid(sourceGrid);
[isMapped, sourceIndex] = ismember(foldedGrid.coordV, referenceCoord, 'rows');
assert(all(isMapped), 'Folded permeability grid could not map to the level set.');
assert(numel(levelSet) == size(baseCoord_cm, 1), ...
    'Level-set and coordinate counts must match.');
referenceLevelSet = levelSet(:) / lengthScale_cm;

K_reference = computePermeabilityTensor( ...
    foldedGrid, referenceLevelSet(sourceIndex), 'noRefine', varargin{:});
referenceArea = sum(foldedGrid.areaT);
K_cm2 = (K_reference / referenceArea) * lengthScale_cm^2;
if ~any(ismember(varargin, 'MainDirectionOnly'))
    K_cm2 = 0.5 * (K_cm2 + K_cm2.');
end
end
