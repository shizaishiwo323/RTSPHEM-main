function precip_ValidateSplitInletMeshAlignment(leftBoundaryOrGrid, config, epsValue)
% precip_ValidateSplitInletMeshAlignment - Ensure split-inlet endpoints fall on left mesh nodes.
%
% The split inlet assigns CaCl2 and Na2CO3 streams by y-coordinate. Each
% yRange endpoint must coincide with a left-boundary mesh node; otherwise a
% boundary edge can straddle two inlet chemistries and smear the mixing line.

if nargin < 3 || isempty(epsValue)
    epsValue = eps;
end

leftNodeY = extractLeftBoundaryY(leftBoundaryOrGrid, epsValue);
requiredY = requiredSplitEndpointY(config);
meshTolerance = max(epsValue * 10, 1e-10);

for iY = 1:numel(requiredY)
    if ~any(abs(leftNodeY - requiredY(iY)) <= meshTolerance)
        error('RTSPHEM:Precipitate:SplitInletNotMeshAligned', ...
            ['split_left_inlet yRange endpoint %.12g cm is not aligned ', ...
            'with a left-boundary mesh node. Adjust splitInletY/yRange or mesh partitions.'], requiredY(iY));
    end
end
end

function leftNodeY = extractLeftBoundaryY(leftBoundaryOrGrid, epsValue)
if isnumeric(leftBoundaryOrGrid)
    coordinates = leftBoundaryOrGrid;
elseif isstruct(leftBoundaryOrGrid) && isfield(leftBoundaryOrGrid, 'coordV')
    coordinates = leftBoundaryOrGrid.coordV;
elseif isobject(leftBoundaryOrGrid) && isprop(leftBoundaryOrGrid, 'coordV')
    coordinates = leftBoundaryOrGrid.coordV;
else
    error('RTSPHEM:Precipitate:InvalidSplitInletMesh', ...
        'Expected a coordinate matrix or grid object/struct with coordV.');
end

if size(coordinates, 2) < 2
    error('RTSPHEM:Precipitate:InvalidSplitInletMesh', ...
        'Split-inlet mesh coordinates must contain x and y columns.');
end
if size(coordinates, 1) == 2 && size(coordinates, 2) > 2
    coordinates = coordinates.';
end

leftNodeY = coordinates(coordinates(:, 1) < epsValue, 2);
if isempty(leftNodeY)
    error('RTSPHEM:Precipitate:SplitInletNotMeshAligned', ...
        'split_left_inlet requires at least one left-boundary mesh node.');
end
end

function requiredY = requiredSplitEndpointY(config)
inletA = cfgget(config, 'inletA', struct());
inletB = cfgget(config, 'inletB', struct());
splitInletY = cfgget(config, 'splitInletY', 0.5 * cfgget(config, 'lengthYAxis', 1));
aRange = cfgget(inletA, 'yRange', [0, splitInletY]);
bRange = cfgget(inletB, 'yRange', [splitInletY, cfgget(config, 'lengthYAxis', Inf)]);
requiredY = unique([min(aRange), max(aRange), min(bRange), max(bRange)]);
end

function value = cfgget(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end
