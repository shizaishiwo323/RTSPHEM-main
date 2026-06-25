function clusters = BuildReactionClusters(state, geometry, connectivity, options)
%BUILDREACTIONCLUSTERS Build disjoint reaction neighborhoods from geometry.

if nargin < 3 || isempty(connectivity)
    connectivity = struct();
end
if nargin < 4 || isempty(options)
    options = struct();
end
clusters = struct([]);
if ~logical(getFieldOrDefault(options, 'useReactionClusters', true))
    return;
end

numCells = size(state.component_moles, 1);
interfaceArea = geometry.interface_area_cm2(:);
waterVolume = geometry.water_volume_cm3(:);
sourceCells = find(interfaceArea > 0 & waterVolume > 0);
if isempty(sourceCells)
    return;
end

neighbors = getNeighbors(connectivity, numCells);
depthCells = getClusterDepth(options);
candidateMembers = cell(numel(sourceCells), 1);
for iSource = 1:numel(sourceCells)
    cellId = sourceCells(iSource);
    members = neighborhoodCells(cellId, neighbors, depthCells);
    members = members(waterVolume(members) > 0);
    candidateMembers{iSource} = members(:);
end
clusterResult = rtm.chemistry.BuildDisjointReactionClusters( ...
    sourceCells, candidateMembers, numCells);
clusters = clusterResult.clusters;
end

function depthCells = getClusterDepth(options)
depthCells = getFieldOrDefault(options, 'reactionClusterDepthCells', 1);
if ~(isscalar(depthCells) && isfinite(depthCells) && depthCells >= 0 && ...
        depthCells == round(depthCells))
    error('RTSPHEM:Chemistry:InvalidReactionClusterDepth', ...
        'options.reactionClusterDepthCells must be a nonnegative integer scalar.');
end
end

function members = neighborhoodCells(sourceCell, neighbors, depthCells)
members = sourceCell;
frontier = sourceCell;
for iDepth = 1:depthCells
    nextFrontier = zeros(0, 1);
    for iCell = 1:numel(frontier)
        nextFrontier = [nextFrontier; neighbors{frontier(iCell)}(:)]; %#ok<AGROW>
    end
    nextFrontier = setdiff(unique(nextFrontier(:)), members);
    members = unique([members; nextFrontier(:)]);
    frontier = nextFrontier;
    if isempty(frontier)
        break;
    end
end
members = members(:);
end

function neighbors = getNeighbors(connectivity, numCells)
neighbors = repmat({zeros(0, 1)}, numCells, 1);
if ~isstruct(connectivity) || ~isfield(connectivity, 'cell_neighbors') || ...
        isempty(connectivity.cell_neighbors)
    return;
end
neighbors = connectivity.cell_neighbors(:);
if numel(neighbors) ~= numCells
    error('RTSPHEM:Chemistry:ConnectivitySizeMismatch', ...
        'connectivity.cell_neighbors must contain one entry per cell.');
end
for iCell = 1:numCells
    ids = neighbors{iCell}(:);
    if any(ids < 1) || any(ids > numCells) || any(ids ~= round(ids))
        error('RTSPHEM:Chemistry:InvalidConnectivity', ...
            'cell_neighbors contains invalid cell indices.');
    end
    neighbors{iCell} = ids;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
