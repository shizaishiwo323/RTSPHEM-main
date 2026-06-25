function [clusterIds, info] = ReactionClusterDiagnostics(options, numCells)
%REACTIONCLUSTERDIAGNOSTICS Build per-cell cluster IDs for chemistry logs.

clusterIds = nan(numCells, 1);
info = struct('count', 0, 'max_membership', 0, ...
    'overlapping_cell_count', 0);
if nargin < 1 || ~isstruct(options) || ...
        ~isfield(options, 'reactionClusters') || isempty(options.reactionClusters)
    return;
end
if ~(isscalar(numCells) && isfinite(numCells) && numCells >= 0 && ...
        numCells == round(numCells))
    error('RTSPHEM:Chemistry:InvalidReactionClusters', ...
        'numCells must be a nonnegative integer scalar.');
end

clusters = options.reactionClusters;
if ~isstruct(clusters)
    error('RTSPHEM:Chemistry:InvalidReactionClusters', ...
        'options.reactionClusters must be a struct array.');
end

membershipCount = zeros(numCells, 1);
for iCluster = 1:numel(clusters)
    clusterCells = cellsForCluster(clusters(iCluster), numCells);
    membershipCount(clusterCells) = membershipCount(clusterCells) + 1;
    clusterIds(clusterCells) = iCluster;
end

info.count = numel(clusters);
if ~isempty(membershipCount)
    info.max_membership = max(membershipCount);
end
info.overlapping_cell_count = nnz(membershipCount > 1);
if info.overlapping_cell_count > 0
    error('RTSPHEM:Chemistry:OverlappingReactionClusters', ...
        'Reaction cluster membership must be disjoint for PHREEQC chemistry.');
end
end

function cellIds = cellsForCluster(cluster, numCells)
cellIds = [];
if isfield(cluster, 'source_cells') && ~isempty(cluster.source_cells)
    cellIds = [cellIds; cluster.source_cells(:)];
end
if isfield(cluster, 'member_cells') && ~isempty(cluster.member_cells)
    cellIds = [cellIds; cluster.member_cells(:)];
end
cellIds = unique(cellIds(:));
if isempty(cellIds)
    return;
end
if any(~isfinite(cellIds)) || any(cellIds < 1) || any(cellIds > numCells) || ...
        any(cellIds ~= round(cellIds))
    error('RTSPHEM:Chemistry:InvalidReactionClusters', ...
        'Reaction cluster cells must be valid integer cell indices.');
end
end
