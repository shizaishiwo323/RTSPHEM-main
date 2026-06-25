function result = BuildDisjointReactionClusters(sourceCells, candidateMembers, numCells)
%BUILDDISJOINTREACTIONCLUSTERS Merge overlapping reaction neighborhoods.
%
% candidateMembers is a cell array, one entry per source cell. Any two
% candidates sharing at least one member cell are merged into one cluster.

sourceCells = sourceCells(:);
if ~iscell(candidateMembers) || numel(candidateMembers) ~= numel(sourceCells)
    error('RTSPHEM:Chemistry:InvalidClusterInput', ...
        'candidateMembers must be a cell array matching sourceCells.');
end
if ~(isscalar(numCells) && isfinite(numCells) && numCells >= 0 && numCells == round(numCells))
    error('RTSPHEM:Chemistry:InvalidClusterInput', ...
        'numCells must be a nonnegative integer scalar.');
end

numSources = numel(sourceCells);
parents = 1:numSources;
cellOwner = zeros(numCells, 1);
normalizedMembers = cell(numSources, 1);

for iSource = 1:numSources
    members = unique([sourceCells(iSource); candidateMembers{iSource}(:)]);
    validateCellIds([sourceCells(iSource); members], numCells);
    normalizedMembers{iSource} = members(:);

    for iMember = 1:numel(members)
        cellId = members(iMember);
        if cellOwner(cellId) == 0
            cellOwner(cellId) = iSource;
        else
            parents = unionRoots(parents, iSource, cellOwner(cellId));
        end
    end
end

rootIds = zeros(numSources, 1);
for iSource = 1:numSources
    [rootIds(iSource), parents] = findRoot(parents, iSource);
end
uniqueRoots = unique(rootIds, 'stable');

clusters = struct('source_cells', {}, 'member_cells', {});
for iRoot = 1:numel(uniqueRoots)
    rootId = uniqueRoots(iRoot);
    groupSourceIndices = find(rootIds == rootId);
    groupSources = sort(sourceCells(groupSourceIndices));
    groupMembers = [];
    for iGroup = 1:numel(groupSourceIndices)
        groupMembers = [groupMembers; normalizedMembers{groupSourceIndices(iGroup)}(:)]; %#ok<AGROW>
    end
    clusters(end + 1).source_cells = unique(groupSources(:)); %#ok<AGROW>
    clusters(end).member_cells = unique(groupMembers(:));
end

clusters = sortClusters(clusters);
clusterIdByCell = zeros(numCells, 1);
membershipCount = zeros(numCells, 1);
for iCluster = 1:numel(clusters)
    members = clusters(iCluster).member_cells;
    clusterIdByCell(members) = iCluster;
    membershipCount(members) = membershipCount(members) + 1;
end

result = struct();
result.clusters = clusters;
result.cluster_id_by_cell = clusterIdByCell;
result.membership_count = membershipCount;
end

function validateCellIds(cellIds, numCells)
if any(~isfinite(cellIds)) || any(cellIds < 1) || any(cellIds > numCells) || ...
        any(cellIds ~= round(cellIds))
    error('RTSPHEM:Chemistry:InvalidClusterCell', ...
        'Cluster cells must be valid integer cell indices.');
end
end

function parents = unionRoots(parents, firstIndex, secondIndex)
[firstRoot, parents] = findRoot(parents, firstIndex);
[secondRoot, parents] = findRoot(parents, secondIndex);
if firstRoot ~= secondRoot
    keepRoot = min(firstRoot, secondRoot);
    dropRoot = max(firstRoot, secondRoot);
    parents(dropRoot) = keepRoot;
end
end

function [root, parents] = findRoot(parents, index)
root = index;
while parents(root) ~= root
    root = parents(root);
end
while parents(index) ~= index
    nextIndex = parents(index);
    parents(index) = root;
    index = nextIndex;
end
end

function clusters = sortClusters(clusters)
if isempty(clusters)
    return;
end
sortKeys = zeros(numel(clusters), 1);
for iCluster = 1:numel(clusters)
    sortKeys(iCluster) = min(clusters(iCluster).member_cells);
end
[~, order] = sort(sortKeys);
clusters = clusters(order);
end
