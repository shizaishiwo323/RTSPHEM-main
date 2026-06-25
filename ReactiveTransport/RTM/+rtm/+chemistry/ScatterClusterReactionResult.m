function cellResult = ScatterClusterReactionResult(clusterResult, map, numCells)
%SCATTERCLUSTERREACTIONRESULT Scatter clustered PHREEQC output to cells.

cellResult = map.before;
cellResult.calciteDissolvedMoles = zeros(numCells, 1);
componentFields = {'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'na_total_mol_cm3', 'cl_total_mol_cm3', 'alkalinity_mol_cm3'};

for iCluster = 1:numel(map.member_cells)
    memberCells = map.member_cells{iCluster};
    sourceCells = map.source_cells{iCluster};
    memberWeights = normalizedWeights( ...
        map.reaction_water_volume_cm3(memberCells));
    clusterWater = sum(map.reaction_water_volume_cm3(memberCells), 'omitnan');

    for iField = 1:numel(componentFields)
        fieldName = componentFields{iField};
        beforeCluster = weightedClusterConcentration( ...
            map.before.(fieldName)(memberCells), ...
            map.reaction_water_volume_cm3(memberCells));
        afterCluster = resultValue(clusterResult, fieldName, iCluster, ...
            beforeCluster);
        totalDelta = (afterCluster - beforeCluster) .* clusterWater;
        allocatedDelta = totalDelta .* memberWeights;
        active = map.reaction_water_volume_cm3(memberCells) > 0;
        updated = map.before.(fieldName)(memberCells);
        updated(active) = updated(active) + allocatedDelta(active) ./ ...
            map.reaction_water_volume_cm3(memberCells(active));
        cellResult.(fieldName)(memberCells) = updated;
    end

    dissolvedMoles = resultValue(clusterResult, 'calciteDissolvedMoles', ...
        iCluster, sum(map.prescribed_calcite_dissolved_moles(sourceCells), ...
        'omitnan'));
    sourceWeights = sourceMineralWeights(map, sourceCells);
    cellResult.calciteDissolvedMoles(sourceCells) = ...
        cellResult.calciteDissolvedMoles(sourceCells) + ...
        dissolvedMoles .* sourceWeights;
end

cellResult.ca_mol_cm3 = cellResult.ca_total_mol_cm3;
cellResult.c_mol_cm3 = cellResult.c_total_mol_cm3;
cellResult.na_mol_cm3 = cellResult.na_total_mol_cm3;
cellResult.cl_mol_cm3 = cellResult.cl_total_mol_cm3;
cellResult.pH = scatterDiagnostic(clusterResult, map, numCells, 'pH', NaN);
cellResult.calciteSI = scatterDiagnostic(clusterResult, map, numCells, ...
    'calciteSI', NaN);
cellResult.chargeBalance = scatterDiagnostic(clusterResult, map, numCells, ...
    'chargeBalance', 0);
cellResult.failedCells = scatterFailedCells(clusterResult, map, numCells);
cellResult.errorMessage = getFieldOrDefault(clusterResult, 'errorMessage', "");
end

function failedCells = scatterFailedCells(result, map, numCells)
failedClusterIds = getFieldOrDefault(result, 'failedCells', []);
failedClusterIds = failedClusterIds(:);
failedClusterIds = failedClusterIds(isfinite(failedClusterIds));
failedCells = [];
for iFailed = 1:numel(failedClusterIds)
    clusterId = failedClusterIds(iFailed);
    if clusterId >= 1 && clusterId <= numel(map.member_cells) && ...
            clusterId == fix(clusterId)
        failedCells = [failedCells; map.member_cells{clusterId}(:)]; %#ok<AGROW>
    elseif clusterId >= 1 && clusterId == fix(clusterId)
        failedCells = [failedCells; clusterId]; %#ok<AGROW>
    end
end
failedCells = unique(failedCells, 'stable');
end

function value = resultValue(result, fieldName, index, defaultValue)
if isfield(result, fieldName) && ~isempty(result.(fieldName))
    values = result.(fieldName)(:);
    if numel(values) >= index
        value = values(index);
        return;
    end
end
value = defaultValue;
end

function values = scatterDiagnostic(result, map, numCells, fieldName, defaultValue)
values = repmat(defaultValue, numCells, 1);
if ~isfield(result, fieldName) || isempty(result.(fieldName))
    return;
end
sourceValues = result.(fieldName)(:);
for iCluster = 1:numel(map.member_cells)
    if numel(sourceValues) >= iCluster
        values(map.member_cells{iCluster}) = sourceValues(iCluster);
    end
end
end

function weights = normalizedWeights(values)
values = values(:);
total = sum(values, 'omitnan');
if total > 0
    weights = values ./ total;
else
    weights = ones(size(values)) ./ numel(values);
end
end

function weights = sourceMineralWeights(map, sourceCells)
prescribed = map.prescribed_calcite_dissolved_moles(sourceCells);
if sum(prescribed, 'omitnan') > 0
    weights = prescribed ./ sum(prescribed, 'omitnan');
    return;
end
capacity = map.calcite_moles(sourceCells);
if sum(capacity, 'omitnan') > 0
    weights = capacity ./ sum(capacity, 'omitnan');
    return;
end
weights = ones(numel(sourceCells), 1) ./ numel(sourceCells);
end

function value = weightedClusterConcentration(values, waterVolume)
totalWater = sum(waterVolume, 'omitnan');
if totalWater <= 0
    value = 0;
else
    value = sum(values(:) .* waterVolume(:), 'omitnan') ./ totalWater;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
