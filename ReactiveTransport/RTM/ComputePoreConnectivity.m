function [numComponents, largestFraction] = ComputePoreConnectivity(grid, fluidTriangles)
% ComputePoreConnectivity returns fluid-triangle component statistics.

fluidIds = find(fluidTriangles);
numFluid = numel(fluidIds);
if numFluid == 0
    numComponents = 0;
    largestFraction = 0;
    return;
end

fluidMap = zeros(grid.numT, 1);
fluidMap(fluidIds) = 1:numFluid;
internal = all(grid.T0E > 0, 2);
adjacent = grid.T0E(internal, :);
fluidEdge = fluidTriangles(adjacent(:, 1)) & fluidTriangles(adjacent(:, 2));
adjacent = adjacent(fluidEdge, :);
poreGraph = graph(fluidMap(adjacent(:, 1)), fluidMap(adjacent(:, 2)), [], numFluid);
componentId = conncomp(poreGraph).';
componentSizes = accumarray(componentId, 1);
numComponents = numel(componentSizes);
largestFraction = max(componentSizes) / numFluid;
end
