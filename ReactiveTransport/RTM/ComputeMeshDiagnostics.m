function stats = ComputeMeshDiagnostics(coord, triangles, levelSetData)
% ComputeMeshDiagnostics - COMSOL-like summary statistics for the RTM mesh.

if nargin < 3
    levelSetData = [];
end

tri = triangles;
p1 = coord(tri(:, 1), :);
p2 = coord(tri(:, 2), :);
p3 = coord(tri(:, 3), :);

e12 = sqrt(sum((p2 - p1).^2, 2));
e23 = sqrt(sum((p3 - p2).^2, 2));
e31 = sqrt(sum((p1 - p3).^2, 2));
areas = 0.5 * abs((p2(:,1) - p1(:,1)) .* (p3(:,2) - p1(:,2)) - ...
                  (p3(:,1) - p1(:,1)) .* (p2(:,2) - p1(:,2)));
quality = 4 * sqrt(3) * areas ./ max(e12.^2 + e23.^2 + e31.^2, eps);

edges = sort([tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])], 2);
edges = unique(edges, 'rows');
edgeLengths = sqrt(sum((coord(edges(:, 2), :) - coord(edges(:, 1), :)).^2, 2));

stats = struct();
stats.num_nodes = size(coord, 1);
stats.num_triangles = size(tri, 1);
stats.num_edges = size(edges, 1);
stats.triangle_area_min = min(areas);
stats.triangle_area_mean = mean(areas);
stats.triangle_area_max = max(areas);
stats.edge_length_min = min(edgeLengths);
stats.edge_length_mean = mean(edgeLengths);
stats.edge_length_max = max(edgeLengths);
stats.triangle_quality_min = min(quality);
stats.triangle_quality_mean = mean(quality);
stats.triangle_quality_max = max(quality);
stats.total_area_cm2 = sum(areas);

if isempty(levelSetData)
    stats.num_pore_triangles = NaN;
    stats.num_solid_triangles = NaN;
    stats.num_interface_triangles = NaN;
else
    lsTri = levelSetData(tri);
    stats.num_pore_triangles = sum(all(lsTri < 0, 2));
    stats.num_solid_triangles = sum(all(lsTri >= 0, 2));
    stats.num_interface_triangles = sum(~all(lsTri < 0, 2) & ~all(lsTri >= 0, 2));
end
end
