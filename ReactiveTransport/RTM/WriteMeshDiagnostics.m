function WriteMeshDiagnostics(outputDir, coord, triangles, levelSetData, interfaceSegments, lengthXAxis, lengthYAxis)
% WriteMeshDiagnostics - Export mesh statistics and a mesh preview figure.

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

stats = ComputeMeshDiagnostics(coord, triangles, levelSetData);
statsTable = struct2table(stats, 'AsArray', true);
writetable(statsTable, fullfile(outputDir, 'mesh_statistics.csv'));
writetable(statsTable, fullfile(outputDir, 'mesh_statistics.xlsx'));

fid = fopen(fullfile(outputDir, 'mesh_statistics.json'), 'w');
if fid == -1
    error('MATLAB:WriteMeshDiagnostics:OpenFailed', 'Cannot write mesh statistics JSON.');
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(stats, 'PrettyPrint', true));
clear cleaner;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 900]);
triplot(triangles, coord(:, 1), coord(:, 2), 'Color', [0.55 0.55 0.55], 'LineWidth', 0.35);
hold on;
for k = 1:numel(interfaceSegments)
    seg = interfaceSegments{k};
    plot(seg(1, :), seg(2, :), 'r-', 'LineWidth', 1.4);
end
axis equal tight;
xlim([0, lengthXAxis]);
ylim([0, lengthYAxis]);
xlabel('X (cm)');
ylabel('Y (cm)');
title(sprintf('RTM Mesh: %d nodes, %d triangles', stats.num_nodes, stats.num_triangles));
grid on;
print(fig, fullfile(outputDir, 'mesh_plot.png'), '-dpng', '-r300');
savefig(fig, fullfile(outputDir, 'mesh_plot.fig'));
close(fig);
end
