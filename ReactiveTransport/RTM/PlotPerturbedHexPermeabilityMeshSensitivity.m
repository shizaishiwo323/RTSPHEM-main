function results = PlotPerturbedHexPermeabilityMeshSensitivity(geometryFile, outputDir)
% Upstream RTSPHEM Kxx convergence for one fixed geometry.
if nargin < 2 || isempty(outputDir), outputDir = fileparts(geometryFile); end
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
rtmDir = fileparts(mfilename('fullpath')); addpath(genpath(fileparts(rtmDir)));
s = load(geometryFile); resolution = [64 43; 96 64; 128 86; 192 128; 256 171; 320 214];
kxx = nan(size(resolution,1),1);
for i = 1:size(resolution,1)
    mesh = FoldedCartesianGrid(2, [0 s.lengthXAxis 0 s.lengthYAxis], resolution(i,:));
    level = -inf(size(mesh.coordinates,1),1);
    for j = 1:numel(s.circleRadii)
        level = max(level, s.circleRadii(j)-hypot(mesh.coordinates(:,1)-s.circleCenters(j,1), ...
            mesh.coordinates(:,2)-s.circleCenters(j,2)));
    end
    K = computePhysicalPermeabilityTensor(mesh.coordinates, mesh.triangles, level);
    kxx(i) = K(1,1) / 9.869233e-12;
end
results = table(resolution(:,1), resolution(:,2), kxx, ...
    'VariableNames', {'nx','ny','Kxx_mD'});
results.relativeChangeFromFinest = abs(kxx-kxx(end))/kxx(end);
writetable(results, fullfile(outputDir, 'M14_24_Kxx_mesh_sensitivity.csv'));
writetable(results, fullfile(outputDir, 'M14_24_Kxx_mesh_sensitivity.xlsx'));
f = figure('Visible','off','Color','w'); plot(resolution(:,1), kxx, 'o-','LineWidth',1.5);
xlabel('Mesh partitions in x'); ylabel('2-D K_{xx} (mD)'); grid on;
exportgraphics(f, fullfile(outputDir, 'M14_24_Kxx_mesh_sensitivity.png'), 'Resolution', 220); close(f);
end
