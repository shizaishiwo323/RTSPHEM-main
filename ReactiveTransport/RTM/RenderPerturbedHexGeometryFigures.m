function RenderPerturbedHexGeometryFigures(geometryDir)
% Render final vector and raster boundary plots from saved circle geometry.
files = dir(fullfile(geometryDir, '*_perturbed_hex_600x400um.mat'));
for i = 1:numel(files)
    s = load(fullfile(files(i).folder, files(i).name));
    f = figure('Visible','off','Color','w'); hold on; axis equal;
    xlim([0 s.lengthXAxis]); ylim([0 s.lengthYAxis]);
    t = linspace(0, 2*pi, 96);
    for j = 1:numel(s.circleRadii)
        plot(s.circleCenters(j,1)+s.circleRadii(j)*cos(t), ...
            s.circleCenters(j,2)+s.circleRadii(j)*sin(t), 'k-');
    end
    xlabel('x (cm)'); ylabel('y (cm)');
    title(sprintf('%s: RTM phi %.3f, 2-D Kxx %.1f mD, Iv %.3f', ...
        string(s.sampleName), s.rtmPorosity, s.rtmKxx_mD, s.Iv));
    stem = erase(files(i).name, '.mat');
    exportgraphics(f, fullfile(geometryDir, strcat(stem, '.png')), 'Resolution', 300);
    exportgraphics(f, fullfile(geometryDir, strcat(stem, '.svg')), 'ContentType','vector');
    close(f);
end
end
