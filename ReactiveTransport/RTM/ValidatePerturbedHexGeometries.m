function report = ValidatePerturbedHexGeometries(geometryDir, nx, ny)
% Validate fixed-domain circle geometries before RTM loading.
if nargin < 2, nx = 128; end
if nargin < 3, ny = 86; end
files = dir(fullfile(geometryDir, '*_perturbed_hex_600x400um.mat'));
assert(numel(files) == 4, 'RTM:GeometryValidation:FileCount', 'Expected four geometry MAT files.');
rtmDir = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(rtmDir)));
rows = cell(numel(files), 9);
for i = 1:numel(files)
    s = load(fullfile(files(i).folder, files(i).name));
    assert(abs(s.lengthXAxis - 0.060) < eps && abs(s.lengthYAxis - 0.040) < eps, ...
        'RTM:GeometryValidation:Domain', 'Geometry must be 600 x 400 um.');
    c = s.circleCenters; r = s.circleRadii(:);
    dx = c(:,1) - c(:,1).'; dy = c(:,2) - c(:,2).';
    d = hypot(dx, dy); d(1:size(d,1)+1:end) = inf;
    minGap = min(d - r - r.', [], 'all');
    [connected, poreFraction] = hasThroughPore(c, r, s.lengthXAxis, s.lengthYAxis);
    Kxx = permeability(c, r, s.lengthXAxis, s.lengthYAxis, nx, ny);
    cellsPerMinRadius = min(r) / max(s.lengthXAxis/nx, s.lengthYAxis/ny);
    rows(i,:) = {string(s.sampleName), s.finalPorosity, s.Iv, numel(r), minGap, ...
        connected, poreFraction, cellsPerMinRadius, Kxx};
end
report = cell2table(rows, 'VariableNames', {'sampleName','porosity','Iv','particleCount', ...
    'minimumGap_cm','leftRightConnected','poreFraction','cellsPerMinRadius','Kxx_mD'});
writetable(report, fullfile(geometryDir, 'geometry_validation.csv'));
writetable(report, fullfile(geometryDir, 'geometry_validation.xlsx'));
assert(all(report.minimumGap_cm >= -1e-12) && all(report.leftRightConnected) && ...
    all(report.cellsPerMinRadius >= 4), 'RTM:GeometryValidation:Failed', ...
    'Geometry overlap, connectivity, or mesh-resolution check failed.');
end

function [connected, poreFraction] = hasThroughPore(centers, radii, Lx, Ly)
nx = 600; ny = 400;
x = ((0:nx-1)+0.5)*Lx/nx; y = ((0:ny-1)+0.5)*Ly/ny;
[X,Y] = meshgrid(x,y); solid = false(ny,nx);
for i = 1:numel(radii), solid = solid | hypot(X-centers(i,1),Y-centers(i,2)) <= radii(i); end
pore = ~solid; cc = bwconncomp(pore, 4);
connected = any(cellfun(@(idx) any(floor((idx-1)/ny)+1 == 1) && ...
    any(floor((idx-1)/ny)+1 == nx), cc.PixelIdxList));
poreFraction = mean(pore,'all');
end

function Kxx = permeability(centers, radii, Lx, Ly, nx, ny)
base = FoldedCartesianGrid(2, [0 Lx 0 Ly], [nx ny]); coord = base.coordinates;
level = -inf(size(coord,1),1);
for i = 1:numel(radii), level = max(level, radii(i)-hypot(coord(:,1)-centers(i,1),coord(:,2)-centers(i,2))); end
K = computePhysicalPermeabilityTensor(coord, base.triangles, level);
Kxx = K(1,1)/9.869233e-12;
end
