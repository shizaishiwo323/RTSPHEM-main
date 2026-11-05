function result = test_UpstreamPermeabilityForSavedGeometry(geometryPath, outputCsv)
% Independent upstream-RTSPHEM permeability check for a saved PNM geometry.
% Input: random_geometry_config.mat created by PNM_beauty3.
% Output: one CSV only; this function does not call or modify PNM_beauty3.

arguments
    geometryPath (1,1) string {mustBeFile}
    outputCsv (1,1) string = ""
end

rtmDir = fileparts(fileparts(mfilename('fullpath')));
reactiveRoot = fileparts(rtmDir);
addpath(genpath(fullfile(reactiveRoot, 'src')));
addpath(fullfile(reactiveRoot, 'HyPHM'));
addpath(fullfile(reactiveRoot, 'HyPHM', 'tools'));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'classes')));

g = load(geometryPath, 'circleCenters', 'circleRadius', 'lengthXAxis', 'lengthYAxis');
if isfield(g, 'circleRadii')
    radii = g.circleRadii(:);
else
    radii = repmat(g.circleRadius, size(g.circleCenters, 1), 1);
end
assert(numel(radii) == size(g.circleCenters, 1), 'Circle/radius count mismatch.');

nyParts = 128; % Matches the legacy PNM_beauty3 default mesh density.
nxParts = nyParts * max(1, round(g.lengthXAxis / g.lengthYAxis));
baseGrid = FoldedCartesianGrid(2, [0, g.lengthXAxis, 0, g.lengthYAxis], [nxParts, nyParts]);
sourceGrid = Grid(baseGrid.coordinates, baseGrid.triangles);
sourceLevelSet = -inf(sourceGrid.numV, 1); % fluid < 0, solid > 0
for i = 1:size(g.circleCenters, 1)
    r = hypot(sourceGrid.coordV(:, 1) - g.circleCenters(i, 1), ...
        sourceGrid.coordV(:, 2) - g.circleCenters(i, 2));
    sourceLevelSet = max(sourceLevelSet, radii(i) - r);
end

K_cm2 = computePhysicalPermeabilityTensor( ...
    sourceGrid.coordV, sourceGrid.V0T, sourceLevelSet);
mdToCm2 = 9.869233e-12;
K_mD = K_cm2 / mdToCm2;
assert(all(diag(K_mD) > 0), 'Permeability tensor must have positive diagonal entries.');
assert(norm(K_mD - K_mD.', 'fro') < 1e-6 * max(1, norm(K_mD, 'fro')), ...
    'Permeability tensor lost symmetry.');
result = table(g.lengthXAxis, g.lengthYAxis, nxParts, nyParts, size(g.circleCenters, 1), ...
    K_cm2(1,1), K_cm2(1,2), K_cm2(2,1), K_cm2(2,2), ...
    K_mD(1,1), K_mD(1,2), K_mD(2,1), K_mD(2,2), ...
    'VariableNames', {'length_x_cm','length_y_cm','nx_parts','ny_parts','circle_count', ...
    'Kxx_cm2','Kxy_cm2','Kyx_cm2','Kyy_cm2','Kxx_mD','Kxy_mD','Kyx_mD','Kyy_mD'});

if strlength(outputCsv) == 0
    outputCsv = fullfile(fileparts(geometryPath), 'upstream_permeability_tensor.csv');
end
writetable(result, outputCsv);
fprintf('Upstream cell-problem Kxx = %.8g mD; saved to %s\n', K_mD(1,1), outputCsv);
end
