function test_PermeabilityCellAreaNormalization
% The upstream-reference dimensionalization must give K proportional to L^2.

rtmDir = fileparts(fileparts(mfilename('fullpath')));
reactiveRoot = fileparts(rtmDir);
addpath(genpath(fullfile(reactiveRoot, 'src')));
addpath(fullfile(reactiveRoot, 'HyPHM'));
addpath(fullfile(reactiveRoot, 'HyPHM', 'tools'));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'classes')));

K1 = solveCell(1);
K2 = solveCell(2);
assert(norm(K2 - 4 * K1, 'fro') <= 1e-8 * norm(K2, 'fro'), ...
    'Physical permeability did not follow the required L^2 scaling.');
end

function K = solveCell(scale)
lx = 0.06 * scale;
ly = 0.04 * scale;
mesh = FoldedCartesianGrid(2, [0, lx, 0, ly], [48, 32]);
source = Grid(mesh.coordinates, mesh.triangles);
radius = 0.012 * scale;
center = [0.5 * lx, 0.5 * ly];
level = radius - hypot(source.coordV(:, 1) - center(1), ...
    source.coordV(:, 2) - center(2));
K = computePhysicalPermeabilityTensor(source.coordV, source.V0T, level);
end
