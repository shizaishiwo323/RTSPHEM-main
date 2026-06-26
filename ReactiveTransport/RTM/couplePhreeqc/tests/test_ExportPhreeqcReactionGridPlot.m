function tests = test_ExportPhreeqcReactionGridPlot
tests = functiontests(localfunctions);
end

function testWritesReactionGridPlotWithSolidBoundary(testCase)
outputDir = tempname;
mkdir(outputDir);
testCase.addTeardown(@() cleanupOutput(outputDir));

grid = struct();
grid.coordV = [0 0; 0.01 0; 0.01 0.01; 0 0.01];
grid.V0T = [1 2 3; 1 3 4];
grid.numT = 2;

xlin = linspace(0, 0.01, 21);
ylin = linspace(0, 0.01, 21);
[X, Y] = meshgrid(xlin, ylin);
phi = 0.0025 - sqrt((X - 0.005).^2 + (Y - 0.005).^2);

outputPath = fullfile(outputDir, 'phreeqc_reaction_grid_initial.png');
ExportPhreeqcReactionGridPlot(outputPath, grid, X, Y, phi, 0.002, 0.01, 0.01);

verifyTrue(testCase, isfile(outputPath));
fileInfo = dir(outputPath);
verifyGreaterThan(testCase, fileInfo.bytes, 0);
end

function testWritesReactionGridPlotToMultipleOutputPaths(testCase)
outputDir = tempname;
mkdir(outputDir);
testCase.addTeardown(@() cleanupOutput(outputDir));

grid = struct();
grid.coordV = [0 0; 0.01 0; 0.01 0.01; 0 0.01];
grid.V0T = [1 2 3; 1 3 4];
grid.numT = 2;

xlin = linspace(0, 0.01, 21);
ylin = linspace(0, 0.01, 21);
[X, Y] = meshgrid(xlin, ylin);
phi = 0.0025 - sqrt((X - 0.005).^2 + (Y - 0.005).^2);

mainPath = fullfile(outputDir, 'phreeqc_reaction_grid_initial.png');
diagnosticPath = fullfile(outputDir, 'mesh_diagnostics', ...
    'phreeqc_reaction_grid_initial.png');
ExportPhreeqcReactionGridPlot({mainPath, diagnosticPath}, ...
    grid, X, Y, phi, 0.002, 0.01, 0.01);

verifyTrue(testCase, isfile(mainPath));
verifyTrue(testCase, isfile(diagnosticPath));
mainInfo = dir(mainPath);
diagnosticInfo = dir(diagnosticPath);
verifyGreaterThan(testCase, mainInfo.bytes, 0);
verifyGreaterThan(testCase, diagnosticInfo.bytes, 0);
end

function cleanupOutput(outputDir)
pngPath = fullfile(outputDir, 'phreeqc_reaction_grid_initial.png');
if isfile(pngPath)
    delete(pngPath);
end
diagnosticPath = fullfile(outputDir, 'mesh_diagnostics', ...
    'phreeqc_reaction_grid_initial.png');
if isfile(diagnosticPath)
    delete(diagnosticPath);
end
diagnosticDir = fullfile(outputDir, 'mesh_diagnostics');
if exist(diagnosticDir, 'dir') == 7
    rmdir(diagnosticDir);
end
if exist(outputDir, 'dir') == 7
    rmdir(outputDir);
end
end
