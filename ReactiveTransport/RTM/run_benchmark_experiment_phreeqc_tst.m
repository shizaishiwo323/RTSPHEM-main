function summary = run_benchmark_experiment_phreeqc_tst(mode)
% run_benchmark_experiment_phreeqc_tst
%
% Reproduce the Part III experimental-validation benchmark from:
% "Simulation of mineral dissolution at the pore scale with evolving
% fluid-solid interfaces: review of approaches and benchmark problem set".
%
% Outputs:
%   outputs/benchmark_experiment_phreeqc_tst/<timestamp>/
%     runs/phreeqc/
%     runs/tst/
%     comparison/
%       comparison_summary.csv
%       comparison_by_nearest_time.csv
%       global_comparison.png
%       benchmark_manifest.json

if nargin < 1 || isempty(mode)
    mode = "full";
end
mode = lower(string(mode));
if ~ismember(mode, ["full", "smoke"])
    error('RTSPHEM:BenchmarkRunner:InvalidMode', 'Mode must be full or smoke.');
end

clc;

scriptDir = fileparts(mfilename('fullpath'));
rtmDir = scriptDir;
projectRoot = fileparts(fileparts(rtmDir));
coupleDir = fullfile(rtmDir, 'couplePhreeqc');
addpath(rtmDir);
addpath(coupleDir);

stamp = datestr(now, 'yyyymmdd_HHMMSS');
packageDir = fullfile(projectRoot, 'outputs', 'benchmark_experiment_phreeqc_tst', char(mode), stamp);
runsDir = fullfile(packageDir, 'runs');
comparisonDir = fullfile(packageDir, 'comparison');
ensureDir(runsDir);
ensureDir(comparisonDir);

fprintf('Benchmark package: %s\n', packageDir);

baseCfg = buildBenchmarkConfig(projectRoot, rtmDir);
if mode == "smoke"
    baseCfg.endTime = 1;
    baseCfg.initialMacroscaleTimeStepSize = 1;
    baseCfg.maximalStep = 1;
    baseCfg.meshTargetElementSizeCm = 0.012;
    baseCfg.dxfResolutionX = 60;
    baseCfg.dxfResolutionY = 34;
    baseCfg.exportEvery = 1;
    baseCfg.phreeqcExportEvery = 1;
    baseCfg.saveMainPlot = false;
    baseCfg.saveIndividualPlots = false;
    baseCfg.saveInterfaceMask = false;
    baseCfg.exportDXF = false;
    baseCfg.saveFinalPlot = false;
    baseCfg.saveMeshDiagnostics = false;
end

phreeqcCfg = baseCfg;
phreeqcCfg.reactionModel = 'phreeqc';
phreeqcCfg.runName = 'benchmark_partIII_phreeqc';
phreeqcCfg.resultsDir = fullfile(runsDir, 'phreeqc');
fprintf('\n=== Running PHREEQC benchmark ===\n');
phreeqcResult = PNM_beauty3(phreeqcCfg);

tstCfg = baseCfg;
tstCfg.reactionModel = 'tst';
tstCfg.runName = 'benchmark_partIII_tst';
tstCfg.resultsDir = fullfile(runsDir, 'tst');
fprintf('\n=== Running TST benchmark ===\n');
tstResult = PNM_beauty3(tstCfg);

summary = exportBenchmarkComparison(tstResult.resultsDir, phreeqcResult.resultsDir, comparisonDir);
writeBenchmarkManifest(packageDir, comparisonDir, baseCfg, tstResult, phreeqcResult);

fprintf('\nBenchmark export complete:\n');
fprintf('  Package: %s\n', packageDir);
fprintf('  Comparison: %s\n', comparisonDir);
end

function cfg = buildBenchmarkConfig(projectRoot, rtmDir)
cfg = struct();
cfg.outputRoot = fullfile(projectRoot, 'outputs', 'benchmark_experiment_phreeqc_tst');

cfg.layoutType = 'external_dxf';
cfg.useExternalGeometry = false;
cfg.useExternalDxfGeometry = true;
cfg.externalGeometryType = 'dxf';
cfg.externalDxfPath = "C:\Users\imgw\Documents\Codex\SIP模拟\sip模拟\data\dissolution_results-Da_40.4424_Pe_4.1640_L_0.1200_square\validation-domin.dxf";
cfg.externalDxfDomainLayerNames = {'domin', 'DOMAIN'};
cfg.externalDxfSolidLayerNames = {'calcite'};
cfg.externalDxfReferenceLength = 1500;
cfg.externalDxfReferenceLengthCm = 0.150;
cfg.externalDxfImportDirection = 'rotate90_cw';
cfg.externalDxfSmoothingSigmaPixels = 1.5;

% PDF Part III / Table 4 parameters.
cfg.inletVelocity = 0.117;                 % cm/s
cfg.flowDirection = 'left_to_right';
cfg.initialHydrogenConcentration = 1.26e-5; % mol/cm^3
cfg.diffusionCoefficient = 5e-5;           % cm^2/s
cfg.molarVolume = 36.9;                    % cm^3/mol
cfg.rateCoefficientTST = 1e-3;             % solver convention, recorded in manifest
cfg.characteristicLength = 0.150;          % channel width, gives Pe ~= 350

cfg.initialCalciumConcentration = 0;
cfg.initialCarbonConcentration = 0;
cfg.initialSodiumConcentration = 0;
cfg.initialChlorideConcentration = 0;
cfg.inletCalciumConcentration = 0;
cfg.inletCarbonConcentration = 0;
cfg.inletSodiumConcentration = 0;
cfg.inletChlorideConcentration = cfg.initialHydrogenConcentration;

cfg.phreeqcDatabasePath = string(fullfile(rtmDir, 'couplePhreeqc', 'phreeqc-m.dat'));
if exist(char(cfg.phreeqcDatabasePath), 'file') ~= 2
    referencePhreeqcDatabase = "C:\Users\imgw\Downloads\RTSPHEM-P-main (1)\RTSPHEM-P-main\SourceCode\phreeqc-m.dat";
    if exist(char(referencePhreeqcDatabase), 'file') == 2
        cfg.phreeqcDatabasePath = referencePhreeqcDatabase;
    else
        cfg.phreeqcDatabasePath = "C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat";
    end
end
cfg.phreeqcTemperatureC = 25;
cfg.phreeqcKineticsCorrectionFactor = 1;
cfg.phreeqcMaxSpecificSurfaceArea = 10;
cfg.phreeqcBadStepMax = 5000;
cfg.phreeqcKineticsTolerance = 1e-8;
cfg.phreeqcMinHForPHMolL = 1e-7;
cfg.phreeqcMinActiveWaterVolumeFraction = 0;
cfg.phreeqcMinActiveWaterVolumeCm3 = 0;
cfg.phreeqcReactNeutralInterfaceCells = false;
cfg.phreeqcSolutionWaterKg = 1;
cfg.phreeqcWriteSolutionWaterLine = false;
cfg.phreeqcKineticsReservoirMoles = 1;
cfg.phreeqcExportEvery = 20;

cfg.initialMacroscaleTimeStepSize = 1.0;
cfg.maximalStep = 120;
cfg.endTime = 12000;
cfg.timeStepperType = 'expmax';
cfg.targetDissolutionSlices = [];
cfg.enableConcentrationCflLimit = false;
cfg.permeabilityRatioThreshold = 1e7;
cfg.maxExperimentWallSeconds = Inf;

% Keep the published physical parameters, but use a practical 2D grid for
% the local PHREEQC-vs-TST comparison so the 12000 s benchmark completes.
cfg.meshTargetElementSizeCm = 0.010;
cfg.meshNumPartitionsX = [];
cfg.meshNumPartitionsY = [];
cfg.numPartitionsMicroscale = 32;
cfg.dxfResolutionX = 90;
cfg.dxfResolutionY = 50;

cfg.exportEvery = 20;
cfg.exportDXF = false;
cfg.saveMainPlot = true;
cfg.saveIndividualPlots = false;
cfg.saveInterfaceMask = false;
cfg.saveRealtimePlot = false;
cfg.saveFigureFiles = false;
cfg.writeExcel = true;
cfg.saveFinalPlot = true;
cfg.saveMeshDiagnostics = false;
cfg.showDebugFigures = false;

cfg.enableNMRSimulation = false;
cfg.enableNMRSurrogate = false;
cfg.enablePNGSimulation = false;
cfg.nmr_method = 'none';
end

function summary = exportBenchmarkComparison(tstRunDir, phreeqcRunDir, outputDir)
tstGlobal = readGlobalEvolution(tstRunDir);
phreeqcGlobal = readGlobalEvolution(phreeqcRunDir);
phreeqcSummary = readOptionalTable(fullfile(phreeqcRunDir, 'phreeqc_results', 'phreeqc_summary_log.csv'));

summary = buildSummaryTable(tstGlobal, phreeqcGlobal);
writetable(summary, fullfile(outputDir, 'comparison_summary.csv'));
writetable(summary, fullfile(outputDir, 'comparison_summary.xlsx'));

aligned = alignByNearestTime(tstGlobal, phreeqcGlobal);
writetable(aligned, fullfile(outputDir, 'comparison_by_nearest_time.csv'));
writetable(aligned, fullfile(outputDir, 'comparison_by_nearest_time.xlsx'));

plotGlobalComparison(tstGlobal, phreeqcGlobal, phreeqcSummary, ...
    fullfile(outputDir, 'global_comparison.png'));
end

function tbl = readGlobalEvolution(runDir)
path = fullfile(runDir, 'global_evolution_log.csv');
if exist(path, 'file') ~= 2
    error('RTSPHEM:BenchmarkExport:MissingGlobalLog', 'Missing global evolution log: %s', path);
end
tbl = readtable(path, 'VariableNamingRule', 'preserve');
end

function tbl = readOptionalTable(path)
if exist(path, 'file') == 2
    tbl = readtable(path, 'VariableNamingRule', 'preserve');
else
    tbl = table();
end
end

function summary = buildSummaryTable(tst, phreeqc)
thresholds = [0.55, 0.60, 0.80, 0.95, 0.99, 1.00];
metric = strings(0, 1);
tstValue = [];
phreeqcValue = [];
ratioPhreeqcToTst = [];

[metric, tstValue, phreeqcValue, ratioPhreeqcToTst] = addMetric(metric, tstValue, phreeqcValue, ratioPhreeqcToTst, "steps", height(tst), height(phreeqc));
[metric, tstValue, phreeqcValue, ratioPhreeqcToTst] = addMetric(metric, tstValue, phreeqcValue, ratioPhreeqcToTst, "final_time_s", tst.time_s(end), phreeqc.time_s(end));
[metric, tstValue, phreeqcValue, ratioPhreeqcToTst] = addMetric(metric, tstValue, phreeqcValue, ratioPhreeqcToTst, "final_porosity", tst.porosity(end), phreeqc.porosity(end));
[metric, tstValue, phreeqcValue, ratioPhreeqcToTst] = addMetric(metric, tstValue, phreeqcValue, ratioPhreeqcToTst, "final_permeability_mD", tst.permeability_mD(end), phreeqc.permeability_mD(end));
[metric, tstValue, phreeqcValue, ratioPhreeqcToTst] = addMetric(metric, tstValue, phreeqcValue, ratioPhreeqcToTst, "final_k_k0", tst.k_k0(end), phreeqc.k_k0(end));

for iThreshold = 1:numel(thresholds)
    name = "time_to_porosity_" + replace(sprintf('%.2f', thresholds(iThreshold)), '.', 'p') + "_s";
    [metric, tstValue, phreeqcValue, ratioPhreeqcToTst] = addMetric( ...
        metric, tstValue, phreeqcValue, ratioPhreeqcToTst, ...
        name, timeToPorosity(tst, thresholds(iThreshold)), timeToPorosity(phreeqc, thresholds(iThreshold)));
end

summary = table(metric, tstValue, phreeqcValue, ratioPhreeqcToTst);
end

function [metric, tstValue, phreeqcValue, ratio] = addMetric(metric, tstValue, phreeqcValue, ratio, name, tstVal, phreeqcVal)
metric(end + 1, 1) = name;
tstValue(end + 1, 1) = tstVal;
phreeqcValue(end + 1, 1) = phreeqcVal;
if isfinite(tstVal) && abs(tstVal) > eps
    ratio(end + 1, 1) = phreeqcVal / tstVal;
else
    ratio(end + 1, 1) = NaN;
end
end

function t = timeToPorosity(tbl, threshold)
idx = find(tbl.porosity >= threshold, 1, 'first');
if isempty(idx)
    t = NaN;
elseif idx == 1
    t = tbl.time_s(1);
else
    p0 = tbl.porosity(idx - 1);
    p1 = tbl.porosity(idx);
    t0 = tbl.time_s(idx - 1);
    t1 = tbl.time_s(idx);
    if abs(p1 - p0) < eps
        t = t1;
    else
        t = t0 + (threshold - p0) * (t1 - t0) / (p1 - p0);
    end
end
end

function aligned = alignByNearestTime(tst, phreeqc)
phreeqcTime = phreeqc.time_s;
n = numel(phreeqcTime);
tstIndex = zeros(n, 1);
for i = 1:n
    [~, tstIndex(i)] = min(abs(tst.time_s - phreeqcTime(i)));
end

aligned = table();
aligned.time_s = phreeqcTime;
aligned.tst_time_s = tst.time_s(tstIndex);
aligned.tst_porosity = tst.porosity(tstIndex);
aligned.phreeqc_porosity = phreeqc.porosity;
aligned.delta_porosity = phreeqc.porosity - tst.porosity(tstIndex);
aligned.tst_k_k0 = tst.k_k0(tstIndex);
aligned.phreeqc_k_k0 = phreeqc.k_k0;
aligned.delta_k_k0 = phreeqc.k_k0 - tst.k_k0(tstIndex);
aligned.tst_rate = tst.avg_dissolution_rate(tstIndex);
aligned.phreeqc_rate = phreeqc.avg_dissolution_rate;
end

function plotGlobalComparison(tst, phreeqc, phreeqcSummary, outputPath)
fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 900]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(tst.time_s, tst.porosity, 'k-', 'LineWidth', 1.5);
hold on;
plot(phreeqc.time_s, phreeqc.porosity, 'r.-', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Porosity');
legend({'TST', 'PHREEQC'}, 'Location', 'best');
grid on;

nexttile;
semilogy(tst.time_s, max(tst.k_k0, eps), 'k-', 'LineWidth', 1.5);
hold on;
semilogy(phreeqc.time_s, max(phreeqc.k_k0, eps), 'r.-', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('k/k0');
legend({'TST', 'PHREEQC'}, 'Location', 'best');
grid on;

nexttile;
semilogy(tst.time_s, max(tst.avg_dissolution_rate, eps), 'k-', 'LineWidth', 1.5);
hold on;
semilogy(phreeqc.time_s, max(phreeqc.avg_dissolution_rate, eps), 'r.-', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Average dissolution rate');
legend({'TST', 'PHREEQC'}, 'Location', 'best');
grid on;

nexttile;
if ~isempty(phreeqcSummary)
    yyaxis left;
    plot(phreeqcSummary.time_s, phreeqcSummary.pH_mean, 'b.-', 'LineWidth', 1.2);
    ylabel('PHREEQC mean pH');
    yyaxis right;
    semilogy(phreeqcSummary.time_s, max(phreeqcSummary.Ca_mean_mol_cm3, eps), 'm.-', 'LineWidth', 1.2);
    ylabel('Ca mean [mol/cm^3]');
    xlabel('Time [s]');
    grid on;
else
    text(0.5, 0.5, 'No PHREEQC summary log', 'HorizontalAlignment', 'center');
    axis off;
end

exportgraphics(fig, outputPath, 'Resolution', 200);
close(fig);
end

function writeBenchmarkManifest(packageDir, comparisonDir, baseCfg, tstResult, phreeqcResult)
manifest = struct();
manifest.schema_version = "benchmark_experiment_phreeqc_tst_v1";
manifest.created_at = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
manifest.source_pdf = "C:\Users\imgw\Downloads\Simulation of mineral dissolution at the pore scale with evolving fluid-solid interfaces_ review of approaches and benchmark problem set.pdf";
manifest.benchmark = "Part III experimental validation";
manifest.paper_parameters = struct( ...
    'fluid_density_g_cm3', 0.92, ...
    'kinematic_viscosity_cm2_s', 0.0261, ...
    'diffusionCoefficient_cm2_s', 5e-5, ...
    'inletVelocity_cm_s', 0.117, ...
    'length_cm', 0.267, ...
    'width_cm', 0.150, ...
    'height_cm', 0.02, ...
    'specificReactiveArea_cm_inv', 78.5, ...
    'rateConstant_mol_cm2_s', 1e-3, ...
    'activityCoefficient_cm3_mol', 1000, ...
    'molarVolume_cm3_mol', 36.9, ...
    'inletConcentration_mol_cm3', 1.26e-5, ...
    'reynolds', 0.671, ...
    'peclet', 350, ...
    'damkohlerII', 3930, ...
    'endTime_s', 12000);
manifest.runner_config = baseCfg;
manifest.runs = struct( ...
    'tst', string(tstResult.resultsDir), ...
    'phreeqc', string(phreeqcResult.resultsDir));
manifest.comparison_dir = string(comparisonDir);
manifest.outputs = struct( ...
    'summary_csv', string(fullfile(comparisonDir, 'comparison_summary.csv')), ...
    'summary_xlsx', string(fullfile(comparisonDir, 'comparison_summary.xlsx')), ...
    'nearest_time_csv', string(fullfile(comparisonDir, 'comparison_by_nearest_time.csv')), ...
    'nearest_time_xlsx', string(fullfile(comparisonDir, 'comparison_by_nearest_time.xlsx')), ...
    'plot_png', string(fullfile(comparisonDir, 'global_comparison.png')));
manifest.notes = [ ...
    "Both solvers use the same DXF geometry, flow, diffusion, inlet concentration, time span, mesh target size, and export settings."; ...
    "The PNM TST coefficient follows the repository solver convention; the PDF Part III rate constant is recorded separately above."; ...
    "RTSPHEM PNM_beauty3 is a 2D solver, so this export uses the Part III mask/DXF and Table 4 parameters as a 2D local benchmark rather than a full 3D reproduction of the paper."; ...
    "NMR/COMSOL synchronization is disabled for this chemistry benchmark export."];

writeJsonFileLocal(fullfile(comparisonDir, 'benchmark_manifest.json'), manifest);
end

function ensureDir(path)
if exist(path, 'dir') ~= 7
    mkdir(path);
end
end

function writeJsonFileLocal(path, data)
jsonText = jsonencode(data, 'PrettyPrint', true);
fid = fopen(path, 'w');
if fid == -1
    error('RTSPHEM:BenchmarkExport:ManifestOpenFailed', 'Cannot write manifest: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonText);
clear cleanupObj;
end
