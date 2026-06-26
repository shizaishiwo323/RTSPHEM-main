function summary = run_benchmark_experiment_phreeqc_tst(mode)
% run_benchmark_experiment_phreeqc_tst
%
% Reproduce Molins et al. Part III experimental calcite-column benchmark
% using two RTSPHEM chemistry paths:
%   1) cfg.phreeqcRunGroup = 'phreeqc_tst_match'
%   2) legacy RTSPHEM TST
%
% Outputs:
%   outputs/benchmark_molins_partIII_phreeqc_tst/<mode>/<timestamp>/
%     runs/phreeqc_tst_match/
%     runs/tst/
%     comparison/
%       comparison_summary.csv/xlsx
%       comparison_by_nearest_time.csv/xlsx
%       molins_partIII_volume_area.png
%       benchmark_manifest.json

if nargin < 1 || isempty(mode)
    mode = "full";
end
mode = lower(string(mode));
if ~ismember(mode, ["full", "partial", "smoke"])
    error('RTSPHEM:BenchmarkRunner:InvalidMode', 'Mode must be full, partial, or smoke.');
end

clc;

scriptDir = fileparts(mfilename('fullpath'));
rtmDir = scriptDir;
projectRoot = fileparts(fileparts(rtmDir));
coupleDir = fullfile(rtmDir, 'couplePhreeqc');
addpath(rtmDir);
addpath(coupleDir);

stamp = datestr(now, 'yyyymmdd_HHMMSS');
packageDir = fullfile(projectRoot, 'outputs', 'benchmark_molins_partIII_phreeqc_tst', char(mode), stamp);
runsDir = fullfile(packageDir, 'runs');
comparisonDir = fullfile(packageDir, 'comparison');
ensureDir(runsDir);
ensureDir(comparisonDir);

fprintf('Molins Part III benchmark package: %s\n', packageDir);

baseCfg = buildPartIIIBenchmarkConfig(projectRoot, rtmDir, mode);

phreeqcCfg = ConfigurePhreeqcRunGroup(baseCfg, 'phreeqc_tst_match');
phreeqcCfg.runName = 'benchmark_molins_partIII_phreeqc_tstmatch';
phreeqcCfg.resultsDir = fullfile(runsDir, 'phreeqc_tst_match');
fprintf('\n=== Running PHREEQC TST-match Part III benchmark ===\n');
phreeqcResult = PNM_beauty3(phreeqcCfg);

tstCfg = baseCfg;
tstCfg.reactionModel = 'tst';
tstCfg.runName = 'benchmark_molins_partIII_tst';
tstCfg.resultsDir = fullfile(runsDir, 'tst');
fprintf('\n=== Running legacy TST Part III benchmark ===\n');
tstResult = PNM_beauty3(tstCfg);

summary = exportPartIIIComparison(tstResult.resultsDir, phreeqcResult.resultsDir, comparisonDir, baseCfg);
writePartIIIBenchmarkManifest(packageDir, comparisonDir, mode, baseCfg, tstResult, phreeqcResult);

fprintf('\nPart III benchmark export complete:\n');
fprintf('  Package: %s\n', packageDir);
fprintf('  Comparison: %s\n', comparisonDir);
end

function cfg = buildPartIIIBenchmarkConfig(projectRoot, rtmDir, mode)
cfg = struct();
cfg.outputRoot = fullfile(projectRoot, 'outputs', 'benchmark_molins_partIII_phreeqc_tst');

% Molins et al. Part III / Table 4. Units are converted to RTSPHEM cm/s/mol.
cfg.layoutType = 'external_tif';
cfg.useExternalGeometry = true;
cfg.externalGeometryType = 'tif';
cfg.tifPath = "C:\Users\imgw\Downloads\10596_2019_9903_MOESM1_ESM\SupplementaryMaterial\InputsCodes\_Part III-MaskImage\mask.tif";
cfg.externalTifSolidValue = 0;
cfg.externalTifLiquidValue = 255;
cfg.externalTifPixelSizeMicron = 0.2667 / 2524 * 1e4;
cfg.externalTifSmoothingSigmaPixels = 1.5;
cfg.externalTifCropTopRows = 72;
cfg.externalTifCropBottomRows = 38;
cfg.externalTifCropLeftCols = 0;
cfg.externalTifCropRightCols = 0;

cfg.partIIIHeightCm = 0.02;
cfg.thicknessCm = cfg.partIIIHeightCm;
cfg.partIIILengthCm = 0.267;
cfg.partIIIWidthCm = 0.150;
cfg.partIIISpecificReactiveAreaCmInv = 78.5;
cfg.characteristicLength = cfg.partIIIWidthCm;
cfg.targetLengthYAxis = cfg.partIIIWidthCm;
cfg.targetAspectRatio = cfg.partIIILengthCm / cfg.partIIIWidthCm;

cfg.inletVelocity = 0.117;
cfg.flowDirection = 'left_to_right';
cfg.initialHydrogenConcentration = 1.26e-5;
cfg.diffusionCoefficient = 5e-5;
cfg.molarVolume = 36.9;
cfg.rateCoefficientTST = 1e-3;

cfg.initialCalciumConcentration = 0;
cfg.initialCarbonConcentration = 0;
cfg.initialSodiumConcentration = 0;
cfg.initialChlorideConcentration = cfg.initialHydrogenConcentration;
cfg.inletCalciumConcentration = 0;
cfg.inletCarbonConcentration = 0;
cfg.inletSodiumConcentration = 0;
cfg.inletChlorideConcentration = cfg.initialHydrogenConcentration;

cfg.phreeqcDatabasePath = ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc_rates.dat');
cfg.phreeqcTemperatureC = 25;
cfg.phreeqcKineticsCorrectionFactor = 1;
cfg.phreeqcBadStepMax = 5000;
cfg.phreeqcKineticsTolerance = 1e-8;
cfg.phreeqcMinHForPHMolL = 1e-7;
cfg.phreeqcMinActiveWaterVolumeFraction = 0;
cfg.phreeqcMinActiveWaterVolumeCm3 = 0;
cfg.phreeqcReactNeutralInterfaceCells = false;
cfg.phreeqcSolutionWaterKg = 1;
cfg.phreeqcWriteSolutionWaterLine = false;
cfg.phreeqcKineticsReservoirMoles = 1;
cfg.phreeqcRunGroup = 'phreeqc_tst_match';
cfg.phreeqcRateLaw = 'tst_match';
cfg.phreeqcTstRateCoefficient = cfg.rateCoefficientTST;
cfg.phreeqcMaxSpecificSurfaceArea = Inf;
cfg.phreeqcExportEvery = 1;

cfg.endTime = 12000;
cfg.processSliceCount = 20;
cfg.initialMacroscaleTimeStepSize = 600;
cfg.maximalStep = 600;
cfg.timeStepperType = 'linear';
cfg.maxTotalTimeSteps = [];
cfg.targetDissolutionSlices = [];
cfg.enableConcentrationCflLimit = false;
cfg.permeabilityRatioThreshold = 1e7;
cfg.maxExperimentWallSeconds = Inf;

cfg.numPartitionsMicroscale = 2 * 64;
cfg.meshTargetElementSizeCm = [];
cfg.meshNumPartitionsX = 256;
cfg.meshNumPartitionsY = 144;
cfg.dxfResolutionX = 300;
cfg.dxfResolutionY = 168;

if mode == "smoke"
    cfg.endTime = 600;
    cfg.processSliceCount = 1;
    cfg.initialMacroscaleTimeStepSize = 600;
    cfg.maximalStep = 600;
    cfg.meshNumPartitionsX = 64;
    cfg.meshNumPartitionsY = 36;
    cfg.dxfResolutionX = 120;
    cfg.dxfResolutionY = 68;
    cfg.exportDXF = false;
    cfg.saveMainPlot = false;
    cfg.saveIndividualPlots = false;
    cfg.saveInterfaceMask = false;
    cfg.saveFinalPlot = false;
elseif mode == "partial"
    cfg.endTime = 2400;
    cfg.processSliceCount = 4;
    cfg.meshNumPartitionsX = 128;
    cfg.meshNumPartitionsY = 72;
    cfg.dxfResolutionX = 180;
    cfg.dxfResolutionY = 102;
    cfg.exportDXF = false;
    cfg.saveMainPlot = false;
    cfg.saveIndividualPlots = false;
    cfg.saveInterfaceMask = false;
    cfg.saveFinalPlot = false;
else
    cfg.exportDXF = true;
    cfg.saveMainPlot = true;
    cfg.saveIndividualPlots = true;
    cfg.saveInterfaceMask = true;
    cfg.saveFinalPlot = true;
end

cfg.exportEvery = 1;
cfg.saveRealtimePlot = false;
cfg.saveFigureFiles = false;
cfg.writeExcel = true;
cfg.saveMeshDiagnostics = false;
cfg.showDebugFigures = false;

cfg.enableNMRSimulation = false;
cfg.enableNMRSurrogate = false;
cfg.enablePNGSimulation = false;
cfg.nmr_method = 'none';
end

function summary = exportPartIIIComparison(tstRunDir, phreeqcRunDir, outputDir, cfg)
tstGlobal = addPartIIIScaledColumns(readGlobalEvolution(tstRunDir), cfg.partIIIHeightCm);
phreeqcGlobal = addPartIIIScaledColumns(readGlobalEvolution(phreeqcRunDir), cfg.partIIIHeightCm);
phreeqcSummary = readOptionalTable(fullfile(phreeqcRunDir, 'phreeqc_results', 'phreeqc_summary_log.csv'));
molinsRefs = loadMolinsPartIIIReferences();

summary = buildPartIIISummaryTable(tstGlobal, phreeqcGlobal, molinsRefs);
writetable(summary, fullfile(outputDir, 'comparison_summary.csv'));
writetable(summary, fullfile(outputDir, 'comparison_summary.xlsx'));

aligned = alignByNearestTime(tstGlobal, phreeqcGlobal);
writetable(aligned, fullfile(outputDir, 'comparison_by_nearest_time.csv'));
writetable(aligned, fullfile(outputDir, 'comparison_by_nearest_time.xlsx'));

plotPartIIIComparison(tstGlobal, phreeqcGlobal, phreeqcSummary, molinsRefs, ...
    fullfile(outputDir, 'molins_partIII_volume_area.png'));
plotPartIIIComparison(tstGlobal, phreeqcGlobal, phreeqcSummary, molinsRefs, ...
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

function tbl = addPartIIIScaledColumns(tbl, heightCm)
tbl.mineral_volume_cm3 = tbl.grain_volume_cm3 * heightCm;
tbl.reactive_surface_area_cm2_3d = tbl.surface_area_cm2 * heightCm;
if height(tbl) > 0 && tbl.mineral_volume_cm3(1) > 0
    tbl.relative_mineral_volume = tbl.mineral_volume_cm3 / tbl.mineral_volume_cm3(1);
else
    tbl.relative_mineral_volume = NaN(height(tbl), 1);
end
if height(tbl) > 0 && tbl.reactive_surface_area_cm2_3d(1) > 0
    tbl.relative_surface_area = tbl.reactive_surface_area_cm2_3d / tbl.reactive_surface_area_cm2_3d(1);
else
    tbl.relative_surface_area = NaN(height(tbl), 1);
end
end

function summary = buildPartIIISummaryTable(tst, phreeqc, refs)
thresholds = [0.90, 0.75, 0.50, 0.25, 0.10];
metric = strings(0, 1);
tstValue = [];
phreeqcValue = [];
referenceValue = [];

[metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, ...
    "steps", height(tst), height(phreeqc), NaN);
[metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, ...
    "final_time_s", tst.time_s(end), phreeqc.time_s(end), NaN);
[metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, ...
    "initial_mineral_volume_cm3", tst.mineral_volume_cm3(1), phreeqc.mineral_volume_cm3(1), firstReferenceValue(refs.volume, "Chombo-Crunch"));
[metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, ...
    "final_mineral_volume_cm3", tst.mineral_volume_cm3(end), phreeqc.mineral_volume_cm3(end), lastReferenceValue(refs.volume, "Chombo-Crunch"));
[metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, ...
    "initial_surface_area_cm2", tst.reactive_surface_area_cm2_3d(1), phreeqc.reactive_surface_area_cm2_3d(1), firstReferenceValue(refs.surface, "Chombo-Crunch"));
[metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, ...
    "final_surface_area_cm2", tst.reactive_surface_area_cm2_3d(end), phreeqc.reactive_surface_area_cm2_3d(end), lastReferenceValue(refs.surface, "Chombo-Crunch"));

for iThreshold = 1:numel(thresholds)
    label = replace(sprintf('%.2f', thresholds(iThreshold)), '.', 'p');
    [metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, ...
        "time_to_relative_volume_" + label + "_s", ...
        timeToRelativeVolume(tst, thresholds(iThreshold)), ...
        timeToRelativeVolume(phreeqc, thresholds(iThreshold)), ...
        timeToRelativeReference(refs.volume, "Chombo-Crunch", thresholds(iThreshold)));
end

summary = table(metric, tstValue, phreeqcValue, referenceValue);
summary.deltaPhreeqcMinusTst = summary.phreeqcValue - summary.tstValue;
summary.deltaTstMinusReference = summary.tstValue - summary.referenceValue;
summary.deltaPhreeqcMinusReference = summary.phreeqcValue - summary.referenceValue;
end

function [metric, tstValue, phreeqcValue, referenceValue] = addSummaryMetric(metric, tstValue, phreeqcValue, referenceValue, name, tstVal, phreeqcVal, refVal)
metric(end + 1, 1) = name;
tstValue(end + 1, 1) = tstVal;
phreeqcValue(end + 1, 1) = phreeqcVal;
referenceValue(end + 1, 1) = refVal;
end

function t = timeToRelativeVolume(tbl, threshold)
idx = find(tbl.relative_mineral_volume <= threshold, 1, 'first');
if isempty(idx)
    t = NaN;
elseif idx == 1
    t = tbl.time_s(1);
else
    v0 = tbl.relative_mineral_volume(idx - 1);
    v1 = tbl.relative_mineral_volume(idx);
    t0 = tbl.time_s(idx - 1);
    t1 = tbl.time_s(idx);
    if abs(v1 - v0) < eps
        t = t1;
    else
        t = t0 + (threshold - v0) * (t1 - t0) / (v1 - v0);
    end
end
end

function t = timeToRelativeReference(entries, label, threshold)
entry = findReference(entries, label);
if isempty(entry) || isempty(entry.value) || entry.value(1) <= 0
    t = NaN;
    return;
end
rel = entry.value / entry.value(1);
idx = find(rel <= threshold, 1, 'first');
if isempty(idx)
    t = NaN;
elseif idx == 1
    t = entry.time_s(1);
else
    v0 = rel(idx - 1);
    v1 = rel(idx);
    t0 = entry.time_s(idx - 1);
    t1 = entry.time_s(idx);
    if abs(v1 - v0) < eps
        t = t1;
    else
        t = t0 + (threshold - v0) * (t1 - t0) / (v1 - v0);
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
aligned.tst_mineral_volume_cm3 = tst.mineral_volume_cm3(tstIndex);
aligned.phreeqc_mineral_volume_cm3 = phreeqc.mineral_volume_cm3;
aligned.delta_mineral_volume_cm3 = phreeqc.mineral_volume_cm3 - tst.mineral_volume_cm3(tstIndex);
aligned.tst_relative_mineral_volume = tst.relative_mineral_volume(tstIndex);
aligned.phreeqc_relative_mineral_volume = phreeqc.relative_mineral_volume;
aligned.tst_surface_area_cm2 = tst.reactive_surface_area_cm2_3d(tstIndex);
aligned.phreeqc_surface_area_cm2 = phreeqc.reactive_surface_area_cm2_3d;
aligned.delta_surface_area_cm2 = phreeqc.reactive_surface_area_cm2_3d - tst.reactive_surface_area_cm2_3d(tstIndex);
aligned.tst_rate = tst.avg_dissolution_rate(tstIndex);
aligned.phreeqc_rate = phreeqc.avg_dissolution_rate;
end

function plotPartIIIComparison(tst, phreeqc, phreeqcSummary, refs, outputPath)
fig = figure('Visible', 'off', 'Position', [100, 100, 1500, 1100]);
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(tst.time_s, tst.mineral_volume_cm3, 'k-o', 'LineWidth', 1.4, 'MarkerSize', 3, 'DisplayName', 'RTSPHEM TST');
hold on;
plot(phreeqc.time_s, phreeqc.mineral_volume_cm3, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12, 'DisplayName', 'RTSPHEM PHREEQC-TST');
plotReferenceGroup(refs.volume);
xlabel('Time [s]');
ylabel('Mineral volume [cm^3]');
title('Molins Part III mineral volume');
legend('Location', 'best');
grid on;

nexttile;
plot(tst.time_s, tst.relative_mineral_volume, 'k-o', 'LineWidth', 1.4, 'MarkerSize', 3, 'DisplayName', 'RTSPHEM TST');
hold on;
plot(phreeqc.time_s, phreeqc.relative_mineral_volume, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12, 'DisplayName', 'RTSPHEM PHREEQC-TST');
plotReferenceGroup(refs.relativeVolume);
xlabel('Time [s]');
ylabel('V/V_0');
title('Normalized mineral volume');
legend('Location', 'best');
grid on;

nexttile;
plot(tst.time_s, tst.reactive_surface_area_cm2_3d, 'k-o', 'LineWidth', 1.4, 'MarkerSize', 3, 'DisplayName', 'RTSPHEM TST');
hold on;
plot(phreeqc.time_s, phreeqc.reactive_surface_area_cm2_3d, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12, 'DisplayName', 'RTSPHEM PHREEQC-TST');
plotReferenceGroup(refs.surface);
xlabel('Time [s]');
ylabel('Reactive surface area [cm^2]');
title('Reactive surface area');
legend('Location', 'best');
grid on;

nexttile;
plot(tst.time_s, tst.relative_surface_area, 'k-o', 'LineWidth', 1.4, 'MarkerSize', 3, 'DisplayName', 'RTSPHEM TST');
hold on;
plot(phreeqc.time_s, phreeqc.relative_surface_area, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12, 'DisplayName', 'RTSPHEM PHREEQC-TST');
plotReferenceGroup(refs.relativeSurface);
xlabel('Time [s]');
ylabel('A/A_0');
title('Normalized surface area');
legend('Location', 'best');
grid on;

nexttile;
semilogy(tst.time_s, max(tst.avg_dissolution_rate, eps), 'k-o', 'LineWidth', 1.4, 'MarkerSize', 3, 'DisplayName', 'RTSPHEM TST');
hold on;
semilogy(phreeqc.time_s, max(phreeqc.avg_dissolution_rate, eps), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12, 'DisplayName', 'RTSPHEM PHREEQC-TST');
plotReferenceGroup(refs.rate, true);
xlabel('Time [s]');
ylabel('Average rate [mol cm^{-2} s^{-1}]');
title('Average reaction rate');
legend('Location', 'best');
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
    title('PHREEQC chemistry summary');
    grid on;
else
    text(0.5, 0.5, 'No PHREEQC summary log', 'HorizontalAlignment', 'center');
    axis off;
end

exportgraphics(fig, outputPath, 'Resolution', 200);
close(fig);
end

function refs = loadMolinsPartIIIReferences()
rootDir = fullfile('C:\Users\imgw\Downloads\10596_2019_9903_MOESM1_ESM', ...
    'SupplementaryMaterial', 'Results', 'Part-III');
refs = struct('volume', {{}}, 'surface', {{}}, 'rate', {{}}, ...
    'relativeVolume', {{}}, 'relativeSurface', {{}}, ...
    'sourceFiles', strings(0, 1), 'rootDir', string(rootDir));
if exist(rootDir, 'dir') ~= 7
    return;
end

chomboDir = fullfile(rootDir, 'chombo-crunch');
refs = addReferenceFile(refs, fullfile(chomboDir, 'volume.csv'), 'volume', 'Chombo-Crunch', 1, 2, 1, 1);
refs = addReferenceFile(refs, fullfile(chomboDir, 'area.csv'), 'surface', 'Chombo-Crunch', 1, 2, 1, 1);

vortexPath = fullfile(rootDir, 'Vortex', 'vol_vortex.csv');
refs = addReferenceFile(refs, vortexPath, 'volume', 'Vortex', 1, 3, 1, 1);
refs = addReferenceFile(refs, vortexPath, 'surface', 'Vortex', 1, 4, 1, 1);
refs = addReferenceFile(refs, vortexPath, 'rate', 'Vortex', 1, 5, 1, 1);

openFoamDir = fullfile(rootDir, 'openFOAM-DBS', 'Figure4');
chomboVolume0 = firstReferenceValue(refs.volume, "Chombo-Crunch");
chomboSurface0 = firstReferenceValue(refs.surface, "Chombo-Crunch");
refs = addReferenceFile(refs, fullfile(openFoamDir, 'vol_sim.csv'), 'relativeVolume', 'OpenFOAM-DBS simulation', 1, 2, 1, 1);
refs = addReferenceFile(refs, fullfile(openFoamDir, 'vol_expt.csv'), 'relativeVolume', 'Experiment', 1, 2, 1, 1);
refs = addReferenceFile(refs, fullfile(openFoamDir, 'area_sim.csv'), 'relativeSurface', 'OpenFOAM-DBS simulation', 1, 2, 1, 1);
refs = addReferenceFile(refs, fullfile(openFoamDir, 'area_expt.csv'), 'relativeSurface', 'Experiment', 1, 2, 1, 1);
refs = addReferenceFile(refs, fullfile(openFoamDir, 'vol_sim.csv'), 'volume', 'OpenFOAM-DBS simulation scaled', 1, 2, 1, chomboVolume0);
refs = addReferenceFile(refs, fullfile(openFoamDir, 'vol_expt.csv'), 'volume', 'Experiment scaled', 1, 2, 1, chomboVolume0);
refs = addReferenceFile(refs, fullfile(openFoamDir, 'area_sim.csv'), 'surface', 'OpenFOAM-DBS simulation scaled', 1, 2, 1, chomboSurface0);
refs = addReferenceFile(refs, fullfile(openFoamDir, 'area_expt.csv'), 'surface', 'Experiment scaled', 1, 2, 1, chomboSurface0);

dissolFoamPath = fullfile(rootDir, 'dissolFOAM', 'rate.csv');
refs = addReferenceFile(refs, dissolFoamPath, 'relativeSurface', 'dissolFOAM', 1, 3, 1, 1);
refs = addReferenceFile(refs, dissolFoamPath, 'surface', 'dissolFOAM', 1, 4, 1, 0.01);
refs = addReferenceFile(refs, dissolFoamPath, 'rate', 'dissolFOAM', 1, 5, 1, 1e-4);
end

function refs = addReferenceFile(refs, path, groupName, label, timeColumn, valueColumn, timeScale, valueScale)
if nargin < 8
    valueScale = 1;
end
if exist(path, 'file') ~= 2
    return;
end
tbl = readtable(path, 'VariableNamingRule', 'preserve');
if width(tbl) < max(timeColumn, valueColumn)
    return;
end
timeValues = double(tbl{:, timeColumn}) * timeScale;
dataValues = double(tbl{:, valueColumn}) * valueScale;
valid = isfinite(timeValues) & isfinite(dataValues);
entry = struct('label', string(label), 'time_s', timeValues(valid), ...
    'value', dataValues(valid), 'sourceFile', string(path));
refs.(groupName){end + 1} = entry;
if ~any(refs.sourceFiles == string(path))
    refs.sourceFiles(end + 1, 1) = string(path);
end
end

function entry = findReference(entries, label)
entry = [];
for iRef = 1:numel(entries)
    if entries{iRef}.label == string(label)
        entry = entries{iRef};
        return;
    end
end
end

function value = firstReferenceValue(entries, label)
entry = findReference(entries, label);
if isempty(entry) || isempty(entry.value)
    value = NaN;
else
    value = entry.value(1);
end
end

function value = lastReferenceValue(entries, label)
entry = findReference(entries, label);
if isempty(entry) || isempty(entry.value)
    value = NaN;
else
    value = entry.value(end);
end
end

function plotReferenceGroup(entries, logScale)
if nargin < 2
    logScale = false;
end
lineStyles = {'--', ':', '-.', '--', ':'};
colors = lines(max(numel(entries), 1));
for iRef = 1:numel(entries)
    entry = entries{iRef};
    y = entry.value;
    if logScale
        y = max(y, eps);
    end
    plot(entry.time_s, y, lineStyles{mod(iRef - 1, numel(lineStyles)) + 1}, ...
        'Color', colors(iRef, :), 'LineWidth', 1.0, ...
        'DisplayName', char(entry.label));
end
end

function writePartIIIBenchmarkManifest(packageDir, comparisonDir, mode, baseCfg, tstResult, phreeqcResult)
molinsRefs = loadMolinsPartIIIReferences();
manifest = struct();
manifest.schema_version = "benchmark_molins_partIII_phreeqc_tst_v1";
manifest.created_at = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
manifest.source_pdf = "C:\Users\imgw\Downloads\Simulation of mineral dissolution at the pore scale with evolving fluid-solid interfaces_ review of approaches and benchmark problem set.pdf";
manifest.mode = string(mode);
manifest.benchmark = "Part III experimental calcite-column moving-boundary benchmark";
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
    'calciteMolecularWeight_g_mol', 100, ...
    'solidDensity_g_cm3', 2.71, ...
    'inletConcentration_mol_cm3', 1.26e-5, ...
    'reynolds', 0.671, ...
    'peclet', 350, ...
    'damkohlerII', 3930, ...
    'endTime_s', 12000);
manifest.runner_config = baseCfg;
manifest.runs = struct( ...
    'tst', string(tstResult.resultsDir), ...
    'phreeqc_tst_match', string(phreeqcResult.resultsDir));
manifest.comparison_dir = string(comparisonDir);
manifest.outputs = struct( ...
    'summary_csv', string(fullfile(comparisonDir, 'comparison_summary.csv')), ...
    'summary_xlsx', string(fullfile(comparisonDir, 'comparison_summary.xlsx')), ...
    'nearest_time_csv', string(fullfile(comparisonDir, 'comparison_by_nearest_time.csv')), ...
    'nearest_time_xlsx', string(fullfile(comparisonDir, 'comparison_by_nearest_time.xlsx')), ...
    'plot_png', string(fullfile(comparisonDir, 'molins_partIII_volume_area.png')), ...
    'global_comparison_png', string(fullfile(comparisonDir, 'global_comparison.png')));
manifest.molins_partIII_references = struct( ...
    'root_dir', molinsRefs.rootDir, ...
    'source_files', molinsRefs.sourceFiles, ...
    'unit_notes', [ ...
    "RTSPHEM is 2D; mineral area and interface length are multiplied by the Part III column height 0.02 cm for comparison."; ...
    "The mask.tif is cropped by top 72 and bottom 38 rows, following the paper statement that these rows are solid layers outside the simulation."; ...
    "The TIF uses value 0 for black solid grains/walls and 255 for pore space in this local copy."; ...
    "OpenFOAM and experiment normalized Figure 4 curves are also plotted after scaling by the Chombo-Crunch initial value for absolute-volume/area panels."; ...
    "dissolFOAM surface area in mm^2 is converted to cm^2; rates in mol/m^2/s are converted to mol/cm^2/s."]);
manifest.notes = [ ...
    "Both RTSPHEM runs use the same cropped Part III image geometry, flow, diffusion, inlet chemistry, time span, mesh density, and export settings."; ...
    "PHREEQC uses phreeqc_tst_match, where PHREEQC updates aqueous chemistry while the kinetic rate law matches the legacy RTSPHEM TST moving-boundary rate."; ...
    "NMR/COMSOL synchronization is disabled for this chemistry benchmark export."; ...
    "Run modes are progressive: smoke is a one-step coarse check, partial is a 2400 s coarse run, and full is 12000 s with 600 s output spacing."];

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
