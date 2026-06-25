function summary = run_benchmark_molins_partII_phreeqc_tst(mode)
% run_benchmark_molins_partII_phreeqc_tst - Single-circle Part II benchmark.
%
% Runs the 0.1 cm x 0.05 cm centered-calcite-disk case with two chemistry
% paths:
%   1) external TST rate + PHREEQC equilibrium closure.
%   2) legacy RTSPHEM TST.

if nargin < 1 || isempty(mode)
    mode = "full";
end
mode = lower(string(mode));
if ~ismember(mode, ["full", "partial", "smoke", "diagnostic"])
    error('RTSPHEM:PartIIRunner:InvalidMode', ...
        'Mode must be full, partial, smoke, or diagnostic.');
end

scriptDir = fileparts(mfilename('fullpath'));
rtmDir = scriptDir;
projectRoot = fileparts(fileparts(rtmDir));
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));

stamp = datestr(now, 'yyyymmdd_HHMMSS');
packageDir = fullfile(projectRoot, 'outputs', ...
    'benchmark_molins_circle_external_tst_phreeqc', char(mode), stamp);
runsDir = fullfile(packageDir, 'runs');
comparisonDir = fullfile(packageDir, 'comparison');
ensureDir(runsDir);
ensureDir(comparisonDir);

fprintf('Molins Part II benchmark package: %s\n', packageDir);
baseCfg = buildPartIIConfig(projectRoot, rtmDir, mode);

phreeqcCfg = ConfigurePhreeqcRunGroup(baseCfg, 'external_tst_phreeqc');
phreeqcCfg.runName = 'benchmark_molins_partII_external_tst_phreeqc';
phreeqcCfg.resultsDir = fullfile(runsDir, 'phreeqc');
fprintf('\n=== Running external TST + PHREEQC closure Part II benchmark ===\n');
phreeqcResult = PNM_beauty3(phreeqcCfg);

tstCfg = baseCfg;
tstCfg.reactionModel = 'tst';
tstCfg.runName = 'benchmark_molins_partII_tst';
tstCfg.resultsDir = fullfile(runsDir, 'tst');
fprintf('\n=== Running legacy TST Part II benchmark ===\n');
tstResult = PNM_beauty3(tstCfg);

summary = exportPartIIComparison(tstResult.resultsDir, phreeqcResult.resultsDir, comparisonDir);
writeManifest(packageDir, comparisonDir, mode, baseCfg, tstResult, phreeqcResult);

fprintf('\nPart II benchmark export complete:\n');
fprintf('  Package: %s\n', packageDir);
fprintf('  Comparison: %s\n', comparisonDir);
end

function cfg = buildPartIIConfig(projectRoot, rtmDir, mode)
cfg = struct();
cfg.outputRoot = fullfile(projectRoot, 'outputs', 'benchmark_molins_circle_external_tst_phreeqc');
cfg.layoutType = 'single_circle';
cfg.useExternalGeometry = false;
cfg.useExternalDxfGeometry = false;
cfg.externalGeometryType = 'dxf';

cfg.circleRadius = 0.01;
cfg.targetLengthYAxis = 0.05;
cfg.targetAspectRatio = 2;
cfg.singleCircleCenter = [0.05, 0.025];
cfg.circleSpacing = 0;
cfg.circleSpacingXLeft = 0;

cfg.inletVelocity = 0.12;
cfg.flowDirection = 'left_to_right';
cfg.initialHydrogenConcentration = 1e-5;
cfg.diffusionCoefficient = 1e-5;
cfg.molarVolume = 36.9;
cfg.rateCoefficientTST = 8.912509381337459e-5;
cfg.characteristicLength = 0.05;

cfg.initialCalciumConcentration = 0;
cfg.initialCarbonConcentration = 0;
cfg.initialSodiumConcentration = 0;
cfg.initialChlorideConcentration = 0;
cfg.inletCalciumConcentration = 0;
cfg.inletCarbonConcentration = 0;
cfg.inletSodiumConcentration = 0;
cfg.inletChlorideConcentration = cfg.initialHydrogenConcentration;

cfg.phreeqcDatabasePath = ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc-m.dat', ...
    struct('databasePolicy', 'exact_local'));
cfg.phreeqcTemperatureC = 25;
cfg.phreeqcKineticsCorrectionFactor = 1;
cfg.phreeqcMaxSpecificSurfaceArea = Inf;
cfg.phreeqcBadStepMax = 5000;
cfg.phreeqcKineticsTolerance = 1e-8;
cfg.phreeqcMinHForPHMolL = 1e-7;
cfg.phreeqcMinActiveWaterVolumeFraction = 0;
cfg.phreeqcMinActiveWaterVolumeCm3 = 0;
cfg.phreeqcReactionWaterVolumeFloorFraction = 0;
cfg.phreeqcReactionWaterVolumeFloorCm3 = 0;
cfg.phreeqcReactionProjectionMinCells = 6;
cfg.phreeqcReactionProjectionMaxCells = 18;
cfg.phreeqcReactionProjectionMaxRings = 3;
cfg.phreeqcReactionProjectionFlowBias = 0;
cfg.phreeqcReactionProjectionMode = 'interface_cell';
cfg.phreeqcReactionWaterVolumeMode = 'cut_cell_agglomerated';
cfg.phreeqcReactNeutralInterfaceCells = false;
cfg.phreeqcSolutionWaterKg = 1;
cfg.phreeqcWriteSolutionWaterLine = false;
cfg.phreeqcKineticsReservoirMoles = 1;
cfg.phreeqcRunGroup = 'external_tst_phreeqc';
cfg.phreeqcRateLaw = 'tst_match';
cfg.phreeqcTstRateCoefficient = cfg.rateCoefficientTST;
cfg.phreeqcExportEvery = 1;

cfg.maximalStep = 60;
cfg.endTime = 2700;
cfg.processSliceCount = 50;
cfg.initialMacroscaleTimeStepSize = 54;
cfg.timeStepperType = 'linear';
cfg.maxTotalTimeSteps = [];
cfg.targetDissolutionSlices = [];
cfg.enableConcentrationCflLimit = false;
cfg.permeabilityRatioThreshold = 1e7;
cfg.maxExperimentWallSeconds = Inf;

cfg.meshTargetElementSizeCm = [];
cfg.meshNumPartitionsX = [];
cfg.meshNumPartitionsY = [];
cfg.numPartitionsMicroscale = 128;
cfg.dxfResolutionX = 200;
cfg.dxfResolutionY = 100;

cfg.exportEvery = 1;
cfg.exportDXF = true;
cfg.saveMainPlot = true;
cfg.saveIndividualPlots = true;
cfg.saveInterfaceMask = true;
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

if mode == "smoke"
    % Guarded external TST + PHREEQC smoke: keep dt below the advective residence time
    % (L/u ~= 0.83 s) so operator splitting can resolve a downstream H tail.
    cfg.endTime = 2.7;
    cfg.processSliceCount = 5;
    cfg.maximalStep = 0.54;
    cfg.initialMacroscaleTimeStepSize = 0.54;
    cfg.exportDXF = false;
    cfg.saveMainPlot = true;
    cfg.saveIndividualPlots = true;
    cfg.saveInterfaceMask = false;
    cfg.saveFinalPlot = false;
elseif mode == "diagnostic"
    % Short diagnostic Part II run. It keeps the same physical parameters
    % and external TST + PHREEQC coupling, but uses a reduced grid so the per-cell
    % PHREEQC batch can be checked before committing to expensive runs.
    cfg.endTime = 5.4;
    cfg.processSliceCount = 10;
    cfg.maximalStep = 0.54;
    cfg.initialMacroscaleTimeStepSize = 0.54;
    cfg.meshNumPartitionsX = 128;
    cfg.meshNumPartitionsY = 64;
    cfg.numPartitionsMicroscale = 64;
    cfg.dxfResolutionX = 160;
    cfg.dxfResolutionY = 80;
    cfg.exportDXF = false;
    cfg.saveMainPlot = true;
    cfg.saveIndividualPlots = true;
    cfg.saveInterfaceMask = false;
    cfg.saveFinalPlot = true;
elseif mode == "partial"
    cfg.endTime = 270;
    cfg.processSliceCount = 5;
    cfg.exportDXF = false;
    cfg.saveFinalPlot = false;
end
end

function summary = exportPartIIComparison(tstRunDir, phreeqcRunDir, outputDir)
tst = readGlobalEvolution(tstRunDir);
phreeqc = readGlobalEvolution(phreeqcRunDir);
aligned = alignByNearestTime(tst, phreeqc);
writetable(aligned, fullfile(outputDir, 'comparison_by_nearest_time.csv'));
writetable(aligned, fullfile(outputDir, 'comparison_by_nearest_time.xlsx'));

summary = buildSummary(aligned);
writetable(summary, fullfile(outputDir, 'comparison_summary.csv'));
writetable(summary, fullfile(outputDir, 'comparison_summary.xlsx'));

plotComparison(aligned, fullfile(outputDir, 'global_comparison.png'));
end

function tableData = readGlobalEvolution(runDir)
csvPath = fullfile(runDir, 'global_evolution_log.csv');
if exist(csvPath, 'file') ~= 2
    error('RTSPHEM:PartIIRunner:MissingGlobalLog', 'Missing global log: %s', csvPath);
end
tableData = readtable(csvPath);
end

function aligned = alignByNearestTime(tst, phreeqc)
numRows = height(tst);
aligned = table();
aligned.time_s = tst.time_s;
aligned.tst_porosity = tst.porosity;
aligned.tst_k_k0 = tst.k_k0;
aligned.tst_rate = tst.avg_dissolution_rate;
aligned.tst_interface_tst_rate = optionalTableColumn(tst, 'interface_tst_rate_mol_cm2_s', NaN);
aligned.tst_flux_apparent_rate = optionalTableColumn(tst, 'flux_apparent_rate_mol_cm2_s', NaN);
aligned.tst_surface_area_cm2 = tst.surface_area_cm2;
aligned.tst_grain_volume_cm3 = tst.grain_volume_cm3;
aligned.phreeqc_time_s = nan(numRows, 1);
aligned.phreeqc_porosity = nan(numRows, 1);
aligned.phreeqc_k_k0 = nan(numRows, 1);
aligned.phreeqc_rate = nan(numRows, 1);
aligned.phreeqc_interface_tst_rate = nan(numRows, 1);
aligned.phreeqc_flux_apparent_rate = nan(numRows, 1);
aligned.phreeqc_surface_area_cm2 = nan(numRows, 1);
aligned.phreeqc_grain_volume_cm3 = nan(numRows, 1);
for iRow = 1:numRows
    [~, idx] = min(abs(phreeqc.time_s - tst.time_s(iRow)));
    aligned.phreeqc_time_s(iRow) = phreeqc.time_s(idx);
    aligned.phreeqc_porosity(iRow) = phreeqc.porosity(idx);
    aligned.phreeqc_k_k0(iRow) = phreeqc.k_k0(idx);
    aligned.phreeqc_rate(iRow) = phreeqc.avg_dissolution_rate(idx);
    phreeqcInterfaceRate = optionalTableColumn(phreeqc, 'interface_tst_rate_mol_cm2_s', NaN);
    phreeqcFluxRate = optionalTableColumn(phreeqc, 'flux_apparent_rate_mol_cm2_s', NaN);
    aligned.phreeqc_interface_tst_rate(iRow) = phreeqcInterfaceRate(idx);
    aligned.phreeqc_flux_apparent_rate(iRow) = phreeqcFluxRate(idx);
    aligned.phreeqc_surface_area_cm2(iRow) = phreeqc.surface_area_cm2(idx);
    aligned.phreeqc_grain_volume_cm3(iRow) = phreeqc.grain_volume_cm3(idx);
end
aligned.delta_porosity = aligned.phreeqc_porosity - aligned.tst_porosity;
aligned.delta_k_k0 = aligned.phreeqc_k_k0 - aligned.tst_k_k0;
aligned.delta_rate = aligned.phreeqc_rate - aligned.tst_rate;
aligned.delta_interface_tst_rate = aligned.phreeqc_interface_tst_rate - aligned.tst_interface_tst_rate;
aligned.delta_flux_apparent_rate = aligned.phreeqc_flux_apparent_rate - aligned.tst_flux_apparent_rate;
aligned.delta_surface_area_cm2 = aligned.phreeqc_surface_area_cm2 - aligned.tst_surface_area_cm2;
aligned.delta_grain_volume_cm3 = aligned.phreeqc_grain_volume_cm3 - aligned.tst_grain_volume_cm3;
end

function values = optionalTableColumn(tableData, columnName, defaultValue)
if any(strcmp(tableData.Properties.VariableNames, columnName))
    values = tableData.(columnName);
else
    values = repmat(defaultValue, height(tableData), 1);
end
end

function summary = buildSummary(aligned)
metric = [
    "steps"
    "final_time_s"
    "final_porosity"
    "final_k_k0"
    "final_rate"
    "final_interface_tst_rate"
    "final_flux_apparent_rate"
    "final_surface_area_cm2"
    "final_grain_volume_cm3"
    "max_abs_delta_porosity"
    "max_abs_delta_k_k0"
    "max_abs_delta_rate"
    "max_abs_delta_interface_tst_rate"
    "max_abs_delta_flux_apparent_rate"
    "max_abs_delta_surface_area_cm2"
    "max_abs_delta_grain_volume_cm3"
    ];
tstValue = [
    height(aligned)
    aligned.time_s(end)
    aligned.tst_porosity(end)
    aligned.tst_k_k0(end)
    aligned.tst_rate(end)
    aligned.tst_interface_tst_rate(end)
    aligned.tst_flux_apparent_rate(end)
    aligned.tst_surface_area_cm2(end)
    aligned.tst_grain_volume_cm3(end)
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    ];
phreeqcValue = [
    height(aligned)
    aligned.phreeqc_time_s(end)
    aligned.phreeqc_porosity(end)
    aligned.phreeqc_k_k0(end)
    aligned.phreeqc_rate(end)
    aligned.phreeqc_interface_tst_rate(end)
    aligned.phreeqc_flux_apparent_rate(end)
    aligned.phreeqc_surface_area_cm2(end)
    aligned.phreeqc_grain_volume_cm3(end)
    max(abs(aligned.delta_porosity), [], 'omitnan')
    max(abs(aligned.delta_k_k0), [], 'omitnan')
    max(abs(aligned.delta_rate), [], 'omitnan')
    max(abs(aligned.delta_interface_tst_rate), [], 'omitnan')
    max(abs(aligned.delta_flux_apparent_rate), [], 'omitnan')
    max(abs(aligned.delta_surface_area_cm2), [], 'omitnan')
    max(abs(aligned.delta_grain_volume_cm3), [], 'omitnan')
    ];
deltaPhreeqcMinusTst = phreeqcValue - tstValue;
summary = table(metric, tstValue, phreeqcValue, deltaPhreeqcMinusTst);
end

function plotComparison(aligned, outputPath)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 900]);
t = aligned.time_s;
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(t, aligned.tst_porosity, 'k-', 'LineWidth', 1.5); hold on;
plot(t, aligned.phreeqc_porosity, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Porosity'); title('Porosity'); legend('TST', 'External TST + PHREEQC', 'Location', 'best'); grid on;

nexttile;
plot(t, aligned.tst_grain_volume_cm3, 'k-', 'LineWidth', 1.5); hold on;
plot(t, aligned.phreeqc_grain_volume_cm3, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Grain volume (cm^3)'); title('Grain Volume'); legend('TST', 'External TST + PHREEQC', 'Location', 'best'); grid on;

nexttile;
plot(t, aligned.tst_surface_area_cm2, 'k-', 'LineWidth', 1.5); hold on;
plot(t, aligned.phreeqc_surface_area_cm2, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Surface area (cm^2)'); title('Surface Area'); legend('TST', 'External TST + PHREEQC', 'Location', 'best'); grid on;

nexttile;
semilogy(t, abs(aligned.tst_rate), 'k-', 'LineWidth', 1.5); hold on;
semilogy(t, abs(aligned.phreeqc_rate), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('|rate| (mol cm^{-2} s^{-1})'); title('Average Dissolution Rate'); legend('TST', 'External TST + PHREEQC', 'Location', 'best'); grid on;

exportgraphics(fig, outputPath, 'Resolution', 200);
close(fig);
end

function writeManifest(packageDir, comparisonDir, mode, baseCfg, tstResult, phreeqcResult)
manifest = struct();
manifest.schema_version = "benchmark_molins_circle_external_tst_phreeqc_unified_h_v1";
manifest.created_at = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
manifest.mode = string(mode);
manifest.benchmark = "Molins Parts I/II circular calcite grain moving-boundary benchmark";
manifest.runner = string(mfilename('fullpath'));
manifest.runner_config = baseCfg;
manifest.runs = struct('tst', string(tstResult.resultsDir), 'phreeqc', string(phreeqcResult.resultsDir));
manifest.phreeqc_chemistry_mode = "external_tst_phreeqc";
manifest.phreeqc_chemistry_semantics = "explicit external TST rate + PHREEQC equilibrium closure";
manifest.comparison_dir = string(comparisonDir);
manifest.outputs = struct( ...
    'summary_csv', string(fullfile(comparisonDir, 'comparison_summary.csv')), ...
    'summary_xlsx', string(fullfile(comparisonDir, 'comparison_summary.xlsx')), ...
    'nearest_time_csv', string(fullfile(comparisonDir, 'comparison_by_nearest_time.csv')), ...
    'nearest_time_xlsx', string(fullfile(comparisonDir, 'comparison_by_nearest_time.xlsx')), ...
    'global_comparison_png', string(fullfile(comparisonDir, 'global_comparison.png')));
manifest.mode_description = describeBenchmarkMode(mode, baseCfg);
manifest.notes = [
    "External TST + PHREEQC closure uses the unified H+ field: the same H+ transport state drives reaction sink, moving boundary, PHREEQC chemistry, and exports."
    "H/Ca/C/Na/Cl are transported without reaction first; the transported H+ field defines the first-order TST CaCO3 amount prescribed to PHREEQC."
    "PHREEQC speciation uses a cut-cell mixing-volume floor for interface sliver cells, so the first-order rate is not artificially capped while extreme mol/kgw inputs are avoided."
    "The moving boundary uses the same external TST + PHREEQC prescribed reaction-rate field; the existing PNM loop applies that field to the next level-set evolution step."
    "Legacy TST is run with the same grid, geometry, flow, diffusion, and time settings for visual comparison."
    "The legacy stability_diagnostics_log mass_res field is a single-H-component diagnostic and is not a strict PHREEQC multi-component conservation metric."
    "Rate comparison should use matched definitions: interface_tst_rate_mol_cm2_s is the interface-area TST/CaCO3 rate, while flux_apparent_rate_mol_cm2_s is the outlet H+ flux-deficit apparent acid-consumption rate."
    ];
writeJsonFileLocal(fullfile(comparisonDir, 'benchmark_manifest.json'), manifest);
end

function description = describeBenchmarkMode(mode, cfg)
description = struct();
description.mode = string(mode);
description.endTime_s = cfg.endTime;
description.processSliceCount = cfg.processSliceCount;
description.initialMacroscaleTimeStepSize_s = cfg.initialMacroscaleTimeStepSize;
description.maximalStep_s = cfg.maximalStep;
description.meshNumPartitionsX = cfg.meshNumPartitionsX;
description.meshNumPartitionsY = cfg.meshNumPartitionsY;
description.numPartitionsMicroscale = cfg.numPartitionsMicroscale;
description.isDiagnosticReducedGrid = strcmpi(string(mode), "diagnostic");
if description.isDiagnosticReducedGrid
    description.note = [
        "Diagnostic mode keeps the same Part II physical parameters and solver path, " + ...
        "but uses a reduced grid and short small-step horizon to test H+ wake, " + ...
        "PHREEQC prescribed/result closure, and rate trend. It is not a " + ...
        "full-resolution Molins benchmark replacement."
        ];
else
    description.note = "Benchmark runner mode.";
end
end

function ensureDir(pathValue)
if exist(pathValue, 'dir') ~= 7
    mkdir(pathValue);
end
end

function writeJsonFileLocal(pathValue, data)
fid = fopen(pathValue, 'w');
if fid == -1
    error('RTSPHEM:PartIIRunner:JsonOpenFailed', 'Cannot write JSON file: %s', pathValue);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(data, PrettyPrint=true));
clear cleanupObj;
end
