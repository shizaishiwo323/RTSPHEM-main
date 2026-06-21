function summary = ComparePhreeqcRunToPe0p1(phreeqcRunDir, baselineDir, outputDir)
% ComparePhreeqcRunToPe0p1 - Compare a PHREEQC RTM run against the Pe0p1 TST baseline.
%
% Inputs:
%   phreeqcRunDir: results directory from run_single_pnm_mine_phreeqc.m.
%   baselineDir  : optional baseline directory. Defaults to outputs/rtm_runs/Pe0p1.
%   outputDir    : optional output directory. Defaults to <phreeqcRunDir>/comparison_to_Pe0p1.
%
% Outputs:
%   summary: MATLAB table also written to comparison_summary.csv.

if nargin < 1 || isempty(phreeqcRunDir)
    error('RTSPHEM:ComparePhreeqc:MissingRunDir', 'phreeqcRunDir is required.');
end

thisDir = fileparts(mfilename('fullpath'));
rtmDir = fileparts(thisDir);
projectRoot = fileparts(fileparts(rtmDir));

if nargin < 2 || isempty(baselineDir)
    baselineDir = fullfile(projectRoot, 'outputs', 'rtm_runs', 'Pe0p1');
end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(phreeqcRunDir, 'comparison_to_Pe0p1');
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

baselineGlobal = readGlobalEvolution(baselineDir);
phreeqcGlobal = readGlobalEvolution(phreeqcRunDir);
phreeqcSummary = readOptionalTable(fullfile(phreeqcRunDir, 'phreeqc_results', 'phreeqc_summary_log.csv'));

summary = buildSummaryTable(baselineGlobal, phreeqcGlobal);
writetable(summary, fullfile(outputDir, 'comparison_summary.csv'));

aligned = alignByNearestTime(baselineGlobal, phreeqcGlobal);
writetable(aligned, fullfile(outputDir, 'comparison_by_nearest_time.csv'));

plotGlobalComparison(baselineGlobal, phreeqcGlobal, phreeqcSummary, ...
    fullfile(outputDir, 'global_comparison.png'));

fprintf('Comparison written to: %s\n', outputDir);
end

function tbl = readGlobalEvolution(runDir)
path = fullfile(runDir, 'global_evolution_log.csv');
if ~exist(path, 'file')
    error('RTSPHEM:ComparePhreeqc:MissingGlobalLog', 'Missing global evolution log: %s', path);
end
tbl = readtable(path, 'VariableNamingRule', 'preserve');
end

function tbl = readOptionalTable(path)
if exist(path, 'file')
    tbl = readtable(path, 'VariableNamingRule', 'preserve');
else
    tbl = table();
end
end

function summary = buildSummaryTable(baseline, phreeqc)
thresholds = [0.55, 0.60, 0.80, 0.95, 0.99, 1.00];
metric = strings(0, 1);
baselineValue = [];
phreeqcValue = [];
ratioPhreeqcToBaseline = [];

[metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline] = addMetric( ...
    metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline, ...
    "steps", height(baseline), height(phreeqc));
[metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline] = addMetric( ...
    metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline, ...
    "final_time_s", baseline.time_s(end), phreeqc.time_s(end));
[metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline] = addMetric( ...
    metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline, ...
    "final_porosity", baseline.porosity(end), phreeqc.porosity(end));
[metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline] = addMetric( ...
    metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline, ...
    "final_permeability_mD", baseline.permeability_mD(end), phreeqc.permeability_mD(end));
[metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline] = addMetric( ...
    metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline, ...
    "final_k_k0", baseline.k_k0(end), phreeqc.k_k0(end));

for iThreshold = 1:numel(thresholds)
    name = "time_to_porosity_" + replace(sprintf('%.2f', thresholds(iThreshold)), '.', 'p') + "_s";
    baselineTime = timeToPorosity(baseline, thresholds(iThreshold));
    phreeqcTime = timeToPorosity(phreeqc, thresholds(iThreshold));
    [metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline] = addMetric( ...
        metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline, ...
        name, baselineTime, phreeqcTime);
end

summary = table(metric, baselineValue, phreeqcValue, ratioPhreeqcToBaseline);
end

function [metric, baselineValue, phreeqcValue, ratio] = addMetric(metric, baselineValue, phreeqcValue, ratio, name, baseVal, phreeqcVal)
metric(end + 1, 1) = name;
baselineValue(end + 1, 1) = baseVal;
phreeqcValue(end + 1, 1) = phreeqcVal;
if isfinite(baseVal) && abs(baseVal) > eps
    ratio(end + 1, 1) = phreeqcVal / baseVal;
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

function aligned = alignByNearestTime(baseline, phreeqc)
phreeqcTime = phreeqc.time_s;
n = numel(phreeqcTime);
baselineIndex = zeros(n, 1);
for i = 1:n
    [~, baselineIndex(i)] = min(abs(baseline.time_s - phreeqcTime(i)));
end

aligned = table();
aligned.time_s = phreeqcTime;
aligned.baseline_time_s = baseline.time_s(baselineIndex);
aligned.baseline_porosity = baseline.porosity(baselineIndex);
aligned.phreeqc_porosity = phreeqc.porosity;
aligned.delta_porosity = phreeqc.porosity - baseline.porosity(baselineIndex);
aligned.baseline_k_k0 = baseline.k_k0(baselineIndex);
aligned.phreeqc_k_k0 = phreeqc.k_k0;
aligned.delta_k_k0 = phreeqc.k_k0 - baseline.k_k0(baselineIndex);
aligned.baseline_rate = baseline.avg_dissolution_rate(baselineIndex);
aligned.phreeqc_rate = phreeqc.avg_dissolution_rate;
end

function plotGlobalComparison(baseline, phreeqc, phreeqcSummary, outputPath)
fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 900]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(baseline.time_s, baseline.porosity, 'k-', 'LineWidth', 1.5);
hold on;
plot(phreeqc.time_s, phreeqc.porosity, 'r.-', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Porosity');
legend({'Pe0p1 TST', 'PHREEQC'}, 'Location', 'best');
grid on;

nexttile;
semilogy(baseline.time_s, max(baseline.k_k0, eps), 'k-', 'LineWidth', 1.5);
hold on;
semilogy(phreeqc.time_s, max(phreeqc.k_k0, eps), 'r.-', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('k/k0');
legend({'Pe0p1 TST', 'PHREEQC'}, 'Location', 'best');
grid on;

nexttile;
semilogy(baseline.time_s, max(baseline.avg_dissolution_rate, eps), 'k-', 'LineWidth', 1.5);
hold on;
semilogy(phreeqc.time_s, max(phreeqc.avg_dissolution_rate, eps), 'r.-', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Dissolution rate');
legend({'Pe0p1 TST', 'PHREEQC'}, 'Location', 'best');
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
