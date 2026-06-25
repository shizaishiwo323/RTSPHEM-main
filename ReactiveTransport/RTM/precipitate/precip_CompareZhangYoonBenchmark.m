function report = precip_CompareZhangYoonBenchmark(runDir, referenceCsv)
%PRECIP_COMPAREZHANGYOONBENCHMARK Compare RTSPHEM outputs with Zhang/Yoon curves.
%
% The reference CSV must contain real digitized literature values. This
% function validates the schema but does not invent missing reference data.

if nargin < 2
    error('RTSPHEM:Precipitate:InvalidComparisonInputs', ...
        'Usage: precip_CompareZhangYoonBenchmark(runDir, referenceCsv).');
end
runDir = char(string(runDir));
referenceCsv = char(string(referenceCsv));
if ~isfolder(runDir)
    error('RTSPHEM:Precipitate:InvalidComparisonInputs', ...
        'runDir does not exist: %s', runDir);
end
if exist(referenceCsv, 'file') ~= 2
    error('RTSPHEM:Precipitate:InvalidComparisonInputs', ...
        'referenceCsv does not exist: %s', referenceCsv);
end

areaCsv = fullfile(runDir, 'precipitation_area_timeseries.csv');
comparisonCsv = fullfile(runDir, 'benchmark_comparison_times_log.csv');
if exist(areaCsv, 'file') ~= 2
    error('RTSPHEM:Precipitate:MissingComparisonInput', ...
        'Missing precipitation area time-series: %s', areaCsv);
end
if exist(comparisonCsv, 'file') ~= 2
    error('RTSPHEM:Precipitate:MissingComparisonInput', ...
        'Missing benchmark comparison log: %s', comparisonCsv);
end

ensureBenchmarkHelpersOnPath();
simulation = readSimulationAreaCurves(areaCsv);
simulationTable = simulation.table;
reference = precip_LoadReferenceCurves(referenceCsv);
referenceTable = reference.table;
benchmarkTable = readtable(comparisonCsv);
runDiagnostics = readRunDiagnostics(runDir);

combinedTable = [simulationTable; referenceTable];
outputCsv = fullfile(runDir, 'zhang_yoon_area_comparison.csv');
outputPng = fullfile(runDir, 'zhang_yoon_area_comparison.png');
outputReport = fullfile(runDir, 'benchmark_comparison_report.md');
writetable(combinedTable, outputCsv);
writeComparisonPlot(combinedTable, outputPng);
writeComparisonReport(outputReport, runDir, referenceCsv, outputCsv, outputPng, ...
    simulationTable, referenceTable, benchmarkTable, simulation.notes, runDiagnostics);

report = struct();
report.runDir = string(runDir);
report.referenceCsv = string(referenceCsv);
report.outputCsv = string(outputCsv);
report.outputPng = string(outputPng);
report.outputReport = string(outputReport);
report.numSimulationRows = height(simulationTable);
report.numReferenceRows = height(referenceTable);
report.numBenchmarkRows = height(benchmarkTable);
report.regions = unique(combinedTable.region, 'stable');
report.cases = unique(combinedTable.("case"), 'stable');
report.runDiagnostics = runDiagnostics;
end

function ensureBenchmarkHelpersOnPath()
moduleRoot = fileparts(mfilename('fullpath'));
benchmarkDir = fullfile(moduleRoot, 'benchmark');
if exist('precip_LoadReferenceCurves', 'file') ~= 2
    addpath(benchmarkDir);
end
end

function simulation = readSimulationAreaCurves(areaCsv)
areaTable = readtable(areaCsv);
requiredColumns = {'time_s', 'first_pore_net_solid_area_cm2', ...
    'first_three_pores_net_solid_area_cm2'};
requireColumns(areaTable, requiredColumns, ...
    'RTSPHEM:Precipitate:InvalidSimulationAreaCurves', areaCsv);

timeMin = areaTable.time_s(:) ./ 60;
hasEntireDomainArea = ismember('total_net_solid_area_cm2', areaTable.Properties.VariableNames);
firstPoreArea = areaTable.first_pore_net_solid_area_cm2(:);
firstThreeArea = areaTable.first_three_pores_net_solid_area_cm2(:);

if hasEntireDomainArea
    entireDomainArea = areaTable.total_net_solid_area_cm2(:);
    source = repmat("RTSPHEM", 3 * numel(timeMin), 1);
    caseName = repmat("local_benchmark", 3 * numel(timeMin), 1);
    region = [repmat("entire_domain", numel(timeMin), 1); ...
        repmat("first_pore", numel(timeMin), 1); ...
        repmat("first_three_pores", numel(timeMin), 1)];
    timeMinCombined = [timeMin; timeMin; timeMin];
    areaCm2 = [entireDomainArea; firstPoreArea; firstThreeArea];
    simulation.notes = strings(0, 1);
else
    source = repmat("RTSPHEM", 2 * numel(timeMin), 1);
    caseName = repmat("local_benchmark", 2 * numel(timeMin), 1);
    region = [repmat("first_pore", numel(timeMin), 1); ...
        repmat("first_three_pores", numel(timeMin), 1)];
    timeMinCombined = [timeMin; timeMin];
    areaCm2 = [firstPoreArea; firstThreeArea];
    simulation.notes = "Legacy simulation area CSV lacks total_net_solid_area_cm2; RTSPHEM entire_domain simulation rows are omitted.";
end
areaNorm = normalizeWithinRegion(areaCm2, region);
note = repmat("from precipitation_area_timeseries.csv", numel(areaCm2), 1);
dataRole = repmat("simulation", numel(areaCm2), 1);

simulation.table = table(source, caseName, region, timeMinCombined, areaNorm, areaCm2, ...
    note, dataRole, 'VariableNames', {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note', 'data_role'});
end

function requireColumns(tableData, requiredColumns, errorId, path)
missing = setdiff(requiredColumns, tableData.Properties.VariableNames);
if ~isempty(missing)
    error(errorId, 'Missing required columns in %s: %s', path, strjoin(missing, ', '));
end
end

function areaNorm = normalizeWithinRegion(areaCm2, region)
areaNorm = NaN(size(areaCm2));
regions = unique(region, 'stable');
for iRegion = 1:numel(regions)
    mask = region == regions(iRegion);
    regionArea = areaCm2(mask);
    scale = max(abs(regionArea), [], 'omitnan');
    if isempty(scale) || ~isfinite(scale) || scale <= 0
        areaNorm(mask) = regionArea;
    else
        areaNorm(mask) = regionArea ./ scale;
    end
end
end

function writeComparisonPlot(combinedTable, outputPng)
regions = unique(combinedTable.region, 'stable');
figureWidth = max(900, 430 * numel(regions));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, figureWidth, 520]);
cleanupObj = onCleanup(@() closeFigureIfOpen(fig));
tiledlayout(fig, 1, numel(regions), 'TileSpacing', 'compact', 'Padding', 'compact');
for iRegion = 1:numel(regions)
    ax = nexttile;
    hold(ax, 'on');
    regionMask = combinedTable.region == regions(iRegion);
    regionData = combinedTable(regionMask, :);
    groups = unique(regionData(:, {'source', 'case', 'data_role'}), 'rows', 'stable');
    for iGroup = 1:height(groups)
        mask = regionData.source == groups.source(iGroup) & ...
            regionData.("case") == groups.("case")(iGroup) & ...
            regionData.data_role == groups.data_role(iGroup);
        groupData = sortrows(regionData(mask, :), 'time_min');
        plot(ax, groupData.time_min, groupData.precipitated_area_norm, ...
            '-o', 'LineWidth', 1.3, 'DisplayName', sprintf('%s %s', ...
            groups.source(iGroup), groups.("case")(iGroup)));
    end
    title(ax, strrep(regions(iRegion), '_', ' '));
    xlabel(ax, 'time [min]');
    ylabel(ax, 'precipitated area [normalized]');
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    legend(ax, 'Location', 'best', 'Interpreter', 'none');
end
print(fig, outputPng, '-dpng', '-r300');
end

function writeComparisonReport(outputReport, runDir, referenceCsv, outputCsv, outputPng, ...
    simulationTable, referenceTable, benchmarkTable, simulationNotes, runDiagnostics)
fid = fopen(outputReport, 'w');
if fid == -1
    error('RTSPHEM:Precipitate:ComparisonReportWriteFailed', ...
        'Cannot write comparison report: %s', outputReport);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '# Zhang/Yoon CaCO3 Precipitation Comparison Report\n\n');
fprintf(fid, '- Run directory: `%s`\n', runDir);
fprintf(fid, '- Reference CSV: `%s`\n', referenceCsv);
fprintf(fid, '- Combined output CSV: `%s`\n', outputCsv);
fprintf(fid, '- Comparison plot: `%s`\n', outputPng);
fprintf(fid, '- Simulation rows: %d\n', height(simulationTable));
fprintf(fid, '- Reference rows: %d\n', height(referenceTable));
fprintf(fid, '- Benchmark captured rows: %d\n\n', height(benchmarkTable));
fprintf(fid, '## Interpretation Boundary\n\n');
fprintf(fid, ['This report compares available RTSPHEM local benchmark outputs against ', ...
    'digitized Zhang/Yoon reference curves. It does not prove strict reproduction ', ...
    'when the geometry is the approximate local cylindrical-post layout or when ', ...
    'the reference CSV contains only partial digitized points.\n\n']);
fprintf(fid, ['The RTSPHEM normalized area in the plot is a simulation-internal normalization ', ...
    'within each simulation region using that region''s own maximum absolute net area. ', ...
    'If a reference row has finite cm2 area but missing normalized area, missing reference ', ...
    'normalized values are filled from cm2 by source/case/region maximum absolute area and ', ...
    'the row note records that conversion.\n\n']);
zhangSelectedRegions = ["zhang_upgradient", "zhang_middle", "zhang_downgradient"];
if any(ismember(referenceTable.region, zhangSelectedRegions))
    fprintf(fid, ['Zhang selected-pore regions are reference-only in this comparison scaffold. ', ...
        'They are source-defined Upgradient/middle/Downgradient pores from Zhang 2010 SI Figure S3(b), ', ...
        'and currently have no matching RTSPHEM simulation rows unless a future simulation export ', ...
        'defines the same source-specific regions.\n\n']);
end
writeRunDiagnosticsSection(fid, runDiagnostics);
if ~isempty(simulationNotes)
    fprintf(fid, '## Simulation Data Notes\n\n');
    for iNote = 1:numel(simulationNotes)
        fprintf(fid, '- %s\n', simulationNotes(iNote));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '## Data Sources\n\n');
sourceCases = unique(referenceTable(:, {'source', 'case', 'region'}), 'rows', 'stable');
for iRow = 1:height(sourceCases)
    fprintf(fid, '- %s / %s / %s\n', sourceCases.source(iRow), ...
        sourceCases.("case")(iRow), sourceCases.region(iRow));
end
clear cleanupObj;
end

function diagnostics = readRunDiagnostics(runDir)
diagnostics = struct();
diagnostics.hasMetadata = false;
diagnostics.phreeqcTransportMaxFactor = Inf;
diagnostics.hasStabilityDiagnostics = false;
diagnostics.stabilityFlags = strings(0, 1);

metadataFile = fullfile(runDir, 'run_metadata.json');
if exist(metadataFile, 'file') == 2
    try
        metadata = jsondecode(fileread(metadataFile));
        diagnostics.hasMetadata = true;
        if isfield(metadata, 'parameters') && ...
                isfield(metadata.parameters, 'phreeqcTransportMaxFactor') && ...
                ~isempty(metadata.parameters.phreeqcTransportMaxFactor)
            diagnostics.phreeqcTransportMaxFactor = double(metadata.parameters.phreeqcTransportMaxFactor);
        end
    catch
        diagnostics.hasMetadata = false;
    end
end

stabilityFile = fullfile(runDir, 'stability_diagnostics_log.csv');
if exist(stabilityFile, 'file') == 2
    stabilityTable = readtable(stabilityFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    if ismember('diagnostic_flag', stabilityTable.Properties.VariableNames)
        diagnostics.hasStabilityDiagnostics = true;
        diagnostics.stabilityFlags = unique(splitDiagnosticFlags(stabilityTable.diagnostic_flag), 'stable');
    end
end
end

function flags = splitDiagnosticFlags(flagColumn)
flags = strings(0, 1);
for iRow = 1:numel(flagColumn)
    rowFlags = split(string(flagColumn(iRow)), ';');
    rowFlags = strtrim(rowFlags);
    rowFlags = rowFlags(rowFlags ~= "" & rowFlags ~= "none");
    flags = [flags; rowFlags(:)]; %#ok<AGROW>
end
end

function writeRunDiagnosticsSection(fid, diagnostics)
hasFiniteLimiter = isfinite(diagnostics.phreeqcTransportMaxFactor);
hasFlags = diagnostics.hasStabilityDiagnostics && ~isempty(diagnostics.stabilityFlags);
if ~hasFiniteLimiter && ~hasFlags
    return;
end

fprintf(fid, '## Numerical Diagnostics\n\n');
if hasFiniteLimiter
    fprintf(fid, ['- finite PHREEQC transport limiter used: ', ...
        '`phreeqcTransportMaxFactor = %.15g`. This is a numerical guard against ', ...
        'transport overshoot before PHREEQC calls, not a Zhang/Yoon literature ', ...
        'parameter or calibration value.\n'], diagnostics.phreeqcTransportMaxFactor);
else
    fprintf(fid, '- PHREEQC transport limiter was not finite in metadata.\n');
end
if hasFlags
    fprintf(fid, '- Stability diagnostic flags observed: `%s`.\n', ...
        strjoin(diagnostics.stabilityFlags, '`, `'));
    if any(ismember(diagnostics.stabilityFlags, ["advective_cfl_gt_1", "mass_balance_drift", "overshoot_c"]))
        fprintf(fid, ['- These flags mean the run should be treated as diagnostic/output-chain ', ...
            'evidence rather than quantitative Zhang/Yoon validation until transport ', ...
            'stability and mass balance are resolved.\n']);
    end
end
fprintf(fid, '\n');
end

function closeFigureIfOpen(fig)
if ishghandle(fig)
    close(fig);
end
end
