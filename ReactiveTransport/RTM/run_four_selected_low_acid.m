function [summary, batchDir] = run_four_selected_low_acid(action, options)
% Run the four selected 600 x 400 um geometries to porosity 0.50.
%
%   run_four_selected_low_acid('design')
%   run_four_selected_low_acid('production')
%   run_four_selected_low_acid('resume', struct('batchDir', <existing>))
%   run_four_selected_low_acid('resume', struct('batchDir', <existing>, ...
%       'samples', "M14_25", 'writeSharedArtifacts', false))

% Runs are sequential to avoid four simultaneous Stokes/PHREEQC memory peaks.

if nargin < 1 || isempty(action); action = 'production'; end
if nargin < 2 || isempty(options); options = struct(); end
action = lower(string(action));
if ~ismember(action, ["design", "production", "resume"])
    error('RTSPHEM:FourGeometry:UnknownAction', ...
        'action must be design, production, or resume.');
end

rtmDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(rtmDir));
geometryRoot = fullfile('C:\Users\imgw\Documents\Codex\震电模拟', ...
    'Wang_SE_experiment', 'generated_geometries', ...
    'selected_relativeK_600x400_384x256');
fullDesign = buildDesign(geometryRoot);
assert(all(isfile(fullDesign.GeometryFile)), 'A selected geometry MAT is missing.');
design = selectDesign(fullDesign, options);
if isfield(options, 'meshNumPartitionsX')
    design.MeshX(:) = options.meshNumPartitionsX;
end
if isfield(options, 'meshNumPartitionsY')
    design.MeshY(:) = options.meshNumPartitionsY;
end

if action == "design"
    summary = design;
    batchDir = "";
    return;
end

if action == "resume"
    if ~isfield(options, 'batchDir') || ~isfolder(options.batchDir)
        error('RTSPHEM:FourGeometry:MissingBatchDir', ...
            'resume requires options.batchDir pointing to an existing batch.');
    end
    batchDir = char(options.batchDir);
else
    if isfield(options, 'batchDir') && strlength(string(options.batchDir)) > 0
        batchDir = char(options.batchDir);
    else
        stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        batchDir = fullfile(projectRoot, 'outputs', 'rtm_batches', ...
            ['selected_four_low_flow_low_acid_phi050_100_' stamp]);
    end
    if ~isfolder(batchDir); mkdir(batchDir); end
end
writeSharedArtifacts = ~isfield(options, 'writeSharedArtifacts') || ...
    logical(options.writeSharedArtifacts);
if writeSharedArtifacts
    writetable(fullDesign, fullfile(batchDir, 'input_geometry_manifest.csv'));
end

summary = table();
for i = 1:height(design)
    baseDir = fullfile(batchDir, char(design.Sample(i)));
    [complete, completedRow] = completedRun(baseDir, design(i, :));
    if complete
        summary = [summary; completedRow]; %#ok<AGROW>
        if writeSharedArtifacts; writeSummary(summary, batchDir); end
        continue;
    end
    resultsDir = baseDir;
    if isfolder(resultsDir)
        stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        resultsDir = [baseDir '_retry_' stamp];
    end
    mkdir(resultsDir);

    reactiveCflSafety = 0.2;
    if design.Sample(i) == "M14_25"
        reactiveCflSafety = 0.05;
    elseif ismember(design.Sample(i), ["M14_24", "M14_26", "M14_34"])
        reactiveCflSafety = 0.05;
    end
    if isfield(options, 'reactiveCflSafety')
        reactiveCflSafety = options.reactiveCflSafety;
    end
    phreeqcGridHmax = 0.001;
    if isfield(options, 'phreeqcGridHmax')
        phreeqcGridHmax = options.phreeqcGridHmax;
    end
    wallLimitSeconds = Inf;
    if isfield(options, 'maxExperimentWallSeconds')
        wallLimitSeconds = options.maxExperimentWallSeconds;
    end
    adaptiveGrowthFactor = 1.05;
    if isfield(options, 'adaptiveGrowthFactor')
        adaptiveGrowthFactor = options.adaptiveGrowthFactor;
    end
    writeCaseManifest(resultsDir, design(i, :), reactiveCflSafety, ...
        phreeqcGridHmax, wallLimitSeconds, adaptiveGrowthFactor);
    cfg = struct( ...
        'resultsDir', resultsDir, ...
        'loadExistingGeometry', true, ...
        'geometryLoadFile', char(design.GeometryFile(i)), ...
        'inletVelocity', design.InletVelocity_cm_s(i), ...
        'acidConcentration_mol_L', design.AcidConcentration_mol_L(i), ...
        'rateCoefficientTST', design.RateCoefficient_mol_dm2_s(i), ...
        'phreeqcKineticsCorrectionFactor', design.PhreeqcRateCorrection(i), ...
        'targetDissolutionSlices', design.TargetSlices(i), ...
        'targetDissolutionFinalPorosity', design.TargetFinalPorosity(i), ...
        'initialMacroscaleTimeStepSize', design.InitialTimeStep_s(i), ...
        'numPartitionsMicroscale', design.MeshY(i), ...
        'meshNumPartitionsX', design.MeshX(i), ...
        'meshNumPartitionsY', design.MeshY(i), ...
        'phreeqcGridMode', 'separate', ...
        'phreeqcGridHmax', phreeqcGridHmax, ...
        'phreeqcExportEvery', Inf, ...
        'phreeqcEarlyExportSteps', 6, ...
        'phreeqcExportSeismoelectric', true, ...
        'phreeqcWriteSpeciesTable', false, ...
        'permeabilityOnExportOnly', true, ...
        'transportSolverBackend', 'cut_cell_fv', ...
        'transportUpwindMode', 'exp', ...
        'transportInletBoundaryMode', 'dirichlet', ...
        'transportFlowProjectionMode', 'conservative_edge_average', ...
        'enableHardCflLimit', true, ...
        'hardCflLimitInitialGrid', true, ...
        'reactiveCflSafety', reactiveCflSafety, ...
        'adaptiveGrowthFactor', adaptiveGrowthFactor, ...
        'maxExperimentWallSeconds', wallLimitSeconds, ...
        'saveRealtimePlot', false, ...
        'showDebugFigures', false, ...
        'saveIndividualPlots', false);

    fprintf('\n[%d/%d] %s -> %s\n', i, height(design), ...
        design.Sample(i), resultsDir);
    timer = tic;
    try
        result = run_phreeqc_se(cfg);
        runtimeSeconds = toc(timer);
        status = "Success";
        message = "";
        if result.timedOut || result.exportedDissolutionSlices ~= design.TargetSlices(i) || ...
                result.finalPorosity < design.TargetFinalPorosity(i) - 1e-6
            status = "Failed";
            message = "Run did not meet timeout, slice-count, or final-porosity gate.";
        end
        row = resultRow(design(i, :), status, message, resultsDir, ...
            result.exportedDissolutionSlices, result.finalPorosity, ...
            result.finalPermeability, result.timedOut, runtimeSeconds);
    catch ME
        report = getReport(ME, 'extended', 'hyperlinks', 'off');
        fid = fopen(fullfile(resultsDir, 'case_failure.txt'), 'w');
        if fid ~= -1
            fprintf(fid, '%s\n', report);
            fclose(fid);
        end
        fprintf(2, '%s\n', report);
        row = resultRow(design(i, :), "Failed", ...
            string(report), ...
            resultsDir, NaN, NaN, NaN, false, toc(timer));
    end
    summary = [summary; row]; %#ok<AGROW>
    if writeSharedArtifacts; writeSummary(summary, batchDir); end
    close all hidden;
end

if any(summary.Status ~= "Success" & summary.Status ~= "SkippedComplete")
    error('RTSPHEM:FourGeometry:Incomplete', ...
        'At least one geometry failed; see batch_summary.csv.');
end
end

function design = selectDesign(design, options)
if ~isfield(options, 'samples') || isempty(options.samples); return; end
requested = string(options.samples);
requested = unique(requested(:), 'stable');
missing = setdiff(requested, design.Sample);
assert(isempty(missing), 'Unknown requested sample: %s', strjoin(missing, ', '));
design = design(ismember(design.Sample, requested), :);
end

function design = buildDesign(root)
sample = ["M14_24"; "M14_25"; "M14_26"; "M14_34"];
geometryFile = fullfile(root, sample + "_circle_geometry_600x400.mat");
selection = readtable(fullfile(root, 'selected_window_summary.csv'), ...
    'VariableNamingRule', 'preserve');
[found, order] = ismember(sample, string(selection.sample));
connected = strcmpi(string(selection.left_right_connected(order)), "true");
assert(all(found) && all(connected), ...
    'Selected geometry connectivity audit is missing or failed.');
expectedInitialPorosity = selection.rtm_porosity(order);
n = numel(sample);
inletVelocity = repmat(1e-4, n, 1);       % cm/s, Pe=0.01
acidConcentration = repmat(0.020, n, 1);  % mol/L HCl
rateCoefficient = repmat(0.0013550135501355, n, 1);
rateCorrection = rateCoefficient / 8.9e-5;
targetSlices = repmat(100, n, 1);
targetFinalPorosity = repmat(0.50, n, 1);
initialTimeStep = [0.005; repmat(0.00075, n - 1, 1)];
meshX = repmat(384, n, 1);
meshY = repmat(256, n, 1);
design = table(sample, geometryFile, expectedInitialPorosity, inletVelocity, acidConcentration, ...
    rateCoefficient, rateCorrection, targetSlices, targetFinalPorosity, ...
    initialTimeStep, meshX, meshY, 'VariableNames', ...
    {'Sample', 'GeometryFile', 'ExpectedInitialPorosity', 'InletVelocity_cm_s', ...
    'AcidConcentration_mol_L', 'RateCoefficient_mol_dm2_s', ...
    'PhreeqcRateCorrection', 'TargetSlices', 'TargetFinalPorosity', ...
    'InitialTimeStep_s', 'MeshX', 'MeshY'});
end

function [complete, row] = completedRun(resultsDir, design)
complete = false;
row = table();
retryDirs = dir([resultsDir '_retry_*']);
[~, newestFirst] = sort([retryDirs.datenum], 'descend');
candidates = [string(fullfile({retryDirs(newestFirst).folder}, ...
    {retryDirs(newestFirst).name})), string(resultsDir)];
for candidate = candidates
    metadataFile = fullfile(candidate, 'run_metadata.json');
    manifestFile = fullfile(candidate, 'case_input_manifest.json');
    if ~isfile(metadataFile) || ~isfile(manifestFile); continue; end
    try
        metadata = jsondecode(fileread(metadataFile));
        manifest = jsondecode(fileread(manifestFile));
        complete = ~metadata.final.timedOut && ...
            metadata.final.exportedDissolutionSlices == design.TargetSlices && ...
            metadata.final.finalPorosity >= design.TargetFinalPorosity - 1e-6 && ...
            isfinite(metadata.final.finalPermeability_mD) && ...
            metadata.parameters.targetDissolutionFinalPorosity == design.TargetFinalPorosity && ...
            strcmp(manifest.geometryFile, char(design.GeometryFile));
        if complete
            row = resultRow(design, "SkippedComplete", "", candidate, ...
                metadata.final.exportedDissolutionSlices, metadata.final.finalPorosity, ...
                metadata.final.finalPermeability_mD, false, metadata.final.runtimeWallSeconds);
            return;
        end
    catch
        complete = false;
    end
end
end

function writeCaseManifest(resultsDir, design, reactiveCflSafety, ...
        phreeqcGridHmax, wallLimitSeconds, adaptiveGrowthFactor)
manifest = struct('sample', char(design.Sample), ...
    'geometryFile', char(design.GeometryFile), ...
    'geometryBytes', dir(design.GeometryFile).bytes, ...
    'targetSlices', design.TargetSlices, ...
    'targetFinalPorosity', design.TargetFinalPorosity, ...
    'inletVelocity_cm_s', design.InletVelocity_cm_s, ...
    'acidConcentration_mol_L', design.AcidConcentration_mol_L, ...
    'meshNumPartitionsX', design.MeshX, ...
    'meshNumPartitionsY', design.MeshY, ...
    'reactiveCflSafety', reactiveCflSafety, ...
    'phreeqcGridHmax_cm', phreeqcGridHmax, ...
    'maxExperimentWallSeconds', wallLimitSeconds, ...
    'adaptiveGrowthFactor', adaptiveGrowthFactor);
fid = fopen(fullfile(resultsDir, 'case_input_manifest.json'), 'w');
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
end

function row = resultRow(design, status, message, resultsDir, slices, porosity, permeability, timedOut, runtime)
row = design;
row.ExportedSlices = slices;
row.FinalPorosity = porosity;
row.FinalPermeability_mD = permeability;
row.TimedOut = timedOut;
row.RuntimeSeconds = runtime;
row.Status = string(status);
row.Message = string(message);
row.ResultsDir = string(resultsDir);
end

function writeSummary(summary, batchDir)
writetable(summary, fullfile(batchDir, 'batch_summary.csv'));
writetable(summary, fullfile(batchDir, 'batch_summary.xlsx'));
end
