function [summary, batchDir] = run_phreeqc_se_batch(stage, options)
% run_phreeqc_se_batch - Run three literature-guided dissolution regimes.
%
%   run_phreeqc_se_batch()              % 3 cases x 100 slices, full dissolution
%   run_phreeqc_se_batch('smoke')       % 3 cases x 10 target slices
%   run_phreeqc_se_batch('production')  % 3 cases x 100 slices, full dissolution
%   run_phreeqc_se_batch('resume')      % enter an existing production batch
%   design = run_phreeqc_se_batch('design')

if nargin < 1 || isempty(stage)
    stage = 'production';
end
if nargin < 2 || isempty(options)
    options = struct();
end
stage = lower(string(stage));
resumeExistingBatch = stage == "resume";
if ~ismember(stage, ["design", "smoke", "production", "resume"])
    error('RTSPHEM:SEBatch:UnknownStage', ...
        'stage must be design, smoke, production, or resume.');
end
if resumeExistingBatch
    stage = "production";
    if ~isfield(options, 'batchDir') || strlength(string(options.batchDir)) == 0
        options.batchDir = input('Existing batch result directory: ', 's');
    end
    options.batchDir = char(strip(string(options.batchDir), '"'));
    if ~isfolder(options.batchDir)
        error('RTSPHEM:SEBatch:MissingBatchDir', ...
            'Existing batch result directory not found: %s', options.batchDir);
    end
end

rtmDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(rtmDir));
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));

design = buildDesign();
if stage == "design"
    summary = design;
    batchDir = "";
    return;
end

if stage == "smoke"
    targetSlices = design.SmokeTargetSlices(1);
    meshPartitions = 64;
    meshX = meshPartitions;
    meshY = meshPartitions;
    targetFinalPorosity = 0.78;
    wallLimitSeconds = 2 * 3600;
else
    targetSlices = design.ProductionTargetSlices(1);
    meshPartitions = 200;
    meshX = 200;
    meshY = 100;
    targetFinalPorosity = 1.0;
    wallLimitSeconds = Inf;
end
earlyExportSteps = 6;
if isfield(options, 'targetSlices')
    targetSlices = options.targetSlices;
end
if isfield(options, 'meshPartitions')
    meshPartitions = options.meshPartitions;
end
if isfield(options, 'meshNumPartitionsX')
    meshX = options.meshNumPartitionsX;
end
if isfield(options, 'meshNumPartitionsY')
    meshY = options.meshNumPartitionsY;
end
if isfield(options, 'targetFinalPorosity')
    targetFinalPorosity = options.targetFinalPorosity;
end
if isfield(options, 'wallLimitSeconds')
    wallLimitSeconds = options.wallLimitSeconds;
end
validateattributes(targetSlices, {'numeric'}, {'scalar', 'integer', '>=', 2});
validateattributes(meshPartitions, {'numeric'}, {'scalar', 'integer', '>=', 16});
validateattributes(meshX, {'numeric'}, {'scalar', 'integer', '>=', 16});
validateattributes(meshY, {'numeric'}, {'scalar', 'integer', '>=', 16});
validateattributes(targetFinalPorosity, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', '<=', 1});

batchRoot = fullfile(projectRoot, 'outputs', 'rtm_batches');
if isfield(options, 'batchDir') && strlength(string(options.batchDir)) > 0
    batchDir = char(options.batchDir);
else
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    batchDir = fullfile(batchRoot, sprintf('se_modes_%s_%s', stage, stamp));
end
if ~exist(batchDir, 'dir')
    mkdir(batchDir);
end
writetable(design, fullfile(batchDir, 'batch_design.csv'));
writetable(design, fullfile(batchDir, 'batch_design.xlsx'));

summary = table();
for i = 1:height(design)
    caseName = sprintf('%s_%s', stage, design.Mode(i));
    resultsDir = fullfile(batchDir, caseName);
    [complete, existingSlices] = completedRun( ...
        resultsDir, targetSlices, earlyExportSteps);
    if complete
        row = resultRow(design(i, :), stage, targetSlices, meshPartitions, ...
            meshX, meshY, "SkippedComplete", "", resultsDir, ...
            existingSlices, NaN, NaN, 0);
        summary = [summary; row]; %#ok<AGROW>
        writeSummary(summary, batchDir);
        continue;
    elseif exist(resultsDir, 'dir')
        % ponytail: preserve partial outputs; retry only this unfinished case.
        retryStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        resultsDir = sprintf('%s_retry_%s', resultsDir, retryStamp);
    end

    cfg = struct();
    cfg.resultsDir = resultsDir;
    cfg.inletVelocity = design.InletVelocity_cm_s(i);
    cfg.acidConcentration_mol_L = design.AcidConcentration_mol_L(i);
    cfg.rateCoefficientTST = design.RateCoefficient_mol_dm2_s(i);
    cfg.phreeqcKineticsCorrectionFactor = design.PhreeqcRateCorrection(i);
    cfg.targetDissolutionSlices = targetSlices;
    cfg.targetDissolutionFinalPorosity = targetFinalPorosity;
    cfg.initialMacroscaleTimeStepSize = design.InitialTimeStep_s(i);
    cfg.numPartitionsMicroscale = meshY;
    cfg.meshNumPartitionsX = meshX;
    cfg.meshNumPartitionsY = meshY;
    cfg.phreeqcExportEvery = Inf;
    cfg.phreeqcEarlyExportSteps = earlyExportSteps;
    cfg.phreeqcExportSeismoelectric = true;
    cfg.phreeqcWriteSpeciesTable = true;
    cfg.phreeqcTransportMaxFactor = Inf;
    cfg.transportSolverBackend = 'cut_cell_fv';
    cfg.transportUpwindMode = 'exp';
    cfg.transportInletBoundaryMode = 'dirichlet';
    cfg.transportFlowProjectionMode = 'conservative_edge_average';
    cfg.transportAccumulationCoefficient = 1;
    cfg.maxExperimentWallSeconds = wallLimitSeconds;
    cfg.enableHardCflLimit = false;
    cfg.hardCflLimitInitialGrid = false;
    cfg.adaptiveGrowthFactor = 2;
    cfg.concentrationHardOvershootThreshold = 1.2;
    cfg.concentrationShrinkFactor = 0.25;
    cfg.saveRealtimePlot = false;
    cfg.showDebugFigures = false;
    cfg.loadExistingGeometry = true;
    cfg.geometryLoadFile = "C:\Users\imgw\Documents\Codex\震电模拟\Wang_SE_experiment\generated_geometries\M14_24_bottom_left_600x400\M14_24_bottom_left_600x400_literal.mat";
    if isfield(options, 'geometryLoadFile'); cfg.geometryLoadFile = options.geometryLoadFile; end
    if isfield(options, 'phreeqcGridMode')
        cfg.phreeqcGridMode = options.phreeqcGridMode;
    end
    if isfield(options, 'phreeqcGridHmax')
        cfg.phreeqcGridHmax = options.phreeqcGridHmax;
    end
    if isfield(options, 'dxfResolutionX')
        cfg.dxfResolutionX = options.dxfResolutionX;
    end
    if isfield(options, 'dxfResolutionY')
        cfg.dxfResolutionY = options.dxfResolutionY;
    end
    if isfield(options, 'isolateLargeLinearSolves')
        cfg.isolateLargeLinearSolves = options.isolateLargeLinearSolves;
    end
    if isfield(options, 'phreeqcWriteSpeciesTable')
        cfg.phreeqcWriteSpeciesTable = options.phreeqcWriteSpeciesTable;
    end

    fprintf('\n[%s %d/%d] %s: Pe=%.4g, DaII=%.4g, DaI=%.4g\n', ...
        upper(stage), i, height(design), design.Mode(i), ...
        design.Pe(i), design.DaII(i), design.DaI(i));
    timer = tic;
    try
        result = run_phreeqc_se(cfg);
        runtimeSeconds = toc(timer);
        exportedSlices = countChemistryRows(result.seismoelectricChemistryCsv);
        finalPorosity = result.finalPorosity;
        finalPermeability = result.finalPermeability;
        resultsDir = result.resultsDir;
        if result.timedOut
            status = "Failed";
            message = "Simulation timed out: " + string(result.timeoutMessage);
        elseif exportedSlices < targetSlices
            status = "Failed";
            message = sprintf("Only %d of %d target slices were exported.", ...
                exportedSlices, targetSlices);
        elseif finalPorosity < targetFinalPorosity - 1e-6
            status = "Failed";
            message = sprintf("Final porosity %.8f did not reach %.8f.", ...
                finalPorosity, targetFinalPorosity);
        else
            status = "Success";
            message = "";
        end
    catch ME
        runtimeSeconds = toc(timer);
        exportedSlices = countChemistryRows(fullfile(resultsDir, ...
            'phreeqc_results', 'seismoelectric_chemistry.csv'));
        status = "Failed";
        message = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
        finalPorosity = NaN;
        finalPermeability = NaN;
        warning('RTSPHEM:SEBatch:CaseFailed', '%s failed: %s', ...
            design.Mode(i), ME.message);
    end

    row = resultRow(design(i, :), stage, targetSlices, meshPartitions, ...
        meshX, meshY, status, message, resultsDir, exportedSlices, ...
        finalPorosity, finalPermeability, runtimeSeconds);
    summary = [summary; row]; %#ok<AGROW>
    writeSummary(summary, batchDir);
    close all hidden;
end

successful = summary.Status == "Success" | summary.Status == "SkippedComplete";
if any(~successful)
    error('RTSPHEM:SEBatch:Incomplete', ...
        '%d of %d cases failed. See batch_summary.csv.', ...
        sum(~successful), height(summary));
end
if stage == "production" && sum(summary.ExportedSlices) < 100
    error('RTSPHEM:SEBatch:TooFewSlices', ...
        'Production batch exported %d slices; at least 100 are required.', ...
        sum(summary.ExportedSlices));
end
end

function design = buildDesign()
% The solver reports DaII = c0*Vm*k*1000*L/D. Soulaine et al. (2017)
% classify the pore-scale regimes with DaI = DaII/Pe.
mode = ["face"; "uniform"; "ramified_wormhole"];
expectedPattern = ["compact_or_face"; "uniform"; "ramified_multichannel"];
Pe = [0.01; 10; 12];
DaII = [0.10; 0.10; 1];
acidMolL = [0.020; 0.005; 0.100];
diffusion = repmat(1e-5, 3, 1);
lengthCm = repmat(1e-3, 3, 1);
molarVolume = repmat(36.9, 3, 1);
databaseAcidRateMolDm2S = repmat(8.9e-5, 3, 1);

inletVelocity = Pe .* diffusion ./ lengthCm;
acidMolCm3 = acidMolL * 1e-3;
rateCoefficient = DaII .* diffusion ./ ...
    (acidMolCm3 .* molarVolume .* 1000 .* lengthCm);
rateCorrection = rateCoefficient ./ databaseAcidRateMolDm2S;
DaI = DaII ./ Pe;
reactionInitialStep = max(1e-5, min(0.1, 0.1 * 0.01 ./ DaII));
advectionInitialStep = 5e-3 ./ max(Pe, 1);
initialTimeStep = min(reactionInitialStep, advectionInitialStep);
smokeTargetSlices = repmat(10, 3, 1);
productionTargetSlices = repmat(100, 3, 1);

design = table(mode, expectedPattern, Pe, DaII, DaI, inletVelocity, ...
    acidMolL, acidMolCm3, rateCoefficient, rateCorrection, diffusion, ...
    lengthCm, molarVolume, initialTimeStep, smokeTargetSlices, ...
    productionTargetSlices, ...
    'VariableNames', {'Mode', 'ExpectedPattern', 'Pe', 'DaII', 'DaI', ...
    'InletVelocity_cm_s', 'AcidConcentration_mol_L', ...
    'AcidConcentration_mol_cm3', 'RateCoefficient_mol_dm2_s', ...
    'PhreeqcRateCorrection', 'DiffusionCoefficient_cm2_s', ...
    'CharacteristicLength_cm', 'MolarVolume_cm3_mol', ...
    'InitialTimeStep_s', 'SmokeTargetSlices', 'ProductionTargetSlices'});
end

function [complete, nSlices] = completedRun(resultsDir, targetSlices, earlyExportSteps)
chemistryCsv = fullfile(resultsDir, 'phreeqc_results', ...
    'seismoelectric_chemistry.csv');
metadataJson = fullfile(resultsDir, 'run_metadata.json');
nSlices = countChemistryRows(chemistryCsv);
complete = nSlices >= targetSlices && ...
    hasEarlyChemistryRows(chemistryCsv, earlyExportSteps) && ...
    exist(metadataJson, 'file') == 2;
if complete
    metadata = jsondecode(fileread(metadataJson));
    complete = isfield(metadata, 'final') && ...
        isfield(metadata.final, 'timedOut') && ~metadata.final.timedOut && ...
        isfield(metadata.final, 'finalPorosity') && ...
        metadata.final.finalPorosity >= 1 - 1e-6;
end
end

function tf = hasEarlyChemistryRows(csvPath, requiredSteps)
tf = false;
if exist(csvPath, 'file') ~= 2
    return;
end
try
    chemistry = readtable(csvPath, 'VariableNamingRule', 'preserve');
    tf = ismember('TimeStep', chemistry.Properties.VariableNames) && ...
        all(ismember((1:requiredSteps)', chemistry.TimeStep(:)));
catch
    tf = false;
end
end

function nRows = countChemistryRows(csvPath)
nRows = 0;
if exist(csvPath, 'file') ~= 2
    return;
end
try
    nRows = height(readtable(csvPath, 'VariableNamingRule', 'preserve'));
catch
    nRows = 0;
end
end

function row = resultRow(caseDesign, stage, targetSlices, meshPartitions, ...
    meshX, meshY, status, message, resultsDir, exportedSlices, ...
    finalPorosity, finalPermeability, runtimeSeconds)
row = caseDesign;
row.Stage = string(stage);
row.TargetSlices = targetSlices;
row.MeshPartitions = meshPartitions;
row.MeshPartitionsX = meshX;
row.MeshPartitionsY = meshY;
row.ExportedSlices = exportedSlices;
row.FinalPorosity = finalPorosity;
row.FinalPermeability_mD = finalPermeability;
row.RuntimeSeconds = runtimeSeconds;
row.Status = string(status);
row.Message = string(message);
row.ResultsDir = string(resultsDir);
row.MatlabRelease = string(version('-release'));
end

function writeSummary(summary, batchDir)
writetable(summary, fullfile(batchDir, 'batch_summary.csv'));
writetable(summary, fullfile(batchDir, 'batch_summary.xlsx'));
end
