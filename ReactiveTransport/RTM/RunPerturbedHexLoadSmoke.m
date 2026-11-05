function summary = RunPerturbedHexLoadSmoke(geometryDir, resultsRoot)
% Load each fixed-domain geometry through the production PHREEQC entrypoint.
if nargin < 2 || isempty(resultsRoot)
    resultsRoot = fullfile(fileparts(fileparts(fileparts(geometryDir))), ...
        'RTSPHEM-main', 'outputs', 'rtm_runs', 'perturbed_hex_phi_only_load_smoke');
end
rtmDir = fileparts(mfilename('fullpath')); addpath(rtmDir);
files = dir(fullfile(geometryDir, '*_perturbed_hex_600x400um.mat'));
assert(numel(files) == 4, 'RTM:GeometrySmoke:FileCount', 'Expected four geometry MAT files.');
rows = cell(numel(files), 6);
for i = 1:numel(files)
    geometryFile = fullfile(files(i).folder, files(i).name);
    sampleName = erase(string(files(i).name), '_perturbed_hex_600x400um.mat');
    existingLog = fullfile(resultsRoot, char(sampleName), 'global_evolution_log.csv');
    existingMetadata = fullfile(resultsRoot, char(sampleName), 'run_metadata.json');
    if exist(existingLog, 'file') == 2 && exist(existingMetadata, 'file') == 2
        metadata = jsondecode(fileread(existingMetadata));
        if isfield(metadata, 'final') && isfield(metadata.final, 'timedOut') && ...
                ~metadata.final.timedOut
            evolution = readtable(existingLog);
            last = evolution(end, :);
            rows(i,:) = {sampleName, last.porosity, last.permeability_mD, ...
                last.permeability_mD, false, string(fullfile(resultsRoot, char(sampleName)))};
            continue;
        end
    end
    cfg = struct('layoutType','random','loadExistingGeometry',true, ...
        'geometryLoadFile',geometryFile,'resultsDir',fullfile(resultsRoot,char(sampleName)), ...
        'inletVelocity',1e-4,'acidConcentration_mol_L',0.02, ...
        'rateCoefficientTST',0.00135501355013550, ...
        'phreeqcKineticsCorrectionFactor',0.00135501355013550/8.9e-5, ...
        'numPartitionsMicroscale',128,'meshNumPartitionsX',128,'meshNumPartitionsY',86, ...
        'targetDissolutionSlices',[],'endTime',1e-6,'maximalStep',1e-6, ...
        'initialMacroscaleTimeStepSize',1e-6,'adaptiveMaxTimeStep',1e-6, ...
        'maxTotalTimeSteps',2,'phreeqcExportEvery',Inf,'saveRealtimePlot',false);
    result = run_phreeqc_se(cfg);
    rows(i,:) = {sampleName, result.finalPorosity, result.initialPermeability, ...
        result.finalPermeability, result.timedOut, string(result.resultsDir)};
end
summary = cell2table(rows, 'VariableNames', {'sampleName','porosity','initialKxx_mD', ...
    'finalKxx_mD','timedOut','resultsDir'});
writetable(summary, fullfile(resultsRoot, 'load_smoke_summary.csv'));
writetable(summary, fullfile(resultsRoot, 'load_smoke_summary.xlsx'));
assert(~any(summary.timedOut), 'RTM:GeometrySmoke:Timeout', 'At least one geometry smoke run timed out.');
end
