function run = run_molins_convergence_suites(mode, options)
%RUN_MOLINS_CONVERGENCE_SUITES User-facing Molins convergence runner.
%
% Modes:
%   diagnostic  - quick PHREEQC integration smoke.
%   integration - external TST + PHREEQC integration suite.
%   full        - PHREEQC integration suite.
%   full_acceptance - plan acceptance matrix for PHREEQC integration with
%                     54, 5.4, 0.54, 0.054 s steps and Molins grid labels
%                     recorded in the result index.

if nargin < 1 || isempty(mode)
    mode = 'diagnostic';
end
if isstruct(mode)
    options = mode;
    mode = 'diagnostic';
elseif nargin < 2 || isempty(options)
    options = struct();
end
if nargin < 2 || isempty(options)
    options = struct();
end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

[modeName, defaults] = modeDefaults(mode);
merged = mergeOptions(defaults, options);
if ~isfield(merged, 'outputRoot') || isempty(merged.outputRoot)
    merged.outputRoot = defaultOutputRoot(scriptDir, modeName);
end

run = rtm.benchmark.RunMolinsConvergenceSuites(merged);
run.mode = string(modeName);
run.runner = string(mfilename('fullpath'));
run.index_path = rewriteIndexWithMode(run);
end

function [modeName, defaults] = modeDefaults(mode)
modeName = lower(strrep(strtrim(char(mode)), '-', '_'));
defaults = struct();
defaults.errorTolerance = Inf;
defaults.minObservedOrder = -Inf;

switch modeName
    case {'diagnostic', 'smoke'}
        modeName = 'diagnostic';
        defaults.kinds = "integration_phreeqc";
        defaults.refinementScales = [0.5; 0.25];
        defaults.totalTime_s = 0.5;
        defaults.maxDisplacementOverH = 1;
    case {'integration', 'phreeqc', 'integration_phreeqc'}
        modeName = 'integration';
        defaults.kinds = "integration_phreeqc";
        defaults.refinementScales = [5.4; 0.54; 0.054];
        defaults.totalTime_s = 5.4;
    case {'full', 'all'}
        modeName = 'full';
        defaults.kinds = "integration_phreeqc";
        defaults.refinementScales = [5.4; 0.54; 0.054];
        defaults.totalTime_s = 5.4;
    case {'full_acceptance', 'acceptance', 'plan_acceptance'}
        modeName = 'full_acceptance';
        defaults.kinds = "integration_phreeqc";
        defaults.refinementScales = [54; 5.4; 0.54; 0.054];
        defaults.totalTime_s = 54;
        defaults.acceptanceMatrix = molinsAcceptanceMatrix();
        defaults.useAcceptanceGrid = true;
        defaults.writePartialCheckpoint = true;
    otherwise
        error('RTSPHEM:Benchmark:UnknownMolinsConvergenceMode', ...
            'Unknown Molins convergence mode: %s.', char(mode));
end
end

function matrix = molinsAcceptanceMatrix()
matrix = struct();
matrix.schema_version = "molins_plan_acceptance_matrix_v1";
matrix.time_steps_s = [54; 5.4; 0.54; 0.054];
matrix.grid_resolutions = ["128x64"; "256x128"; "512x256"];
matrix.partI_reference.initial_surface_area_cm2 = 0.0628;
matrix.partI_reference.initial_solid_volume_cm3 = 3.14e-4;
matrix.partI_reference.mean_effective_rate_mol_cm2_s = 4.32e-8;
matrix.partII_observation_times_min = [15; 30; 45];
matrix.completion_requires_real_molins_geometry = true;
matrix.completion_requires_phreeqc_integration_benchmark = true;
end

function merged = mergeOptions(defaults, options)
merged = defaults;
if ~isstruct(options)
    error('RTSPHEM:Benchmark:InvalidMolinsConvergenceOptions', ...
        'options must be a struct.');
end
fields = fieldnames(options);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    merged.(fieldName) = options.(fieldName);
end
end

function outputRoot = defaultOutputRoot(scriptDir, modeName)
stamp = datestr(now, 'yyyymmdd_HHMMSS');
outputRoot = fullfile(scriptDir, 'outputs', ...
    ['molins_convergence_' char(modeName) '_' stamp]);
end

function path = rewriteIndexWithMode(run)
path = run.index_path;
if strlength(path) == 0
    return;
end
fid = fopen(path, 'w');
if fid == -1
    error('RTSPHEM:Benchmark:IndexOpenFailed', ...
        'Cannot rewrite Molins convergence index: %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(makeJsonSafe(slimRunIndex(run))));
clear cleanup;
end

function index = slimRunIndex(run)
index = run;
if ~isfield(index, 'suites')
    return;
end
for iSuite = 1:numel(index.suites)
    if isfield(index.suites(iSuite), 'runs')
        index.suites(iSuite).runs = slimRuns(index.suites(iSuite).runs);
    end
end
end

function runs = slimRuns(runs)
for iRun = 1:numel(runs)
    if isfield(runs(iRun), 'summary')
        runs(iRun).summary = slimSummary(runs(iRun).summary);
    end
end
end

function summary = slimSummary(summary)
if ~isstruct(summary) || isempty(summary)
    return;
end
keepFields = {'time_s', 'accepted_steps', 'rejected_steps', 'run_name', ...
    'accepted', 'failure_message', 'initial_porosity', ...
    'initial_surface_area_cm2', 'initial_solid_volume_cm3', ...
    'initial_mineral_moles', 'final_mineral_moles', ...
    'mineral_dissolved_moles', 'solid_volume_change_cm3', ...
    'final_solid_volume_cm3', 'final_surface_area_cm2', ...
    'final_porosity', 'mean_effective_rate_mol_cm2_s', ...
    'max_component_mass_residual_moles', 'max_displacement_over_h', ...
    'benchmark_mesh', 'acceptance_matrix'};
fieldNames = fieldnames(summary);
for iField = 1:numel(fieldNames)
    if ~ismember(fieldNames{iField}, keepFields)
        summary = rmfield(summary, fieldNames{iField});
    end
end
end

function value = makeJsonSafe(value)
if isa(value, 'function_handle')
    value = string(sprintf('<function_handle:%s>', func2str(value)));
elseif isstruct(value)
    fieldNames = fieldnames(value);
    for iValue = 1:numel(value)
        for iField = 1:numel(fieldNames)
            fieldName = fieldNames{iField};
            value(iValue).(fieldName) = makeJsonSafe(value(iValue).(fieldName));
        end
    end
elseif iscell(value)
    for iCell = 1:numel(value)
        value{iCell} = makeJsonSafe(value{iCell});
    end
end
end
