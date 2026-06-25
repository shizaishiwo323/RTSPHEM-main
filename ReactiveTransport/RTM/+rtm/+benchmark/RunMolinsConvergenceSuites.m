function run = RunMolinsConvergenceSuites(options)
%RUNMOLINSCONVERGENCESUITES Run selected Molins convergence suites.

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isfield(options, 'outputRoot') || isempty(options.outputRoot)
    error('RTSPHEM:Benchmark:MissingOutputRoot', ...
        'options.outputRoot is required.');
end
outputRoot = char(options.outputRoot);
if exist(outputRoot, 'dir') ~= 7
    mkdir(outputRoot);
end

kinds = string(getOption(options, 'kinds', ...
    ["partI_strict"; "partII_strict"; "integration_phreeqc"]));
kinds = kinds(:);
suites = struct([]);
for iKind = 1:numel(kinds)
    suiteDir = fullfile(outputRoot, char(kinds(iKind)));
    if exist(suiteDir, 'dir') ~= 7
        mkdir(suiteDir);
    end
    kindOptions = optionsForKind(options, kinds(iKind));
    suiteOptions = rtm.benchmark.CreateMolinsDriverCaseOptions( ...
        kinds(iKind), optionsWithOutputDir(kindOptions, suiteDir));
    suite = rtm.benchmark.RunConvergenceSuite(suiteOptions);
    if isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix)
        suite.acceptance_matrix = options.acceptanceMatrix;
    end
    suites = appendStructWithUnifiedFields(suites, suite);
end

run = struct();
run.schema_version = "molins_convergence_suites_v1";
run.created_at = string(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
run.output_root = string(outputRoot);
run.kinds = kinds;
run.suites = suites;
if isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix)
    run.acceptance_matrix = options.acceptanceMatrix;
end
run.accepted = all(arrayfun(@(suite) suite.report.accepted, suites));
run.index_path = writeIndex(outputRoot, run);
end

function kindOptions = optionsForKind(options, kind)
kindOptions = options;
if ~(isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix))
    return;
end
if isfield(options, 'observableName') && ~isempty(options.observableName)
    return;
end
normalizedKind = lower(strrep(char(kind), '-', '_'));
switch normalizedKind
    case {'parti_strict', 'part_i_strict', 'part1_strict'}
        kindOptions.observableName = 'mean_effective_rate_mol_cm2_s';
        reference = referenceValueFromMatrix(options.acceptanceMatrix, ...
            {'partI_reference', 'mean_effective_rate_mol_cm2_s'}, NaN);
        if isfinite(reference)
            kindOptions.referenceTargetValue = reference;
            kindOptions.referenceRelativeTolerance = 0.05;
            kindOptions.referenceTargetRunNamePattern = 'dt_0p054';
        end
    case {'partii_strict', 'part_ii_strict', 'part2_strict'}
        kindOptions.observableName = 'final_solid_volume_cm3';
    case {'integration_phreeqc', 'molins_geometry_phreeqc'}
        kindOptions.observableName = 'max_component_mass_residual_moles';
        kindOptions.observableMaximumTolerance = 1e-8;
end
end

function value = referenceValueFromMatrix(matrix, path, defaultValue)
value = matrix;
for iPath = 1:numel(path)
    if ~isstruct(value) || ~isfield(value, path{iPath}) || isempty(value.(path{iPath}))
        value = defaultValue;
        return;
    end
    value = value.(path{iPath});
end
if ~(isscalar(value) && isfinite(value))
    value = defaultValue;
end
end

function updated = optionsWithOutputDir(options, suiteDir)
updated = options;
updated.outputDir = suiteDir;
if isfield(updated, 'outputRoot')
    updated = rmfield(updated, 'outputRoot');
end
end

function path = writeIndex(outputRoot, run)
path = string(fullfile(outputRoot, 'molins_convergence_index.json'));
index = makeJsonSafe(slimRunIndex(run));
index.index_path = path;
fid = fopen(path, 'w');
if fid == -1
    error('RTSPHEM:Benchmark:IndexOpenFailed', ...
        'Cannot write Molins convergence index: %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(index));
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
    'accepted', 'failure_message', 'time_step_s', 'grid_spacing_cm', ...
    'grid_resolution', 'initial_porosity', ...
    'initial_surface_area_cm2', 'initial_solid_volume_cm3', ...
    'initial_mineral_moles', 'final_mineral_moles', ...
    'mineral_dissolved_moles', 'reaction_realized_moles', ...
    'solid_volume_change_cm3', ...
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

function structs = appendStructWithUnifiedFields(structs, value)
if isempty(structs)
    structs = value;
    return;
end
allFields = union(fieldnames(structs), fieldnames(value));
structs = addMissingFields(structs, allFields);
value = addMissingFields(value, allFields);
structs(end + 1, 1) = value;
end

function structs = addMissingFields(structs, fields)
for iField = 1:numel(fields)
    fieldName = fields{iField};
    if ~isfield(structs, fieldName)
        [structs.(fieldName)] = deal([]);
    end
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
