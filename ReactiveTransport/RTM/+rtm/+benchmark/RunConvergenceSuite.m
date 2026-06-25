function suite = RunConvergenceSuite(options)
%RUNCONVERGENCESUITE Execute a refinement suite and build convergence evidence.

if nargin < 1 || isempty(options) || ~isstruct(options)
    options = struct();
end
if ~isfield(options, 'runFunction') || isempty(options.runFunction)
    error('RTSPHEM:Benchmark:MissingRunFunction', ...
        'options.runFunction is required.');
end
if ~isfield(options, 'observableFunction') || isempty(options.observableFunction)
    error('RTSPHEM:Benchmark:MissingObservableFunction', ...
        'options.observableFunction is required.');
end

suiteName = string(getOption(options, 'suiteName', 'convergence_suite'));
refinementScales = requireRefinementScales(options);
runNames = getRunNames(options, numel(refinementScales));
runs = repmat(emptyRunRecord(), numel(refinementScales), 1);
writePartialCheckpoint = logical(getOption(options, ...
    'writePartialCheckpoint', false));

for iRun = 1:numel(refinementScales)
    runInfo = struct('name', runNames(iRun), ...
        'index', iRun, ...
        'refinement_scale', refinementScales(iRun), ...
        'suite_name', suiteName);
    runs(iRun).name = char(runNames(iRun));
    runs(iRun).refinement_scale = refinementScales(iRun);
    try
        summary = options.runFunction(refinementScales(iRun), runInfo);
        runs(iRun).summary = summary;
        try
            observableValue = options.observableFunction(summary);
            runs(iRun).observable_value = validateObservableValue(observableValue);
            runs(iRun).accepted = logical(getFieldOrDefault(summary, 'accepted', true));
            runs(iRun).failure_message = "";
        catch observableErr
            runs(iRun).observable_value = NaN;
            runs(iRun).accepted = false;
            runs(iRun).failure_message = string(observableErr.message);
        end
    catch err
        runs(iRun).observable_value = NaN;
        runs(iRun).accepted = false;
        runs(iRun).summary = struct();
        runs(iRun).failure_message = string(err.message);
    end
    if writePartialCheckpoint && isfield(options, 'outputDir') && ...
            ~isempty(options.outputDir)
        writePartialRunCheckpoint(options.outputDir, suiteName, runs, iRun);
    end
end

reportRuns = runs;
for iRun = 1:numel(reportRuns)
    if ~isfinite(reportRuns(iRun).observable_value)
        reportRuns(iRun).observable_value = fallbackObservableValue(reportRuns, iRun);
    end
end
reportOptions = struct();
reportOptions.observableName = getOption(options, 'observableName', 'observable');
reportOptions.errorTolerance = getOption(options, 'errorTolerance', Inf);
reportOptions.minObservedOrder = getOption(options, 'minObservedOrder', -Inf);
reportOptions.referenceTargetValue = getOption(options, ...
    'referenceTargetValue', NaN);
reportOptions.referenceRelativeTolerance = getOption(options, ...
    'referenceRelativeTolerance', Inf);
reportOptions.referenceTargetRunNamePattern = getOption(options, ...
    'referenceTargetRunNamePattern', "");
reportOptions.observableMaximumTolerance = getOption(options, ...
    'observableMaximumTolerance', Inf);
report = rtm.benchmark.BuildConvergenceReport(reportRuns, reportOptions);
if isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix)
    report.acceptance_matrix = options.acceptanceMatrix;
end

suite = struct();
suite.suite_name = suiteName;
suite.runs = runs;
suite.report = report;
if isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix)
    suite.acceptance_matrix = options.acceptanceMatrix;
end
suite.report_path = "";
if isfield(options, 'outputDir') && ~isempty(options.outputDir)
    suite.report_path = rtm.benchmark.WriteConvergenceReport(options.outputDir, report);
end
end

function run = emptyRunRecord()
run = struct();
run.name = '';
run.refinement_scale = NaN;
run.observable_value = NaN;
run.accepted = false;
run.summary = struct();
run.failure_message = "";
end

function scales = requireRefinementScales(options)
if ~isfield(options, 'refinementScales') || isempty(options.refinementScales)
    error('RTSPHEM:Benchmark:MissingRefinementScales', ...
        'options.refinementScales is required.');
end
scales = options.refinementScales(:);
if numel(scales) < 2 || any(~isfinite(scales)) || any(scales <= 0)
    error('RTSPHEM:Benchmark:InvalidRefinementScale', ...
        'options.refinementScales must contain at least two positive finite values.');
end
end

function names = getRunNames(options, numRuns)
if isfield(options, 'runNames') && ~isempty(options.runNames)
    names = string(options.runNames(:));
    if numel(names) ~= numRuns
        error('RTSPHEM:Benchmark:InvalidRunNames', ...
            'options.runNames must match refinementScales length.');
    end
else
    names = strings(numRuns, 1);
    for iRun = 1:numRuns
        names(iRun) = "run" + iRun;
    end
end
end

function value = validateObservableValue(value)
if ~(isscalar(value) && isfinite(value))
    error('RTSPHEM:Benchmark:InvalidObservableValue', ...
        'Observable function must return a finite scalar.');
end
end

function value = fallbackObservableValue(runs, index)
finiteValues = [runs.observable_value];
finiteValues = finiteValues(isfinite(finiteValues));
if isempty(finiteValues)
    value = 0;
elseif index > 1 && isfinite(runs(index - 1).observable_value)
    value = runs(index - 1).observable_value;
else
    value = finiteValues(end);
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function path = writePartialRunCheckpoint(outputDir, suiteName, runs, completedRunCount)
path = fullfile(char(outputDir), 'benchmark_convergence_partial.json');
payload = struct();
payload.schema_version = "benchmark_convergence_partial_v1";
payload.suite_name = suiteName;
payload.created_at = string(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
payload.completed_run_count = completedRunCount;
payload.run_count = numel(runs);
payload.runs = makeCheckpointRuns(runs, completedRunCount);
fid = fopen(path, 'w');
if fid == -1
    error('RTSPHEM:Benchmark:PartialCheckpointOpenFailed', ...
        'Cannot write benchmark partial checkpoint: %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(makeJsonSafe(payload)));
clear cleanup;
end

function records = makeCheckpointRuns(runs, completedRunCount)
records = repmat(emptyRunRecord(), completedRunCount, 1);
for iRun = 1:completedRunCount
    records(iRun) = runs(iRun);
    records(iRun).summary = slimSummary(records(iRun).summary);
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
    'benchmark_mesh', 'observation_times_s', 'observation_records'};
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
