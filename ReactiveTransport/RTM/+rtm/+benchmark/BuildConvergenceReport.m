function report = BuildConvergenceReport(runs, options)
%BUILDCONVERGENCEREPORT Summarize grid/time-step benchmark convergence.

if nargin < 2 || isempty(options)
    options = struct();
end
validateRuns(runs);

observableName = string(getOption(options, 'observableName', 'observable'));
errorTolerance = getScalarOption(options, 'errorTolerance', Inf, false);
minObservedOrder = getScalarOption(options, 'minObservedOrder', -Inf, false);
referenceTargetValue = getOptionalFiniteScalar(options, ...
    'referenceTargetValue', NaN);
referenceRelativeTolerance = getScalarOption(options, ...
    'referenceRelativeTolerance', Inf, false);
referenceTargetRunNamePattern = string(getOption(options, ...
    'referenceTargetRunNamePattern', ""));
observableMaximumTolerance = getScalarOption(options, ...
    'observableMaximumTolerance', Inf, false);
[refinementDimension, refinementField] = resolveRefinementDimension(options);

names = strings(numel(runs), 1);
scales = zeros(numel(runs), 1);
legacyScales = zeros(numel(runs), 1);
timeSteps = nan(numel(runs), 1);
gridSpacings = nan(numel(runs), 1);
gridResolutions = strings(numel(runs), 1);
values = zeros(numel(runs), 1);
acceptedRuns = false(numel(runs), 1);
for iRun = 1:numel(runs)
    names(iRun) = string(getFieldOrDefault(runs(iRun), 'name', "run" + iRun));
    legacyScales(iRun) = scalarRunField(runs(iRun), 'refinement_scale');
    scales(iRun) = scalarRunField(runs(iRun), refinementField);
    timeSteps(iRun) = optionalScalarRunField(runs(iRun), 'time_step_s');
    gridSpacings(iRun) = optionalScalarRunField(runs(iRun), 'grid_spacing_cm');
    gridResolutions(iRun) = optionalStringRunField(runs(iRun), ...
        'grid_resolution');
    values(iRun) = scalarRunField(runs(iRun), 'observable_value');
    acceptedRuns(iRun) = logical(getFieldOrDefault(runs(iRun), 'accepted', true));
end

refinementGroups = refinementGroupLabels(refinementDimension, gridResolutions, ...
    timeSteps);
validateRefinementScales(scales, refinementGroups, refinementField);

referenceValue = values(end);
errorsToReference = abs(values - referenceValue);
relativeErrorsToTarget = nan(size(values));
maxRelativeErrorToTarget = NaN;
referenceTargetApplies = true(size(values));
hasReferenceTarget = isfinite(referenceTargetValue);
if hasReferenceTarget
    if strlength(referenceTargetRunNamePattern) > 0
        referenceTargetApplies = contains(names, referenceTargetRunNamePattern);
    end
    if ~any(referenceTargetApplies)
        error('RTSPHEM:Benchmark:InvalidConvergenceOption', ...
            'referenceTargetRunNamePattern did not match any run names.');
    end
    denominator = max(abs(referenceTargetValue), eps);
    relativeErrorsToTarget(referenceTargetApplies) = ...
        abs(values(referenceTargetApplies) - referenceTargetValue) ./ denominator;
    maxRelativeErrorToTarget = max( ...
        relativeErrorsToTarget(referenceTargetApplies), [], 'omitnan');
else
    referenceTargetApplies(:) = false;
end
orders = observedOrders(scales, values, refinementGroups);
failureReasons = strings(0, 1);
if any(~acceptedRuns)
    failureReasons(end + 1, 1) = "benchmark run rejected";
end
if max(errorsToReference(1:end-1), [], 'omitnan') > errorTolerance
    failureReasons(end + 1, 1) = "error exceeds tolerance";
end
if ~isempty(orders) && min(orders, [], 'omitnan') < minObservedOrder
    failureReasons(end + 1, 1) = "observed order below tolerance";
end
if hasReferenceTarget && maxRelativeErrorToTarget > referenceRelativeTolerance
    failureReasons(end + 1, 1) = ...
        "reference target error exceeds tolerance";
end
maxObservableValue = max(values, [], 'omitnan');
if maxObservableValue > observableMaximumTolerance
    failureReasons(end + 1, 1) = ...
        "observable exceeds maximum tolerance";
end

report = struct();
report.observable_name = observableName;
report.run_names = names;
report.refinement_dimension = refinementDimension;
report.refinement_field = string(refinementField);
report.refinement_scales = scales;
report.legacy_refinement_scales = legacyScales;
report.time_steps_s = timeSteps;
report.grid_spacings_cm = gridSpacings;
report.grid_resolutions = gridResolutions;
report.refinement_groups = refinementGroups;
report.refinement_matrix = buildRefinementMatrix(names, gridResolutions, ...
    gridSpacings, timeSteps, values, acceptedRuns);
report.joint_acceptance_matrix = report.refinement_matrix;
report.time_convergence_report = buildTimeConvergenceReports( ...
    report.refinement_matrix);
report.grid_convergence_report = buildGridConvergenceReports( ...
    report.refinement_matrix);
report.observable_values = values;
report.reference_value = referenceValue;
report.errors_to_reference = errorsToReference;
report.reference_target_value = referenceTargetValue;
report.reference_relative_tolerance = referenceRelativeTolerance;
report.relative_errors_to_target = relativeErrorsToTarget;
report.reference_target_run_name_pattern = referenceTargetRunNamePattern;
report.reference_target_applies = referenceTargetApplies;
report.max_relative_error_to_target = maxRelativeErrorToTarget;
report.observable_maximum_tolerance = observableMaximumTolerance;
report.max_observable_value = maxObservableValue;
report.max_error_to_reference = max(errorsToReference(1:end-1), [], 'omitnan');
if isempty(report.max_error_to_reference)
    report.max_error_to_reference = 0;
end
report.observed_orders = orders;
report.error_tolerance = errorTolerance;
report.min_observed_order = minObservedOrder;
report.accepted_runs = acceptedRuns;
report.accepted = isempty(failureReasons);
report.failure_reasons = failureReasons;
end

function reports = buildTimeConvergenceReports(matrix)
numGrids = numel(matrix.grid_resolutions);
reports = repmat(struct('grid_resolution', "", ...
    'grid_spacing_cm', NaN, ...
    'time_steps_s', [], ...
    'observable_values', [], ...
    'accepted_runs', [], ...
    'run_names', strings(1, 0)), numGrids, 1);
for iGrid = 1:numGrids
    reports(iGrid).grid_resolution = matrix.grid_resolutions(iGrid);
    reports(iGrid).grid_spacing_cm = matrix.grid_spacings_cm(iGrid);
    reports(iGrid).time_steps_s = matrix.time_steps_s;
    reports(iGrid).observable_values = matrix.observable_values(iGrid, :);
    reports(iGrid).accepted_runs = matrix.accepted_runs(iGrid, :);
    reports(iGrid).run_names = matrix.run_names(iGrid, :);
end
end

function reports = buildGridConvergenceReports(matrix)
numTimes = numel(matrix.time_steps_s);
reports = repmat(struct('time_step_s', NaN, ...
    'grid_resolutions', strings(0, 1), ...
    'grid_spacings_cm', [], ...
    'observable_values', [], ...
    'accepted_runs', [], ...
    'run_names', strings(0, 1)), numTimes, 1);
for iTime = 1:numTimes
    reports(iTime).time_step_s = matrix.time_steps_s(iTime);
    reports(iTime).grid_resolutions = matrix.grid_resolutions;
    reports(iTime).grid_spacings_cm = matrix.grid_spacings_cm;
    reports(iTime).observable_values = matrix.observable_values(:, iTime);
    reports(iTime).accepted_runs = matrix.accepted_runs(:, iTime);
    reports(iTime).run_names = matrix.run_names(:, iTime);
end
end

function matrix = buildRefinementMatrix(names, gridResolutions, gridSpacings, ...
        timeSteps, values, acceptedRuns)
gridLabels = gridResolutions(:);
for iRun = 1:numel(gridLabels)
    if strlength(gridLabels(iRun)) == 0
        if isfinite(gridSpacings(iRun))
            gridLabels(iRun) = "h_" + string(gridSpacings(iRun));
        else
            gridLabels(iRun) = "grid_" + string(iRun);
        end
    end
end

[uniqueGridLabels, gridIndex] = uniqueStable(gridLabels);
uniqueGridSpacings = nan(numel(uniqueGridLabels), 1);
for iGrid = 1:numel(uniqueGridLabels)
    matching = find(gridIndex == iGrid & isfinite(gridSpacings), 1);
    if ~isempty(matching)
        uniqueGridSpacings(iGrid) = gridSpacings(matching);
    end
end

finiteTimes = timeSteps(isfinite(timeSteps));
if isempty(finiteTimes)
    uniqueTimeSteps = NaN;
else
    uniqueTimeSteps = unique(finiteTimes(:).', 'stable');
end
observableMatrix = nan(numel(uniqueGridLabels), numel(uniqueTimeSteps));
acceptedMatrix = false(numel(uniqueGridLabels), numel(uniqueTimeSteps));
runNameMatrix = strings(numel(uniqueGridLabels), numel(uniqueTimeSteps));

for iRun = 1:numel(values)
    iGrid = gridIndex(iRun);
    if isfinite(timeSteps(iRun))
        iTime = find(uniqueTimeSteps == timeSteps(iRun), 1);
    else
        iTime = 1;
    end
    if isempty(iTime)
        continue;
    end
    observableMatrix(iGrid, iTime) = values(iRun);
    acceptedMatrix(iGrid, iTime) = acceptedRuns(iRun);
    runNameMatrix(iGrid, iTime) = names(iRun);
end

matrix = struct();
matrix.grid_resolutions = uniqueGridLabels;
matrix.grid_spacings_cm = uniqueGridSpacings;
matrix.time_steps_s = uniqueTimeSteps;
matrix.observable_values = observableMatrix;
matrix.accepted_runs = acceptedMatrix;
matrix.run_names = runNameMatrix;
end

function [uniqueValues, index] = uniqueStable(values)
uniqueValues = strings(0, 1);
index = zeros(numel(values), 1);
for iValue = 1:numel(values)
    existing = find(uniqueValues == values(iValue), 1);
    if isempty(existing)
        uniqueValues(end + 1, 1) = values(iValue); %#ok<AGROW>
        existing = numel(uniqueValues);
    end
    index(iValue) = existing;
end
end

function [dimension, fieldName] = resolveRefinementDimension(options)
dimension = lower(string(getOption(options, 'refinementDimension', ...
    'refinement_scale')));
switch dimension
    case {"refinement_scale", "legacy", "scale"}
        dimension = "refinement_scale";
        fieldName = 'refinement_scale';
    case {"time", "time_step", "time_step_s", "dt"}
        dimension = "time";
        fieldName = 'time_step_s';
    case {"grid", "grid_spacing", "grid_spacing_cm", "space", "spatial"}
        dimension = "grid";
        fieldName = 'grid_spacing_cm';
    otherwise
        error('RTSPHEM:Benchmark:InvalidConvergenceOption', ...
            'Unsupported refinementDimension: %s.', char(dimension));
end
end

function groups = refinementGroupLabels(refinementDimension, gridResolutions, ...
        timeSteps)
groups = strings(size(gridResolutions));
switch refinementDimension
    case "time"
        if any(strlength(gridResolutions) > 0)
            groups = gridResolutions;
        end
    case "grid"
        if any(isfinite(timeSteps))
            groups = "dt_" + string(timeSteps);
        end
end
end

function validateRefinementScales(scales, groups, refinementField)
if any(~isfinite(scales)) || any(scales <= 0)
    error('RTSPHEM:Benchmark:InvalidRefinementScale', ...
        'runs.%s must contain positive finite values.', refinementField);
end
uniqueGroups = unique(groups, 'stable');
for iGroup = 1:numel(uniqueGroups)
    mask = groups == uniqueGroups(iGroup);
    groupScales = scales(mask);
    if numel(groupScales) > 1 && any(diff(groupScales) >= 0)
        error('RTSPHEM:Benchmark:InvalidRefinementScale', ...
            'runs.%s must be strictly decreasing positive values within each refinement group.', ...
            refinementField);
    end
end
end

function validateRuns(runs)
if ~isstruct(runs) || numel(runs) < 2
    error('RTSPHEM:Benchmark:InvalidConvergenceRuns', ...
        'runs must be a struct array with at least two entries.');
end
requiredFields = {'refinement_scale', 'observable_value'};
for iField = 1:numel(requiredFields)
    if ~isfield(runs, requiredFields{iField})
        error('RTSPHEM:Benchmark:InvalidConvergenceRuns', ...
            'runs.%s is required.', requiredFields{iField});
    end
end
end

function orders = observedOrders(scales, values, groups)
orders = zeros(0, 1);
uniqueGroups = unique(groups, 'stable');
for iGroup = 1:numel(uniqueGroups)
    mask = find(groups == uniqueGroups(iGroup));
    if numel(mask) < 3
        continue;
    end
    groupOrders = zeros(numel(mask) - 2, 1);
    for iOrder = 1:numel(groupOrders)
        i0 = mask(iOrder);
        i1 = mask(iOrder + 1);
        i2 = mask(iOrder + 2);
        numerator = max(abs(values(i0) - values(i1)), eps);
        denominator = max(abs(values(i1) - values(i2)), eps);
        refinementRatio = scales(i0) ./ scales(i1);
        groupOrders(iOrder) = log(numerator ./ denominator) ./ log(refinementRatio);
    end
    orders = [orders; groupOrders]; %#ok<AGROW>
end
end

function value = scalarRunField(run, fieldName)
if ~isfield(run, fieldName)
    error('RTSPHEM:Benchmark:InvalidConvergenceRuns', ...
        'runs.%s is required.', fieldName);
end
value = run.(fieldName);
if ~(isscalar(value) && isfinite(value))
    error('RTSPHEM:Benchmark:InvalidConvergenceRuns', ...
        'runs.%s must be a finite scalar.', fieldName);
end
end

function value = optionalScalarRunField(run, fieldName)
value = NaN;
if ~isfield(run, fieldName) || isempty(run.(fieldName))
    return;
end
candidate = run.(fieldName);
if ~(isscalar(candidate) && isfinite(candidate))
    error('RTSPHEM:Benchmark:InvalidConvergenceRuns', ...
        'runs.%s must be a finite scalar when present.', fieldName);
end
value = candidate;
end

function value = optionalStringRunField(run, fieldName)
value = "";
if ~isfield(run, fieldName) || isempty(run.(fieldName))
    return;
end
value = string(run.(fieldName));
if ~isscalar(value)
    error('RTSPHEM:Benchmark:InvalidConvergenceRuns', ...
        'runs.%s must be scalar text when present.', fieldName);
end
end

function value = getScalarOption(options, fieldName, defaultValue, requirePositive)
value = getOption(options, fieldName, defaultValue);
if ~(isscalar(value) && (isfinite(value) || isinf(value)))
    error('RTSPHEM:Benchmark:InvalidConvergenceOption', ...
        'options.%s must be a finite scalar.', fieldName);
end
if requirePositive && value <= 0
    error('RTSPHEM:Benchmark:InvalidConvergenceOption', ...
        'options.%s must be positive.', fieldName);
end
end

function value = getOptionalFiniteScalar(options, fieldName, defaultValue)
value = getOption(options, fieldName, defaultValue);
if isnan(value)
    return;
end
if ~(isscalar(value) && isfinite(value))
    error('RTSPHEM:Benchmark:InvalidConvergenceOption', ...
        'options.%s must be a finite scalar.', fieldName);
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
