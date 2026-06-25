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

names = strings(numel(runs), 1);
scales = zeros(numel(runs), 1);
values = zeros(numel(runs), 1);
acceptedRuns = false(numel(runs), 1);
for iRun = 1:numel(runs)
    names(iRun) = string(getFieldOrDefault(runs(iRun), 'name', "run" + iRun));
    scales(iRun) = scalarRunField(runs(iRun), 'refinement_scale');
    values(iRun) = scalarRunField(runs(iRun), 'observable_value');
    acceptedRuns(iRun) = logical(getFieldOrDefault(runs(iRun), 'accepted', true));
end

if any(~isfinite(scales)) || any(scales <= 0) || any(diff(scales) >= 0)
    error('RTSPHEM:Benchmark:InvalidRefinementScale', ...
        'runs.refinement_scale must be strictly decreasing positive values.');
end

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
orders = observedOrders(scales, values);
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
report.refinement_scales = scales;
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

function orders = observedOrders(scales, values)
if numel(values) < 3
    orders = zeros(0, 1);
    return;
end
orders = zeros(numel(values) - 2, 1);
for iOrder = 1:numel(orders)
    numerator = max(abs(values(iOrder) - values(iOrder + 1)), eps);
    denominator = max(abs(values(iOrder + 1) - values(iOrder + 2)), eps);
    refinementRatio = scales(iOrder) ./ scales(iOrder + 1);
    orders(iOrder) = log(numerator ./ denominator) ./ log(refinementRatio);
end
end

function value = scalarRunField(run, fieldName)
value = run.(fieldName);
if ~(isscalar(value) && isfinite(value))
    error('RTSPHEM:Benchmark:InvalidConvergenceRuns', ...
        'runs.%s must be a finite scalar.', fieldName);
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
