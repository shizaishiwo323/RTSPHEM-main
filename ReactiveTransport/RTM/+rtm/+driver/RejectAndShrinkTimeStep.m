function updated = RejectAndShrinkTimeStep(solverState, reason, options)
%REJECTANDSHRINKTIMESTEP Record a rejected step and shrink dt for retry.

if nargin < 3 || isempty(options)
    options = struct();
end
if nargin < 2 || isempty(reason)
    reason = "";
end

dtSeconds = requirePositiveScalar(solverState, 'dt_s');
shrinkFactor = getOption(options, 'shrinkFactor', 0.5);
minDtSeconds = getOption(options, 'minDt_s', 1e-8);
maxRetries = getOption(options, 'maxRetries', 12);
if ~(isscalar(shrinkFactor) && isfinite(shrinkFactor) && shrinkFactor > 0 && shrinkFactor < 1)
    error('RTSPHEM:Driver:InvalidRetryOption', ...
        'shrinkFactor must be between 0 and 1.');
end
if ~(isscalar(minDtSeconds) && isfinite(minDtSeconds) && minDtSeconds > 0)
    error('RTSPHEM:Driver:InvalidRetryOption', ...
        'minDt_s must be a positive finite scalar.');
end
if ~(isscalar(maxRetries) && isfinite(maxRetries) && maxRetries >= 0 && maxRetries == round(maxRetries))
    error('RTSPHEM:Driver:InvalidRetryOption', ...
        'maxRetries must be a nonnegative integer scalar.');
end

oldRejectedSteps = getFieldOrDefault(solverState, 'rejected_steps', 0);
newDtSeconds = dtSeconds .* shrinkFactor;
updated = solverState;
updated.dt_s = newDtSeconds;
updated.rejected_steps = oldRejectedSteps + 1;
updated.abort = false;
updated.abort_reason = "";

entry = struct();
entry.index = updated.rejected_steps;
entry.reason = string(reason);
entry.old_dt_s = dtSeconds;
entry.new_dt_s = newDtSeconds;
entry.time_s = getFieldOrDefault(solverState, 'time_s', NaN);
updated.rejection_log = appendLog(getFieldOrDefault(solverState, 'rejection_log', []), entry);

if oldRejectedSteps >= maxRetries
    updated.abort = true;
    updated.abort_reason = "maximum retries exceeded";
elseif newDtSeconds < minDtSeconds
    updated.abort = true;
    updated.abort_reason = "minimum dt reached";
end
end

function value = requirePositiveScalar(structValue, fieldName)
if ~isfield(structValue, fieldName)
    error('RTSPHEM:Driver:MissingTimeStep', ...
        'solverState.%s is required.', fieldName);
end
value = structValue.(fieldName);
if ~(isscalar(value) && isfinite(value) && value > 0)
    error('RTSPHEM:Driver:InvalidTimeStep', ...
        'solverState.%s must be a positive finite scalar.', fieldName);
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function logValue = appendLog(logValue, entry)
if isempty(logValue)
    logValue = entry;
else
    logValue(end + 1) = entry;
end
end
