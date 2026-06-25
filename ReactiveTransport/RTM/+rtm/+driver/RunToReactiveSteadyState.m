function info = RunToReactiveSteadyState(driver, options)
%RUNTOREACTIVESTEADYSTATE Advance RT windows until reaction rate is steady.
%
% This helper implements the Molins Part I quasi-steady initializer contract:
% transport and chemistry advance through the conservative driver, while the
% initializer only decides when the window-averaged realized mineral rate has
% stabilized.

if nargin < 2 || isempty(options)
    options = struct();
end
validateDriver(driver);

windowSeconds = getScalarOption(options, 'window_s', 1.0, true);
maxWindows = getIntegerOption(options, 'maxWindows', 50);
absoluteTolerance = getScalarOption(options, ...
    'rateAbsoluteTolerance_mol_s', 1e-15, false);
relativeTolerance = getScalarOption(options, ...
    'rateRelativeTolerance', 1e-6, false);
requiredConsecutivePasses = getIntegerOption(options, ...
    'requiredConsecutivePasses', 1);

windowRates = zeros(0, 1);
summaries = struct([]);
consecutivePasses = 0;
converged = false;
reason = "max_windows";

for iWindow = 1:maxWindows
    summary = driver.runRtSubcycles(windowSeconds);
    summaries = appendSummary(summaries, summary);

    windowRate = realizedRateFromSummary(summary);
    windowRates(end + 1, 1) = windowRate; %#ok<AGROW>

    if abs(windowRate) <= absoluteTolerance
        converged = true;
        reason = "absolute_rate_tolerance";
        break;
    end

    if numel(windowRates) >= 2
        previousRate = windowRates(end - 1);
        relativeChange = abs(windowRate - previousRate) ./ ...
            max(abs(previousRate), eps);
        if relativeChange <= relativeTolerance
            consecutivePasses = consecutivePasses + 1;
        else
            consecutivePasses = 0;
        end
        if consecutivePasses >= requiredConsecutivePasses
            converged = true;
            reason = "relative_rate_tolerance";
            break;
        end
    end
end

info = struct();
info.converged = converged;
info.reason = reason;
info.windows = numel(windowRates);
info.time_s = totalTimeFromSummaries(summaries);
info.steady_rate_mol_s = lastOrDefault(windowRates, NaN);
info.window_rates_mol_s = windowRates;
info.summaries = summaries;
info.state = lastSummaryField(summaries, 'state', struct());
info.geometry = lastSummaryField(summaries, 'geometry', struct());
info.solver_state = lastSummaryField(summaries, 'solver_state', struct());
end

function validateDriver(driver)
if isempty(driver) || ~ismethod(driver, 'runRtSubcycles')
    error('RTSPHEM:Driver:InvalidSteadyStateDriver', ...
        'driver must provide a runRtSubcycles(totalTimeSeconds) method.');
end
end

function value = getScalarOption(options, fieldName, defaultValue, requirePositive)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && isfinite(value))
    error('RTSPHEM:Driver:InvalidSteadyStateOption', ...
        'options.%s must be a finite scalar.', fieldName);
end
if requirePositive
    valid = value > 0;
else
    valid = value >= 0;
end
if ~valid
    error('RTSPHEM:Driver:InvalidSteadyStateOption', ...
        'options.%s is outside the allowed range.', fieldName);
end
end

function value = getIntegerOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && isfinite(value) && value >= 1 && value == round(value))
    error('RTSPHEM:Driver:InvalidSteadyStateOption', ...
        'options.%s must be a positive integer.', fieldName);
end
end

function summaries = appendSummary(summaries, summary)
if isempty(summaries)
    summaries = summary;
else
    summaries(end + 1) = summary;
end
end

function rate = realizedRateFromSummary(summary)
acceptedTime = 0;
realizedMoles = 0;
for iStep = 1:numel(summary.step_results)
    stepResult = summary.step_results(iStep);
    if isfield(stepResult, 'diagnostics') && stepResult.diagnostics.accepted
        acceptedTime = acceptedTime + stepResult.transport_ledger.dt_s;
        realizedMoles = realizedMoles + sum( ...
            stepResult.reaction_result.realized_interface_moles(:), 'omitnan');
    end
end
if acceptedTime > 0
    rate = realizedMoles ./ acceptedTime;
else
    rate = 0;
end
end

function value = totalTimeFromSummaries(summaries)
value = 0;
for iSummary = 1:numel(summaries)
    value = value + summaries(iSummary).time_s;
end
end

function value = lastOrDefault(values, defaultValue)
if isempty(values)
    value = defaultValue;
else
    value = values(end);
end
end

function value = lastSummaryField(summaries, fieldName, defaultValue)
if isempty(summaries) || ~isfield(summaries(end), fieldName)
    value = defaultValue;
else
    value = summaries(end).(fieldName);
end
end
