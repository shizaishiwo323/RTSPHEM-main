function comparison = precip_RunYoonCase5ShortComparison(spec, options)
% precip_RunYoonCase5ShortComparison - Compare short Case 1 and Case 5 dissolution.
%
% The precipitation phase uses the base spec. The dissolution phase starts
% from the same precipitated state and compares dissolutionFactor = 1 and 300.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

precipitationEndTimeS = getFieldOrDefault(options, 'precipitationEndTime_s', 13 * 60);
dissolutionEndTimeS = getFieldOrDefault(options, 'dissolutionEndTime_s', 30 * 60);
dtS = getFieldOrDefault(options, 'dt_s', 30);

precipitationRun = precip_RunYoonCase1Short(spec, struct( ...
    'endTime_s', precipitationEndTimeS, 'dt_s', dtS, ...
    'refreshComponentsEachStep', getFieldOrDefault(options, ...
    'refreshComponentsEachStep', true), ...
    'coupleDiffusionFeedback', getFieldOrDefault(options, ...
    'coupleDiffusionFeedback', false)));

caseNames = ["case_1"; "case_5"];
dissolutionFactor = [1; 300];
finalPrecipitateMoles = zeros(2, 1);
finalPrecipitatedArea_cm2 = zeros(2, 1);
dissolutionOnlyFactor = true(2, 1);

for iCase = 1:2
    caseSpec = spec;
    caseSpec.dissolutionFactor = dissolutionFactor(iCase);
    state = precipitationRun.state;
    state = runDissolutionPhase(state, caseSpec, dissolutionEndTimeS, dtS);
    metrics = precip_ComputeYoonAreaMetrics(state, caseSpec, struct());
    finalPrecipitateMoles(iCase) = sum(state.precipitateMoles(:));
    finalPrecipitatedArea_cm2(iCase) = metrics.totalPrecipitatedArea_cm2;
end

comparison = table(caseNames, dissolutionFactor, finalPrecipitateMoles, ...
    finalPrecipitatedArea_cm2, dissolutionOnlyFactor);
end

function state = runDissolutionPhase(state, spec, endTimeS, dtS)
numSteps = ceil(endTimeS / dtS);
for iStep = 1:numSteps
    localDt = min(dtS, endTimeS - (iStep - 1) * dtS);
    state = setInitialWaterComponents(state, spec);
    state = precip_AdvanceReaction(state, spec, localDt);
    [state, ~] = precip_UpdateFlowMask(state, spec);
end
end

function state = setInitialWaterComponents(state, spec)
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    state.components.(fieldName) = ones(size(state.Vm)) .* spec.initial.(fieldName);
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
