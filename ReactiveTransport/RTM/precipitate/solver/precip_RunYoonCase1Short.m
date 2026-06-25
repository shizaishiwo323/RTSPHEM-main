function result = precip_RunYoonCase1Short(spec, options)
% precip_RunYoonCase1Short - Fixed-geometry Yoon Case 1 short smoke driver.
%
% Inputs:
%   spec    - benchmark spec from precip_ZhangYoonBenchmarkSpec.
%   options - endTime_s, dt_s, and optional initial state controls.
%
% Output:
%   result  - final state, area timeseries, and mass diagnostics.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

endTimeS = getFieldOrDefault(options, 'endTime_s', 13 * 60);
dtS = getFieldOrDefault(options, 'dt_s', 30);
if endTimeS <= 0 || dtS <= 0
    error('RTSPHEM:Precipitate:InvalidYoonCaseStep', ...
        'endTime_s and dt_s must be positive.');
end

state = precip_YoonMicrocontinuumSolver('initialize', spec);
state = seedSplitInletComponents(state, spec, options);
initialComponents = state.components;
initialState = state;
numSteps = ceil(endTimeS / dtS);
flowSolverFcn = getFieldOrDefault(options, 'flowSolverFcn', []);
transportMode = getFieldOrDefault(options, 'transportMode', 'analytic');
reactionSubcycling = getFieldOrDefault(options, 'reactionSubcycling', false);
transportOptions = struct();
transportOptions.advectiveCfl = getFieldOrDefault(options, ...
    'transportAdvectiveCfl', 0.5);
transportOptions.diffusiveCfl = getFieldOrDefault(options, ...
    'transportDiffusiveCfl', 0.25);
regionMasks = precip_BuildYoonRegionMasks(spec);
flowField = precip_RecomputeYoonFlowField(state, spec, flowSolverFcn);
numFlowRecomputationsTotal = double(~logical(getFieldOrDefault(flowField, ...
    'isProxy', false)));

timeS = zeros(numSteps, 1);
areaTotal = zeros(numSteps, 1);
areaFirstPore = zeros(numSteps, 1);
areaFirstThreePores = zeros(numSteps, 1);
naRelative = zeros(numSteps, 1);
clRelative = zeros(numSteps, 1);
accepted = false(numSteps, 1);
numBlockedCells = zeros(numSteps, 1);
numNewBlockedCells = zeros(numSteps, 1);
topologyChanged = false(numSteps, 1);
numFlowRecomputations = zeros(numSteps, 1);
relativePermeability = zeros(numSteps, 1);
pressureDropRelative = zeros(numSteps, 1);
flowRateRelative = zeros(numSteps, 1);
flowSolver = strings(numSteps, 1);
flowIsProxy = false(numSteps, 1);
flowIsStokes = false(numSteps, 1);
flowLinearResidualRelative = NaN(numSteps, 1);
transportAcceptedDt = zeros(numSteps, 1);
transportRequestedDt = zeros(numSteps, 1);
transportTotalAdvancedDt = zeros(numSteps, 1);
transportAcceptedSubstepCount = zeros(numSteps, 1);
transportRejectedSteps = zeros(numSteps, 1);
transportBoundaryClosure = zeros(numSteps, 1);
transportCandidateName = strings(numSteps, 1);
transportStableDt = zeros(numSteps, 1);
transportStabilityLimiter = strings(numSteps, 1);
reactionTotalCaRelativeError = zeros(numSteps, 1);
reactionTotalCRelativeError = zeros(numSteps, 1);
reactionTotalAlkalinityRelativeError = zeros(numSteps, 1);
reactionNaRelativeError = zeros(numSteps, 1);
reactionClRelativeError = zeros(numSteps, 1);
reactionMassAccepted = false(numSteps, 1);
reactionRequestedDt = zeros(numSteps, 1);
reactionAcceptedDt = zeros(numSteps, 1);
reactionTotalAdvancedDt = zeros(numSteps, 1);
reactionAcceptedSubstepCount = zeros(numSteps, 1);
reactionMaxAcceptedVmChange = zeros(numSteps, 1);
reactionStableDt = zeros(numSteps, 1);
reactionStabilityLimiter = strings(numSteps, 1);
targetTimesS = getFieldOrDefault(options, 'targetTimes_s', []);
targetTimesS = targetTimesS(:);
snapshotRows = initializeSnapshotRows(numel(targetTimesS));

for iStep = 1:numSteps
    localDt = min(dtS, endTimeS - (iStep - 1) * dtS);
    if getFieldOrDefault(options, 'refreshComponentsEachStep', false)
        state = seedSplitInletComponents(state, spec, options);
    end
    [state, transportStep] = advanceTransportIfRequested(state, spec, ...
        localDt, transportMode, transportOptions, flowField);
    [state, reactionStep] = advanceReaction(state, spec, localDt, ...
        reactionSubcycling);
    [state, flowStep] = precip_UpdateFlowMask(state, spec);
    if flowStep.topologyChanged
        flowField = precip_RecomputeYoonFlowField(state, spec, flowSolverFcn);
        numFlowRecomputationsTotal = numFlowRecomputationsTotal + 1;
    end
    metrics = precip_ComputeYoonAreaMetrics(state, spec, regionMasks);
    ledger = precip_ComputeComponentMassLedger(initialComponents, state.components, spec);
    reactionLedger = precip_ComputeYoonReactionMassLedger(initialState, ...
        state, spec);

    timeS(iStep) = iStep * dtS;
    areaTotal(iStep) = metrics.totalPrecipitatedArea_cm2;
    areaFirstPore(iStep) = metrics.firstPorePrecipitatedArea_cm2;
    areaFirstThreePores(iStep) = metrics.firstThreePoresPrecipitatedArea_cm2;
    naRelative(iStep) = ledger.relative.Na_total;
    clRelative(iStep) = ledger.relative.Cl_total;
    numBlockedCells(iStep) = flowStep.numBlockedCells;
    numNewBlockedCells(iStep) = flowStep.numNewBlockedCells;
    topologyChanged(iStep) = flowStep.topologyChanged;
    numFlowRecomputations(iStep) = numFlowRecomputationsTotal;
    relativePermeability(iStep) = flowField.relativePermeability;
    pressureDropRelative(iStep) = flowField.pressureDropRelative;
    flowRateRelative(iStep) = flowField.flowRateRelative;
    flowSolver(iStep) = string(flowField.solver);
    flowIsProxy(iStep) = logical(getFieldOrDefault(flowField, 'isProxy', false));
    flowIsStokes(iStep) = logical(getFieldOrDefault(flowField, 'isStokes', false));
    flowLinearResidualRelative(iStep) = getFieldOrDefault(flowField, ...
        'linearResidualRelative', NaN);
    transportRequestedDt(iStep) = transportStep.requestedDt_s;
    transportAcceptedDt(iStep) = transportStep.acceptedDt_s;
    transportTotalAdvancedDt(iStep) = transportStep.totalAdvancedDt_s;
    transportAcceptedSubstepCount(iStep) = transportStep.acceptedSubstepCount;
    transportRejectedSteps(iStep) = transportStep.rejectedStepCount;
    transportBoundaryClosure(iStep) = transportStep.maxBoundaryClosureError;
    transportCandidateName(iStep) = transportStep.candidateName;
    transportStableDt(iStep) = transportStep.stableDt_s;
    transportStabilityLimiter(iStep) = transportStep.stabilityLimiter;
    reactionTotalCaRelativeError(iStep) = reactionLedger.relative.totalCa_mol;
    reactionTotalCRelativeError(iStep) = reactionLedger.relative.totalC_mol;
    reactionTotalAlkalinityRelativeError(iStep) = ...
        reactionLedger.relative.totalAlkalinityEq;
    reactionNaRelativeError(iStep) = reactionLedger.relative.Na_total_mol;
    reactionClRelativeError(iStep) = reactionLedger.relative.Cl_total_mol;
    reactionMassAccepted(iStep) = reactionLedger.accepted;
    reactionRequestedDt(iStep) = reactionStep.requestedDt_s;
    reactionAcceptedDt(iStep) = reactionStep.acceptedDt_s;
    reactionTotalAdvancedDt(iStep) = reactionStep.totalAdvancedDt_s;
    reactionAcceptedSubstepCount(iStep) = reactionStep.acceptedSubstepCount;
    reactionMaxAcceptedVmChange(iStep) = reactionStep.maxAcceptedVmChange;
    reactionStableDt(iStep) = reactionStep.stableDt_s;
    reactionStabilityLimiter(iStep) = reactionStep.stabilityLimiter;
    accepted(iStep) = all(isfinite(state.Vm(:))) && ...
        min(state.components.Ca_total(:)) >= -1e-12 && ...
        min(state.components.C_total(:)) >= -1e-12 && ...
        min(state.components.Na_total(:)) >= -1e-12 && ...
        min(state.components.Cl_total(:)) >= -1e-12;

    snapshotRows = captureDueSnapshots(snapshotRows, targetTimesS, ...
        timeS(iStep), state, metrics);
end

result = struct();
result.caseName = 'yoon_case1_short_fixed_geometry';
result.state = state;
result.areaTimeseries = table(timeS, areaTotal, areaFirstPore, ...
    areaFirstThreePores, 'VariableNames', {'time_s', ...
    'totalPrecipitatedArea_cm2', 'firstPorePrecipitatedArea_cm2', ...
    'firstThreePoresPrecipitatedArea_cm2'});
result.massLedger = struct();
result.massLedger.relative = struct('Na_total', naRelative, 'Cl_total', clRelative);
result.massLedger.accepted = accepted;
result.reactionMassLedger = table(timeS, reactionTotalCaRelativeError, ...
    reactionTotalCRelativeError, reactionTotalAlkalinityRelativeError, ...
    reactionNaRelativeError, reactionClRelativeError, reactionMassAccepted, ...
    'VariableNames', {'time_s', 'totalCaRelativeError', ...
    'totalCRelativeError', 'totalAlkalinityRelativeError', ...
    'naRelativeError', 'clRelativeError', 'accepted'});
result.reactionDiagnostics = table(timeS, reactionRequestedDt, ...
    reactionAcceptedDt, reactionTotalAdvancedDt, reactionAcceptedSubstepCount, ...
    reactionMaxAcceptedVmChange, reactionStableDt, reactionStabilityLimiter, ...
    'VariableNames', {'time_s', 'requestedDt_s', 'acceptedDt_s', ...
    'totalAdvancedDt_s', 'acceptedSubstepCount', 'maxAcceptedVmChange', ...
    'stableDt_s', 'stabilityLimiter'});
result.flowDiagnostics = table(timeS, numBlockedCells, numNewBlockedCells, ...
    topologyChanged, numFlowRecomputations, relativePermeability, ...
    pressureDropRelative, flowRateRelative, flowSolver, flowIsProxy, ...
    flowIsStokes, flowLinearResidualRelative, 'VariableNames', {'time_s', ...
    'numBlockedCells', 'numNewBlockedCells', 'topologyChanged', ...
    'numFlowRecomputations', 'relativePermeability', ...
    'pressureDropRelative', 'flowRateRelative', 'flowSolver', ...
    'flowIsProxy', 'flowIsStokes', 'flowLinearResidualRelative'});
result.flowField = flowField;
result.regionMaskSource = regionMasks.source;
result.regionMaskFile = getFieldOrDefault(regionMasks, 'filePath', '');
result.transportMode = transportMode;
result.transportDiagnostics = table(timeS, transportCandidateName, ...
    transportRequestedDt, transportAcceptedDt, transportTotalAdvancedDt, ...
    transportAcceptedSubstepCount, transportRejectedSteps, ...
    transportBoundaryClosure, transportStableDt, transportStabilityLimiter, ...
    'VariableNames', {'time_s', 'candidateName', 'requestedDt_s', ...
    'acceptedDt_s', 'totalAdvancedDt_s', 'acceptedSubstepCount', ...
    'rejectedStepCount', 'maxBoundaryClosureError', 'stableDt_s', ...
    'stabilityLimiter'});
result.snapshots = snapshotRowsToTable(snapshotRows);
end

function state = seedSplitInletComponents(state, spec, options)
if isfield(options, 'initialState') && ~isempty(options.initialState)
    state = options.initialState;
    return;
end

fractionA = analyticSplitFraction(spec, state.grid.xCenters_cm, ...
    state.grid.yCenters_cm, state, options);
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    state.components.(fieldName) = fractionA .* spec.inletA.(fieldName) + ...
        (1 - fractionA) .* spec.inletB.(fieldName);
end
end

function [state, reactionStep] = advanceReaction(state, spec, dtS, reactionSubcycling)
if reactionSubcycling
    [state, reactionStep] = precip_RunYoonReactionSubcycle(state, spec, dtS);
    return;
end
beforeVm = state.Vm;
state = precip_AdvanceReaction(state, spec, dtS);
reactionStep = struct();
reactionStep.requestedDt_s = dtS;
reactionStep.acceptedDt_s = dtS;
reactionStep.totalAdvancedDt_s = dtS;
reactionStep.acceptedSubstepCount = 1;
reactionStep.maxAcceptedVmChange = max(abs(state.Vm(:) - beforeVm(:)));
reactionStep.stableDt_s = dtS;
reactionStep.stabilityLimiter = "none";
end

function fractionA = analyticSplitFraction(spec, xCenters, yCenters, state, options)
[xGrid, yGrid] = meshgrid(xCenters, yCenters);
u = max(spec.darcyVelocity_cm_s, eps);
d = diffusionField(spec, state, options);
timeSinceInlet = max(xGrid, spec.dx_cm) ./ u;
sigma = sqrt(2 .* d .* timeSinceInlet);
argument = (yGrid - spec.splitInletY_cm) ./ max(sqrt(2) .* sigma, eps);
fractionA = 0.5 .* erfc(argument);
fractionA = min(max(fractionA, 0), 1);
end

function d = diffusionField(spec, state, options)
if getFieldOrDefault(options, 'coupleDiffusionFeedback', false) && ...
        isfield(state, 'effectiveDiffusivity_cm2_s') && ...
        ~isempty(state.effectiveDiffusivity_cm2_s)
    columnD = mean(max(state.effectiveDiffusivity_cm2_s, 0), 1, 'omitnan');
    d = repmat(max(columnD, eps), spec.numY, 1);
else
    d = max(spec.diffusionCoefficient_cm2_s, eps);
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function [state, transportStep] = advanceTransportIfRequested(state, spec, dtS, transportMode, transportOptions, flowField)
transportStep = struct();
transportStep.candidateName = string(transportMode);
transportStep.acceptedDt_s = 0;
transportStep.requestedDt_s = dtS;
transportStep.totalAdvancedDt_s = 0;
transportStep.acceptedSubstepCount = 0;
transportStep.rejectedStepCount = 0;
transportStep.maxBoundaryClosureError = 0;
transportStep.stableDt_s = dtS;
transportStep.stabilityLimiter = "none";

switch transportMode
    case 'analytic'
        transportStep.candidateName = "analytic_split_inlet";
        transportStep.acceptedDt_s = dtS;
        transportStep.totalAdvancedDt_s = dtS;
        transportStep.acceptedSubstepCount = 1;
    case 'finite_volume'
        [state.components, transportStep] = runFiniteVolumeTransportInterval( ...
            state.components, spec, dtS, transportOptions, flowField);
    otherwise
        error('RTSPHEM:Precipitate:InvalidTransportMode', ...
            'Unsupported transportMode: %s.', transportMode);
end
end

function [components, transportStep] = runFiniteVolumeTransportInterval(components, spec, totalDtS, transportOptions, flowField)
remainingDt = totalDtS;
advancedDt = 0;
acceptedSubstepCount = 0;
rejectedStepCount = 0;
maxBoundaryClosureError = 0;
lastAcceptedDt = 0;
if isstruct(flowField) && isfield(flowField, 'velocityX_cm_s') && ...
        ~isempty(flowField.velocityX_cm_s)
    transportOptions.flowField = flowField;
end
[stableDt, stabilityDiagnostics] = precip_ComputeTransportStableDt(spec, ...
    transportOptions);
transportStepOptions = struct('boundaryMode', 'split_inlet');
if isfield(transportOptions, 'flowField')
    transportStepOptions.flowField = transportOptions.flowField;
end
candidateFcn = @(candidateComponents, dt) ...
    precip_AdvanceConservativeTransport2D(candidateComponents, spec, dt, ...
    transportStepOptions);

while remainingDt > max(10 * eps(totalDtS), 1e-15)
    componentsBeforeSubstep = components;
    attemptedDt = min(remainingDt, stableDt);
    if ~isfinite(attemptedDt)
        attemptedDt = remainingDt;
    end
    subcycle = precip_RunTransportSubcycle(componentsBeforeSubstep, spec, ...
        attemptedDt, candidateFcn, struct('massTolerance', Inf));
    components = subcycle.components;
    [~, boundaryLedger] = precip_AdvanceConservativeTransport2D( ...
        componentsBeforeSubstep, spec, subcycle.acceptedDt_s, ...
        transportStepOptions);
    acceptedSubstepCount = acceptedSubstepCount + 1;
    rejectedStepCount = rejectedStepCount + subcycle.rejectedStepCount;
    lastAcceptedDt = subcycle.acceptedDt_s;
    advancedDt = advancedDt + subcycle.acceptedDt_s;
    remainingDt = max(totalDtS - advancedDt, 0);
    maxBoundaryClosureError = max(maxBoundaryClosureError, ...
        boundaryLedger.maxBoundaryClosureError);
end

transportStep = struct();
transportStep.candidateName = "finite_volume_split_inlet";
transportStep.requestedDt_s = totalDtS;
transportStep.acceptedDt_s = lastAcceptedDt;
transportStep.totalAdvancedDt_s = advancedDt;
transportStep.acceptedSubstepCount = acceptedSubstepCount;
transportStep.rejectedStepCount = rejectedStepCount;
transportStep.maxBoundaryClosureError = maxBoundaryClosureError;
transportStep.stableDt_s = stableDt;
transportStep.stabilityLimiter = stabilityDiagnostics.limiter;
end

function rows = initializeSnapshotRows(numTargets)
rows = repmat(struct('targetTime_s', NaN, 'capturedTime_s', NaN, ...
    'Vm', [], 'totalPrecipitatedArea_cm2', NaN, ...
    'firstPorePrecipitatedArea_cm2', NaN, ...
    'firstThreePoresPrecipitatedArea_cm2', NaN, 'captured', false), ...
    numTargets, 1);
end

function rows = captureDueSnapshots(rows, targetTimesS, currentTimeS, state, metrics)
for iTarget = 1:numel(targetTimesS)
    if rows(iTarget).captured
        continue;
    end
    if currentTimeS + 10 * eps(currentTimeS) >= targetTimesS(iTarget)
        rows(iTarget).targetTime_s = targetTimesS(iTarget);
        rows(iTarget).capturedTime_s = currentTimeS;
        rows(iTarget).Vm = state.Vm;
        rows(iTarget).totalPrecipitatedArea_cm2 = metrics.totalPrecipitatedArea_cm2;
        rows(iTarget).firstPorePrecipitatedArea_cm2 = ...
            metrics.firstPorePrecipitatedArea_cm2;
        rows(iTarget).firstThreePoresPrecipitatedArea_cm2 = ...
            metrics.firstThreePoresPrecipitatedArea_cm2;
        rows(iTarget).captured = true;
    end
end
end

function snapshotTable = snapshotRowsToTable(rows)
if isempty(rows)
    snapshotTable = table([], [], {}, [], [], [], 'VariableNames', ...
        {'targetTime_s', 'capturedTime_s', 'Vm', ...
        'totalPrecipitatedArea_cm2', 'firstPorePrecipitatedArea_cm2', ...
        'firstThreePoresPrecipitatedArea_cm2'});
    return;
end
capturedRows = rows([rows.captured]);
if isempty(capturedRows)
    snapshotTable = table([], [], {}, [], [], [], 'VariableNames', ...
        {'targetTime_s', 'capturedTime_s', 'Vm', ...
        'totalPrecipitatedArea_cm2', 'firstPorePrecipitatedArea_cm2', ...
        'firstThreePoresPrecipitatedArea_cm2'});
    return;
end
snapshotTable = table([capturedRows.targetTime_s]', ...
    [capturedRows.capturedTime_s]', {capturedRows.Vm}', ...
    [capturedRows.totalPrecipitatedArea_cm2]', ...
    [capturedRows.firstPorePrecipitatedArea_cm2]', ...
    [capturedRows.firstThreePoresPrecipitatedArea_cm2]', ...
    'VariableNames', {'targetTime_s', 'capturedTime_s', 'Vm', ...
    'totalPrecipitatedArea_cm2', 'firstPorePrecipitatedArea_cm2', ...
    'firstThreePoresPrecipitatedArea_cm2'});
end
