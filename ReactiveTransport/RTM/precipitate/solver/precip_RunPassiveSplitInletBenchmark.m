function result = precip_RunPassiveSplitInletBenchmark(spec, options)
% precip_RunPassiveSplitInletBenchmark - No-reaction split-inlet transport smoke.
%
% Inputs:
%   spec    - struct from precip_ZhangYoonBenchmarkSpec.
%   options - optional controls; currently records numSteps for metadata.
%
% Output:
%   result  - conservative component fields and diagnostic mass/mixing checks.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

state = precip_YoonMicrocontinuumSolver('initialize', spec);
numSteps = getFieldOrDefault(options, 'numSteps', 1);
transportMode = getFieldOrDefault(options, 'transportMode', 'finite_volume');
dtS = getFieldOrDefault(options, 'dt_s', 1);
initialComponentSource = getFieldOrDefault(options, 'initialComponentSource', ...
    defaultInitialComponentSource(transportMode));
fractionA = splitInletFractionA(spec, state.grid.xCenters_cm, state.grid.yCenters_cm);
targetComponents = passiveSplitInletComponents(spec, fractionA);
initialComponents = initialPassiveComponents(spec, state, targetComponents, options);
[subcycle, components, massLedger, boundaryLedger] = runTransportCandidate( ...
    transportMode, initialComponents, spec, state, numSteps, dtS, options);

result = struct();
result.componentNames = spec.componentNames;
result.components = components;
result.fractionInletA = fractionA;
result.usedFiniteConcentrationLimiter = false;
result.numSteps = numSteps;
result.dt_s = dtS;
result.transportMode = transportMode;
result.initialComponentSource = initialComponentSource;
result.massLedger = massLedger;
result.mixingWidth = precip_TestMixingWidth(fractionA, spec);
result.transportSubcycle = subcycle;
if ~isempty(boundaryLedger)
    result.transportBoundaryLedger = boundaryLedger;
end
acceptance = buildAcceptanceManifest(result, options);
result.passiveTransportAccepted = acceptance.passiveTransportAccepted;
result.passiveTransportAcceptance = acceptance;
outputManifestPath = string(getFieldOrDefault(options, ...
    'outputManifestPath', ""));
if strlength(outputManifestPath) > 0
    writeJsonManifest(char(outputManifestPath), acceptance);
    result.outputManifestPath = outputManifestPath;
else
    result.outputManifestPath = "";
end
end

function components = passiveSplitInletCandidate(~, ~, spec, state)
fractionA = splitInletFractionA(spec, state.grid.xCenters_cm, state.grid.yCenters_cm);
components = passiveSplitInletComponents(spec, fractionA);
end

function [subcycle, components, massLedger, boundaryLedger] = ...
    runTransportCandidate(transportMode, initialComponents, spec, state, ...
    numSteps, dtS, options)
switch transportMode
    case 'analytic'
        candidateFcn = @(components, dt) passiveSplitInletCandidate( ...
            components, dt, spec, state);
        attemptedDt = numSteps;
        defaultMassTolerance = 1e-12;
        subcycle = precip_RunTransportSubcycle(initialComponents, spec, ...
            attemptedDt, candidateFcn, struct('massTolerance', ...
            getFieldOrDefault(options, 'massTolerance', defaultMassTolerance)));
        subcycle.candidateName = 'analytic_split_inlet';
        components = subcycle.components;
        massLedger = precip_ComputeComponentMassLedger(initialComponents, ...
            components, spec);
        boundaryLedger = [];
    case 'finite_volume'
        state.components = initialComponents;
        state = precip_RefreshYoonComponentMolesFromAqueous(state, spec);
        initialState = state;
        [state, boundaryLedger, transportDiagnostics] = ...
            runFiniteVolumePassiveSteps(state, spec, numSteps, dtS);
        components = state.components;
        subcycle = struct();
        subcycle.components = components;
        subcycle.acceptedDt_s = transportDiagnostics.totalAdvancedDt_s;
        subcycle.totalAdvancedDt_s = transportDiagnostics.totalAdvancedDt_s;
        subcycle.acceptedSubstepCount = ...
            transportDiagnostics.acceptedSubstepCount;
        subcycle.rejectedStepCount = transportDiagnostics.rejectedStepCount;
        subcycle.diagnostics = struct('dt_s', ...
            transportDiagnostics.totalAdvancedDt_s, ...
            'accepted', true, 'reason', "accepted", ...
            'maxRelativeMassError', boundaryLedger.maxBoundaryClosureError);
        subcycle.candidateName = 'finite_volume_split_inlet';
        massLedger = computeStateInventoryLedger(initialState, state, spec);
    otherwise
        error('RTSPHEM:Precipitate:InvalidTransportMode', ...
            'Unsupported transportMode: %s.', transportMode);
end
end

function [state, boundaryLedger, diagnostics] = runFiniteVolumePassiveSteps( ...
    state, spec, numSteps, dtS)
boundaryLedger = initializeBoundaryLedger(spec);
diagnostics = struct();
diagnostics.totalAdvancedDt_s = 0;
diagnostics.acceptedSubstepCount = 0;
diagnostics.rejectedStepCount = 0;
for iStep = 1:numSteps
    remainingDt = dtS;
    while remainingDt > max(10 * eps(dtS), 1e-15)
        [stableDt, ~] = precip_ComputeMicrocontinuumStableDt(state, spec, ...
            stateTransportOptions(state, struct()));
        attemptedDt = min(remainingDt, stableDt);
        if ~isfinite(attemptedDt) || attemptedDt <= 0
            attemptedDt = remainingDt;
        end
        [state, stepLedger] = precip_AdvanceMicrocontinuumTransport2D( ...
            state, spec, attemptedDt, ...
            stateTransportOptions(state, struct('boundaryMode', 'split_inlet')));
        boundaryLedger = accumulateBoundaryLedger(boundaryLedger, ...
            stepLedger);
        diagnostics.totalAdvancedDt_s = diagnostics.totalAdvancedDt_s + ...
            attemptedDt;
        diagnostics.acceptedSubstepCount = ...
            diagnostics.acceptedSubstepCount + 1;
        remainingDt = max(remainingDt - attemptedDt, 0);
    end
end
boundaryLedger = finalizeBoundaryLedger(boundaryLedger);
end

function ledger = initializeBoundaryLedger(spec)
ledger = struct();
ledger.massChange = struct();
ledger.boundaryIntegral = struct();
ledger.boundaryClosureError = struct();
ledger.relativeClosureError = struct();
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    ledger.massChange.(fieldName) = 0;
    ledger.boundaryIntegral.(fieldName) = 0;
    ledger.boundaryClosureError.(fieldName) = 0;
    ledger.relativeClosureError.(fieldName) = NaN;
end
ledger.maxBoundaryClosureError = 0;
ledger.maxBoundaryClosureRelativeError = NaN;
end

function ledger = accumulateBoundaryLedger(ledger, stepLedger)
fields = fieldnames(ledger.massChange);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    ledger.massChange.(fieldName) = ledger.massChange.(fieldName) + ...
        stepLedger.massChange.(fieldName);
    ledger.boundaryIntegral.(fieldName) = ...
        ledger.boundaryIntegral.(fieldName) + ...
        (stepLedger.massChange.(fieldName) - ...
        stepLedger.boundaryClosureError.(fieldName));
    ledger.boundaryClosureError.(fieldName) = ...
        ledger.boundaryClosureError.(fieldName) + ...
        stepLedger.boundaryClosureError.(fieldName);
end
end

function ledger = finalizeBoundaryLedger(ledger)
fields = fieldnames(ledger.massChange);
maxAbsClosure = 0;
maxRelativeClosure = 0;
for iField = 1:numel(fields)
    fieldName = fields{iField};
    closure = ledger.boundaryClosureError.(fieldName);
    denominator = max([abs(ledger.massChange.(fieldName)), ...
        abs(ledger.boundaryIntegral.(fieldName)), eps]);
    ledger.relativeClosureError.(fieldName) = abs(closure) / denominator;
    maxAbsClosure = max(maxAbsClosure, abs(closure));
    maxRelativeClosure = max(maxRelativeClosure, ...
        ledger.relativeClosureError.(fieldName));
end
ledger.maxBoundaryClosureError = maxAbsClosure;
ledger.maxBoundaryClosureRelativeError = maxRelativeClosure;
end

function ledger = computeStateInventoryLedger(initialState, currentState, spec)
ledger = struct();
ledger.before = struct();
ledger.after = struct();
ledger.delta = struct();
ledger.relative = struct();
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    before = sum(initialState.componentMoles.(fieldName)(:));
    after = sum(currentState.componentMoles.(fieldName)(:));
    ledger.before.(fieldName) = before;
    ledger.after.(fieldName) = after;
    ledger.delta.(fieldName) = after - before;
    ledger.relative.(fieldName) = ledger.delta.(fieldName) / ...
        max(abs(before), eps);
end
end

function components = initialPassiveComponents(spec, state, targetComponents, options)
transportMode = getFieldOrDefault(options, 'transportMode', 'finite_volume');
source = getFieldOrDefault(options, 'initialComponentSource', ...
    defaultInitialComponentSource(transportMode));
switch source
    case 'target'
        components = targetComponents;
    case 'initial'
        components = state.components;
    otherwise
        error('RTSPHEM:Precipitate:InvalidInitialComponentSource', ...
            'Unsupported initialComponentSource: %s.', source);
end
end

function source = defaultInitialComponentSource(transportMode)
switch transportMode
    case 'finite_volume'
        source = 'initial';
    otherwise
        source = 'target';
end
end

function components = passiveSplitInletComponents(spec, fractionA)
components = struct();
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    components.(fieldName) = fractionA .* spec.inletA.(fieldName) + ...
        (1 - fractionA) .* spec.inletB.(fieldName);
end
end

function fractionA = splitInletFractionA(spec, xCenters, yCenters)
[xGrid, yGrid] = meshgrid(xCenters, yCenters);
u = max(spec.darcyVelocity_cm_s, eps);
d = max(spec.diffusionCoefficient_cm2_s, eps);
timeSinceInlet = max(xGrid, spec.dx_cm) ./ u;
sigma = sqrt(2 .* d .* timeSinceInlet);
argument = (yGrid - spec.splitInletY_cm) ./ max(sqrt(2) .* sigma, eps);
fractionA = 0.5 .* erfc(argument);
fractionA = min(max(fractionA, 0), 1);
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function options = stateTransportOptions(state, options)
if isfield(state, 'substrateMask') && ~isempty(state.substrateMask)
    options.substrateMask = state.substrateMask;
end
if isfield(state, 'blockedMask') && ~isempty(state.blockedMask)
    options.blockedMask = state.blockedMask;
end
if isfield(state, 'effectiveDiffusivity_cm2_s') && ...
        ~isempty(state.effectiveDiffusivity_cm2_s)
    options.effectiveDiffusivity_cm2_s = state.effectiveDiffusivity_cm2_s;
end
end

function manifest = buildAcceptanceManifest(result, options)
maxAcceptedMassRelativeError = getFieldOrDefault(options, ...
    'maxAcceptedMassRelativeError', 1e-4);
maxAcceptedMixingSymmetryError = getFieldOrDefault(options, ...
    'maxAcceptedMixingSymmetryError', 0.05);
maxAcceptedRejectedStepCount = getFieldOrDefault(options, ...
    'maxAcceptedRejectedStepCount', 0);
minConcentration = minimumNonnegativeComponent(result.components);
maxAbsInertRelativeMassError = max(abs([ ...
    result.massLedger.relative.Na_total, ...
    result.massLedger.relative.Cl_total]));

manifest = struct();
manifest.runner = 'precip_RunPassiveSplitInletBenchmark';
manifest.transportMode = result.transportMode;
manifest.initialComponentSource = result.initialComponentSource;
manifest.componentNames = result.componentNames;
manifest.numSteps = result.numSteps;
manifest.dt_s = result.dt_s;
manifest.usedFiniteConcentrationLimiter = ...
    result.usedFiniteConcentrationLimiter;
manifest.minConcentration = minConcentration;
manifest.maxAbsInertRelativeMassError = maxAbsInertRelativeMassError;
manifest.maxAcceptedMassRelativeError = maxAcceptedMassRelativeError;
manifest.maxAbsInertBoundaryClosureRelativeError = ...
    maxInertBoundaryClosureRelativeError(result);
manifest.maxAcceptedBoundaryClosureRelativeError = getFieldOrDefault( ...
    options, 'maxAcceptedBoundaryClosureRelativeError', 1e-8);
manifest.mixingSymmetryError = result.mixingWidth.symmetryError;
manifest.maxAcceptedMixingSymmetryError = ...
    maxAcceptedMixingSymmetryError;
manifest.rejectedStepCount = result.transportSubcycle.rejectedStepCount;
manifest.maxAcceptedRejectedStepCount = maxAcceptedRejectedStepCount;
manifest.passiveTransportAccepted = ...
    ~manifest.usedFiniteConcentrationLimiter && ...
    manifest.minConcentration >= -1e-15 && ...
    manifest.maxAbsInertBoundaryClosureRelativeError <= ...
        manifest.maxAcceptedBoundaryClosureRelativeError && ...
    manifest.mixingSymmetryError <= maxAcceptedMixingSymmetryError && ...
    manifest.rejectedStepCount <= maxAcceptedRejectedStepCount;
end

function value = maxInertBoundaryClosureRelativeError(result)
if ~isfield(result, 'transportBoundaryLedger') || ...
        ~isfield(result.transportBoundaryLedger, 'relativeClosureError')
    value = max(abs([result.massLedger.relative.Na_total, ...
        result.massLedger.relative.Cl_total]));
    return;
end
value = max(abs([ ...
    result.transportBoundaryLedger.relativeClosureError.Na_total, ...
    result.transportBoundaryLedger.relativeClosureError.Cl_total]));
end

function minValue = minimumNonnegativeComponent(components)
fields = {'Ca_total', 'C_total', 'Na_total', 'Cl_total'};
minValue = Inf;
for iField = 1:numel(fields)
    if isfield(components, fields{iField})
        minValue = min(minValue, min(components.(fields{iField})(:)));
    end
end
if isinf(minValue)
    minValue = NaN;
end
end

function writeJsonManifest(path, manifest)
parentDir = fileparts(path);
if ~isempty(parentDir) && ~isfolder(parentDir)
    mkdir(parentDir);
end
fid = fopen(path, 'w');
if fid < 0
    error('RTSPHEM:Precipitate:ManifestWriteFailed', ...
        'Could not open manifest for writing: %s.', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
clear cleanupObj;
end
