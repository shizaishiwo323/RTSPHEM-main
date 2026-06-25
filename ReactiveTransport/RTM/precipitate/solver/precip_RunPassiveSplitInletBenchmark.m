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
transportMode = getFieldOrDefault(options, 'transportMode', 'analytic');
dtS = getFieldOrDefault(options, 'dt_s', 1);
fractionA = splitInletFractionA(spec, state.grid.xCenters_cm, state.grid.yCenters_cm);
targetComponents = passiveSplitInletComponents(spec, fractionA);
initialComponents = initialPassiveComponents(spec, state, targetComponents, options);
[candidateFcn, candidateName, attemptedDt, defaultMassTolerance] = ...
    transportCandidate(transportMode, spec, state, numSteps, dtS);
subcycle = precip_RunTransportSubcycle(initialComponents, spec, attemptedDt, ...
    candidateFcn, struct('massTolerance', getFieldOrDefault(options, ...
    'massTolerance', defaultMassTolerance)));
subcycle.candidateName = candidateName;
components = subcycle.components;

result = struct();
result.componentNames = spec.componentNames;
result.components = components;
result.fractionInletA = fractionA;
result.usedFiniteConcentrationLimiter = false;
result.numSteps = numSteps;
result.dt_s = dtS;
result.transportMode = transportMode;
result.massLedger = precip_ComputeComponentMassLedger(components, components, spec);
result.mixingWidth = precip_TestMixingWidth(fractionA, spec);
result.transportSubcycle = subcycle;
if strcmp(transportMode, 'finite_volume')
    [~, result.transportBoundaryLedger] = precip_AdvanceConservativeTransport2D( ...
        initialComponents, spec, subcycle.acceptedDt_s, ...
        struct('boundaryMode', 'split_inlet'));
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

function components = finiteVolumeSplitInletCandidate(components, dt, spec)
components = precip_AdvanceConservativeTransport2D(components, spec, dt, ...
    struct('boundaryMode', 'split_inlet'));
end

function [candidateFcn, candidateName, attemptedDt, defaultMassTolerance] = ...
    transportCandidate(transportMode, spec, state, numSteps, dtS)
switch transportMode
    case 'analytic'
        candidateFcn = @(components, dt) passiveSplitInletCandidate( ...
            components, dt, spec, state);
        candidateName = 'analytic_split_inlet';
        attemptedDt = numSteps;
        defaultMassTolerance = 1e-12;
    case 'finite_volume'
        candidateFcn = @(components, dt) finiteVolumeSplitInletCandidate( ...
            components, dt, spec);
        candidateName = 'finite_volume_split_inlet';
        attemptedDt = numSteps * dtS;
        defaultMassTolerance = Inf;
    otherwise
        error('RTSPHEM:Precipitate:InvalidTransportMode', ...
            'Unsupported transportMode: %s.', transportMode);
end
end

function components = initialPassiveComponents(spec, state, targetComponents, options)
source = getFieldOrDefault(options, 'initialComponentSource', 'target');
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
manifest.componentNames = result.componentNames;
manifest.numSteps = result.numSteps;
manifest.dt_s = result.dt_s;
manifest.usedFiniteConcentrationLimiter = ...
    result.usedFiniteConcentrationLimiter;
manifest.minConcentration = minConcentration;
manifest.maxAbsInertRelativeMassError = maxAbsInertRelativeMassError;
manifest.maxAcceptedMassRelativeError = maxAcceptedMassRelativeError;
manifest.mixingSymmetryError = result.mixingWidth.symmetryError;
manifest.maxAcceptedMixingSymmetryError = ...
    maxAcceptedMixingSymmetryError;
manifest.rejectedStepCount = result.transportSubcycle.rejectedStepCount;
manifest.maxAcceptedRejectedStepCount = maxAcceptedRejectedStepCount;
manifest.passiveTransportAccepted = ...
    ~manifest.usedFiniteConcentrationLimiter && ...
    manifest.minConcentration >= -1e-15 && ...
    manifest.maxAbsInertRelativeMassError <= ...
        maxAcceptedMassRelativeError && ...
    manifest.mixingSymmetryError <= maxAcceptedMixingSymmetryError && ...
    manifest.rejectedStepCount <= maxAcceptedRejectedStepCount;
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
