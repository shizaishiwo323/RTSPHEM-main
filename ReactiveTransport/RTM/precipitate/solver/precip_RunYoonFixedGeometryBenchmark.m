function result = precip_RunYoonFixedGeometryBenchmark(spec, options)
% precip_RunYoonFixedGeometryBenchmark - Diagnostic fixed-geometry Yoon runner.
%
% Inputs:
%   spec    - benchmark spec from precip_ZhangYoonBenchmarkSpec.
%   options - outputRoot, endTime_s, dt_s, targetTimes_s, and driver options.
%
% Output:
%   result  - run result, export manifest, and manifest path.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

[spec, geometryPackage] = applyGeometryPackageIfConfigured(spec);
targetTimesS = getFieldOrDefault(options, 'targetTimes_s', [13, 18, 118] * 60);
endTimeS = getFieldOrDefault(options, 'endTime_s', max(targetTimesS));
dtS = getFieldOrDefault(options, 'dt_s', 30);
outputRoot = getFieldOrDefault(options, 'outputRoot', defaultOutputRoot());

runOptions = options;
runOptions.endTime_s = endTimeS;
runOptions.dt_s = dtS;
runOptions.targetTimes_s = targetTimesS;
run = precip_RunYoonCase1Short(spec, runOptions);
snapshotManifest = precip_ExportYoonSnapshots(outputRoot, run.snapshots, spec);
flowDiagnosticsCsv = fullfile(outputRoot, 'yoon_flow_diagnostics.csv');
writeFlowDiagnostics(flowDiagnosticsCsv, run.flowDiagnostics);
transportDiagnosticsCsv = fullfile(outputRoot, ...
    'yoon_transport_diagnostics.csv');
writeFlowDiagnostics(transportDiagnosticsCsv, run.transportDiagnostics);
reactionMassLedgerCsv = fullfile(outputRoot, 'yoon_reaction_mass_ledger.csv');
writeFlowDiagnostics(reactionMassLedgerCsv, run.reactionMassLedger);
reactionDiagnosticsCsv = fullfile(outputRoot, 'yoon_reaction_diagnostics.csv');
writeFlowDiagnostics(reactionDiagnosticsCsv, run.reactionDiagnostics);

manifest = struct();
manifest.runner = 'precip_RunYoonFixedGeometryBenchmark';
manifest.modelFamily = spec.modelFamily;
manifest.caseName = run.caseName;
manifest.outputRoot = outputRoot;
manifest.endTime_s = endTimeS;
manifest.dt_s = dtS;
manifest.targetTimes_s = targetTimesS(:)';
manifest.numSnapshots = snapshotManifest.numSnapshots;
manifest.capturedTargetTimes_s = snapshotManifest.capturedTargetTimes_s;
manifest.missingTargetTimes_s = setdiff(targetTimesS(:)', ...
    manifest.capturedTargetTimes_s);
manifest.targetSnapshotsComplete = ...
    isempty(manifest.missingTargetTimes_s) && ...
    manifest.numSnapshots >= numel(targetTimesS);
manifest.snapshotDir = snapshotManifest.snapshotDir;
manifest.areaCsv = snapshotManifest.areaCsv;
manifest.flowDiagnosticsCsv = flowDiagnosticsCsv;
manifest.transportDiagnosticsCsv = transportDiagnosticsCsv;
manifest.reactionMassLedgerCsv = reactionMassLedgerCsv;
manifest.reactionDiagnosticsCsv = reactionDiagnosticsCsv;
manifest.substrateGeometrySource = run.state.substrateGeometrySource;
manifest.substrateMaskFile = run.state.substrateMaskFile;
manifest.regionMaskSource = run.regionMaskSource;
manifest.regionMaskFile = run.regionMaskFile;
manifest.geometryPackageDir = '';
manifest.geometryPackageIsQuantitative = false;
manifest.geometryPackageName = '';
manifest.geometryPackageNote = '';
if ~isempty(geometryPackage)
    manifest.geometryPackageDir = char(geometryPackage.packageDir);
    manifest.geometryPackageIsQuantitative = ...
        geometryPackage.isQuantitativeGeometry;
    manifest.geometryPackageName = char(geometryPackage.packageName);
    manifest.geometryPackageNote = geometryPackage.note;
    manifest.geometryPackageAssetFilesVerified = ...
        geometryPackage.assetFilesVerified;
    manifest.geometryPackageNumSubstrateCells = ...
        geometryPackage.numSubstrateCells;
    manifest.geometryPackageNumFirstPoreCells = ...
        geometryPackage.numFirstPoreCells;
    manifest.geometryPackageNumFirstThreePoresCells = ...
        geometryPackage.numFirstThreePoresCells;
end
manifest.finalNumFlowRecomputations = ...
    run.flowDiagnostics.numFlowRecomputations(end);
manifest.finalRelativePermeability = ...
    run.flowDiagnostics.relativePermeability(end);
manifest.finalPressureDropRelative = ...
    run.flowDiagnostics.pressureDropRelative(end);
manifest.finalFlowRateRelative = ...
    run.flowDiagnostics.flowRateRelative(end);
manifest.finalFlowSolver = char(run.flowDiagnostics.flowSolver(end));
manifest.finalFlowIsProxy = run.flowDiagnostics.flowIsProxy(end);
manifest.finalFlowIsStokes = run.flowDiagnostics.flowIsStokes(end);
manifest.finalFlowLinearResidualRelative = ...
    run.flowDiagnostics.flowLinearResidualRelative(end);
manifest.maxAcceptedFlowLinearResidualRelative = getFieldOrDefault( ...
    options, 'maxAcceptedFlowLinearResidualRelative', 1e-8);
manifest.flowTopologyChangedAny = any(run.flowDiagnostics.topologyChanged);
manifest.totalFlowTopologyChangedSteps = ...
    nnz(run.flowDiagnostics.topologyChanged);
manifest.flowRecomputedAfterTopologyChange = ...
    manifest.flowTopologyChangedAny && ...
    manifest.finalNumFlowRecomputations >= ...
        manifest.totalFlowTopologyChangedSteps;
manifest.flowBackendProductionValidated = getFieldOrDefault(options, ...
    'flowBackendProductionValidated', false);
manifest.flowFeedbackAccepted = manifest.flowTopologyChangedAny && ...
    manifest.flowRecomputedAfterTopologyChange && ...
    manifest.finalFlowIsStokes && ~manifest.finalFlowIsProxy && ...
    manifest.flowBackendProductionValidated;
manifest.finalTransportCandidate = ...
    char(run.transportDiagnostics.candidateName(end));
manifest.finalTransportRejectedStepCount = ...
    run.transportDiagnostics.rejectedStepCount(end);
manifest.finalTransportBoundaryClosureError = ...
    run.transportDiagnostics.maxBoundaryClosureError(end);
manifest.finalTransportStableDt_s = ...
    run.transportDiagnostics.stableDt_s(end);
manifest.finalTransportStabilityLimiter = ...
    char(run.transportDiagnostics.stabilityLimiter(end));
manifest.finalReactionTotalCaRelativeError = ...
    run.reactionMassLedger.totalCaRelativeError(end);
manifest.finalReactionTotalCRelativeError = ...
    run.reactionMassLedger.totalCRelativeError(end);
manifest.finalReactionTotalAlkalinityRelativeError = ...
    run.reactionMassLedger.totalAlkalinityRelativeError(end);
manifest.finalReactionMassAccepted = run.reactionMassLedger.accepted(end);
manifest.finalReactionAcceptedSubstepCount = ...
    run.reactionDiagnostics.acceptedSubstepCount(end);
manifest.finalReactionMaxAcceptedVmChange = ...
    run.reactionDiagnostics.maxAcceptedVmChange(end);
manifest.finalReactionStabilityLimiter = ...
    char(run.reactionDiagnostics.stabilityLimiter(end));
morphology = computeCenterBandMorphology(run.state, spec);
manifest.centerBandMorphologyAccepted = ...
    morphology.centerBandMorphologyAccepted;
manifest.finalPrecipitateCentroidY_cm = ...
    morphology.finalPrecipitateCentroidY_cm;
manifest.finalPrecipitateCentroidDistanceFromSplit_cm = ...
    morphology.finalPrecipitateCentroidDistanceFromSplit_cm;
manifest.finalActivePrecipitateArea_cm2 = ...
    morphology.finalActivePrecipitateArea_cm2;
manifest.centerBandMorphologyTolerance_cm = ...
    morphology.centerBandMorphologyTolerance_cm;
manifest.yoonRateSource = char(string(getFieldOrDefault(spec.yoonRate, ...
    'source', '')));
manifest.yoonRateK1 = getFieldOrDefault(spec.yoonRate, 'k1', NaN);
manifest.yoonRateK2 = getFieldOrDefault(spec.yoonRate, 'k2', NaN);
manifest.yoonRateK3 = getFieldOrDefault(spec.yoonRate, 'k3', NaN);
manifest.yoonRateUnits = char(string(getFieldOrDefault(spec.yoonRate, ...
    'units', '')));
manifest.yoonRateSourceDoi = char(string(getFieldOrDefault(spec.yoonRate, ...
    'sourceDoi', '')));
manifest.yoonRateSourceEquation = char(string(getFieldOrDefault(spec.yoonRate, ...
    'sourceEquation', '')));
manifest.yoonRateMineralPhase = char(string(getFieldOrDefault(spec.yoonRate, ...
    'mineralPhase', '')));
manifest.yoonRateSourceValuesVerified = getFieldOrDefault(spec.yoonRate, ...
    'sourceValuesVerified', false);
manifest.yoonRateConstantsLocked = isLockedYoonRateSource( ...
    manifest.yoonRateSource);
manifest.usesFiniteConcentrationLimiter = false;
manifest.concentrationLimiterNote = ['Yoon micro-continuum diagnostic path ' ...
    'uses conservative no-clipping transport/reaction operators; the legacy ' ...
    'PHREEQC transport upper-bound limiter is not part of this runner.'];
manifest.matFiles = cellstr(snapshotManifest.matFiles);
manifest.outputEvidenceFilesVerified = verifyOutputEvidenceFiles(manifest);
manifest.isQuantitativeBenchmark = false;
manifest.note = ['Fixed-geometry diagnostic Yoon-path output; not a ' ...
    'grid-converged Zhang/Yoon reproduction.'];

manifestPath = fullfile(outputRoot, 'yoon_fixed_geometry_manifest.json');
writeJsonManifest(manifestPath, manifest);

result = struct();
result.run = run;
result.manifest = manifest;
result.manifestPath = manifestPath;
end

function [spec, geometryPackage] = applyGeometryPackageIfConfigured(spec)
geometryPackage = [];
if ~isfield(spec, 'geometryPackageDir') || isempty(spec.geometryPackageDir)
    return;
end
geometryPackage = precip_LoadZhangYoonGeometryPackage( ...
    spec.geometryPackageDir, spec);
if ~geometryPackage.isQuantitativeGeometry
    error('RTSPHEM:Precipitate:IncompleteGeometryPackage', ...
        'Configured geometryPackageDir is incomplete: %s', geometryPackage.note);
end
spec.substrateMaskFile = char(geometryPackage.substrateMaskFile);
spec.regionMaskFile = char(geometryPackage.regionMaskFile);
end

function writeFlowDiagnostics(path, flowDiagnostics)
parentDir = fileparts(path);
if ~isfolder(parentDir)
    mkdir(parentDir);
end
writetable(flowDiagnostics, path);
end

function outputRoot = defaultOutputRoot()
moduleRoot = fileparts(fileparts(mfilename('fullpath')));
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
outputRoot = fullfile(moduleRoot, 'outputs', 'yoon_fixed_geometry', timestamp);
end

function tf = verifyOutputEvidenceFiles(manifest)
tf = isfolder(manifest.outputRoot) && isfolder(manifest.snapshotDir) && ...
    isfile(manifest.areaCsv) && isfile(manifest.flowDiagnosticsCsv) && ...
    isfile(manifest.transportDiagnosticsCsv) && ...
    isfile(manifest.reactionMassLedgerCsv) && ...
    isfile(manifest.reactionDiagnosticsCsv) && ...
    ~isempty(manifest.matFiles) && all(isfile(string(manifest.matFiles)));
end

function writeJsonManifest(path, manifest)
parentDir = fileparts(path);
if ~isfolder(parentDir)
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

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function tf = isLockedYoonRateSource(rateSource)
source = lower(string(rateSource));
tf = strlength(source) > 0 && ~contains(source, "smoke") && ...
    ~contains(source, "pending");
end

function morphology = computeCenterBandMorphology(state, spec)
activeMask = state.Vm > spec.areaVmThreshold;
if isfield(state, 'substrateMask') && ~isempty(state.substrateMask)
    activeMask = activeMask & ~state.substrateMask;
end
cellAreaCm2 = spec.dx_cm * spec.dy_cm;
activeAreaCm2 = nnz(activeMask) * cellAreaCm2;
if activeAreaCm2 > 0
    yCenters = ((1:spec.numY)' - 0.5) .* spec.dy_cm;
    yGrid = repmat(yCenters, 1, spec.numX);
    centroidY = mean(yGrid(activeMask));
    distanceFromSplit = abs(centroidY - spec.splitInletY_cm);
else
    centroidY = NaN;
    distanceFromSplit = NaN;
end
toleranceCm = 0.25 * spec.lengthYAxis_cm;
morphology = struct();
morphology.finalActivePrecipitateArea_cm2 = activeAreaCm2;
morphology.finalPrecipitateCentroidY_cm = centroidY;
morphology.finalPrecipitateCentroidDistanceFromSplit_cm = distanceFromSplit;
morphology.centerBandMorphologyTolerance_cm = toleranceCm;
morphology.centerBandMorphologyAccepted = activeAreaCm2 > 0 && ...
    isfinite(distanceFromSplit) && distanceFromSplit <= toleranceCm;
end
