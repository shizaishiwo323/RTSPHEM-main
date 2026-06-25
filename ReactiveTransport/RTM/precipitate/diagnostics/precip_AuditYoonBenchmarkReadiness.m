function audit = precip_AuditYoonBenchmarkReadiness(evidence)
% precip_AuditYoonBenchmarkReadiness - Gate quantitative Zhang/Yoon claims.
%
% The audit consumes already-collected evidence structs and reports whether
% each hard requirement for a quantitative benchmark claim is satisfied. It is
% intentionally conservative: diagnostic smoke manifests should fail this gate.

if nargin < 1 || ~isstruct(evidence)
    error('RTSPHEM:Precipitate:InvalidReadinessEvidence', ...
        'evidence must be a struct.');
end

fixedManifest = getStructField(evidence, 'fixedGeometryManifest');
gridManifest = getStructField(evidence, 'gridConvergenceManifest');
reference = getStructField(evidence, 'reference');
chemistryManifest = getStructField(evidence, 'chemistryManifest');
passiveTransportManifest = getStructField(evidence, ...
    'passiveTransportManifest');
caseComparisonManifest = getStructField(evidence, ...
    'caseComparisonManifest');
diffusionFeedbackManifest = getStructField(evidence, ...
    'diffusionFeedbackManifest');

rows = {
    'target_times', ...
    hasRequiredTargetTimes(fixedManifest), ...
    'Fixed-geometry run includes 13, 18, and 118 min targets.'
    'geometry_package', ...
    hasAcceptedGeometryPackage(fixedManifest), ...
    'Calibrated geometry package and exact region masks are quantitative.'
    'reference_package', ...
    hasAcceptedReferencePackage(reference), ...
    'Reference curves have complete source/digitization provenance.'
    'grid_convergence', ...
    hasAcceptedGridConvergence(gridManifest), ...
    'Production 10/5/2.5 um grid-convergence evidence is quantitative.'
    'production_stokes', ...
    hasAcceptedProductionStokes(fixedManifest), ...
    'Flow feedback uses a production-validated Stokes backend.'
    'flow_feedback', ...
    hasAcceptedFlowFeedback(fixedManifest), ...
    'Vm >= 0.6 topology changes trigger accepted Stokes flow recomputation.'
    'reaction_mass', ...
    hasAcceptedReactionMass(fixedManifest), ...
    'Final reaction mass ledger is accepted.'
    'center_band_morphology', ...
    hasAcceptedCenterBandMorphology(fixedManifest), ...
    'Precipitation morphology is centered on the split-inlet mixing band.'
    'rate_constants_locked', ...
    hasLockedRateConstants(fixedManifest), ...
    'Vaterite rate constants are literature-locked for benchmark use.'
    'chemistry_speciation', ...
    hasAcceptedChemistrySpeciation(chemistryManifest), ...
    'PHREEQC/Yoon speciation cross-validation is accepted.'
    'passive_transport', ...
    hasAcceptedPassiveTransport(passiveTransportManifest), ...
    'No-reaction split-inlet conservative transport checks are accepted.'
    'case1_case5', ...
    hasAcceptedCaseComparison(caseComparisonManifest), ...
    'Case 1/Case 5 comparison evidence is accepted.'
    'diffusion_feedback', ...
    hasAcceptedDiffusionFeedback(diffusionFeedbackManifest), ...
    'Diffusion-feedback Case 2/1/3 sensitivity evidence is accepted.'
    'no_finite_clipping', ...
    hasNoFiniteClippingEvidence(fixedManifest), ...
    'Benchmark evidence explicitly avoids finite concentration clipping.'
    'quantitative_manifest', ...
    hasAcceptedQuantitativeManifest(fixedManifest), ...
    'Benchmark manifest is quantitative and carries required output evidence.'
    };

requirementId = strings(size(rows, 1), 1);
ready = false(size(rows, 1), 1);
detail = strings(size(rows, 1), 1);
for iRow = 1:size(rows, 1)
    requirementId(iRow) = string(rows{iRow, 1});
    ready(iRow) = logical(rows{iRow, 2});
    detail(iRow) = string(rows{iRow, 3});
end

audit = struct();
audit.requirements = table(requirementId, ready, detail);
audit.allReady = all(ready);
audit.isQuantitativeBenchmarkAllowed = audit.allReady;
if audit.allReady
    audit.note = 'All audited Zhang/Yoon quantitative-readiness gates pass.';
else
    audit.note = 'One or more Zhang/Yoon quantitative-readiness gates fail.';
end
end

function s = getStructField(parent, fieldName)
if isfield(parent, fieldName) && isstruct(parent.(fieldName))
    s = parent.(fieldName);
else
    s = struct();
end
end

function value = getLogicalField(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = logical(s.(fieldName));
else
    value = defaultValue;
end
end

function tf = hasRequiredTargetTimes(manifest)
requiredTimes = [13, 18, 118] * 60;
if ~isfield(manifest, 'targetTimes_s') || isempty(manifest.targetTimes_s)
    tf = false;
    return;
end
targetTimes = double(manifest.targetTimes_s(:)');
tf = all(ismembertol(requiredTimes, targetTimes, 1e-9));
if tf && isfield(manifest, 'targetSnapshotsComplete') && ...
        ~isempty(manifest.targetSnapshotsComplete)
    tf = logical(manifest.targetSnapshotsComplete);
else
    tf = false;
end
if ~tf || ~isfield(manifest, 'capturedTargetTimes_s') || ...
        isempty(manifest.capturedTargetTimes_s)
    tf = false;
    return;
end

capturedTimes = double(manifest.capturedTargetTimes_s(:)');
tf = all(ismembertol(requiredTimes, capturedTimes, 1e-9));
if ~tf
    return;
end

if isfield(manifest, 'missingTargetTimes_s') && ...
        ~isempty(manifest.missingTargetTimes_s)
    tf = false;
    return;
end

[numSnapshots, hasNumSnapshots] = getFiniteScalarDouble(manifest, ...
    'numSnapshots');
tf = hasNumSnapshots && numSnapshots >= numel(requiredTimes);
end

function tf = hasTargetGridSequence(manifest)
requiredSpacing = [10, 5, 2.5];
if ~isfield(manifest, 'requestedTargetGridSpacing_um') || ...
        isempty(manifest.requestedTargetGridSpacing_um)
    tf = false;
    return;
end
spacing = double(manifest.requestedTargetGridSpacing_um(:)');
tf = all(ismembertol(requiredSpacing, spacing, 1e-9));
end

function tf = hasAcceptedGridConvergence(manifest)
tf = hasTargetGridSequence(manifest) && ...
    getLogicalField(manifest, 'gridConvergenceAccepted', false);
if ~tf
    return;
end

[maxRelativeDifference, hasMaxDifference] = getFiniteScalarDouble( ...
    manifest, 'maxRelativeTotalAreaDifferenceFromFinest');
[tolerance, hasTolerance] = getFiniteScalarDouble( ...
    manifest, 'gridConvergenceTolerance');

tf = getLogicalField(manifest, 'targetGridSpacingSequenceComplete', false) && ...
    getLogicalField(manifest, 'actualGridSpacingSequenceComplete', false) && ...
    getLogicalField(manifest, 'gridConvergenceWithinTolerance', false) && ...
    getLogicalField(manifest, 'productionGridConvergenceValidated', false) && ...
    hasMaxDifference && hasTolerance && tolerance >= 0 && ...
    maxRelativeDifference <= tolerance && hasActualGridEvidence(manifest);
end

function tf = hasActualGridEvidence(manifest)
requiredSpacing = [10, 5, 2.5];
[dx, hasDx] = getFiniteNumericVector(manifest, 'actualDx_um', 3);
[dy, hasDy] = getFiniteNumericVector(manifest, 'actualDy_um', 3);
[numX, hasNumX] = getFiniteNumericVector(manifest, 'actualNumX', 3);
[numY, hasNumY] = getFiniteNumericVector(manifest, 'actualNumY', 3);
tf = hasDx && hasDy && hasNumX && hasNumY && ...
    all(dx > 0) && all(dy > 0) && all(numX > 0) && all(numY > 0) && ...
    hasRequiredSpacingValues(dx, requiredSpacing) && ...
    hasRequiredSpacingValues(dy, requiredSpacing) && ...
    gridCountsIncreaseAsSpacingDecreases(dx, numX) && ...
    gridCountsIncreaseAsSpacingDecreases(dy, numY);
end

function tf = hasRequiredSpacingValues(spacing, requiredSpacing)
tf = all(ismembertol(double(requiredSpacing(:)'), double(spacing(:)'), 1e-9));
end

function tf = gridCountsIncreaseAsSpacingDecreases(spacing, counts)
[sortedSpacing, order] = sort(double(spacing(:)), 'descend');
sortedCounts = double(counts(order));
tf = all(diff(sortedSpacing) < 0) && all(diff(sortedCounts) > 0);
end

function tf = hasAcceptedGeometryPackage(manifest)
tf = getLogicalField(manifest, 'geometryPackageIsQuantitative', false) && ...
    hasNonemptyTextField(manifest, 'geometryPackageDir') && ...
    hasNonemptyTextField(manifest, 'geometryPackageName') && ...
    hasNonemptyTextField(manifest, 'geometryPackageNote') && ...
    hasNonemptyTextField(manifest, 'substrateMaskFile') && ...
    hasNonemptyTextField(manifest, 'regionMaskFile');
if ~tf
    return;
end

note = lower(string(manifest.geometryPackageNote));
tf = ~contains(note, "missing") && ~contains(note, "incomplete");
if ~tf
    return;
end

[numSubstrateCells, hasSubstrateCells] = getFiniteScalarDouble( ...
    manifest, 'geometryPackageNumSubstrateCells');
[numFirstPoreCells, hasFirstPoreCells] = getFiniteScalarDouble( ...
    manifest, 'geometryPackageNumFirstPoreCells');
[numFirstThreePoresCells, hasFirstThreePoresCells] = getFiniteScalarDouble( ...
    manifest, 'geometryPackageNumFirstThreePoresCells');
tf = getLogicalField(manifest, ...
    'geometryPackageAssetFilesVerified', false) && ...
    hasSubstrateCells && hasFirstPoreCells && hasFirstThreePoresCells && ...
    numSubstrateCells > 0 && numFirstPoreCells > 0 && ...
    numFirstThreePoresCells >= numFirstPoreCells;
end

function tf = hasAcceptedReferencePackage(reference)
tf = getLogicalField(reference, 'isQuantitativeBenchmark', false) && ...
    isfield(reference, 'digitizationPackage') && ...
    isstruct(reference.digitizationPackage);
if ~tf
    return;
end

package = reference.digitizationPackage;
requiredFields = {'packageDir', 'manifestPath', 'packageName', ...
    'sourceFigure', 'screenshotFile', 'webPlotDigitizerProjectFile', ...
    'calibrationCsv', 'rawExportCsv', 'conversionScript', ...
    'uncertaintyCsv', 'convertedReferenceCsv', 'note'};
for iField = 1:numel(requiredFields)
    if ~hasNonemptyTextField(package, requiredFields{iField})
        tf = false;
        return;
    end
end

if isfield(package, 'isQuantitativeBenchmark') && ...
        ~isempty(package.isQuantitativeBenchmark)
    tf = tf && logical(package.isQuantitativeBenchmark);
else
    tf = false;
end
tf = tf && isfield(package, 'missingAssets') && ...
    isempty(package.missingAssets);
[numReferenceRows, hasReferenceRows] = getFiniteScalarDouble(package, ...
    'numReferenceRows');
tf = tf && getLogicalField(package, 'assetFilesVerified', false) && ...
    hasReferenceRows && numReferenceRows > 0;
note = lower(string(package.note));
tf = tf && ~contains(note, "missing") && ~contains(note, "incomplete");
end

function tf = hasAcceptedFlowFeedback(manifest)
tf = getLogicalField(manifest, 'flowFeedbackAccepted', false) && ...
    getLogicalField(manifest, 'flowTopologyChangedAny', false) && ...
    getLogicalField(manifest, 'flowRecomputedAfterTopologyChange', false) && ...
    hasAcceptedProductionStokes(manifest);
if ~tf
    return;
end

[topologyChangedSteps, hasTopologyCount] = getFiniteScalarDouble( ...
    manifest, 'totalFlowTopologyChangedSteps');
[flowRecomputations, hasRecomputeCount] = getFiniteScalarDouble( ...
    manifest, 'finalNumFlowRecomputations');
[relativePermeability, hasRelativePermeability] = getFiniteScalarDouble( ...
    manifest, 'finalRelativePermeability');
[pressureDropRelative, hasPressureDropRelative] = getFiniteScalarDouble( ...
    manifest, 'finalPressureDropRelative');
[flowRateRelative, hasFlowRateRelative] = getFiniteScalarDouble( ...
    manifest, 'finalFlowRateRelative');

tf = hasTopologyCount && hasRecomputeCount && ...
    topologyChangedSteps >= 1 && ...
    flowRecomputations >= topologyChangedSteps && ...
    hasRelativePermeability && hasPressureDropRelative && ...
    hasFlowRateRelative && relativePermeability >= 0 && ...
    relativePermeability <= 1 && pressureDropRelative >= 1 && ...
    flowRateRelative >= 0;
end

function tf = hasAcceptedProductionStokes(manifest)
tf = getLogicalField(manifest, 'finalFlowIsStokes', false) && ...
    ~getLogicalField(manifest, 'finalFlowIsProxy', true) && ...
    getLogicalField(manifest, 'flowBackendProductionValidated', false) && ...
    hasNonemptyTextField(manifest, 'finalFlowSolver');
if ~tf
    return;
end

solverName = lower(string(manifest.finalFlowSolver));
tf = contains(solverName, "stokes") && ...
    ~contains(solverName, "proxy") && ...
    ~contains(solverName, "darcy");
if ~tf
    return;
end

[linearResidual, hasLinearResidual] = getFiniteScalarDouble( ...
    manifest, 'finalFlowLinearResidualRelative');
[maxAcceptedResidual, hasResidualThreshold] = getFiniteScalarDouble( ...
    manifest, 'maxAcceptedFlowLinearResidualRelative');
tf = hasLinearResidual && hasResidualThreshold && ...
    maxAcceptedResidual >= 0 && linearResidual <= maxAcceptedResidual;
end

function tf = hasAcceptedReactionMass(manifest)
tf = getLogicalField(manifest, 'finalReactionMassAccepted', false);
if ~tf
    return;
end

[caError, hasCaError] = getFiniteScalarDouble( ...
    manifest, 'finalReactionTotalCaRelativeError');
[cError, hasCError] = getFiniteScalarDouble( ...
    manifest, 'finalReactionTotalCRelativeError');
[alkalinityError, hasAlkalinityError] = getFiniteScalarDouble( ...
    manifest, 'finalReactionTotalAlkalinityRelativeError');

maxAcceptedRelativeError = 1e-10;
tf = hasCaError && hasCError && hasAlkalinityError && ...
    max(abs([caError, cError, alkalinityError])) <= ...
        maxAcceptedRelativeError;
end

function tf = hasAcceptedCenterBandMorphology(manifest)
tf = getLogicalField(manifest, 'centerBandMorphologyAccepted', false);
if ~tf
    return;
end

[activeArea, hasActiveArea] = getFiniteScalarDouble( ...
    manifest, 'finalActivePrecipitateArea_cm2');
[centroidY, hasCentroidY] = getFiniteScalarDouble( ...
    manifest, 'finalPrecipitateCentroidY_cm');
[distanceFromSplit, hasDistance] = getFiniteScalarDouble( ...
    manifest, 'finalPrecipitateCentroidDistanceFromSplit_cm');
[tolerance, hasTolerance] = getFiniteScalarDouble( ...
    manifest, 'centerBandMorphologyTolerance_cm');

tf = hasActiveArea && hasCentroidY && hasDistance && hasTolerance && ...
    activeArea > 0 && tolerance >= 0 && distanceFromSplit <= tolerance;
end

function tf = hasAcceptedCaseComparison(manifest)
tf = getLogicalField(manifest, 'case1Case5Accepted', false) && ...
    getLogicalField(manifest, 'caseTargetTimesComplete', false) && ...
    getLogicalField(manifest, 'caseReactionMassAcceptedAll', false) && ...
    getLogicalField(manifest, ...
        'case5DissolutionFactorAppliedOnlyToDissolution', false) && ...
    getLogicalField(manifest, ...
        'case5FinalPrecipitateMolesLessThanCase1', false) && ...
    getLogicalField(manifest, 'case5FinalAreaNoGreaterThanCase1', false) && ...
    getLogicalField(manifest, 'productionComparisonValidated', false);
if ~tf
    return;
end

[case1Moles, hasCase1Moles] = getFiniteScalarDouble( ...
    manifest, 'case1FinalPrecipitateMoles');
[case5Moles, hasCase5Moles] = getFiniteScalarDouble( ...
    manifest, 'case5FinalPrecipitateMoles');
[case1Area, hasCase1Area] = getFiniteScalarDouble( ...
    manifest, 'case1FinalTotalPrecipitatedArea_cm2');
[case5Area, hasCase5Area] = getFiniteScalarDouble( ...
    manifest, 'case5FinalTotalPrecipitatedArea_cm2');

tf = hasCase1Moles && hasCase5Moles && hasCase1Area && hasCase5Area && ...
    case1Moles >= 0 && case5Moles >= 0 && ...
    case1Area >= 0 && case5Area >= 0 && ...
    case5Moles < case1Moles && case5Area <= case1Area;
end

function tf = hasAcceptedDiffusionFeedback(manifest)
tf = getLogicalField(manifest, 'diffusionFeedbackAccepted', false) && ...
    hasDiffusionExponentSequence(manifest) && ...
    getLogicalField(manifest, 'areaTrendNonincreasing', false) && ...
    getLogicalField(manifest, 'massTrendDecreasing', false);

if ~tf
    return;
end

[areaValues, hasAreaValues] = getFiniteNumericVector( ...
    manifest, 'totalPrecipitatedArea_cm2', 3);
[moleValues, hasMoleValues] = getFiniteNumericVector( ...
    manifest, 'totalPrecipitateMoles', 3);
[diffusivityValues, hasDiffusivityValues] = getFiniteNumericVector( ...
    manifest, 'finalMeanEffectiveDiffusivity_cm2_s', 3);
[feedbackValues, hasFeedbackValues] = getLogicalVector( ...
    manifest, 'feedbackCoupled', 3);

tf = hasAreaValues && hasMoleValues && hasDiffusivityValues && ...
    hasFeedbackValues && all(areaValues >= 0) && all(moleValues >= 0) && ...
    all(diffusivityValues >= 0) && all(feedbackValues) && ...
    isNonincreasingNumeric(areaValues) && isStrictlyDecreasingNumeric(moleValues);
end

function tf = hasAcceptedChemistrySpeciation(manifest)
tf = getLogicalField(manifest, 'chemistrySpeciationAccepted', false) && ...
    getLogicalField(manifest, 'hasIphreeqc', false) && ...
    getLogicalField(manifest, 'isQuantitativeAcceptance', false);
if ~tf
    return;
end

[maxPhDifference, hasPhDifference] = getFiniteScalarDouble( ...
    manifest, 'maxAbsPhDifference');
[maxSiDifference, hasSiDifference] = getFiniteScalarDouble( ...
    manifest, 'maxAbsSiVateriteDifference');
[maxAcceptedPhDifference, hasPhThreshold] = getFiniteScalarDouble( ...
    manifest, 'maxAcceptedPhDifference');
[maxAcceptedSiDifference, hasSiThreshold] = getFiniteScalarDouble( ...
    manifest, 'maxAcceptedSiVateriteDifference');
[expectedInletAPh, hasExpectedInletAPh] = getFiniteScalarDouble( ...
    manifest, 'expectedInletAPh');
[expectedInletBPh, hasExpectedInletBPh] = getFiniteScalarDouble( ...
    manifest, 'expectedInletBPh');
[inletAPhIphreeqc, hasInletAPh] = getFiniteScalarDouble( ...
    manifest, 'inletAPhIphreeqc');
[inletBPhIphreeqc, hasInletBPh] = getFiniteScalarDouble( ...
    manifest, 'inletBPhIphreeqc');
[maxInletPhError, hasInletPhError] = getFiniteScalarDouble( ...
    manifest, 'maxAbsInletPhError');
[maxAcceptedInletPhError, hasInletPhThreshold] = getFiniteScalarDouble( ...
    manifest, 'maxAcceptedInletPhError');
[yoonPeakFraction, hasYoonPeak] = getFiniteScalarDouble( ...
    manifest, 'yoonPeakFractionInletA');
[iphreeqcPeakFraction, hasIphreeqcPeak] = getFiniteScalarDouble( ...
    manifest, 'iphreeqcPeakFractionInletA');

tf = hasPhDifference && hasSiDifference && ...
    hasPhThreshold && hasSiThreshold && hasExpectedInletAPh && ...
    hasExpectedInletBPh && hasInletAPh && hasInletBPh && ...
    hasInletPhError && hasInletPhThreshold && hasYoonPeak && ...
    hasIphreeqcPeak && maxAcceptedPhDifference >= 0 && ...
    maxAcceptedSiDifference >= 0 && maxAcceptedInletPhError >= 0 && ...
    maxPhDifference <= maxAcceptedPhDifference && ...
    maxSiDifference <= maxAcceptedSiDifference && ...
    abs(inletAPhIphreeqc - expectedInletAPh) <= maxAcceptedInletPhError && ...
    abs(inletBPhIphreeqc - expectedInletBPh) <= maxAcceptedInletPhError && ...
    maxInletPhError <= maxAcceptedInletPhError && ...
    getLogicalField(manifest, 'inletPhAccepted', false) && ...
    yoonPeakFraction > 0 && yoonPeakFraction < 1 && ...
    iphreeqcPeakFraction > 0 && iphreeqcPeakFraction < 1 && ...
    getLogicalField(manifest, 'yoonPeakFractionInternal', false) && ...
    getLogicalField(manifest, 'iphreeqcPeakFractionInternal', false);
end

function tf = hasAcceptedPassiveTransport(manifest)
tf = getLogicalField(manifest, 'passiveTransportAccepted', false);
if ~tf || ~isfield(manifest, 'usedFiniteConcentrationLimiter') || ...
        isempty(manifest.usedFiniteConcentrationLimiter) || ...
        logical(manifest.usedFiniteConcentrationLimiter)
    tf = false;
    return;
end

[minConcentration, hasMinConcentration] = getFiniteScalarDouble( ...
    manifest, 'minConcentration');
[maxMassError, hasMassError] = getFiniteScalarDouble( ...
    manifest, 'maxAbsInertRelativeMassError');
[maxAcceptedMassError, hasMassThreshold] = getFiniteScalarDouble( ...
    manifest, 'maxAcceptedMassRelativeError');
[mixingError, hasMixingError] = getFiniteScalarDouble( ...
    manifest, 'mixingSymmetryError');
[maxAcceptedMixingError, hasMixingThreshold] = getFiniteScalarDouble( ...
    manifest, 'maxAcceptedMixingSymmetryError');
[rejectedStepCount, hasRejectedCount] = getFiniteScalarDouble( ...
    manifest, 'rejectedStepCount');
[maxAcceptedRejectedCount, hasRejectedThreshold] = getFiniteScalarDouble( ...
    manifest, 'maxAcceptedRejectedStepCount');

tf = hasMinConcentration && hasMassError && hasMassThreshold && ...
    hasMixingError && hasMixingThreshold && hasRejectedCount && ...
    hasRejectedThreshold && minConcentration >= -1e-15 && ...
    maxMassError <= maxAcceptedMassError && ...
    mixingError <= maxAcceptedMixingError && ...
    rejectedStepCount <= maxAcceptedRejectedCount && ...
    hasRequiredConservativeComponents(manifest);
end

function tf = hasRequiredConservativeComponents(manifest)
requiredComponents = ["Ca_total", "C_total", "Na_total", "Cl_total", ...
    "Alkalinity"];
forbiddenComponents = ["H_total", "H+", "hTransport", "free_H"];
tf = isfield(manifest, 'componentNames') && ~isempty(manifest.componentNames);
if ~tf
    return;
end
componentNames = string(manifest.componentNames(:));
tf = all(ismember(requiredComponents, componentNames)) && ...
    ~any(ismember(forbiddenComponents, componentNames));
end

function tf = hasAcceptedQuantitativeManifest(manifest)
tf = getLogicalField(manifest, 'isQuantitativeBenchmark', false);
if ~tf
    return;
end

tf = getLogicalField(manifest, 'outputEvidenceFilesVerified', false);
if ~tf
    return;
end

requiredPathFields = {'outputRoot', 'snapshotDir', 'areaCsv', ...
    'flowDiagnosticsCsv', 'transportDiagnosticsCsv', ...
    'reactionMassLedgerCsv', 'reactionDiagnosticsCsv'};
for iField = 1:numel(requiredPathFields)
    if ~hasNonemptyTextField(manifest, requiredPathFields{iField})
        tf = false;
        return;
    end
end

tf = hasNonemptyTextListField(manifest, 'matFiles');
end

function tf = hasDiffusionExponentSequence(manifest)
requiredExponent = [0, 2, 3];
if ~isfield(manifest, 'diffusionExponent') || ...
        isempty(manifest.diffusionExponent)
    tf = false;
    return;
end
exponent = double(manifest.diffusionExponent(:)');
tf = numel(exponent) == numel(requiredExponent) && ...
    all(abs(exponent - requiredExponent) <= 1e-12);
end

function tf = hasNoFiniteClippingEvidence(manifest)
if isfield(manifest, 'usesFiniteConcentrationLimiter') && ...
        ~isempty(manifest.usesFiniteConcentrationLimiter)
    tf = ~logical(manifest.usesFiniteConcentrationLimiter);
else
    tf = false;
end

if isfield(manifest, 'phreeqcTransportMaxFactor') && ...
        ~isempty(manifest.phreeqcTransportMaxFactor)
    tf = tf && isNonFiniteNumeric(manifest.phreeqcTransportMaxFactor);
end
end

function tf = hasLockedRateConstants(manifest)
tf = getLogicalField(manifest, 'yoonRateConstantsLocked', false);
if ~tf || ~isfield(manifest, 'yoonRateSource') || ...
        isempty(manifest.yoonRateSource)
    tf = false;
    return;
end
rateSource = lower(string(manifest.yoonRateSource));
tf = ~contains(rateSource, "smoke") && ~contains(rateSource, "pending");
if ~tf
    return;
end

[k1, hasK1] = getFiniteScalarDouble(manifest, 'yoonRateK1');
[k2, hasK2] = getFiniteScalarDouble(manifest, 'yoonRateK2');
[k3, hasK3] = getFiniteScalarDouble(manifest, 'yoonRateK3');
if isfield(manifest, 'yoonRateUnits') && ~isempty(manifest.yoonRateUnits)
    units = string(manifest.yoonRateUnits);
else
    units = "";
end

if isfield(manifest, 'yoonRateMineralPhase') && ...
        ~isempty(manifest.yoonRateMineralPhase)
    mineralPhase = string(manifest.yoonRateMineralPhase);
else
    mineralPhase = "";
end

tf = hasK1 && hasK2 && hasK3 && ...
    all([k1, k2, k3] >= 0) && ...
    strlength(units) > 0 && units == "mol_cm-2_s-1" && ...
    hasNonemptyTextField(manifest, 'yoonRateSourceDoi') && ...
    hasNonemptyTextField(manifest, 'yoonRateSourceEquation') && ...
    getLogicalField(manifest, 'yoonRateSourceValuesVerified', false) && ...
    mineralPhase == "Vaterite";
end

function tf = isNonFiniteNumeric(value)
if isnumeric(value)
    tf = all(~isfinite(double(value(:))));
else
    tf = false;
end
end

function [value, tf] = getFiniteScalarDouble(s, fieldName)
value = NaN;
tf = false;
if isfield(s, fieldName) && isnumeric(s.(fieldName)) && ...
        isscalar(s.(fieldName))
    value = double(s.(fieldName));
    tf = isfinite(value);
end
end

function [values, tf] = getFiniteNumericVector(s, fieldName, expectedLength)
values = [];
tf = false;
if isfield(s, fieldName) && isnumeric(s.(fieldName)) && ...
        numel(s.(fieldName)) == expectedLength
    values = double(s.(fieldName)(:));
    tf = all(isfinite(values));
end
end

function [values, tf] = getLogicalVector(s, fieldName, expectedLength)
values = [];
tf = false;
if isfield(s, fieldName) && numel(s.(fieldName)) == expectedLength
    rawValues = s.(fieldName)(:);
    if islogical(rawValues) || isnumeric(rawValues)
        values = logical(rawValues);
        tf = all(double(rawValues) == 0 | double(rawValues) == 1);
    end
end
end

function tf = isNonincreasingNumeric(values)
values = double(values(:));
tolerance = 10 * eps(max(1, max(abs(values))));
tf = all(diff(values) <= tolerance);
end

function tf = isStrictlyDecreasingNumeric(values)
values = double(values(:));
tolerance = 10 * eps(max(1, max(abs(values))));
tf = all(diff(values) < -tolerance);
end

function tf = hasNonemptyTextField(s, fieldName)
tf = isfield(s, fieldName) && ~isempty(s.(fieldName)) && ...
    strlength(strtrim(string(s.(fieldName)))) > 0;
end

function tf = hasNonemptyTextListField(s, fieldName)
tf = false;
if ~isfield(s, fieldName) || isempty(s.(fieldName))
    return;
end
values = string(s.(fieldName));
tf = ~isempty(values) && all(strlength(strtrim(values(:))) > 0);
end
