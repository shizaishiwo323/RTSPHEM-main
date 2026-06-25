function result = precip_RunYoonSpeciationMixingSeries(spec, options)
% precip_RunYoonSpeciationMixingSeries - Compare zero-dimensional mixing chemistry.
%
% Builds the Yoon split-inlet mixing series, evaluates the local carbonate
% equilibrium model, and optionally evaluates an injected PHREEQC speciation
% backend for side-by-side diagnostics.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2 || isempty(options)
    options = struct();
end

fractions = getOption(options, 'fractions', (0:0.01:1)');
fractions = fractions(:);
maxAcceptedPhDifference = getOption(options, 'maxAcceptedPhDifference', 0);
maxAcceptedSiVateriteDifference = getOption(options, ...
    'maxAcceptedSiVateriteDifference', 0);
maxAcceptedInletPhError = getOption(options, ...
    'maxAcceptedInletPhError', 0.05);
samples = precip_BuildYoonMixingSeries(spec, fractions);
yoonChem = precip_YoonCarbonateEquilibrium(samples, spec);

comparison = table();
comparison.fractionInletA = fractions;
comparison.pH_yoon = yoonChem.pH(:);
comparison.aH_yoon = yoonChem.aH(:);
comparison.aCa_yoon = yoonChem.aCa(:);
comparison.aCO3_yoon = yoonChem.aCO3(:);
comparison.aH2CO3_yoon = yoonChem.aH2CO3(:);
comparison.omegaVaterite_yoon = yoonChem.omegaVaterite(:);
comparison.siVaterite_yoon = yoonChem.saturationIndexVaterite(:);

hasIphreeqc = isfield(spec, 'iphreeqcSpeciationFcn') && ...
    ~isempty(spec.iphreeqcSpeciationFcn);
if hasIphreeqc
    [iphreeqcChem, iphreeqcMetadata] = spec.iphreeqcSpeciationFcn(samples, spec);
    comparison.pH_iphreeqc = iphreeqcChem.pH(:);
    comparison.aH_iphreeqc = iphreeqcChem.aH(:);
    comparison.aCa_iphreeqc = iphreeqcChem.aCa(:);
    comparison.aCO3_iphreeqc = iphreeqcChem.aCO3(:);
    comparison.aH2CO3_iphreeqc = iphreeqcChem.aH2CO3(:);
    comparison.omegaVaterite_iphreeqc = iphreeqcChem.omegaVaterite(:);
    comparison.siVaterite_iphreeqc = iphreeqcChem.saturationIndexVaterite(:);
    comparison.pHAbsDiff = abs(comparison.pH_iphreeqc - comparison.pH_yoon);
    comparison.siVateriteAbsDiff = abs(comparison.siVaterite_iphreeqc - ...
        comparison.siVaterite_yoon);
else
    iphreeqcMetadata = struct();
    comparison.pH_iphreeqc = nan(size(fractions));
    comparison.aH_iphreeqc = nan(size(fractions));
    comparison.aCa_iphreeqc = nan(size(fractions));
    comparison.aCO3_iphreeqc = nan(size(fractions));
    comparison.aH2CO3_iphreeqc = nan(size(fractions));
    comparison.omegaVaterite_iphreeqc = nan(size(fractions));
    comparison.siVaterite_iphreeqc = nan(size(fractions));
    comparison.pHAbsDiff = nan(size(fractions));
    comparison.siVateriteAbsDiff = nan(size(fractions));
end

summary = buildSummary(comparison, hasIphreeqc, ...
    maxAcceptedPhDifference, maxAcceptedSiVateriteDifference, ...
    maxAcceptedInletPhError, spec);

outputCsv = string(getOption(options, 'outputCsv', ""));
if strlength(outputCsv) > 0
    writetable(comparison, char(outputCsv));
end
outputManifestPath = string(getOption(options, 'outputManifestPath', ""));
if strlength(outputManifestPath) > 0
    writeJsonManifest(char(outputManifestPath), buildManifest(summary, ...
        outputCsv, fractions, hasIphreeqc));
end

result = struct();
result.table = comparison;
result.summary = summary;
result.samples = samples;
result.yoonMetadata = struct( ...
    'backend', 'yoon_equilibrium', ...
    'role', 'aqueous_speciation', ...
    'source', 'local_yoon_carbonate_equilibrium');
result.iphreeqcMetadata = iphreeqcMetadata;
result.outputCsv = outputCsv;
result.outputManifestPath = outputManifestPath;
end

function summary = buildSummary(comparison, hasIphreeqc, ...
    maxAcceptedPhDifference, maxAcceptedSiVateriteDifference, ...
    maxAcceptedInletPhError, spec)
summary = struct();
summary.numSamples = height(comparison);
summary.mixingFractionStep = inferUniformFractionStep( ...
    comparison.fractionInletA);
summary.mixingFractionsCover101 = hasDefault101FractionCoverage( ...
    comparison.fractionInletA);
[summary.yoonPeakOmegaVaterite, yoonPeakIdx] = ...
    max(comparison.omegaVaterite_yoon);
summary.yoonPeakFractionInletA = comparison.fractionInletA(yoonPeakIdx);
summary.maxAbsPhDifference = max(comparison.pHAbsDiff, [], 'omitnan');
summary.maxAbsSiVateriteDifference = max(comparison.siVateriteAbsDiff, ...
    [], 'omitnan');
summary.hasIphreeqc = hasIphreeqc;
summary.maxAcceptedPhDifference = maxAcceptedPhDifference;
summary.maxAcceptedSiVateriteDifference = maxAcceptedSiVateriteDifference;
summary.expectedInletAPh = spec.inletA.pH;
summary.expectedInletBPh = spec.inletB.pH;
summary.maxAcceptedInletPhError = maxAcceptedInletPhError;
summary.inletAPhYoon = getPhAtFraction(comparison, 'pH_yoon', 1);
summary.inletBPhYoon = getPhAtFraction(comparison, 'pH_yoon', 0);
if hasIphreeqc
    [summary.iphreeqcPeakOmegaVaterite, iphreeqcPeakIdx] = ...
        max(comparison.omegaVaterite_iphreeqc);
    summary.iphreeqcPeakFractionInletA = ...
        comparison.fractionInletA(iphreeqcPeakIdx);
    summary.inletAPhIphreeqc = getPhAtFraction(comparison, 'pH_iphreeqc', 1);
    summary.inletBPhIphreeqc = getPhAtFraction(comparison, 'pH_iphreeqc', 0);
else
    summary.iphreeqcPeakOmegaVaterite = NaN;
    summary.iphreeqcPeakFractionInletA = NaN;
    summary.inletAPhIphreeqc = NaN;
    summary.inletBPhIphreeqc = NaN;
end
inletPhErrors = abs([summary.inletAPhYoon, summary.inletBPhYoon, ...
    summary.inletAPhIphreeqc, summary.inletBPhIphreeqc] - ...
    [summary.expectedInletAPh, summary.expectedInletBPh, ...
    summary.expectedInletAPh, summary.expectedInletBPh]);
summary.maxAbsInletPhError = max(inletPhErrors, [], 'omitnan');
summary.inletPhAccepted = all(isfinite(inletPhErrors)) && ...
    summary.maxAbsInletPhError <= maxAcceptedInletPhError;
summary.yoonPeakFractionInternal = ...
    isInternalFraction(summary.yoonPeakFractionInletA);
summary.iphreeqcPeakFractionInternal = ...
    isInternalFraction(summary.iphreeqcPeakFractionInletA);
summary.isQuantitativeAcceptance = hasIphreeqc && ...
    isfinite(summary.maxAbsPhDifference) && ...
    isfinite(summary.maxAbsSiVateriteDifference) && ...
    summary.maxAbsPhDifference <= maxAcceptedPhDifference && ...
    summary.maxAbsSiVateriteDifference <= maxAcceptedSiVateriteDifference && ...
    summary.inletPhAccepted && summary.yoonPeakFractionInternal && ...
    summary.iphreeqcPeakFractionInternal;
summary.chemistrySpeciationAccepted = summary.isQuantitativeAcceptance;
end

function value = getPhAtFraction(comparison, phField, fraction)
value = NaN;
matches = find(abs(comparison.fractionInletA - fraction) <= 10 * eps, 1);
if ~isempty(matches) && ismember(phField, comparison.Properties.VariableNames)
    values = comparison.(phField);
    value = values(matches);
end
end

function tf = isInternalFraction(value)
tf = isfinite(value) && value > 0 && value < 1;
end

function step = inferUniformFractionStep(fractions)
if numel(fractions) < 2
    step = NaN;
    return;
end
d = diff(fractions(:));
if all(abs(d - d(1)) <= 1e-12)
    step = d(1);
else
    step = NaN;
end
end

function tf = hasDefault101FractionCoverage(fractions)
fractions = fractions(:);
tf = numel(fractions) == 101 && ...
    abs(fractions(1)) <= 1e-12 && ...
    abs(fractions(end) - 1) <= 1e-12 && ...
    all(abs(diff(fractions) - 0.01) <= 1e-12);
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function manifest = buildManifest(summary, outputCsv, fractions, hasIphreeqc)
manifest = struct();
manifest.runner = 'precip_RunYoonSpeciationMixingSeries';
manifest.numSamples = summary.numSamples;
manifest.fractions = fractions(:)';
manifest.mixingFractionStep = summary.mixingFractionStep;
manifest.mixingFractionsCover101 = summary.mixingFractionsCover101;
manifest.hasIphreeqc = hasIphreeqc;
manifest.outputCsv = char(outputCsv);
manifest.maxAbsPhDifference = summary.maxAbsPhDifference;
manifest.maxAbsSiVateriteDifference = summary.maxAbsSiVateriteDifference;
manifest.maxAcceptedPhDifference = summary.maxAcceptedPhDifference;
manifest.maxAcceptedSiVateriteDifference = ...
    summary.maxAcceptedSiVateriteDifference;
manifest.expectedInletAPh = summary.expectedInletAPh;
manifest.expectedInletBPh = summary.expectedInletBPh;
manifest.inletAPhYoon = summary.inletAPhYoon;
manifest.inletBPhYoon = summary.inletBPhYoon;
manifest.inletAPhIphreeqc = summary.inletAPhIphreeqc;
manifest.inletBPhIphreeqc = summary.inletBPhIphreeqc;
manifest.maxAbsInletPhError = summary.maxAbsInletPhError;
manifest.maxAcceptedInletPhError = summary.maxAcceptedInletPhError;
manifest.inletPhAccepted = summary.inletPhAccepted;
manifest.chemistrySpeciationAccepted = ...
    summary.chemistrySpeciationAccepted;
manifest.isQuantitativeAcceptance = summary.isQuantitativeAcceptance;
manifest.yoonPeakFractionInletA = summary.yoonPeakFractionInletA;
manifest.iphreeqcPeakFractionInletA = summary.iphreeqcPeakFractionInletA;
manifest.yoonPeakFractionInternal = summary.yoonPeakFractionInternal;
manifest.iphreeqcPeakFractionInternal = summary.iphreeqcPeakFractionInternal;
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
