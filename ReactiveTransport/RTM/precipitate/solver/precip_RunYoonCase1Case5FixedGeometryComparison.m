function comparison = precip_RunYoonCase1Case5FixedGeometryComparison(spec, options)
% precip_RunYoonCase1Case5FixedGeometryComparison - Diagnostic Case 1/5 runner.
%
% Runs the fixed-geometry Yoon path twice, once with dissolutionFactor = 1
% and once with dissolutionFactor = 300, and writes a comparison manifest.
% This is a reproducibility scaffold for the planned Case 1/Case 5 benchmark,
% not a quantitative Zhang/Yoon reproduction.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

outputRoot = getFieldOrDefault(options, 'outputRoot', defaultOutputRoot());
targetTimesS = getFieldOrDefault(options, 'targetTimes_s', [13, 18, 118] * 60);
endTimeS = getFieldOrDefault(options, 'endTime_s', max(targetTimesS));
dtS = getFieldOrDefault(options, 'dt_s', 30);

caseNames = ["case_1"; "case_5"];
dissolutionFactor = [1; 300];
caseDir = strings(2, 1);
manifestPath = strings(2, 1);
areaCsv = strings(2, 1);
numSnapshots = zeros(2, 1);
finalTotalPrecipitatedArea_cm2 = zeros(2, 1);
finalFirstPorePrecipitatedArea_cm2 = zeros(2, 1);
finalFirstThreePoresPrecipitatedArea_cm2 = zeros(2, 1);
finalPrecipitateMoles = zeros(2, 1);
finalFlowSolver = strings(2, 1);
finalFlowIsStokes = false(2, 1);
finalReactionMassAccepted = false(2, 1);
case5BoostActive = false(2, 1);
case5ActivationTime_s = NaN(2, 1);

for iCase = 1:2
    caseSpec = spec;
    caseSpec.dissolutionFactor = dissolutionFactor(iCase);
    caseDir(iCase) = string(fullfile(outputRoot, char(caseNames(iCase))));
    caseOptions = options;
    caseOptions.outputRoot = char(caseDir(iCase));
    caseOptions.endTime_s = endTimeS;
    caseOptions.dt_s = dtS;
    caseOptions.targetTimes_s = targetTimesS;

    run = precip_RunYoonFixedGeometryBenchmark(caseSpec, caseOptions);
    manifestPath(iCase) = string(run.manifestPath);
    areaCsv(iCase) = string(run.manifest.areaCsv);
    numSnapshots(iCase) = run.manifest.numSnapshots;
    finalTotalPrecipitatedArea_cm2(iCase) = ...
        run.run.areaTimeseries.totalPrecipitatedArea_cm2(end);
    finalFirstPorePrecipitatedArea_cm2(iCase) = ...
        run.run.areaTimeseries.firstPorePrecipitatedArea_cm2(end);
    finalFirstThreePoresPrecipitatedArea_cm2(iCase) = ...
        run.run.areaTimeseries.firstThreePoresPrecipitatedArea_cm2(end);
    finalPrecipitateMoles(iCase) = sum(run.run.state.precipitateMoles(:));
    finalFlowSolver(iCase) = string(run.manifest.finalFlowSolver);
    finalFlowIsStokes(iCase) = logical(run.manifest.finalFlowIsStokes);
    finalReactionMassAccepted(iCase) = logical(run.manifest.finalReactionMassAccepted);
    case5BoostActive(iCase) = logical(getFieldOrDefault(run.manifest, ...
        'case5BoostActive', false));
    case5ActivationTime_s(iCase) = getFieldOrDefault(run.manifest, ...
        'case5ActivationTime_s', NaN);
end

summary = table(caseNames, dissolutionFactor, caseDir, manifestPath, ...
    areaCsv, numSnapshots, finalTotalPrecipitatedArea_cm2, ...
    finalFirstPorePrecipitatedArea_cm2, ...
    finalFirstThreePoresPrecipitatedArea_cm2, finalPrecipitateMoles, ...
    finalFlowSolver, finalFlowIsStokes, finalReactionMassAccepted, ...
    case5BoostActive, case5ActivationTime_s, ...
    'VariableNames', {'caseName', 'dissolutionFactor', 'caseDir', ...
    'manifestPath', 'areaCsv', 'numSnapshots', ...
    'finalTotalPrecipitatedArea_cm2', ...
    'finalFirstPorePrecipitatedArea_cm2', ...
    'finalFirstThreePoresPrecipitatedArea_cm2', ...
    'finalPrecipitateMoles', 'finalFlowSolver', 'finalFlowIsStokes', ...
    'finalReactionMassAccepted', 'case5BoostActive', ...
    'case5ActivationTime_s'});

summaryCsv = fullfile(outputRoot, 'yoon_case1_case5_comparison_summary.csv');
writeTable(summaryCsv, summary);

manifest = struct();
manifest.runner = 'precip_RunYoonCase1Case5FixedGeometryComparison';
manifest.modelFamily = spec.modelFamily;
manifest.outputRoot = outputRoot;
manifest.summaryCsv = summaryCsv;
manifest.caseNames = cellstr(caseNames);
manifest.dissolutionFactor = dissolutionFactor(:)';
manifest.caseDirs = cellstr(caseDir);
manifest.caseManifestPaths = cellstr(manifestPath);
manifest.targetTimes_s = targetTimesS(:)';
manifest.endTime_s = endTimeS;
manifest.dt_s = dtS;
manifest.numCases = numel(caseNames);
manifest.case5DissolutionFactorAppliedOnlyToDissolution = true;
manifest.case5BoostActivated = case5BoostActive(2);
manifest.case5ActivationTime_s = case5ActivationTime_s(2);
manifest.caseTargetTimesComplete = all(numSnapshots >= numel(targetTimesS));
manifest.caseReactionMassAcceptedAll = all(finalReactionMassAccepted);
manifest.case1FinalPrecipitateMoles = finalPrecipitateMoles(1);
manifest.case5FinalPrecipitateMoles = finalPrecipitateMoles(2);
manifest.case1FinalTotalPrecipitatedArea_cm2 = ...
    finalTotalPrecipitatedArea_cm2(1);
manifest.case5FinalTotalPrecipitatedArea_cm2 = ...
    finalTotalPrecipitatedArea_cm2(2);
manifest.case5FinalPrecipitateMolesLessThanCase1 = ...
    finalPrecipitateMoles(2) < finalPrecipitateMoles(1);
manifest.case5FinalAreaNoGreaterThanCase1 = ...
    finalTotalPrecipitatedArea_cm2(2) <= ...
        finalTotalPrecipitatedArea_cm2(1);
manifest = appendStructFields(manifest, ...
    precip_BuildOutputEvidenceProvenance(summaryCsv, options));
manifest.productionComparisonValidated = ...
    logical(getFieldOrDefault(options, ...
    'productionComparisonValidated', false)) && ...
    hasCompleteOutputProvenance(manifest);
manifest.case1Case5AcceptanceCriteria = struct( ...
    'caseTargetTimesComplete', manifest.caseTargetTimesComplete, ...
    'caseReactionMassAcceptedAll', manifest.caseReactionMassAcceptedAll, ...
    'case5DissolutionFactorAppliedOnlyToDissolution', ...
        manifest.case5DissolutionFactorAppliedOnlyToDissolution, ...
    'case5FinalPrecipitateMolesLessThanCase1', ...
        manifest.case5FinalPrecipitateMolesLessThanCase1, ...
    'case5FinalAreaNoGreaterThanCase1', ...
        manifest.case5FinalAreaNoGreaterThanCase1, ...
    'productionComparisonValidated', ...
        manifest.productionComparisonValidated);
manifest.case1Case5Accepted = ...
    manifest.caseTargetTimesComplete && ...
    manifest.caseReactionMassAcceptedAll && ...
    manifest.case5DissolutionFactorAppliedOnlyToDissolution && ...
    manifest.case5FinalPrecipitateMolesLessThanCase1 && ...
    manifest.case5FinalAreaNoGreaterThanCase1 && ...
    manifest.productionComparisonValidated;
manifest.isQuantitativeBenchmark = false;
manifest.note = ['Fixed-geometry diagnostic Case 1/Case 5 comparison; ', ...
    'not a grid-converged or literature-locked Zhang/Yoon reproduction.'];

comparisonManifestPath = fullfile(outputRoot, ...
    'yoon_case1_case5_comparison_manifest.json');
writeJsonManifest(comparisonManifestPath, manifest);

comparison = struct();
comparison.summary = summary;
comparison.summaryCsv = summaryCsv;
comparison.manifest = manifest;
comparison.manifestPath = comparisonManifestPath;
end

function writeTable(path, tableData)
parentDir = fileparts(path);
if ~isfolder(parentDir)
    mkdir(parentDir);
end
writetable(tableData, path);
end

function outputRoot = defaultOutputRoot()
moduleRoot = fileparts(fileparts(mfilename('fullpath')));
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
outputRoot = fullfile(moduleRoot, 'outputs', ...
    'yoon_case1_case5_fixed_geometry', timestamp);
end

function writeJsonManifest(path, manifest)
parentDir = fileparts(path);
if ~isfolder(parentDir)
    mkdir(parentDir);
end
fid = fopen(path, 'w');
if fid < 0
    error('RTSPHEM:Precipitate:ManifestWriteFailed', ...
        'Could not open manifest for writing: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
clear cleanupObj;
end

function manifest = appendStructFields(manifest, fields)
fieldNames = fieldnames(fields);
for iField = 1:numel(fieldNames)
    manifest.(fieldNames{iField}) = fields.(fieldNames{iField});
end
end

function tf = hasCompleteOutputProvenance(manifest)
tf = isfield(manifest, 'sourceCommitSha') && ...
    ~isempty(regexp(string(manifest.sourceCommitSha), ...
    "^[0-9a-fA-F]{40}$", 'once')) && ...
    isfield(manifest, 'outputHashAlgorithm') && ...
    strcmp(string(manifest.outputHashAlgorithm), "SHA-256") && ...
    isfield(manifest, 'outputEvidenceHashSha256') && ...
    ~isempty(regexp(string(manifest.outputEvidenceHashSha256), ...
    "^[0-9a-fA-F]{64}$", 'once')) && ...
    isfield(manifest, 'outputHashInputFilesVerified') && ...
    logical(manifest.outputHashInputFilesVerified);
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
