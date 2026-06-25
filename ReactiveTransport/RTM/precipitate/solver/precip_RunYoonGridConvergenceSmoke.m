function result = precip_RunYoonGridConvergenceSmoke(spec, options)
% precip_RunYoonGridConvergenceSmoke - Multi-grid diagnostic wrapper.
%
% Runs the fixed-geometry Yoon diagnostic on several grid sizes and writes a
% machine-readable summary. This is a smoke convergence audit, not a
% quantitative Zhang/Yoon reproduction.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

if isfield(options, 'gridCases') && ~isempty(options.gridCases)
    gridCases = options.gridCases;
else
    targetGridSpacingUm = getFieldOrDefault(options, ...
        'targetGridSpacing_um', []);
    if isempty(targetGridSpacingUm)
        gridCases = defaultGridCases(spec);
    else
        gridCases = gridCasesFromTargetSpacing(spec, targetGridSpacingUm);
    end
end
outputRoot = getFieldOrDefault(options, 'outputRoot', defaultOutputRoot());
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

numCases = numel(gridCases);
labels = strings(numCases, 1);
numX = zeros(numCases, 1);
numY = zeros(numCases, 1);
targetGridSpacingUm = nan(numCases, 1);
dxUm = zeros(numCases, 1);
dyUm = zeros(numCases, 1);
endTimeS = zeros(numCases, 1);
finalTotalArea = zeros(numCases, 1);
finalFirstPoreArea = zeros(numCases, 1);
finalFirstThreePoresArea = zeros(numCases, 1);
finalPrecipitateMoles = zeros(numCases, 1);
caseManifestPath = strings(numCases, 1);

for iCase = 1:numCases
    caseSpec = applyGridCase(spec, gridCases(iCase));
    labels(iCase) = string(gridCases(iCase).label);
    numX(iCase) = caseSpec.numX;
    numY(iCase) = caseSpec.numY;
    targetGridSpacingUm(iCase) = getGridCaseField(gridCases(iCase), ...
        'targetGridSpacing_um', nan);
    dxUm(iCase) = caseSpec.dx_cm * 1e4;
    dyUm(iCase) = caseSpec.dy_cm * 1e4;

    caseOptions = options;
    caseOptions.outputRoot = fullfile(outputRoot, char(labels(iCase)));
    if isfield(caseOptions, 'gridCases')
        caseOptions = rmfield(caseOptions, 'gridCases');
    end
    runResult = precip_RunYoonFixedGeometryBenchmark(caseSpec, caseOptions);

    endTimeS(iCase) = runResult.manifest.endTime_s;
    areaTimeseries = runResult.run.areaTimeseries;
    finalTotalArea(iCase) = areaTimeseries.totalPrecipitatedArea_cm2(end);
    finalFirstPoreArea(iCase) = ...
        areaTimeseries.firstPorePrecipitatedArea_cm2(end);
    finalFirstThreePoresArea(iCase) = ...
        areaTimeseries.firstThreePoresPrecipitatedArea_cm2(end);
    finalPrecipitateMoles(iCase) = ...
        sum(runResult.run.state.precipitateMoles(:));
    caseManifestPath(iCase) = string(runResult.manifestPath);
end

relativeDifference = computeRelativeDifferenceFromFinest( ...
    finalTotalArea, dxUm, dyUm);
isQuantitativeBenchmark = false(numCases, 1);

summary = table(labels, numX, numY, targetGridSpacingUm, dxUm, dyUm, endTimeS, ...
    finalTotalArea, finalFirstPoreArea, finalFirstThreePoresArea, ...
    finalPrecipitateMoles, relativeDifference, caseManifestPath, ...
    isQuantitativeBenchmark, 'VariableNames', {'label', 'numX', 'numY', ...
    'targetGridSpacing_um', 'dx_um', 'dy_um', 'endTime_s', ...
    'finalTotalPrecipitatedArea_cm2', 'finalFirstPorePrecipitatedArea_cm2', ...
    'finalFirstThreePoresPrecipitatedArea_cm2', ...
    'finalPrecipitateMoles', ...
    'relativeTotalAreaDifferenceFromFinest', 'caseManifestPath', ...
    'isQuantitativeBenchmark'});

summaryCsv = fullfile(outputRoot, 'yoon_grid_convergence_summary.csv');
writetable(summary, summaryCsv);

manifest = struct();
manifest.runner = 'precip_RunYoonGridConvergenceSmoke';
manifest.outputRoot = outputRoot;
manifest.summaryCsv = summaryCsv;
manifest.numGridCases = numCases;
manifest.caseLabels = cellstr(labels');
manifest.caseManifestPaths = cellstr(caseManifestPath');
manifest.caseManifestFilesVerified = ...
    isfile(summaryCsv) && all(isfile(caseManifestPath));
manifest.requestedTargetGridSpacing_um = targetGridSpacingUm(:)';
manifest.actualDx_um = dxUm(:)';
manifest.actualDy_um = dyUm(:)';
manifest.actualNumX = numX(:)';
manifest.actualNumY = numY(:)';
manifest.maxRelativeTotalAreaDifferenceFromFinest = ...
    max(relativeDifference);
manifest.requiredTargetGridSpacing_um = [10, 5, 2.5];
manifest.targetGridSpacingSequenceComplete = ...
    hasTargetGridSequence(targetGridSpacingUm, ...
        manifest.requiredTargetGridSpacing_um);
manifest.actualGridSpacingSequenceComplete = ...
    hasTargetGridSequence(dxUm, manifest.requiredTargetGridSpacing_um) && ...
    hasTargetGridSequence(dyUm, manifest.requiredTargetGridSpacing_um);
manifest.gridConvergenceTolerance = getFieldOrDefault(options, ...
    'gridConvergenceTolerance', 0.05);
manifest.gridConvergenceWithinTolerance = ...
    manifest.maxRelativeTotalAreaDifferenceFromFinest <= ...
        manifest.gridConvergenceTolerance;
manifest.productionGridConvergenceValidated = getFieldOrDefault(options, ...
    'productionGridConvergenceValidated', false);
manifest.gridConvergenceAcceptanceCriteria = struct( ...
    'targetGridSpacingSequenceComplete', ...
        manifest.targetGridSpacingSequenceComplete, ...
    'gridConvergenceWithinTolerance', ...
        manifest.gridConvergenceWithinTolerance, ...
    'productionGridConvergenceValidated', ...
        manifest.productionGridConvergenceValidated);
manifest.gridConvergenceAccepted = ...
    manifest.targetGridSpacingSequenceComplete && ...
    manifest.gridConvergenceWithinTolerance && ...
    manifest.productionGridConvergenceValidated;
manifest.isQuantitativeBenchmark = manifest.gridConvergenceAccepted;
manifest.note = ['Grid diagnostic smoke wrapper; outputs are not a ' ...
    'grid-converged Zhang/Yoon benchmark claim.'];

manifestPath = fullfile(outputRoot, 'yoon_grid_convergence_manifest.json');
writeJsonManifest(manifestPath, manifest);

result = struct();
result.summary = summary;
result.manifest = manifest;
result.manifestPath = manifestPath;
end

function gridCases = defaultGridCases(spec)
baseNumX = min(max(spec.numX, 8), 30);
baseNumY = min(max(spec.numY, 8), 30);
gridCases = struct( ...
    'label', {'coarse', 'medium', 'fine'}, ...
    'numX', {max(4, ceil(baseNumX / 2)), baseNumX, baseNumX * 2}, ...
    'numY', {max(4, ceil(baseNumY / 2)), baseNumY, baseNumY * 2});
end

function gridCases = gridCasesFromTargetSpacing(spec, targetSpacingUm)
targetSpacingUm = double(targetSpacingUm(:)');
if isempty(targetSpacingUm) || any(~isfinite(targetSpacingUm)) || ...
        any(targetSpacingUm <= 0)
    error('RTSPHEM:Precipitate:InvalidTargetGridSpacing', ...
        'targetGridSpacing_um must contain positive finite values.');
end

gridCases = repmat(struct('label', '', 'numX', 0, 'numY', 0, ...
    'targetGridSpacing_um', nan), 1, numel(targetSpacingUm));
for iCase = 1:numel(targetSpacingUm)
    spacingCm = targetSpacingUm(iCase) * 1e-4;
    gridCases(iCase).label = spacingLabel(targetSpacingUm(iCase));
    gridCases(iCase).numX = max(1, round(spec.lengthXAxis_cm / spacingCm));
    gridCases(iCase).numY = max(1, round(spec.lengthYAxis_cm / spacingCm));
    gridCases(iCase).targetGridSpacing_um = targetSpacingUm(iCase);
end
end

function label = spacingLabel(spacingUm)
spacingText = regexprep(sprintf('%.12g', spacingUm), '\.', 'p');
label = ['dx_', spacingText, 'um'];
end

function caseSpec = applyGridCase(spec, gridCase)
caseSpec = spec;
caseSpec.numX = gridCase.numX;
caseSpec.numY = gridCase.numY;
caseSpec.dx_cm = caseSpec.lengthXAxis_cm / caseSpec.numX;
caseSpec.dy_cm = caseSpec.lengthYAxis_cm / caseSpec.numY;
caseSpec.cellVolume_cm3 = ...
    caseSpec.dx_cm * caseSpec.dy_cm * caseSpec.thickness_cm;
if isfield(caseSpec, 'splitInletY_cm') && isempty(caseSpec.splitInletY_cm)
    caseSpec.splitInletY_cm = 0.5 * caseSpec.lengthYAxis_cm;
end
caseSpec.inletA.yRange_cm = [0, caseSpec.splitInletY_cm];
caseSpec.inletB.yRange_cm = ...
    [caseSpec.splitInletY_cm, caseSpec.lengthYAxis_cm];
end

function relativeDifference = computeRelativeDifferenceFromFinest( ...
    finalTotalArea, dxUm, dyUm)
cellSize = max(dxUm, dyUm);
[~, finestIndex] = min(cellSize);
referenceArea = finalTotalArea(finestIndex);
denominator = max(abs(referenceArea), eps);
relativeDifference = abs(finalTotalArea - referenceArea) ./ denominator;
end

function tf = hasTargetGridSequence(spacingUm, requiredSpacingUm)
spacing = double(spacingUm(:)');
required = double(requiredSpacingUm(:)');
tf = all(isfinite(spacing)) && ...
    all(ismembertol(required, spacing, 1e-9));
end

function outputRoot = defaultOutputRoot()
moduleRoot = fileparts(fileparts(mfilename('fullpath')));
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
outputRoot = fullfile(moduleRoot, 'outputs', 'yoon_grid_convergence', timestamp);
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

function value = getGridCaseField(gridCase, fieldName, defaultValue)
if isfield(gridCase, fieldName) && ~isempty(gridCase.(fieldName))
    value = gridCase.(fieldName);
else
    value = defaultValue;
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
