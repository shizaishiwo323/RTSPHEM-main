function result = precip_RunYoonProductionFlowValidation(spec, options)
% precip_RunYoonProductionFlowValidation - Validate production Stokes flow evidence.
%
% Runs a small set of fail-closed flow-validation cases through the
% `hyphm_stokes` bridge. This does not make a benchmark quantitative by itself;
% it creates an auditable manifest that readiness gates can require alongside
% fixed-geometry Case 1 evidence.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

outputRoot = getFieldOrDefault(options, 'outputRoot', defaultOutputRoot());
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

validationSpec = spec;
validationSpec.yoonFlowSolver = 'hyphm_stokes';
if isfield(options, 'hyphmStokesSolverFcn') && ...
        ~isempty(options.hyphmStokesSolverFcn)
    validationSpec.hyphmStokesSolverFcn = options.hyphmStokesSolverFcn;
end

caseNames = ["empty_channel"; "single_obstacle"; "blocked_column"];
numCases = numel(caseNames);
linearResidual = NaN(numCases, 1);
divergenceResidual = NaN(numCases, 1);
boundaryClosure = NaN(numCases, 1);
relativePermeability = NaN(numCases, 1);
inletFlux = NaN(numCases, 1);
outletFlux = NaN(numCases, 1);
numBlockedCells = zeros(numCases, 1);

for iCase = 1:numCases
    state = makeValidationState(validationSpec, caseNames(iCase));
    flow = precip_RecomputeYoonFlowField(state, validationSpec);
    linearResidual(iCase) = getFieldOrDefault(flow, ...
        'linearResidualRelative', NaN);
    divergenceResidual(iCase) = getFieldOrDefault(flow, ...
        'maxDivergenceResidual_s_inv', NaN);
    boundaryClosure(iCase) = getFieldOrDefault(flow, ...
        'boundaryFluxClosureRelativeError', NaN);
    relativePermeability(iCase) = getFieldOrDefault(flow, ...
        'relativePermeability', NaN);
    inletFlux(iCase) = getFieldOrDefault(flow, 'inletFlux_cm3_s', NaN);
    outletFlux(iCase) = getFieldOrDefault(flow, 'outletFlux_cm3_s', NaN);
    numBlockedCells(iCase) = getFieldOrDefault(flow, 'numBlockedCells', 0);
end

summary = table(caseNames, linearResidual, divergenceResidual, ...
    boundaryClosure, relativePermeability, inletFlux, outletFlux, ...
    numBlockedCells, 'VariableNames', {'caseName', ...
    'linearResidualRelative', 'divergenceResidual_s_inv', ...
    'boundaryClosureRelativeError', 'relativePermeability', ...
    'inletFlux_cm3_s', 'outletFlux_cm3_s', 'numBlockedCells'});
summaryCsv = fullfile(outputRoot, 'yoon_production_flow_validation.csv');
writetable(summary, summaryCsv);

maxAcceptedLinear = getFieldOrDefault(options, ...
    'maxAcceptedLinearResidualRelative', 1e-8);
maxAcceptedDivergence = getFieldOrDefault(options, ...
    'maxAcceptedDivergenceResidual_s_inv', 1e-8);
maxAcceptedClosure = getFieldOrDefault(options, ...
    'maxAcceptedBoundaryClosureRelativeError', 1e-8);

manifest = struct();
manifest.runner = 'precip_RunYoonProductionFlowValidation';
manifest.solver = 'hyphm_stokes';
manifest.isProxy = false;
manifest.isStokes = true;
manifest.outputRoot = outputRoot;
manifest.summaryCsv = summaryCsv;
manifest.validationSummaryCsvVerified = isfile(summaryCsv);
manifest.requiredCaseNames = cellstr(caseNames);
manifest.validatedCaseNames = cellstr(summary.caseName(:)');
manifest.requiredCasesComplete = all(ismember(caseNames, summary.caseName));
manifest.maxLinearResidualRelative = max(linearResidual);
manifest.maxAcceptedLinearResidualRelative = maxAcceptedLinear;
manifest.maxDivergenceResidual_s_inv = max(divergenceResidual);
manifest.maxAcceptedDivergenceResidual_s_inv = maxAcceptedDivergence;
manifest.maxBoundaryClosureRelativeError = max(boundaryClosure);
manifest.maxAcceptedBoundaryClosureRelativeError = maxAcceptedClosure;
manifest.minRelativePermeability = min(relativePermeability);
manifest.maxRelativePermeability = max(relativePermeability);
manifest = appendStructFields(manifest, ...
    precip_BuildOutputEvidenceProvenance(summaryCsv, options));
manifest.productionFlowValidationAccepted = ...
    manifest.validationSummaryCsvVerified && ...
    hasCompleteOutputProvenance(manifest) && ...
    manifest.requiredCasesComplete && ...
    all(isfinite(linearResidual)) && ...
    all(isfinite(divergenceResidual)) && ...
    all(isfinite(boundaryClosure)) && ...
    all(isfinite(relativePermeability)) && ...
    manifest.maxLinearResidualRelative <= maxAcceptedLinear && ...
    manifest.maxDivergenceResidual_s_inv <= maxAcceptedDivergence && ...
    manifest.maxBoundaryClosureRelativeError <= maxAcceptedClosure && ...
    manifest.minRelativePermeability >= 0 && ...
    manifest.maxRelativePermeability <= 1;

manifestPath = fullfile(outputRoot, ...
    'yoon_production_flow_validation_manifest.json');
writeJsonManifest(manifestPath, manifest);

result = struct();
result.summary = summary;
result.summaryCsv = summaryCsv;
result.manifest = manifest;
result.manifestPath = manifestPath;
end

function state = makeValidationState(spec, caseName)
state = precip_YoonMicrocontinuumSolver('initialize', spec, ...
    struct('substrateMask', false(spec.numY, spec.numX)));
state.Vm(:) = 0;
state.blockedMask(:) = false;
switch string(caseName)
    case "empty_channel"
    case "single_obstacle"
        iy = max(1, round(0.5 * spec.numY));
        ix = max(1, round(0.5 * spec.numX));
        state.substrateMask(iy, ix) = true;
    case "blocked_column"
        ix = max(2, min(spec.numX - 1, round(0.5 * spec.numX)));
        yRange = 2:max(2, spec.numY - 1);
        state.blockedMask(yRange, ix) = true;
        state.Vm(yRange, ix) = max(getFieldOrDefault(spec, ...
            'blockedVmThreshold', 0.6), 0.6);
    otherwise
        error('RTSPHEM:Precipitate:UnknownFlowValidationCase', ...
            'Unsupported flow validation case: %s.', caseName);
end
end

function outputRoot = defaultOutputRoot()
moduleRoot = fileparts(fileparts(mfilename('fullpath')));
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
outputRoot = fullfile(moduleRoot, 'outputs', ...
    'yoon_production_flow_validation', timestamp);
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
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
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
