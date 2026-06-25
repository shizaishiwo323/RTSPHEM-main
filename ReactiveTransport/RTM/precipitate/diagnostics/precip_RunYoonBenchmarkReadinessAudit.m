function report = precip_RunYoonBenchmarkReadinessAudit(options)
% precip_RunYoonBenchmarkReadinessAudit - Write readiness audit CSV/JSON.
%
% Inputs:
%   options.fixedGeometryManifestPath
%   options.gridConvergenceManifestPath
%   options.chemistryManifestPath (optional)
%   options.passiveTransportManifestPath (optional)
%   options.caseComparisonManifestPath (optional)
%   options.diffusionFeedbackManifestPath (optional)
%   options.flowValidationManifestPath (optional)
%   options.referenceCsv
%   options.referenceDigitizationPackageDir (optional)
%   options.outputRoot (optional)

if nargin < 1 || ~isstruct(options)
    error('RTSPHEM:Precipitate:InvalidReadinessAuditOptions', ...
        'options must be a struct.');
end

outputRoot = getFieldOrDefault(options, 'outputRoot', pwd);
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

evidence = struct();
evidence.fixedGeometryManifest = loadJsonStruct( ...
    requireOption(options, 'fixedGeometryManifestPath'));
evidence.gridConvergenceManifest = loadJsonStruct( ...
    requireOption(options, 'gridConvergenceManifestPath'));
chemistryManifestPath = getFieldOrDefault(options, ...
    'chemistryManifestPath', '');
if ~isempty(chemistryManifestPath)
    evidence.chemistryManifest = loadJsonStruct(chemistryManifestPath);
else
    evidence.chemistryManifest = struct();
end
passiveTransportManifestPath = getFieldOrDefault(options, ...
    'passiveTransportManifestPath', '');
if ~isempty(passiveTransportManifestPath)
    evidence.passiveTransportManifest = loadJsonStruct( ...
        passiveTransportManifestPath);
else
    evidence.passiveTransportManifest = struct();
end
caseComparisonManifestPath = getFieldOrDefault(options, ...
    'caseComparisonManifestPath', '');
if ~isempty(caseComparisonManifestPath)
    evidence.caseComparisonManifest = loadJsonStruct( ...
        caseComparisonManifestPath);
else
    evidence.caseComparisonManifest = struct();
end
diffusionFeedbackManifestPath = getFieldOrDefault(options, ...
    'diffusionFeedbackManifestPath', '');
if ~isempty(diffusionFeedbackManifestPath)
    evidence.diffusionFeedbackManifest = loadJsonStruct( ...
        diffusionFeedbackManifestPath);
else
    evidence.diffusionFeedbackManifest = struct();
end
flowValidationManifestPath = getFieldOrDefault(options, ...
    'flowValidationManifestPath', '');
if ~isempty(flowValidationManifestPath)
    evidence.flowValidationManifest = loadJsonStruct( ...
        flowValidationManifestPath);
else
    evidence.flowValidationManifest = struct();
end

referenceCsv = requireOption(options, 'referenceCsv');
referencePackageDir = getFieldOrDefault(options, ...
    'referenceDigitizationPackageDir', []);
if isempty(referencePackageDir)
    evidence.reference = precip_LoadReferenceCurves(referenceCsv);
else
    evidence.reference = precip_LoadReferenceCurves(referenceCsv, ...
        referencePackageDir);
end

audit = precip_AuditYoonBenchmarkReadiness(evidence);
ladder = precip_AuditYoonBenchmarkLadder(evidence);

requirementsCsv = fullfile(outputRoot, ...
    'yoon_benchmark_readiness_requirements.csv');
writetable(audit.requirements, requirementsCsv);
ladderCsv = fullfile(outputRoot, 'yoon_benchmark_ladder.csv');
writetable(formatLadderStagesForCsv(ladder.stages), ladderCsv);

manifest = struct();
manifest.runner = 'precip_RunYoonBenchmarkReadinessAudit';
manifest.fixedGeometryManifestPath = char(string( ...
    options.fixedGeometryManifestPath));
manifest.gridConvergenceManifestPath = char(string( ...
    options.gridConvergenceManifestPath));
manifest.chemistryManifestPath = char(string(chemistryManifestPath));
manifest.passiveTransportManifestPath = char(string( ...
    passiveTransportManifestPath));
manifest.caseComparisonManifestPath = char(string( ...
    caseComparisonManifestPath));
manifest.diffusionFeedbackManifestPath = char(string( ...
    diffusionFeedbackManifestPath));
manifest.flowValidationManifestPath = char(string( ...
    flowValidationManifestPath));
manifest.referenceCsv = char(string(referenceCsv));
manifest.requirementsCsv = requirementsCsv;
manifest.ladderCsv = ladderCsv;
manifest.allReady = audit.allReady;
manifest.isQuantitativeBenchmarkAllowed = ...
    audit.isQuantitativeBenchmarkAllowed;
manifest.allStagesReady = ladder.allStagesReady;
manifest.firstFailedStageId = char(ladder.firstFailedStageId);
manifest.failedStageIds = cellstr(ladder.stages.stageId( ...
    ladder.stages.entered & ~ladder.stages.ready));
manifest.failedRequirementIds = cellstr( ...
    audit.requirements.requirementId(~audit.requirements.ready));
manifest.note = audit.note;

manifestPath = fullfile(outputRoot, ...
    'yoon_benchmark_readiness_manifest.json');
writeJsonManifest(manifestPath, manifest);

report = struct();
report.audit = audit;
report.ladder = ladder;
report.requirementsCsv = requirementsCsv;
report.ladderCsv = ladderCsv;
report.manifest = manifest;
report.manifestPath = manifestPath;
end

function csvTable = formatLadderStagesForCsv(stages)
csvTable = stages;
csvTable.requirementIds = cellfun(@joinTextList, ...
    csvTable.requirementIds, 'UniformOutput', false);
csvTable.blockingRequirementIds = cellfun(@joinTextList, ...
    csvTable.blockingRequirementIds, 'UniformOutput', false);
end

function value = joinTextList(values)
value = char(strjoin(string(values(:)'), ';'));
end

function data = loadJsonStruct(path)
path = char(string(path));
if exist(path, 'file') ~= 2
    error('RTSPHEM:Precipitate:ReadinessEvidenceNotFound', ...
        'Readiness evidence file does not exist: %s', path);
end
data = jsondecode(fileread(path));
end

function value = requireOption(options, fieldName)
if ~isfield(options, fieldName) || isempty(options.(fieldName))
    error('RTSPHEM:Precipitate:MissingReadinessAuditOption', ...
        'Missing required readiness audit option: %s', fieldName);
end
value = options.(fieldName);
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
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
