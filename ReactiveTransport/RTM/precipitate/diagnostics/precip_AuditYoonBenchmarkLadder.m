function ladder = precip_AuditYoonBenchmarkLadder(evidence)
% precip_AuditYoonBenchmarkLadder - Enforce ordered B0-B8 benchmark stages.
%
% This audit maps the readiness requirements onto the formal Zhang/Yoon
% benchmark ladder. A failed stage stops entry into every later stage.

if nargin < 1 || ~isstruct(evidence)
    error('RTSPHEM:Precipitate:InvalidBenchmarkLadderEvidence', ...
        'evidence must be a struct.');
end

readiness = precip_AuditYoonBenchmarkReadiness(evidence);
definitions = stageDefinitions();

numStages = numel(definitions);
stageId = strings(numStages, 1);
stageName = strings(numStages, 1);
requirementIds = cell(numStages, 1);
entered = false(numStages, 1);
ready = false(numStages, 1);
blockingRequirementIds = cell(numStages, 1);

canEnter = true;
firstFailedStageId = "";
for iStage = 1:numStages
    stageId(iStage) = definitions(iStage).stageId;
    stageName(iStage) = definitions(iStage).stageName;
    requirementIds{iStage} = cellstr(definitions(iStage).requirementIds);
    entered(iStage) = canEnter;
    if ~canEnter
        ready(iStage) = false;
        blockingRequirementIds{iStage} = {};
        continue;
    end

    [stageReady, failedRequirements] = requirementsReady( ...
        readiness, definitions(iStage).requirementIds);
    ready(iStage) = stageReady;
    blockingRequirementIds{iStage} = cellstr(failedRequirements);
    if ~stageReady
        canEnter = false;
        if strlength(firstFailedStageId) == 0
            firstFailedStageId = stageId(iStage);
        end
    end
end

ladder = struct();
ladder.stages = table(stageId, stageName, requirementIds, entered, ready, ...
    blockingRequirementIds);
ladder.allStagesReady = all(entered & ready);
ladder.firstFailedStageId = firstFailedStageId;
ladder.readinessAudit = readiness;
if ladder.allStagesReady
    ladder.note = 'All ordered B0-B8 Zhang/Yoon benchmark stages pass.';
else
    ladder.note = ['One or more ordered Zhang/Yoon benchmark stages fail; ', ...
        'later stages are not entered.'];
end
end

function definitions = stageDefinitions()
definitions = struct( ...
    'stageId', {'B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8'}, ...
    'stageName', { ...
    'No-reaction split-inlet mixing', ...
    'Zero-dimensional PHREEQC/Yoon chemistry', ...
    'Fixed-Vm steady reaction transport', ...
    'Case 1 n=2 fixed-flow benchmark', ...
    'Case 2/1/3 n=0/2/3 diffusion feedback', ...
    'Vm >= 0.6 flow-feedback update', ...
    'Case 5 df=300 delayed activation', ...
    '10/5/2.5 um grid convergence', ...
    '13/18/118 min image and area comparison'}, ...
    'requirementIds', { ...
    ["passive_transport"], ...
    ["chemistry_speciation"], ...
    ["reaction_mass", "no_finite_clipping"], ...
    ["rate_constants_locked", "center_band_morphology", ...
        "quantitative_manifest"], ...
    ["diffusion_feedback"], ...
    ["production_stokes", "flow_feedback"], ...
    ["case1_case5"], ...
    ["grid_convergence"], ...
    ["target_times", "geometry_package", "reference_package"]});
end

function [tf, failedRequirements] = requirementsReady(readiness, requiredIds)
failedRequirements = strings(0, 1);
tf = true;
for iRequirement = 1:numel(requiredIds)
    requirementId = string(requiredIds(iRequirement));
    row = readiness.requirements.requirementId == requirementId;
    if ~any(row) || ~all(readiness.requirements.ready(row))
        tf = false;
        failedRequirements(end + 1, 1) = requirementId; %#ok<AGROW>
    end
end
end
