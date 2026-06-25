function [cellResult, info] = RunPhreeqcReactionClusters( ...
        state, geometry, hMolCm3, prescribedMoles, clusters, options)
%RUNPHREEQCREACTIONCLUSTERS Run PHREEQC once on disjoint reaction clusters.
%
% The input state remains cell-based. This helper aggregates extensive
% quantities over each disjoint cluster, executes the configured PHREEQC batch
% runner on the cluster state, then scatters component and mineral deltas back
% to the original cells.

if nargin < 6 || isempty(options)
    options = struct();
end

[clusterState, map] = rtm.chemistry.AggregateReactionClusterState( ...
    state, geometry, hMolCm3, prescribedMoles, clusters, options);
if logical(getOption(options, 'omitPrescribedCalciteReaction', false)) && ...
        isfield(clusterState, 'prescribed_calcite_dissolved_moles')
    clusterState = rmfield(clusterState, 'prescribed_calcite_dissolved_moles');
end

runner = getRunBatchFunction(options);
rawBatchResult = runner(clusterState, options);
cellResult = rtm.chemistry.ScatterClusterReactionResult( ...
    rawBatchResult, map, size(state.component_moles, 1));

info = struct();
info.cluster_count = numel(clusters);
info.membership_count = map.membership_count;
info.cluster_state = clusterState;
info.map = map;
info.raw_batch_result = rawBatchResult;
end

function runner = getRunBatchFunction(options)
if isfield(options, 'runBatchFunction') && ~isempty(options.runBatchFunction)
    runner = options.runBatchFunction;
else
    runner = @RunPhreeqcCalciteBatch;
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
