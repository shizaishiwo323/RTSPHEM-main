function reactionInput = BuildStrictMolinsReactionInput(state, geometry, connectivity, options)
%BUILDSTRICTMOLINSREACTIONINPUT Prepare interface-based strict Molins options.

if nargin < 3 || isempty(connectivity)
    connectivity = struct();
end
if nargin < 4 || isempty(options)
    options = struct();
end

reactionInput = struct();
reactionInput.rate_constant_cm_s = getFieldOrDefault(options, 'rate_constant_cm_s', 0);
reactionInput.maxReactantFraction = getFieldOrDefault(options, ...
    'maxReactantFraction', Inf);
reactionInput.maxMineralFraction = getFieldOrDefault(options, ...
    'maxMineralFraction', Inf);
reactionInput.reaction_time_integration = getFieldOrDefault(options, ...
    'reaction_time_integration', 'explicit_euler');
reactionInput.interfaceState = rtm.chemistry.ReconstructFluidSideState( ...
    state, geometry, connectivity, options);
reactionInput.interfaceQuadrature = rtm.geometry.BuildInterfaceQuadrature(geometry);
reactionInput.reactionClusters = rtm.chemistry.BuildReactionClusters( ...
    state, geometry, connectivity, options);
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
