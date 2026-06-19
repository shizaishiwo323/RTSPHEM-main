function settings = ResolveAdaptivePorosityStepSettings(config)
% ResolveAdaptivePorosityStepSettings centralizes adaptive porosity defaults.

if nargin < 1
    config = struct();
end

settings = struct();
settings.porosityStepUpperFactor = cfggetLocal(config, 'porosityStepUpperFactor', 1.0);
if isTargetSliceMode(config)
    defaultGrowthFactor = 4.0;
else
    defaultGrowthFactor = 2.0;
end
settings.adaptiveGrowthFactor = cfggetLocal(config, 'adaptiveGrowthFactor', defaultGrowthFactor);
end

function value = cfggetLocal(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end

function tf = isTargetSliceMode(config)
tf = isstruct(config) && isfield(config, 'targetDissolutionSlices') && ...
    ~isempty(config.targetDissolutionSlices) && ...
    isfinite(double(config.targetDissolutionSlices(1))) && ...
    double(config.targetDissolutionSlices(1)) > 0;
end
