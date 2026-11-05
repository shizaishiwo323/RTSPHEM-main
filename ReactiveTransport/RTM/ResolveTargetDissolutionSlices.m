function settings = ResolveTargetDissolutionSlices(targetSliceCount, initialPorosity, explicitMaxTotalTimeSteps, fallbackPorosityStepTarget, safetyFactor, minExtraSteps, finalPorosity)
% ResolveTargetDissolutionSlices converts a desired full-dissolution slice
% count into adaptive-porosity time stepping settings.

if nargin < 5 || isempty(safetyFactor)
    safetyFactor = 2.0;
end
if nargin < 6 || isempty(minExtraSteps)
    minExtraSteps = 20;
end
if nargin < 7 || isempty(finalPorosity)
    finalPorosity = 1;
end
validateattributes(finalPorosity, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', '<=', 1});

settings = struct( ...
    'enabled', false, ...
    'targetSliceCount', [], ...
    'finalPorosity', finalPorosity, ...
    'porosityStepTarget', fallbackPorosityStepTarget, ...
    'maxTotalTimeSteps', explicitMaxTotalTimeSteps);

if isempty(targetSliceCount)
    return;
end

targetSliceCount = double(targetSliceCount(1));
if ~isfinite(targetSliceCount) || targetSliceCount < 2
    error('RTM:InvalidTargetDissolutionSlices', ...
        'targetDissolutionSlices must be empty or a finite value >= 2.');
end

targetSliceCount = round(targetSliceCount);
if nargin < 2 || isempty(initialPorosity) || ~isfinite(initialPorosity)
    dissolvablePorosityRange = finalPorosity;
else
    initialPorosity = min(max(double(initialPorosity(1)), 0), 1);
    if finalPorosity <= initialPorosity
        error('RTM:InvalidTargetDissolutionFinalPorosity', ...
            'targetDissolutionFinalPorosity must exceed the initial porosity.');
    end
    dissolvablePorosityRange = finalPorosity - initialPorosity;
end

settings.enabled = true;
settings.targetSliceCount = targetSliceCount;
settings.porosityStepTarget = dissolvablePorosityRange / (targetSliceCount - 1);

if isempty(explicitMaxTotalTimeSteps)
    safetyFactor = max(1.0, double(safetyFactor(1)));
    minExtraSteps = max(0, round(double(minExtraSteps(1))));
    settings.maxTotalTimeSteps = max( ...
        targetSliceCount + minExtraSteps, ...
        ceil(targetSliceCount * safetyFactor));
else
    settings.maxTotalTimeSteps = max(2, round(double(explicitMaxTotalTimeSteps(1))));
end
end
