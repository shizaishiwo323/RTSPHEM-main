function [limitedStepSize, wasLimited, targetBoundaryPorosity] = LimitTargetDissolutionStepSize( ...
    proposedStepSize, currentStepSize, porosityDelta, currentPorosity, ...
    initialPorosity, lastExportIndex, targetSliceCount, minStepSize, maxStepSize, boundaryFactor)
% LimitTargetDissolutionStepSize caps dt so target-slice mode does not jump
% beyond the second upcoming dissolution-progress export bucket.

limitedStepSize = proposedStepSize;
wasLimited = false;
targetBoundaryPorosity = NaN;

if nargin < 10 || isempty(boundaryFactor)
    boundaryFactor = 0.95;
end

if ~all(isfinite([proposedStepSize, currentStepSize, porosityDelta, ...
        currentPorosity, initialPorosity, lastExportIndex, targetSliceCount]))
    return;
end
if proposedStepSize <= 0 || currentStepSize <= 0 || porosityDelta <= eps
    return;
end

targetSliceCount = max(2, round(double(targetSliceCount(1))));
lastExportIndex = max(0, min(targetSliceCount, round(double(lastExportIndex(1)))));
if lastExportIndex >= targetSliceCount - 1
    return;
end

initialPorosity = min(max(double(initialPorosity(1)), 0), 1 - eps);
currentPorosity = min(max(double(currentPorosity(1)), initialPorosity), 1);
nextBoundaryIndex = min(targetSliceCount, lastExportIndex + 2);
targetBoundaryPorosity = initialPorosity + ...
    (1 - initialPorosity) * double(nextBoundaryIndex) / double(targetSliceCount);
remainingPorosity = targetBoundaryPorosity - currentPorosity;
if remainingPorosity <= eps
    return;
end

porositySlope = porosityDelta / currentStepSize;
if porositySlope <= eps
    return;
end

boundaryFactor = min(1.5, max(0.5, double(boundaryFactor(1))));
boundaryStepSize = boundaryFactor * remainingPorosity / porositySlope;
boundaryStepSize = min(maxStepSize, max(minStepSize, boundaryStepSize));

if boundaryStepSize < limitedStepSize
    limitedStepSize = boundaryStepSize;
    wasLimited = true;
end
end
