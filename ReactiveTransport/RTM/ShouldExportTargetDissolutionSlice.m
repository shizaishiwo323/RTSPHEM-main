function [doExport, nextExportIndex] = ShouldExportTargetDissolutionSlice( ...
    stepIndex, currentPorosity, initialPorosity, targetSliceCount, lastExportIndex, forceExport, finalPorosity)
% ShouldExportTargetDissolutionSlice decides whether a target dissolution
% progress slice should be exported for the current RTM step.

targetSliceCount = max(2, round(double(targetSliceCount(1))));
if nargin < 7 || isempty(finalPorosity)
    finalPorosity = 1;
end
lastExportIndex = max(0, round(double(lastExportIndex(1))));
nextExportIndex = min(lastExportIndex, targetSliceCount);
doExport = false;

if forceExport
    doExport = true;
    nextExportIndex = targetSliceCount;
    return;
end

if stepIndex <= 1 && lastExportIndex < 1
    doExport = true;
    nextExportIndex = 1;
    return;
end

if ~isfinite(currentPorosity) || ~isfinite(initialPorosity)
    return;
end

initialPorosity = min(max(double(initialPorosity(1)), 0), 1 - eps);
finalPorosity = min(max(double(finalPorosity(1)), initialPorosity + eps), 1);
currentPorosity = min(max(double(currentPorosity(1)), initialPorosity), finalPorosity);
progress = (currentPorosity - initialPorosity) / ...
    max(eps, finalPorosity - initialPorosity);
bucketTolerance = 0.05;
progressIndex = min(targetSliceCount, max(1, ...
    floor(progress * (targetSliceCount - 1) + bucketTolerance) + 1));
if finalPorosity >= 1 - eps
    % Full-dissolution runs reserve the last slice for the grain-volume stop.
    progressIndex = min(progressIndex, targetSliceCount - 1);
end

if progressIndex > lastExportIndex
    doExport = true;
    nextExportIndex = progressIndex;
end
end
