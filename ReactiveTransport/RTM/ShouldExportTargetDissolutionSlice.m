function [doExport, nextExportIndex] = ShouldExportTargetDissolutionSlice( ...
    stepIndex, currentPorosity, initialPorosity, targetSliceCount, lastExportIndex, forceExport)
% ShouldExportTargetDissolutionSlice decides whether a target dissolution
% progress slice should be exported for the current RTM step.

targetSliceCount = max(2, round(double(targetSliceCount(1))));
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
currentPorosity = min(max(double(currentPorosity(1)), initialPorosity), 1);
progress = (currentPorosity - initialPorosity) / max(eps, 1 - initialPorosity);
bucketTolerance = 0.05;
progressIndex = min(targetSliceCount, max(1, floor(progress * targetSliceCount + bucketTolerance)));

if progressIndex > lastExportIndex
    doExport = true;
    nextExportIndex = progressIndex;
end
end
