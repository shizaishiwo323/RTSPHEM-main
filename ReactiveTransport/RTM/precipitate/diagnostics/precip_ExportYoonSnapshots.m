function manifest = precip_ExportYoonSnapshots(outputRoot, snapshotTable, spec)
% precip_ExportYoonSnapshots - Export Yoon Vm snapshots and area metrics.
%
% Inputs:
%   outputRoot    - output directory.
%   snapshotTable - table from precip_RunYoonCase1Short result.snapshots.
%   spec          - benchmark spec.
%
% Output:
%   manifest      - exported file paths and count.

if nargin < 3
    spec = struct();
end
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end
snapshotDir = fullfile(outputRoot, 'yoon_snapshots');
if ~isfolder(snapshotDir)
    mkdir(snapshotDir);
end

numSnapshots = height(snapshotTable);
matFiles = strings(numSnapshots, 1);
for iSnapshot = 1:numSnapshots
    targetTime_s = snapshotTable.targetTime_s(iSnapshot);
    capturedTime_s = snapshotTable.capturedTime_s(iSnapshot);
    Vm = snapshotTable.Vm{iSnapshot};
    filename = snapshotFilename(targetTime_s);
    outputPath = fullfile(snapshotDir, filename);
    matFiles(iSnapshot) = string(outputPath);
    save(outputPath, 'Vm', 'targetTime_s', 'capturedTime_s', 'spec');
end

areaTable = snapshotTable(:, {'targetTime_s', 'capturedTime_s', ...
    'totalPrecipitatedArea_cm2', 'firstPorePrecipitatedArea_cm2', ...
    'firstThreePoresPrecipitatedArea_cm2'});
areaCsv = fullfile(snapshotDir, 'yoon_snapshot_area_metrics.csv');
writetable(areaTable, areaCsv);

manifest = struct();
manifest.snapshotDir = snapshotDir;
manifest.numSnapshots = numSnapshots;
manifest.capturedTargetTimes_s = snapshotTable.targetTime_s(:)';
manifest.capturedTimes_s = snapshotTable.capturedTime_s(:)';
manifest.matFiles = matFiles;
manifest.areaCsv = areaCsv;
end

function filename = snapshotFilename(targetTime_s)
targetMin = targetTime_s / 60;
if abs(targetMin - round(targetMin)) < 1e-9
    filename = sprintf('yoon_vm_snapshot_%03dmin.mat', round(targetMin));
else
    filename = sprintf('yoon_vm_snapshot_%06gs.mat', targetTime_s);
end
end
