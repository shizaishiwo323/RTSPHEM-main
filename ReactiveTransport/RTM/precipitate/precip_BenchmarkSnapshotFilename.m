function filename = precip_BenchmarkSnapshotFilename(targetSeconds)
%PRECIP_BENCHMARKSNAPSHOTFILENAME Stable filename for benchmark snapshots.

if ~isscalar(targetSeconds) || ~isnumeric(targetSeconds) || ~isfinite(targetSeconds) || targetSeconds < 0
    error('RTSPHEM:Precipitate:InvalidSnapshotTime', ...
        'targetSeconds must be a finite nonnegative numeric scalar.');
end

minutes = targetSeconds / 60;
if targetSeconds >= 60 && abs(minutes - round(minutes)) < 1e-9
    filename = sprintf('benchmark_snapshot_%03dmin.png', round(minutes));
elseif abs(targetSeconds - round(targetSeconds)) < 1e-9
    filename = sprintf('benchmark_snapshot_%03ds.png', round(targetSeconds));
else
    secondsLabel = regexprep(sprintf('%.6g', targetSeconds), '\.', 'p');
    filename = sprintf('benchmark_snapshot_%ss.png', secondsLabel);
end
end
