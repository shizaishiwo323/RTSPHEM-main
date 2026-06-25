function path = WriteConvergenceReport(outputDir, report)
%WRITECONVERGENCEREPORT Persist benchmark convergence evidence as JSON.

outputDir = char(outputDir);
if exist(outputDir, 'dir') ~= 7
    error('RTSPHEM:Benchmark:MissingOutputDirectory', ...
        'Output directory does not exist: %s.', outputDir);
end
if nargin < 2 || isempty(report) || ~isstruct(report)
    error('RTSPHEM:Benchmark:InvalidConvergenceReport', ...
        'report must be a nonempty struct.');
end

path = string(fullfile(outputDir, 'benchmark_convergence_report.json'));
fid = fopen(path, 'w');
if fid == -1
    error('RTSPHEM:Benchmark:ReportOpenFailed', ...
        'Cannot write benchmark convergence report: %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(report));
clear cleanup;
end
