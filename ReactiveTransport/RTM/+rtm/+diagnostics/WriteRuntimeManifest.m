function path = WriteRuntimeManifest(outputDir, manifest)
%WRITERUNTIMEMANIFEST Persist RTM provenance manifest as JSON.

outputDir = char(outputDir);
if exist(outputDir, 'dir') ~= 7
    error('RTSPHEM:Diagnostics:MissingOutputDirectory', ...
        'Output directory does not exist: %s.', outputDir);
end
if nargin < 2 || isempty(manifest) || ~isstruct(manifest)
    error('RTSPHEM:Diagnostics:InvalidRuntimeManifest', ...
        'manifest must be a nonempty struct.');
end

path = string(fullfile(outputDir, 'runtime_manifest.json'));
fid = fopen(path, 'w');
if fid == -1
    error('RTSPHEM:Diagnostics:ManifestOpenFailed', ...
        'Cannot write runtime manifest: %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(manifest));
clear cleanup;
end
