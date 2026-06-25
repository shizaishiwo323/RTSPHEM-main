function provenance = precip_BuildOutputEvidenceProvenance(outputFiles, options)
% precip_BuildOutputEvidenceProvenance - Hash output files and record code source.
%
% Inputs:
%   outputFiles - char/string/cellstr file paths that support a validation claim.
%   options     - optional sourceCommitSha override for tests or archived runs.
%
% Output:
%   provenance  - struct with sourceCommitSha and SHA-256 output hash evidence.

if nargin < 2
    options = struct();
end

files = string(outputFiles);
files = files(:);
sourceCommitSha = getFieldOrDefault(options, 'sourceCommitSha', "");
if strlength(sourceCommitSha) == 0
    sourceCommitSha = queryGitCommitSha();
end

fileHashes = strings(numel(files), 1);
filesVerified = false(numel(files), 1);
for iFile = 1:numel(files)
    filesVerified(iFile) = isfile(files(iFile));
    if filesVerified(iFile)
        fileHashes(iFile) = computeFileSha256(files(iFile));
    end
end

provenance = struct();
provenance.sourceCommitSha = char(sourceCommitSha);
provenance.outputHashAlgorithm = 'SHA-256';
provenance.outputHashInputFiles = cellstr(files);
provenance.outputHashInputFilesVerified = all(filesVerified);
provenance.outputFileSha256 = cellstr(fileHashes);
if provenance.outputHashInputFilesVerified
    provenance.outputEvidenceHashSha256 = char(computeTextSha256( ...
        strjoin(files + "=" + fileHashes, newline)));
else
    provenance.outputEvidenceHashSha256 = '';
end
end

function sha = queryGitCommitSha()
sha = "";
moduleRoot = fileparts(fileparts(mfilename('fullpath')));
[status, output] = system(sprintf('git -C "%s" rev-parse HEAD', moduleRoot));
if status == 0
    candidate = string(strtrim(output));
    if strlength(candidate) > 0
        sha = candidate;
    end
end
end

function sha = computeFileSha256(path)
fid = fopen(path, 'r');
if fid < 0
    error('RTSPHEM:Precipitate:HashReadFailed', ...
        'Could not open file for hashing: %s.', path);
end
cleanupObj = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
clear cleanupObj;
sha = computeBytesSha256(bytes);
end

function sha = computeTextSha256(text)
sha = computeBytesSha256(uint8(char(text)));
end

function sha = computeBytesSha256(bytes)
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes(:));
rawHash = typecast(int8(digest.digest()), 'uint8');
sha = lower(string(reshape(dec2hex(rawHash, 2).', 1, [])));
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
