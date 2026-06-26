function runtime = ResolveRuntime(rtmDir, config)
%RESOLVERUNTIME Resolve PHREEQC database and engine provenance for RTM.

if nargin < 1 || isempty(rtmDir)
    rtmDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2 || isempty(config)
    config = struct();
end

phreeqcConfig = getFieldOrDefault(config, 'phreeqc', config);
engineType = string(getFieldOrDefault(phreeqcConfig, 'engine', ...
    getFieldOrDefault(phreeqcConfig, 'engineType', 'iphreeqc_com')));
databasePolicy = string(getFieldOrDefault(phreeqcConfig, ...
    'databasePolicy', 'allow_fallback'));
databaseName = string(getFieldOrDefault(phreeqcConfig, ...
    'databaseName', 'phreeqc_rates.dat'));
configuredDatabasePath = string(getFieldOrDefault(phreeqcConfig, ...
    'databasePath', ""));
comProgId = string(getFieldOrDefault(phreeqcConfig, ...
    'comProgId', 'IPhreeqcCOM.Object'));

runtime = struct();
runtime.databasePath = "";
runtime.databaseName = databaseName;
runtime.databasePolicy = databasePolicy;
runtime.databaseSha256 = "";
runtime.databaseSizeBytes = 0;
runtime.engineType = engineType;
runtime.engineVersion = string(getFieldOrDefault(phreeqcConfig, ...
    'engineVersion', ""));
runtime.comProgId = comProgId;
runtime.isAvailable = false;

if strlength(configuredDatabasePath) > 0
    if exist(char(configuredDatabasePath), 'file') ~= 2
        error('RTSPHEM:Phreeqc:MissingConfiguredDatabase', ...
            'Configured PHREEQC database does not exist: %s.', ...
            configuredDatabasePath);
    end
    runtime.databasePath = configuredDatabasePath;
elseif ~strcmpi(databasePolicy, "not_used")
    runtime.databasePath = resolveDatabasePath(rtmDir, databaseName, ...
        databasePolicy);
end
[runtime.databaseSha256, runtime.databaseSizeBytes] = ...
    databaseManifest(runtime.databasePath);

runtime.isAvailable = engineAvailable(engineType, phreeqcConfig);
end

function databasePath = resolveDatabasePath(rtmDir, databaseName, databasePolicy)
normalizedPolicy = lower(strrep(strtrim(char(databasePolicy)), '-', '_'));
runtimeDatabaseDir = fullfile(char(rtmDir), 'phreeqc', 'database');
exactLocalPath = fullfile(runtimeDatabaseDir, char(databaseName));

switch normalizedPolicy
    case 'exact_local'
        if exist(exactLocalPath, 'file') == 2
            databasePath = string(exactLocalPath);
            return;
        end
        error('RTSPHEM:Phreeqc:MissingExactLocalDatabase', ...
            ['Cannot locate the exact local PHREEQC database required by ', ...
            'this run: %s. exact_local mode does not allow fallback ', ...
            'databases.'], exactLocalPath);
    case 'allow_fallback'
        databasePath = firstExistingDatabase({
            exactLocalPath
            fullfile(runtimeDatabaseDir, 'phreeqc_rates.dat')
            fullfile(runtimeDatabaseDir, 'phreeqc-m.dat')
            fullfile(runtimeDatabaseDir, 'phreeqc.dat')
            fullfile(char(rtmDir), 'couplePhreeqc', 'phreeqc-m.dat')
            char(databaseName)
            'C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat'
            });
    otherwise
        error('RTSPHEM:Phreeqc:UnknownDatabasePolicy', ...
            'Unknown PHREEQC database policy: %s.', normalizedPolicy);
end
end

function databasePath = firstExistingDatabase(candidates)
seen = strings(0, 1);
for iCandidate = 1:numel(candidates)
    candidate = string(candidates{iCandidate});
    if any(strcmpi(seen, candidate))
        continue;
    end
    seen(end + 1, 1) = candidate; %#ok<AGROW>
    if exist(char(candidate), 'file') == 2
        databasePath = candidate;
        return;
    end
end
error('RTSPHEM:Phreeqc:MissingDatabase', ...
    'Cannot locate a PHREEQC database for this runtime configuration.');
end

function [sha256, sizeBytes] = databaseManifest(databasePath)
sha256 = "";
sizeBytes = 0;
if strlength(databasePath) == 0 || exist(char(databasePath), 'file') ~= 2
    return;
end
fileInfo = dir(char(databasePath));
sizeBytes = fileInfo.bytes;
sha256 = string(computeSha256(char(databasePath)));
end

function tf = engineAvailable(engineType, phreeqcConfig)
normalized = lower(strrep(strtrim(char(engineType)), '-', '_'));
switch normalized
    case {'none', 'mock'}
        tf = true;
    case 'iphreeqc_com'
        tf = false;
        if logical(getFieldOrDefault(phreeqcConfig, ...
                'checkEngineAvailability', false))
            progId = char(getFieldOrDefault(phreeqcConfig, ...
                'comProgId', 'IPhreeqcCOM.Object'));
            try
                engine = actxserver(progId);
                cleanup = onCleanup(@() releaseComObject(engine));
                tf = ~isempty(engine);
                clear cleanup;
            catch
                tf = false;
            end
        end
    case 'phreeqcrm'
        tf = false;
    otherwise
        tf = false;
end
end

function releaseComObject(engine)
try
    delete(engine);
catch
end
end

function hash = computeSha256(filePath)
messageDigest = java.security.MessageDigest.getInstance('SHA-256');
inputStream = java.io.FileInputStream(java.io.File(filePath));
cleanup = onCleanup(@() inputStream.close());
buffer = zeros(8192, 1, 'int8');
while true
    bytesRead = inputStream.read(buffer, 0, numel(buffer));
    if bytesRead == -1
        break;
    end
    messageDigest.update(buffer, 0, bytesRead);
end
hashBytes = typecast(messageDigest.digest(), 'uint8');
hash = lower(reshape(dec2hex(hashBytes).', 1, []));
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
