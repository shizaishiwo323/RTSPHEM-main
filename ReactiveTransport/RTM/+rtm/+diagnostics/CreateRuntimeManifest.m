function manifest = CreateRuntimeManifest(config, overrides)
%CREATERUNTIMEMANIFEST Build a provenance manifest structure for RTM runs.

if nargin < 2 || isempty(overrides)
    overrides = struct();
end

manifest = struct();
manifest.generated_at = string(datetime('now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
manifest.solver_architecture = string(getFieldOrDefault(config, 'solverArchitecture', 'legacy'));
manifest.chemistry_mode = string(getNestedField(config, {'chemistry', 'mode'}, ...
    getFieldOrDefault(config, 'chemistryMode', 'unknown')));
manifest.chemistry = normalizeManifestValue(getFieldOrDefault(config, ...
    'chemistry', struct()));
manifest.transport_backend = string(getNestedField(config, {'transport', 'backend'}, 'legacy_hyphm'));
manifest.operator_order = string(getNestedField(config, {'time', 'mode'}, 'transient_snia'));
manifest.mass_tolerances = getFieldOrDefault(config, 'mass', struct());
manifest.geometry = normalizeManifestValue(getFieldOrDefault(config, 'geometry', struct()));
manifest.time = normalizeManifestValue(getFieldOrDefault(config, 'time', struct()));
manifest.failure = normalizeManifestValue(getFieldOrDefault(config, 'failure', struct()));

manifest.phreeqc = phreeqcManifest(config);
manifest.git = gitManifest(overrides);
manifest.matlab_version = string(version);
manifest.platform = string(computer);
manifest.units = unitsManifest();
end

function value = phreeqcManifest(config)
databasePath = string(getNestedField(config, {'phreeqc', 'databasePath'}, ""));
value = struct();
value.engine = string(getNestedField(config, {'phreeqc', 'engine'}, ...
    getFieldOrDefault(config, 'phreeqcEngine', 'unknown')));
value.engine_version = string(getNestedField(config, {'phreeqc', 'engineVersion'}, ""));
value.com_progid = string(getNestedField(config, {'phreeqc', 'comProgId'}, ""));
value.persist_session = logical(getNestedField(config, {'phreeqc', 'persistSession'}, false));
value.use_run_string = logical(getNestedField(config, {'phreeqc', 'useRunString'}, false));
value.database_policy = string(getNestedField(config, {'phreeqc', 'databasePolicy'}, ...
    getFieldOrDefault(config, 'phreeqcDatabasePolicy', 'unknown')));
value.database_name = string(getNestedField(config, {'phreeqc', 'databaseName'}, ""));
value.database_path = databasePath;
value.database_sha256 = "";
value.database_size_bytes = 0;

if strlength(databasePath) > 0 && exist(char(databasePath), 'file') == 2
    fileInfo = dir(char(databasePath));
    value.database_size_bytes = fileInfo.bytes;
    value.database_sha256 = string(computeSha256(char(databasePath)));
end
end

function value = gitManifest(overrides)
value = struct();
value.commit = string(getFieldOrDefault(overrides, 'gitCommit', ""));
value.branch = string(getFieldOrDefault(overrides, 'gitBranch', ""));
if value.commit == ""
    [status, text] = system('git rev-parse HEAD');
    if status == 0
        value.commit = string(strtrim(text));
    end
end
if value.branch == ""
    [status, text] = system('git branch --show-current');
    if status == 0
        value.branch = string(strtrim(text));
    end
end
end

function value = unitsManifest()
value = struct();
value.length = "cm";
value.area = "cm^2";
value.volume = "cm^3";
value.time = "s";
value.component_state = "mol";
value.concentration = "mol/cm^3";
value.interface_rate = "mol/cm^2/s";
value.cell_rate = "mol/s";
end

function value = normalizeManifestValue(value)
if isa(value, 'function_handle')
    value = string(sprintf('<function_handle:%s>', func2str(value)));
elseif ischar(value)
    value = string(value);
elseif iscell(value)
    for iCell = 1:numel(value)
        value{iCell} = normalizeManifestValue(value{iCell});
    end
elseif isstruct(value)
    fieldNames = fieldnames(value);
    for iValue = 1:numel(value)
        for iField = 1:numel(fieldNames)
            fieldName = fieldNames{iField};
            value(iValue).(fieldName) = normalizeManifestValue(value(iValue).(fieldName));
        end
    end
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

function value = getNestedField(structValue, path, defaultValue)
value = structValue;
for iPath = 1:numel(path)
    if ~isstruct(value) || ~isfield(value, path{iPath})
        value = defaultValue;
        return;
    end
    value = value.(path{iPath});
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end
