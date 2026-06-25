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
manifest.flow = normalizeManifestValue(getFieldOrDefault(config, 'flow', struct()));
manifest.time = normalizeManifestValue(getFieldOrDefault(config, 'time', struct()));
manifest.failure = normalizeManifestValue(getFieldOrDefault(config, 'failure', struct()));
manifest.acceptance_tolerances = acceptanceTolerancesManifest(config);

manifest.phreeqc = phreeqcManifest(config);
manifest.git = gitManifest(overrides);
manifest.matlab_version = string(version);
manifest.platform = string(computer);
manifest.units = unitsManifest(config);
end

function value = acceptanceTolerancesManifest(config)
value = struct();
value.mass_absolute_tolerance_mol = getNestedField(config, ...
    {'mass', 'absoluteTolerance_mol'}, 1e-14);
value.mass_relative_tolerance = getNestedField(config, ...
    {'mass', 'relativeTolerance'}, 1e-8);
value.mass_global_relative_tolerance = getNestedField(config, ...
    {'mass', 'globalRelativeTolerance'}, 1e-6);
value.solid_absolute_tolerance_cm3 = getNestedField(config, ...
    {'geometry', 'solidAbsoluteTolerance_cm3'}, 1e-14);
value.solid_relative_tolerance = getNestedField(config, ...
    {'geometry', 'solidRelativeTolerance'}, 1e-8);
value.charge_absolute_tolerance_eq = getNestedField(config, ...
    {'chemistry', 'chargeAbsoluteTolerance_eq'}, Inf);
value.calcite_stoichiometry_absolute_tolerance_mol = getNestedField(config, ...
    {'chemistry', 'calciteStoichiometryAbsoluteTolerance_mol'}, 1e-14);
value.calcite_stoichiometry_relative_tolerance = getNestedField(config, ...
    {'chemistry', 'calciteStoichiometryRelativeTolerance'}, 1e-8);
value.max_displacement_over_h = getNestedField(config, ...
    {'geometry', 'maxDisplacementOverH'}, 0.25);
value.flow_absolute_tolerance_cm3_s = getNestedField(config, ...
    {'flow', 'absoluteTolerance_cm3_s'}, 1e-12);
value.flow_relative_tolerance = getNestedField(config, ...
    {'flow', 'relativeTolerance'}, Inf);
end

function value = phreeqcManifest(config)
rtmDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
runtime = rtm.phreeqc.ResolveRuntime(rtmDir, config);
value = struct();
value.engine = runtime.engineType;
value.engine_type = runtime.engineType;
value.engine_version = runtime.engineVersion;
value.com_progid = runtime.comProgId;
value.is_available = logical(runtime.isAvailable);
value.persist_session = logical(getNestedField(config, {'phreeqc', 'persistSession'}, false));
value.use_run_string = logical(getNestedField(config, {'phreeqc', 'useRunString'}, false));
value.database_policy = runtime.databasePolicy;
value.database_name = runtime.databaseName;
value.database_path = runtime.databasePath;
value.database_sha256 = runtime.databaseSha256;
value.database_size_bytes = runtime.databaseSizeBytes;
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

function value = unitsManifest(config)
value = rtm.units.NormalizeRtmUnits(getFieldOrDefault(config, ...
    'units', struct()));
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
