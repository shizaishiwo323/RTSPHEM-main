function config = ValidateReactiveTransportConfig(config)
%VALIDATEREACTIVETRANSPORTCONFIG Fill defaults and reject invalid RTM modes.

if nargin < 1 || isempty(config)
    config = struct();
end

config.solverArchitecture = char(getFieldOrDefault(config, ...
    'solverArchitecture', 'conservative_v2'));
config.strictMassConservation = logical(getFieldOrDefault(config, ...
    'strictMassConservation', strcmp(config.solverArchitecture, 'conservative_v2')));

config.chemistry = ensureStructField(config, 'chemistry');
config.chemistry.mode = canonicalChoice(getFieldOrDefault(config.chemistry, ...
    'mode', 'strict_molins'));
config.chemistry.chargeAbsoluteTolerance_eq = getFieldOrDefault( ...
    config.chemistry, 'chargeAbsoluteTolerance_eq', 1e-8);

config.transport = ensureStructField(config, 'transport');
config.transport.backend = canonicalChoice(getFieldOrDefault(config.transport, ...
    'backend', 'cut_cell_fv'));

config.time = ensureStructField(config, 'time');
config.time.mode = canonicalChoice(getFieldOrDefault(config.time, ...
    'mode', 'quasi_steady_geometry'));

config.cutCell = ensureStructField(config, 'cutCell');
config.cutCell.minFluidFraction = getFieldOrDefault(config.cutCell, ...
    'minFluidFraction', 0.1);
config.cutCell.smallCellMethod = char(getFieldOrDefault(config.cutCell, ...
    'smallCellMethod', 'disjoint_agglomeration'));

config.failure = ensureStructField(config, 'failure');
config.failure.maxRetries = getFieldOrDefault(config.failure, 'maxRetries', 12);
config.failure.shrinkFactor = getFieldOrDefault(config.failure, 'shrinkFactor', 0.5);
config.failure.minDt_s = getFieldOrDefault(config.failure, 'minDt_s', 1e-8);

config.mass = ensureStructField(config, 'mass');
config.mass.absoluteTolerance_mol = getFieldOrDefault(config.mass, ...
    'absoluteTolerance_mol', 1e-14);
config.mass.relativeTolerance = getFieldOrDefault(config.mass, ...
    'relativeTolerance', 1e-8);
config.mass.globalRelativeTolerance = getFieldOrDefault(config.mass, ...
    'globalRelativeTolerance', 1e-6);

config.phreeqc = ensureStructField(config, 'phreeqc');
config.phreeqc.engine = canonicalChoice(getFieldOrDefault(config.phreeqc, ...
    'engine', defaultPhreeqcEngine(config.chemistry.mode)));
config.phreeqc.databasePolicy = canonicalChoice(getFieldOrDefault(config.phreeqc, ...
    'databasePolicy', defaultDatabasePolicy(config.chemistry.mode)));
config.phreeqc.databaseName = char(getFieldOrDefault(config.phreeqc, ...
    'databaseName', 'phreeqc-m.dat'));
config.phreeqc.persistSession = logical(getFieldOrDefault(config.phreeqc, ...
    'persistSession', true));
config.phreeqc.useRunString = logical(getFieldOrDefault(config.phreeqc, ...
    'useRunString', true));

config.benchmark = ensureStructField(config, 'benchmark');
config.benchmark.enabled = logical(getFieldOrDefault(config.benchmark, 'enabled', false));

validateChoices(config);
validateCrossFieldRules(config);
end

function validateChoices(config)
mustBeOneOf(config.solverArchitecture, {'legacy', 'conservative_v2'}, ...
    'RTSPHEM:Config:UnknownSolverArchitecture');
mustBeOneOf(config.chemistry.mode, ...
    {'strict_molins', 'external_tst_phreeqc', 'phreeqc_kinetics'}, ...
    'RTSPHEM:Config:UnknownChemistryMode');
mustBeOneOf(config.transport.backend, {'cut_cell_fv', 'legacy_hyphm'}, ...
    'RTSPHEM:Config:UnknownTransportBackend');
mustBeOneOf(config.time.mode, {'quasi_steady_geometry', 'transient_snia'}, ...
    'RTSPHEM:Config:UnknownTimeMode');
mustBeOneOf(config.phreeqc.databasePolicy, ...
    {'exact_local', 'allow_fallback', 'not_used'}, ...
    'RTSPHEM:Config:UnknownDatabasePolicy');
mustBeOneOf(config.phreeqc.engine, ...
    {'none', 'iphreeqc_com', 'phreeqcrm', 'mock'}, ...
    'RTSPHEM:Config:UnknownPhreeqcEngine');
end

function validateCrossFieldRules(config)
if strcmp(config.chemistry.mode, 'strict_molins') && ~strcmp(config.phreeqc.engine, 'none')
    error('RTSPHEM:Config:StrictMolinsUsesNoPhreeqc', ...
        'strict_molins must not configure a PHREEQC engine.');
end
if config.benchmark.enabled && ~strcmp(config.phreeqc.engine, 'none') && ...
        ~strcmp(config.phreeqc.databasePolicy, 'exact_local')
    error('RTSPHEM:Config:BenchmarkRequiresExactLocalDatabase', ...
        'Benchmark runs with PHREEQC must use exact_local database policy.');
end
if strcmp(config.transport.backend, 'legacy_hyphm') && config.strictMassConservation
    error('RTSPHEM:Config:LegacyTransportNotStrict', ...
        'legacy_hyphm transport cannot claim strict mass conservation.');
end
if strcmp(config.solverArchitecture, 'conservative_v2') && ...
        strcmp(config.transport.backend, 'legacy_hyphm')
    error('RTSPHEM:Config:ConservativeRequiresCutCellTransport', ...
        'conservative_v2 requires cut_cell_fv transport.');
end
end

function value = defaultPhreeqcEngine(chemistryMode)
if strcmp(chemistryMode, 'strict_molins')
    value = 'none';
else
    value = 'iphreeqc_com';
end
end

function value = defaultDatabasePolicy(chemistryMode)
if strcmp(chemistryMode, 'strict_molins')
    value = 'not_used';
else
    value = 'exact_local';
end
end

function value = ensureStructField(config, fieldName)
if ~isfield(config, fieldName) || isempty(config.(fieldName))
    value = struct();
else
    value = config.(fieldName);
    if ~isstruct(value)
        error('RTSPHEM:Config:InvalidSection', ...
            'config.%s must be a struct.', fieldName);
    end
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function value = canonicalChoice(value)
value = lower(strrep(strtrim(char(value)), '-', '_'));
end

function mustBeOneOf(value, choices, errorId)
if ~any(strcmp(value, choices))
    error(errorId, 'Unsupported configuration value: %s.', value);
end
end
