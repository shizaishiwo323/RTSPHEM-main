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
    'mode', 'external_tst_phreeqc'));
config.chemistry.chargeAbsoluteTolerance_eq = nonnegativeScalarField( ...
    config.chemistry, 'chargeAbsoluteTolerance_eq', 1e-8, ...
    'RTSPHEM:Config:InvalidChemistryTolerance');
config.chemistry.calciteStoichiometryAbsoluteTolerance_mol = ...
    nonnegativeScalarField(config.chemistry, ...
    'calciteStoichiometryAbsoluteTolerance_mol', 1e-14, ...
    'RTSPHEM:Config:InvalidChemistryTolerance');
config.chemistry.calciteStoichiometryRelativeTolerance = ...
    nonnegativeScalarField(config.chemistry, ...
    'calciteStoichiometryRelativeTolerance', 1e-8, ...
    'RTSPHEM:Config:InvalidChemistryTolerance');
config.chemistry.semantics = string(getFieldOrDefault(config.chemistry, ...
    'semantics', chemistrySemantics(config.chemistry.mode)));
rawBasis = getFieldOrDefault(config.chemistry, 'basis', []);
if isempty(rawBasis) || isStaleStrictBasisForPhreeqc(rawBasis, config.chemistry.mode)
    rawBasis = defaultChemistryBasis(config.chemistry.mode);
end
config.chemistry.basis = normalizeNameList(rawBasis);
config.chemistry.derived = normalizeNameList(getFieldOrDefault( ...
    config.chemistry, 'derived', defaultChemistryDerived(config.chemistry.mode)));

config.transport = ensureStructField(config, 'transport');
config.transport.backend = canonicalChoice(getFieldOrDefault(config.transport, ...
    'backend', 'cut_cell_fv'));

config.time = ensureStructField(config, 'time');
config.time.mode = canonicalChoice(getFieldOrDefault(config.time, ...
    'mode', 'quasi_steady_geometry'));
config.time.rt = ensureStructField(config.time, 'rt');
config.time.rt.initialDt_s = positiveScalarField(config.time.rt, ...
    'initialDt_s', 0.01, 'RTSPHEM:Config:InvalidTimeConfig');
config.time.rt.maxDt_s = positiveScalarField(config.time.rt, ...
    'maxDt_s', 0.1, 'RTSPHEM:Config:InvalidTimeConfig', true);
config.time.rt.advectiveCfl = nonnegativeScalarField(config.time.rt, ...
    'advectiveCfl', 0.5, 'RTSPHEM:Config:InvalidTimeConfig');
config.time.rt.diffusiveNumber = nonnegativeScalarField(config.time.rt, ...
    'diffusiveNumber', 0.5, 'RTSPHEM:Config:InvalidTimeConfig');
config.time.rt.maxReactantFraction = nonnegativeScalarField(config.time.rt, ...
    'maxReactantFraction', 0.1, 'RTSPHEM:Config:InvalidTimeConfig', true);
config.time.geometry = ensureStructField(config.time, 'geometry');
config.time.geometry.boundaryCfl = nonnegativeScalarField( ...
    config.time.geometry, 'boundaryCfl', 0.25, ...
    'RTSPHEM:Config:InvalidTimeConfig');
config.time.geometry.maxMineralFraction = nonnegativeScalarField( ...
    config.time.geometry, 'maxMineralFraction', 0.2, ...
    'RTSPHEM:Config:InvalidTimeConfig', true);
config.time.geometry.maxDt_s = positiveScalarField(config.time.geometry, ...
    'maxDt_s', 60, 'RTSPHEM:Config:InvalidTimeConfig', true);

config.cutCell = ensureStructField(config, 'cutCell');
config.cutCell.minFluidFraction = closedUnitScalarField(config.cutCell, ...
    'minFluidFraction', 0.1, 'RTSPHEM:Config:InvalidCutCellConfig');
config.cutCell.smallCellMethod = canonicalChoice(getFieldOrDefault( ...
    config.cutCell, 'smallCellMethod', 'disjoint_agglomeration'));

config.failure = ensureStructField(config, 'failure');
config.failure.maxRetries = nonnegativeIntegerField(config.failure, ...
    'maxRetries', 12, 'RTSPHEM:Config:InvalidFailureRetry');
config.failure.shrinkFactor = openUnitScalarField(config.failure, ...
    'shrinkFactor', 0.5, 'RTSPHEM:Config:InvalidFailureRetry');
config.failure.minDt_s = positiveScalarField(config.failure, ...
    'minDt_s', 1e-8, 'RTSPHEM:Config:InvalidFailureRetry');

config.mass = ensureStructField(config, 'mass');
config.mass.absoluteTolerance_mol = nonnegativeScalarField(config.mass, ...
    'absoluteTolerance_mol', 1e-14, ...
    'RTSPHEM:Config:InvalidMassTolerance');
config.mass.relativeTolerance = nonnegativeScalarField(config.mass, ...
    'relativeTolerance', 1e-8, ...
    'RTSPHEM:Config:InvalidMassTolerance');
config.mass.globalRelativeTolerance = nonnegativeScalarField(config.mass, ...
    'globalRelativeTolerance', 1e-6, ...
    'RTSPHEM:Config:InvalidMassTolerance');

config.geometry = ensureStructField(config, 'geometry');
config.geometry.maxDisplacementOverH = nonnegativeScalarField( ...
    config.geometry, 'maxDisplacementOverH', 0.25, ...
    'RTSPHEM:Config:InvalidGeometryTolerance', true);
config.geometry.solidAbsoluteTolerance_cm3 = nonnegativeScalarField( ...
    config.geometry, 'solidAbsoluteTolerance_cm3', 1e-14, ...
    'RTSPHEM:Config:InvalidGeometryTolerance');
config.geometry.solidRelativeTolerance = nonnegativeScalarField( ...
    config.geometry, 'solidRelativeTolerance', 1e-8, ...
    'RTSPHEM:Config:InvalidGeometryTolerance');

config.flow = ensureStructField(config, 'flow');
config.flow.absoluteTolerance_cm3_s = nonnegativeScalarField( ...
    config.flow, 'absoluteTolerance_cm3_s', 1e-12, ...
    'RTSPHEM:Config:InvalidFlowTolerance');
config.flow.relativeTolerance = nonnegativeScalarField( ...
    config.flow, 'relativeTolerance', 1e-8, ...
    'RTSPHEM:Config:InvalidFlowTolerance');

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
config.units = rtm.units.NormalizeRtmUnits(getFieldOrDefault(config, ...
    'units', struct()));

validateChoices(config);
validateChemistryBasis(config);
validateCrossFieldRules(config);
end

function validateChoices(config)
mustBeOneOf(config.solverArchitecture, {'legacy', 'conservative_v2'}, ...
    'RTSPHEM:Config:UnknownSolverArchitecture');
mustBeOneOf(config.chemistry.mode, ...
    {'external_tst_phreeqc', 'phreeqc_kinetics'}, ...
    'RTSPHEM:Config:UnknownChemistryMode');
mustBeOneOf(config.transport.backend, {'cut_cell_fv', 'legacy_hyphm'}, ...
    'RTSPHEM:Config:UnknownTransportBackend');
mustBeOneOf(config.time.mode, ...
    {'quasi_steady_geometry', 'fixed_geometry_steady_rt', 'transient_snia'}, ...
    'RTSPHEM:Config:UnknownTimeMode');
mustBeOneOf(config.cutCell.smallCellMethod, {'disjoint_agglomeration'}, ...
    'RTSPHEM:Config:UnknownSmallCellMethod');
mustBeOneOf(config.phreeqc.databasePolicy, ...
    {'exact_local', 'allow_fallback', 'not_used'}, ...
    'RTSPHEM:Config:UnknownDatabasePolicy');
mustBeOneOf(config.phreeqc.engine, ...
    {'none', 'iphreeqc_com', 'phreeqcrm', 'mock'}, ...
    'RTSPHEM:Config:UnknownPhreeqcEngine');
end

function validateCrossFieldRules(config)
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

function validateChemistryBasis(config)
requiredBasis = {'Ca', 'C', 'Na', 'Cl'};
missing = setdiff(requiredBasis, config.chemistry.basis);
if ~isempty(missing)
    error('RTSPHEM:Config:InvalidChemistryBasis', ...
        'PHREEQC chemistry basis must include Ca, C, Na, and Cl totals.');
end
forbidden = lower(strrep(strrep({'H+', 'pH', 'HCO3-', 'CO3-2'}, ...
    ' ', ''), '_', ''));
canonicalBasis = lower(strrep(strrep(config.chemistry.basis, ' ', ''), '_', ''));
if any(ismember(canonicalBasis, forbidden))
    error('RTSPHEM:Config:InvalidChemistryBasis', ...
        'PHREEQC chemistry basis cannot include derived H/pH/carbonate species.');
end
end

function names = defaultChemistryBasis(chemistryMode)
switch chemistryMode
    case {'external_tst_phreeqc', 'phreeqc_kinetics'}
        names = {'Ca', 'C', 'Na', 'Cl', 'Alkalinity'};
    otherwise
        names = {};
end
end

function names = defaultChemistryDerived(chemistryMode)
switch chemistryMode
    case {'external_tst_phreeqc', 'phreeqc_kinetics'}
        names = {'pH', 'H+', 'HCO3-', 'CO3-2', 'SI_Calcite'};
    otherwise
        names = {};
end
end

function names = normalizeNameList(value)
if isstring(value)
    names = cellstr(value(:).');
elseif ischar(value)
    names = {value};
elseif iscell(value)
    names = value;
else
    error('RTSPHEM:Config:InvalidChemistryBasis', ...
        'chemistry basis and derived lists must be char, string, or cell arrays.');
end
names = cellfun(@char, names(:).', 'UniformOutput', false);
names = cellfun(@strtrim, names, 'UniformOutput', false);
if any(cellfun(@isempty, names))
    error('RTSPHEM:Config:InvalidChemistryBasis', ...
        'chemistry basis and derived lists cannot contain empty names.');
end
end

function tf = isStaleStrictBasisForPhreeqc(value, chemistryMode)
tf = false;
try
    names = normalizeNameList(value);
catch
    return;
end
tf = isequal(names, {'H_reactant'});
end

function value = defaultPhreeqcEngine(chemistryMode)
value = 'iphreeqc_com';
end

function value = defaultDatabasePolicy(chemistryMode)
value = 'exact_local';
end

function value = chemistrySemantics(chemistryMode)
switch chemistryMode
    case 'external_tst_phreeqc'
        value = "explicit external TST rate + PHREEQC equilibrium closure";
    case 'phreeqc_kinetics'
        value = "PHREEQC native KINETICS/RATES";
    otherwise
        value = "unknown";
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

function value = nonnegativeScalarField(structValue, fieldName, defaultValue, errorId, allowInf)
if nargin < 5
    allowInf = false;
end
value = getFieldOrDefault(structValue, fieldName, defaultValue);
if allowInf
    isValid = isscalar(value) && (isfinite(value) || isinf(value)) && value >= 0;
else
    isValid = isscalar(value) && isfinite(value) && value >= 0;
end
if ~isValid
    error(errorId, 'config.%s must be a nonnegative finite scalar.', fieldName);
end
end

function value = positiveScalarField(structValue, fieldName, defaultValue, errorId, allowInf)
if nargin < 5
    allowInf = false;
end
value = getFieldOrDefault(structValue, fieldName, defaultValue);
if allowInf
    isValid = isscalar(value) && (isfinite(value) || isinf(value)) && value > 0;
else
    isValid = isscalar(value) && isfinite(value) && value > 0;
end
if ~isValid
    error(errorId, 'config.%s must be a positive finite scalar.', fieldName);
end
end

function value = openUnitScalarField(structValue, fieldName, defaultValue, errorId)
value = getFieldOrDefault(structValue, fieldName, defaultValue);
if ~(isscalar(value) && isfinite(value) && value > 0 && value < 1)
    error(errorId, 'config.%s must be between 0 and 1.', fieldName);
end
end

function value = closedUnitScalarField(structValue, fieldName, defaultValue, errorId)
value = getFieldOrDefault(structValue, fieldName, defaultValue);
if ~(isscalar(value) && isfinite(value) && value >= 0 && value <= 1)
    error(errorId, 'config.%s must be between 0 and 1.', fieldName);
end
end

function value = nonnegativeIntegerField(structValue, fieldName, defaultValue, errorId)
value = getFieldOrDefault(structValue, fieldName, defaultValue);
if ~(isscalar(value) && isfinite(value) && value >= 0 && value == round(value))
    error(errorId, 'config.%s must be a nonnegative integer scalar.', fieldName);
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
