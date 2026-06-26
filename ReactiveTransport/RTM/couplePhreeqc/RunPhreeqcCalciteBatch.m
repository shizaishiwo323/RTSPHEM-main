function result = RunPhreeqcCalciteBatch(state, options)
% RunPhreeqcCalciteBatch - Run one PHREEQC calcite chemistry batch.

arguments
    state struct
    options struct
end

numCells = numel(state.h_mol_cm3);
activeCells = selectActivePhreeqcCells(state, options);
result = initializeInactiveResult(state, options);

workDir = char(getOption(options, 'workDir', tempdir));
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

stepIndex = getOption(options, 'timeStepIndex', 0);
inputPath = fullfile(workDir, sprintf('phreeqc_calcite_step_%04d.phr', stepIndex));
result.inputPath = string(inputPath);
databasePath = char(getOption(options, 'databasePath', 'phreeqc.dat'));
result.databasePath = string(databasePath);
hasSession = isfield(options, 'phreeqcSession') && ~isempty(options.phreeqcSession);
writeInputFiles = logical(getOption(options, 'writeInputFiles', ~hasSession));
result.inputWritten = false;
result.phreeqcSessionReused = false;
result.phreeqcRunMethod = "RunFile";

if ~any(activeCells)
    result = addReactionMassLedger(result, state);
    return;
end

activeState = subsetState(state, activeCells);
inputText = BuildCalcitePhreeqcInput(activeState, options);
if writeInputFiles
    fid = fopen(inputPath, 'w');
    if fid == -1
        error('RTSPHEM:Phreeqc:InputOpenFailed', 'Cannot write PHREEQC input: %s', inputPath);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', inputText);
    clear cleaner;
    result.inputWritten = true;
end

if hasSession
    session = options.phreeqcSession;
    session.loadDatabaseExact(databasePath);
    rawOutput = session.runString(inputText);
    manifest = session.getDatabaseManifest();
    result.databaseSha256 = manifest.databaseSha256;
    result.phreeqcSessionReused = true;
    result.phreeqcRunMethod = "RunString";
else
    if ~writeInputFiles
        error('RTSPHEM:Phreeqc:InputFileRequired', ...
            'RunFile PHREEQC mode requires writeInputFiles=true.');
    end
    [iphreeqc, ownsEngine] = createIPhreeqcEngine(options);
    if ownsEngine
        cleanupPhreeqc = onCleanup(@() releaseComObject(iphreeqc));
    else
        cleanupPhreeqc = onCleanup(@() []);
    end
    loadStatus = invokePhreeqcCommand(iphreeqc, 'LoadDatabase', ...
        {databasePath}, 'RTSPHEM:Phreeqc:LoadDatabaseFailed');
    assertPhreeqcStatus(loadStatus, iphreeqc, ...
        'RTSPHEM:Phreeqc:LoadDatabaseFailed', ...
        sprintf('PHREEQC LoadDatabase failed for %s', databasePath));
    runStatus = invokePhreeqcCommand(iphreeqc, 'RunFile', {inputPath}, ...
        'RTSPHEM:Phreeqc:RunFileFailed');
    assertPhreeqcStatus(runStatus, iphreeqc, ...
        'RTSPHEM:Phreeqc:RunFileFailed', ...
        sprintf('PHREEQC RunFile failed for %s', inputPath));
    rawOutput = iphreeqc.GetSelectedOutputArray;
    clear cleanupPhreeqc;
end

activeResult = ParsePhreeqcSelectedOutput(rawOutput, nnz(activeCells));
activeResult = applyHydrogenActivityForTstMatch(activeResult, options);
activeResult = scaleKineticDissolutionToCellInventory(activeResult, activeState, options);
activeResult = DiagnosticInferCalciteDissolutionFromTotals( ...
    activeResult, activeState, getOption(options, 'timeStepSize', 1));
activeResult = applyPrescribedCalciteDissolution(activeResult, activeState, options);
result = mergeActiveResult(result, activeResult, activeCells);
result.inputPath = string(inputPath);
result.databasePath = string(databasePath);
result = addReactionMassLedger(result, state);
end

function result = applyHydrogenActivityForTstMatch(result, options)
defaultRateLaw = getOption(options, 'phreeqcRateLaw', 'database_calcite');
rateLaw = getOption(options, 'rateLaw', defaultRateLaw);
rateLaw = lower(strrep(strtrim(char(rateLaw)), '-', '_'));
if ~ismember(rateLaw, {'tst_match', 'calcite_tst_match', 'phreeqc_tst_match', ...
        'external_tst_phreeqc', 'legacy_phreeqc_tst_match'})
    return;
end
if isfield(result, 'pH') && ~isempty(result.pH)
    result.h_activity_mol_cm3 = 10 .^ (-result.pH(:)) ./ 1000;
else
    result.h_activity_mol_cm3 = result.h_mol_cm3;
end
end

function result = scaleKineticDissolutionToCellInventory(result, state, options)
numCells = numel(result.h_mol_cm3);
timeStepSize = getOption(options, 'timeStepSize', 1);
if ~isfield(result, 'calciteDeltaMoles') || isempty(result.calciteDeltaMoles)
    return;
end
if ~isfield(result, 'calciteKinDeltaRate_mol_s') || ...
        isempty(result.calciteKinDeltaRate_mol_s)
    result.calciteKinDeltaRate_mol_s = result.calciteRate_mol_s;
end
if ~isfield(result, 'calciteRawKineticDeltaMoles') || ...
        isempty(result.calciteRawKineticDeltaMoles)
    result.calciteRawKineticDeltaMoles = result.calciteDeltaMoles;
end
result.calciteRawKinDeltaRate_mol_s = max(-result.calciteRawKineticDeltaMoles(:), 0) ./ ...
    max(timeStepSize, eps);

waterVolume = optionalColumn(state, 'water_volume_cm3', numCells, 0);
calciteMoles = optionalColumn(state, 'calcite_moles', numCells, Inf);
solutionWaterKg = getOption(options, 'solutionWaterKg', 1);

scale = max(waterVolume(:), 0) * 1e-3 ./ max(solutionWaterKg, eps);
rawDissolvedMoles = max(-result.calciteDeltaMoles(:), 0);
dissolvedMoles = rawDissolvedMoles .* scale;
dissolvedMoles = min(dissolvedMoles, max(calciteMoles(:), 0));
dissolvedMoles(~isfinite(dissolvedMoles)) = 0;

result.calciteDissolvedMoles = dissolvedMoles;
result.calciteDeltaMoles = -dissolvedMoles;
result.calciteRate_mol_s = dissolvedMoles ./ max(timeStepSize, eps);
result.calciteKinDeltaRate_mol_s = result.calciteRate_mol_s;
result.calcite_cell_rate_mol_s = result.calciteRate_mol_s;
end

function result = applyPrescribedCalciteDissolution(result, state, options)
numCells = numel(result.h_mol_cm3);
if ~isfield(state, 'prescribed_calcite_dissolved_moles') || ...
        isempty(state.prescribed_calcite_dissolved_moles)
    return;
end

timeStepSize = getOption(options, 'timeStepSize', 1);
prescribed = optionalColumn(state, 'prescribed_calcite_dissolved_moles', numCells, 0);
dissolvedMoles = max(prescribed(:), 0);
dissolvedMoles(~isfinite(dissolvedMoles)) = 0;

result.calciteDissolvedMoles = dissolvedMoles;
result.calciteDeltaMoles = -dissolvedMoles;
result.calciteRate_mol_s = dissolvedMoles ./ max(timeStepSize, eps);
result.calciteKinDeltaRate_mol_s = result.calciteRate_mol_s;
result.calcite_cell_rate_mol_s = result.calciteRate_mol_s;
end

function activeCells = selectActivePhreeqcCells(state, options)
numCells = numel(state.h_mol_cm3);
interfaceArea = optionalColumn(state, 'interface_area_cm2', numCells, 0);
waterVolume = optionalColumn(state, 'water_volume_cm3', numCells, 0);
calciteMoles = optionalColumn(state, 'calcite_moles', numCells, 0);
ca = optionalColumn(state, 'ca_mol_cm3', numCells, 0);
c = optionalColumn(state, 'c_mol_cm3', numCells, 0);
na = optionalColumn(state, 'na_mol_cm3', numCells, 0);
cl = optionalColumn(state, 'cl_mol_cm3', numCells, 0);
h = optionalColumn(state, 'h_mol_cm3', numCells, 0);
prescribed = optionalColumn(state, 'prescribed_calcite_dissolved_moles', numCells, 0);
concentrationTol = getOption(options, 'activeConcentrationToleranceMolCm3', 1e-30);
minActiveWaterVolume = getOption(options, 'minActiveWaterVolumeCm3', 0);
minActiveWaterVolumeFraction = getOption(options, 'minActiveWaterVolumeFraction', 0);
reactNeutralInterfaceCells = logical(getOption(options, 'reactNeutralInterfaceCells', false));
positiveWaterVolume = waterVolume(waterVolume > 0);
if minActiveWaterVolumeFraction > 0 && ~isempty(positiveWaterVolume)
    minActiveWaterVolume = max(minActiveWaterVolume, max(positiveWaterVolume) * minActiveWaterVolumeFraction);
end
hasWater = waterVolume(:) > minActiveWaterVolume;
hasInterface = interfaceArea(:) > 0 & calciteMoles(:) > 0;
hasChemistry = any([h(:), ca(:), c(:), na(:), cl(:)] > concentrationTol, 2);
hasPrescribedReaction = prescribed(:) > 0;
isExternalTstClosure = usesExternalTstClosure(options);
includeExternalTstBackgroundCells = logical(getOption( ...
    options, 'activeExternalTstChemistryOnlyCells', false));
if isExternalTstClosure && ~includeExternalTstBackgroundCells
    activeCells = hasWater & (hasPrescribedReaction ...
        | (reactNeutralInterfaceCells & hasInterface));
else
    activeCells = hasWater & (hasChemistry | hasPrescribedReaction ...
        | (reactNeutralInterfaceCells & hasInterface));
end
end

function tf = usesExternalTstClosure(options)
defaultRateLaw = getOption(options, 'phreeqcRateLaw', 'database_calcite');
rateLaw = getOption(options, 'rateLaw', defaultRateLaw);
rateLaw = lower(strrep(strtrim(char(rateLaw)), '-', '_'));
tf = ismember(rateLaw, {'tst_match', 'calcite_tst_match', ...
    'phreeqc_tst_match', 'external_tst_phreeqc', ...
    'legacy_phreeqc_tst_match'});
end

function activeState = subsetState(state, activeCells)
activeState = struct();
fields = fieldnames(state);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    value = state.(fieldName);
    if isnumeric(value) && isvector(value) && numel(value) == numel(activeCells)
        activeState.(fieldName) = value(activeCells);
    else
        activeState.(fieldName) = value;
    end
end
end

function result = initializeInactiveResult(state, options)
numCells = numel(state.h_mol_cm3);
h = optionalColumn(state, 'h_mol_cm3', numCells, 0);
ca = optionalColumn(state, 'ca_mol_cm3', numCells, 0);
c = optionalColumn(state, 'c_mol_cm3', numCells, 0);
na = optionalColumn(state, 'na_mol_cm3', numCells, 0);
cl = optionalColumn(state, 'cl_mol_cm3', numCells, 0);
minHForPH = getOption(options, 'minHForPHMolL', 1e-7);

result = struct();
result.pH = -log10(max(h(:) * 1000, minHForPH));
result.chargeBalance = zeros(numCells, 1);
result.ca_total_mol_cm3 = max(ca(:), 0);
result.c_total_mol_cm3 = max(c(:), 0);
result.na_total_mol_cm3 = max(na(:), 0);
result.cl_total_mol_cm3 = max(cl(:), 0);
result.h_mol_cm3 = max(h(:), 0);
result.ca_mol_cm3 = result.ca_total_mol_cm3;
result.hco3_mol_cm3 = zeros(numCells, 1);
result.co3_mol_cm3 = zeros(numCells, 1);
result.cl_mol_cm3 = result.cl_total_mol_cm3;
result.na_mol_cm3 = result.na_total_mol_cm3;
result.calciteSI = NaN(numCells, 1);
result.calciteDeltaMoles = zeros(numCells, 1);
result.calciteDissolvedMoles = zeros(numCells, 1);
result.calciteKinDeltaRate_mol_s = zeros(numCells, 1);
result.calciteRawKinDeltaRate_mol_s = zeros(numCells, 1);
result.calciteRawKineticDeltaMoles = zeros(numCells, 1);
result.calciteKineticReactantMoles = NaN(numCells, 1);
result.calciteRate_mol_s = zeros(numCells, 1);
result.calcite_cell_rate_mol_s = zeros(numCells, 1);
result.solutionNumber = (1:numCells)';
end

function result = mergeActiveResult(result, activeResult, activeCells)
fields = fieldnames(activeResult);
numCells = numel(activeCells);
numActive = nnz(activeCells);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    value = activeResult.(fieldName);
    if isnumeric(value) && isvector(value) && numel(value) == numActive
        if isfield(result, fieldName)
            merged = result.(fieldName);
        else
            merged = NaN(numCells, 1);
        end
        merged(activeCells) = value(:);
        result.(fieldName) = merged;
    else
        result.(fieldName) = value;
    end
end
result.solutionNumber = (1:numCells)';
end

function result = addReactionMassLedger(result, state)
numCells = numel(result.h_mol_cm3);
waterVolume = optionalColumn(state, 'water_volume_cm3', numCells, 0);
componentSpecs = {
    'h', 'h_mol_cm3', 'h_mol_cm3'
    'ca', 'ca_total_mol_cm3', 'ca_mol_cm3'
    'c', 'c_total_mol_cm3', 'c_mol_cm3'
    'na', 'na_total_mol_cm3', 'na_mol_cm3'
    'cl', 'cl_total_mol_cm3', 'cl_mol_cm3'
    };

for iSpec = 1:size(componentSpecs, 1)
    name = componentSpecs{iSpec, 1};
    resultField = componentSpecs{iSpec, 2};
    stateFallbackField = componentSpecs{iSpec, 3};
    delta = computeComponentDeltaMoles( ...
        result, state, waterVolume, resultField, stateFallbackField);
    result.(sprintf('water_phase_%s_delta_moles', name)) = delta;
    result.(sprintf('water_phase_%s_delta_moles_total', name)) = ...
        sum(delta, 'omitnan');
end

dissolved = optionalColumn(result, 'calciteDissolvedMoles', numCells, 0);
caDelta = result.water_phase_ca_delta_moles(:);
cDelta = result.water_phase_c_delta_moles(:);
result.calcite_ca_stoich_residual_moles = caDelta - dissolved(:);
result.calcite_c_stoich_residual_moles = cDelta - dissolved(:);
result.calcite_ca_stoich_residual_moles_total = ...
    sum(result.calcite_ca_stoich_residual_moles, 'omitnan');
result.calcite_c_stoich_residual_moles_total = ...
    sum(result.calcite_c_stoich_residual_moles, 'omitnan');
result.calcite_stoich_max_abs_residual_moles = max(abs([ ...
    result.calcite_ca_stoich_residual_moles(:); ...
    result.calcite_c_stoich_residual_moles(:)]), [], 'omitnan');
if isempty(result.calcite_stoich_max_abs_residual_moles)
    result.calcite_stoich_max_abs_residual_moles = 0;
end
end

function delta = computeComponentDeltaMoles( ...
    result, state, waterVolume, resultField, stateFallbackField)
numCells = numel(waterVolume);
if isfield(result, resultField) && ~isempty(result.(resultField))
    after = result.(resultField)(:);
else
    after = zeros(numCells, 1);
end
if isfield(state, resultField) && ~isempty(state.(resultField))
    before = state.(resultField)(:);
elseif isfield(state, stateFallbackField) && ~isempty(state.(stateFallbackField))
    before = state.(stateFallbackField)(:);
else
    before = zeros(numCells, 1);
end
if numel(after) ~= numCells || numel(before) ~= numCells
    delta = NaN(numCells, 1);
    return;
end
delta = (after - before) .* waterVolume(:);
delta(~isfinite(delta)) = NaN;
end

function values = optionalColumn(state, fieldName, numCells, defaultValue)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    values = state.(fieldName)(:);
else
    values = repmat(defaultValue, numCells, 1);
end
if numel(values) ~= numCells
    error('RTSPHEM:Phreeqc:StateSizeMismatch', ...
        'State field %s has %d values, expected %d.', fieldName, numel(values), numCells);
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function releaseComObject(obj)
try
    delete(obj);
catch
end
end

function [engine, ownsEngine] = createIPhreeqcEngine(options)
factory = getOption(options, 'engineFactory', []);
if ~isempty(factory)
    engine = factory();
    ownsEngine = false;
else
    engine = actxserver('IPhreeqcCOM.Object');
    ownsEngine = true;
end
end

function status = invokePhreeqcCommand(engine, methodName, args, errorId)
try
    status = engine.(methodName)(args{:});
catch ME
    error(errorId, '%s threw an error: %s', methodName, ME.message);
end
if isempty(status)
    status = 0;
end
end

function assertPhreeqcStatus(status, engine, errorId, message)
if double(status) == 0
    return;
end
errorString = "";
try
    errorString = string(engine.GetErrorString());
catch
end
if strlength(errorString) > 0
    error(errorId, '%s: %s', message, char(errorString));
else
    error(errorId, '%s with status %g.', message, double(status));
end
end
