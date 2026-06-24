function result = precip_RunPhreeqcCalciteBatchSigned(state, options)
% precip_RunPhreeqcCalciteBatchSigned - Run PHREEQC and keep signed calcite.

arguments
    state struct
    options struct
end

activeCells = selectActivePhreeqcCells(state, options);
result = initializeInactiveResult(state, options);

workDir = char(getOption(options, 'workDir', tempdir));
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

stepIndex = getOption(options, 'timeStepIndex', 0);
inputPath = fullfile(workDir, sprintf('phreeqc_calcite_signed_step_%04d.phr', stepIndex));
result.inputPath = string(inputPath);
result.databasePath = string(char(getOption(options, 'databasePath', 'phreeqc.dat')));

if ~any(activeCells)
    return;
end

activeState = subsetState(state, activeCells);
inputText = precip_BuildCalcitePhreeqcInputSigned(activeState, options);
fid = fopen(inputPath, 'w');
if fid == -1
    error('RTSPHEM:Precipitate:InputOpenFailed', 'Cannot write PHREEQC input: %s', inputPath);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s', inputText);
clear cleaner;

databasePath = char(getOption(options, 'databasePath', 'phreeqc.dat'));
iphreeqc = actxserver('IPhreeqcCOM.Object');
cleanupPhreeqc = onCleanup(@() releaseComObject(iphreeqc));
iphreeqc.LoadDatabase(databasePath);
iphreeqc.RunFile(inputPath);
rawOutput = iphreeqc.GetSelectedOutputArray;

dt = getOption(options, 'timeStepSize', 1);
activeResult = precip_ParsePhreeqcSelectedOutputSigned(rawOutput, nnz(activeCells), dt);
activeResult = precip_ScaleSignedCalciteDeltaToCellInventory(activeResult, activeState, options);
activeResult = precip_ApplyPrescribedCalciteDissolutionSigned(activeResult, activeState, options);
result = mergeActiveResult(result, activeResult, activeCells);
result.inputPath = string(inputPath);
result.databasePath = string(databasePath);
clear cleanupPhreeqc;
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
activeCells = hasWater & (hasChemistry | (reactNeutralInterfaceCells & hasInterface));
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
result.calcitePrecipitatedMoles = zeros(numCells, 1);
result.calciteDissolvedMoles = zeros(numCells, 1);
result.calciteSignedRate_mol_s = zeros(numCells, 1);
result.calciteRate_mol_s = zeros(numCells, 1);
result.calciteRate_mol_dm2_s = zeros(numCells, 1);
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

function values = optionalColumn(state, fieldName, numCells, defaultValue)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    values = state.(fieldName)(:);
else
    values = repmat(defaultValue, numCells, 1);
end
if numel(values) ~= numCells
    error('RTSPHEM:Precipitate:StateSizeMismatch', ...
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
