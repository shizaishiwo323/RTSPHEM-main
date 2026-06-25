function [chem, metadata] = precip_IPhreeqcSpeciation(samples, spec)
% precip_IPhreeqcSpeciation - Zero-dimensional PHREEQC aqueous speciation.
%
% Inputs:
%   samples - Yoon conservative component samples in mol/cm^3.
%   spec    - benchmark spec with phreeqcDatabasePath and optional
%             iphreeqcSession fields.
%
% Output:
%   chem    - pH, activities, and Vaterite saturation state.
%   metadata - PHREEQC session/database provenance.

if nargin < 2 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end

databasePath = resolveDatabasePath(spec);
session = getFieldOrDefault(spec, 'iphreeqcSession', []);
ownsSession = isempty(session);
if ownsSession
    runtimeConfig = getFieldOrDefault(spec, 'iphreeqcRuntimeConfig', struct());
    session = rtm.phreeqc.PhreeqcSession(runtimeConfig);
    cleanupSession = onCleanup(@() session.close());
else
    cleanupSession = onCleanup(@() []);
end

inputText = buildSpeciationInput(samples, spec);
session.loadDatabaseExact(databasePath);
rawOutput = session.runString(inputText);
chem = parseSelectedOutput(rawOutput, numel(samples.Ca_total));
manifest = session.getDatabaseManifest();

metadata = struct();
metadata.backend = 'iphreeqc_speciation';
metadata.role = 'aqueous_speciation';
metadata.source = ['iphreeqc_session:' char(manifest.engineType)];
metadata.databasePath = manifest.databasePath;
metadata.databaseSha256 = manifest.databaseSha256;
metadata.numSolutions = numel(samples.Ca_total);
metadata.inputText = string(inputText);
clear cleanupSession;
end

function inputText = buildSpeciationInput(samples, spec)
numSolutions = numel(samples.Ca_total);
temperatureC = getFieldOrDefault(spec, 'temperature_C', 25);
equilibriumPH = precip_YoonCarbonateEquilibrium(samples, spec).pH;

lines = strings(0, 1);
lines(end + 1) = "TITLE RTSPHEM Yoon Vaterite aqueous speciation";
lines(end + 1) = "PHASES";
lines(end + 1) = "Vaterite";
lines(end + 1) = "CaCO3 = Ca+2 + CO3-2";
lines(end + 1) = sprintf('log_k %.15g', log10(spec.vateriteKsp));
for iSolution = 1:numSolutions
    fixedPH = NaN;
    if isfield(samples, 'fixedPH') && numel(samples.fixedPH) >= iSolution
        fixedPH = samples.fixedPH(iSolution);
    end

    lines(end + 1) = sprintf('SOLUTION %d', iSolution);
    lines(end + 1) = sprintf('temp %.15g', temperatureC);
    if isfinite(fixedPH)
        lines(end + 1) = sprintf('pH %.15g', fixedPH);
    else
        lines(end + 1) = sprintf('pH %.15g', equilibriumPH(iSolution));
    end
    lines(end + 1) = "units mol/l";
    lines(end + 1) = sprintf('Ca %.15g', max(samples.Ca_total(iSolution), 0) * 1000);
    lines(end + 1) = sprintf('C %.15g', max(samples.C_total(iSolution), 0) * 1000);
    lines(end + 1) = sprintf('Na %.15g', max(samples.Na_total(iSolution), 0) * 1000);
    lines(end + 1) = sprintf('Cl %.15g', max(samples.Cl_total(iSolution), 0) * 1000);
    lines(end + 1) = sprintf('Alkalinity %.15g', samples.Alkalinity(iSolution) * 1000);
end

lines(end + 1) = "SELECTED_OUTPUT";
lines(end + 1) = "-reset false";
lines(end + 1) = "-simulation true";
lines(end + 1) = "-state true";
lines(end + 1) = "-solution true";
lines(end + 1) = "-pH true";
lines(end + 1) = "-charge_balance true";
lines(end + 1) = "-totals Ca C Na Cl";
lines(end + 1) = "-molalities H+ Ca+2 CO3-2 H2CO3";
lines(end + 1) = "-saturation_indices Vaterite";
lines(end + 1) = "USER_PUNCH";
lines(end + 1) = "-headings la_H+ la_Ca+2 la_CO3-2 la_H2CO3";
lines(end + 1) = "-start";
lines(end + 1) = "10 PUNCH LA(""H+""), LA(""Ca+2""), LA(""CO3-2""), LA(""H2CO3"")";
lines(end + 1) = "-end";
lines(end + 1) = "END";
lines(end + 1) = "RUN_CELLS";
lines(end + 1) = sprintf('-cells 1-%d', numSolutions);
lines(end + 1) = "END";

inputText = char(strjoin(lines, newline));
end

function chem = parseSelectedOutput(rawOutput, expectedRows)
if isempty(rawOutput) || size(rawOutput, 1) < 2
    error('RTSPHEM:Precipitate:InvalidPhreeqcSpeciationOutput', ...
        'PHREEQC selected output is empty.');
end

headings = string(rawOutput(1, :));
data = rawOutput(2:end, :);
if size(data, 1) < expectedRows
    error('RTSPHEM:Precipitate:InvalidPhreeqcSpeciationOutput', ...
        'PHREEQC selected output has %d rows, expected %d.', ...
        size(data, 1), expectedRows);
end
data = data((end - expectedRows + 1):end, :);

chem = struct();
chem.pH = requireNumericColumn(data, headings, "pH");
chem.aH = activityColumn(data, headings, "la_H+", ["m_H+", "H+"]);
chem.aCa = activityColumn(data, headings, "la_Ca+2", ["m_Ca+2", "Ca+2"]);
chem.aCO3 = activityColumn(data, headings, "la_CO3-2", ["m_CO3-2", "CO3-2"]);
chem.aH2CO3 = activityColumn(data, headings, "la_H2CO3", ["m_H2CO3", "H2CO3"]);
chem.saturationIndexVaterite = requireNumericColumn(data, headings, ...
    ["si_Vaterite", "si_vaterite", "Vaterite"]);
chem.omegaVaterite = 10 .^ chem.saturationIndexVaterite;
end

function values = activityColumn(data, headings, logActivityNames, molalityNames)
logActivityIndex = findHeading(headings, logActivityNames);
if ~isnan(logActivityIndex)
    values = 10 .^ numericColumnByIndex(data, logActivityIndex);
    return;
end
values = requireNumericColumn(data, headings, molalityNames);
end

function values = requireNumericColumn(data, headings, candidates)
idx = findHeading(headings, candidates);
if isnan(idx)
    error('RTSPHEM:Precipitate:InvalidPhreeqcSpeciationOutput', ...
        'PHREEQC selected output is missing required heading: %s.', ...
        strjoin(string(candidates), ', '));
end
values = numericColumnByIndex(data, idx);
end

function values = numericColumnByIndex(data, idx)
values = zeros(size(data, 1), 1);
for iRow = 1:size(data, 1)
    item = data{iRow, idx};
    if isnumeric(item)
        values(iRow) = item;
    elseif ischar(item) || isstring(item)
        values(iRow) = str2double(item);
    else
        values(iRow) = NaN;
    end
end
if any(~isfinite(values))
    error('RTSPHEM:Precipitate:InvalidPhreeqcSpeciationOutput', ...
        'PHREEQC selected output contains non-finite numeric values.');
end
end

function idx = findHeading(headings, candidates)
candidates = string(candidates);
normalizedHeadings = normalizeHeading(headings);
for iCandidate = 1:numel(candidates)
    normalizedCandidate = normalizeHeading(candidates(iCandidate));
    matches = find(normalizedHeadings == normalizedCandidate, 1, 'first');
    if ~isempty(matches)
        idx = matches;
        return;
    end
end
idx = NaN;
end

function value = normalizeHeading(value)
value = lower(regexprep(string(value), '\s+', ''));
value = regexprep(value, '\(.*?\)', '');
end

function databasePath = resolveDatabasePath(spec)
databasePath = getFieldOrDefault(spec, 'phreeqcDatabasePath', "");
if strlength(string(databasePath)) > 0
    databasePath = char(databasePath);
    return;
end
databasePath = getFieldOrDefault(spec, 'databasePath', "");
if strlength(string(databasePath)) > 0
    databasePath = char(databasePath);
    return;
end

chemistryDir = fileparts(mfilename('fullpath'));
moduleDir = fileparts(chemistryDir);
rtmDir = fileparts(moduleDir);
if exist('ResolvePhreeqcDatabasePath', 'file') == 2
    databasePath = char(ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc.dat'));
    return;
end

candidates = {
    fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc.dat')
    fullfile(rtmDir, 'phreeqc', 'database', 'phreeqc-m.dat')
    'C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat'
    };
for iCandidate = 1:numel(candidates)
    if exist(candidates{iCandidate}, 'file') == 2
        databasePath = candidates{iCandidate};
        return;
    end
end

error('RTSPHEM:Precipitate:MissingPhreeqcDatabase', ...
    'Cannot locate a PHREEQC database for Yoon speciation.');
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
