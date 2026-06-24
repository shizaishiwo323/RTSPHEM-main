function result = precip_ParsePhreeqcSelectedOutputSigned(rawOutput, expectedRows, timeStepSize)
% precip_ParsePhreeqcSelectedOutputSigned - Parse PHREEQC output preserving signs.

if isempty(rawOutput) || size(rawOutput, 1) < 2
    error('RTSPHEM:Precipitate:EmptySelectedOutput', 'PHREEQC selected output is empty.');
end
if nargin < 2 || isempty(expectedRows)
    expectedRows = size(rawOutput, 1) - 1;
end
if nargin < 3
    timeStepSize = [];
end

headings = string(rawOutput(1, :));
data = rawOutput(2:end, :);
if size(data, 1) < expectedRows
    error('RTSPHEM:Precipitate:SelectedOutputTooShort', ...
        'PHREEQC selected output has %d rows, expected at least %d.', size(data, 1), expectedRows);
end
data = data((end - expectedRows + 1):end, :);

result = struct();
result.solutionNumber = numericColumn(data, headings, ["soln", "solution"], NaN);
result.pH = numericColumn(data, headings, "pH", NaN);
result.chargeBalance = numericColumn(data, headings, ["charge_balance", "charge", "charge(eq)"], NaN);
result.ca_total_mol_cm3 = numericColumn(data, headings, "Ca", 0) / 1000;
result.c_total_mol_cm3 = numericColumn(data, headings, "C", 0) / 1000;
result.na_total_mol_cm3 = numericColumn(data, headings, "Na", 0) / 1000;
result.cl_total_mol_cm3 = numericColumn(data, headings, "Cl", 0) / 1000;
result.h_mol_cm3 = numericColumn(data, headings, ["m_H+", "H+"], 0) / 1000;
result.ca_mol_cm3 = numericColumn(data, headings, ["m_Ca+2", "Ca+2"], 0) / 1000;
result.hco3_mol_cm3 = numericColumn(data, headings, ["m_HCO3-", "HCO3-"], 0) / 1000;
result.co3_mol_cm3 = numericColumn(data, headings, ["m_CO3-2", "CO3-2"], 0) / 1000;
result.cl_mol_cm3 = numericColumn(data, headings, ["m_Cl-", "Cl-"], 0) / 1000;
result.na_mol_cm3 = numericColumn(data, headings, ["m_Na+", "Na+"], 0) / 1000;
result.calciteSI = numericColumn(data, headings, ["si_Calcite", "si_calcite", "Calcite"], NaN);
result.calciteDeltaMoles = requiredNumericColumn(data, headings, ...
    ["KIN_DELTA_Calcite", "kin_delta_Calcite"], 'KIN_DELTA_Calcite');
result.calcitePrecipitatedMoles = max(result.calciteDeltaMoles, 0);
result.calciteDissolvedMoles = max(-result.calciteDeltaMoles, 0);

rawRate = requiredNumericColumn(data, headings, ...
    ["RATE_Calcite", "rate_Calcite"], 'RATE_Calcite');
dissolutionPositiveRate = max(rawRate, 0);
if ~isempty(timeStepSize) && isfinite(timeStepSize) && timeStepSize > 0
    result.calciteSignedRate_mol_s = result.calciteDeltaMoles ./ timeStepSize;
else
    result.calciteSignedRate_mol_s = -rawRate;
end
result.calciteDissolutionRate_mol_s = dissolutionPositiveRate;
result.calciteRate_mol_s = dissolutionPositiveRate;
result.calciteRate_mol_dm2_s = result.calciteRate_mol_s;
end

function values = requiredNumericColumn(data, headings, candidates, label)
idx = findHeading(headings, candidates);
if isnan(idx)
    error('RTSPHEM:Precipitate:MissingSelectedOutputHeading', ...
        'PHREEQC selected output is missing required heading: %s.', label);
end

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
end

function values = numericColumn(data, headings, candidates, defaultValue)
idx = findHeading(headings, candidates);
if isnan(idx)
    values = repmat(defaultValue, size(data, 1), 1);
    return;
end

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
end

function idx = findHeading(headings, candidates)
candidates = string(candidates);
normalizedHeadings = lower(regexprep(headings, '\s+', ''));
normalizedHeadings = regexprep(normalizedHeadings, '\(.*?\)', '');
for iCandidate = 1:numel(candidates)
    normalizedCandidate = lower(regexprep(candidates(iCandidate), '\s+', ''));
    normalizedCandidate = regexprep(normalizedCandidate, '\(.*?\)', '');
    matches = find(normalizedHeadings == normalizedCandidate, 1, 'first');
    if ~isempty(matches)
        idx = matches;
        return;
    end
end
idx = NaN;
end
