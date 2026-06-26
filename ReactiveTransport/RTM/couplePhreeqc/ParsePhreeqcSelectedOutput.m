function result = ParsePhreeqcSelectedOutput(rawOutput, expectedRows)
% ParsePhreeqcSelectedOutput - Convert IPhreeqc selected output to RTM fields.

if isempty(rawOutput) || size(rawOutput, 1) < 2
    error('RTSPHEM:Phreeqc:EmptySelectedOutput', 'PHREEQC selected output is empty.');
end
if nargin < 2 || isempty(expectedRows)
    expectedRows = size(rawOutput, 1) - 1;
end

headings = string(rawOutput(1, :));
data = rawOutput(2:end, :);
if size(data, 1) < expectedRows
    error('RTSPHEM:Phreeqc:SelectedOutputTooShort', ...
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
result.alkalinity_mol_cm3 = numericColumn(data, headings, ...
    ["Alkalinity", "alk"], 0) / 1000;
result.h_mol_cm3 = numericColumn(data, headings, ["m_H+", "H+"], 0) / 1000;
result.ca_mol_cm3 = numericColumn(data, headings, ["m_Ca+2", "Ca+2"], 0) / 1000;
result.hco3_mol_cm3 = numericColumn(data, headings, ["m_HCO3-", "HCO3-"], 0) / 1000;
result.co3_mol_cm3 = numericColumn(data, headings, ["m_CO3-2", "CO3-2"], 0) / 1000;
result.cl_mol_cm3 = numericColumn(data, headings, ["m_Cl-", "Cl-"], 0) / 1000;
result.na_mol_cm3 = numericColumn(data, headings, ["m_Na+", "Na+"], 0) / 1000;
result.calciteSI = numericColumn(data, headings, ["si_Calcite", "si_calcite", "Calcite"], NaN);
legacyDelta = numericColumn(data, headings, ["KIN_DELTA_Calcite", "kin_delta_Calcite"], NaN);
kineticDelta = numericColumn(data, headings, ...
    ["dk_Calcite", "d_Calcite", "delta_Calcite", "dk_Calcite(mol)"], NaN);
result.calciteKineticReactantMoles = numericColumn(data, headings, ...
    ["k_Calcite", "kin_Calcite", "Calcite"], NaN);
result.calciteRawKineticDeltaMoles = kineticDelta;
result.calciteRawUserPunchDeltaMoles = legacyDelta;
result.calciteDeltaMoles = preferFiniteColumn(kineticDelta, legacyDelta, 0);
result.calciteDissolvedMoles = max(-result.calciteDeltaMoles, 0);
result.calciteRawUserPunchRate_mol_s = numericColumn(data, headings, ...
    ["RATE_Calcite", "rate_Calcite"], NaN);
result.calciteRate_mol_s = preferFiniteColumn(result.calciteRawUserPunchRate_mol_s, [], 0);
result.calciteKinDeltaRate_mol_s = result.calciteRate_mol_s;
result.calcite_cell_rate_mol_s = result.calciteRate_mol_s;
end

function values = preferFiniteColumn(primary, secondary, defaultValue)
if nargin < 2 || isempty(secondary)
    secondary = [];
end
values = primary(:);
if isempty(values)
    values = zeros(size(secondary(:)));
end
if ~isempty(secondary)
    secondary = secondary(:);
    if numel(values) ~= numel(secondary)
        values = repmat(defaultValue, numel(secondary), 1);
    end
    replace = ~isfinite(values);
    values(replace) = secondary(replace);
end
values(~isfinite(values)) = defaultValue;
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
