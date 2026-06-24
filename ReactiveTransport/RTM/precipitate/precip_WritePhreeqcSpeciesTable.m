function precip_WritePhreeqcSpeciesTable(outputDir, stepIndex, timeSeconds, grid, speciesData)
% precip_WritePhreeqcSpeciesTable - Save triangle-wise signed PHREEQC output.

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

numRows = grid.numT;
timestep = repmat(stepIndex, numRows, 1);
time_s = repmat(timeSeconds, numRows, 1);
triangle = (1:numRows)';
x_cm = grid.baryT(:, 1);
y_cm = grid.baryT(:, 2);

T = table(timestep, time_s, triangle, x_cm, y_cm);
fields = {'h_mol_cm3', 'h_activity_mol_cm3', 'ca_mol_cm3', 'hco3_mol_cm3', 'co3_mol_cm3', ...
    'cl_mol_cm3', 'na_mol_cm3', 'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'cl_total_mol_cm3', 'na_total_mol_cm3', 'pH', 'chargeBalance', ...
    'calciteSI', 'calciteRate_mol_s', 'calciteDissolutionRate_mol_s', ...
    'calciteRatePerArea_mol_cm2_s', 'calciteSignedRate_mol_s', ...
    'calciteDeltaMoles', 'calcitePrecipitatedMoles', 'calciteDissolvedMoles'};
for iField = 1:numel(fields)
    fieldName = fields{iField};
    if isfield(speciesData, fieldName)
        T.(fieldName) = speciesData.(fieldName)(:);
    end
end

writetable(T, fullfile(outputDir, sprintf('phreeqc_species_%04d.csv', stepIndex)));
end
