function ExportPhreeqcSpeciesPlots(speciesDir, stepIndex, timeSeconds, grid, levelSetData, speciesData, plotOptions)
% ExportPhreeqcSpeciesPlots - Save per-species concentration maps for PHREEQC mode.

if ~exist(speciesDir, 'dir')
    mkdir(speciesDir);
end

species = {
    'H', 'h_mol_cm3', 'H+ concentration';
    'Ca', 'ca_mol_cm3', 'Ca2+ concentration';
    'HCO3', 'hco3_mol_cm3', 'HCO3- concentration';
    'CO3', 'co3_mol_cm3', 'CO3^{2-} concentration';
    'Cl', 'cl_mol_cm3', 'Cl- concentration';
    'Na', 'na_mol_cm3', 'Na+ concentration';
    'pH', 'pH', 'pH';
    'calcite_rate', 'calciteRatePerArea_mol_cm2_s', 'Calcite reaction rate [mol cm^{-2} s^{-1}]'
};

fontName = getOption(plotOptions, 'fontName', 'Arial');
fontSize = getOption(plotOptions, 'fontSize', 12);
titleFontSize = getOption(plotOptions, 'titleFontSize', 13);

for iSpecies = 1:size(species, 1)
    label = species{iSpecies, 1};
    fieldName = species{iSpecies, 2};
    titleText = species{iSpecies, 3};
    if ~isfield(speciesData, fieldName)
        continue;
    end

    outDir = fullfile(speciesDir, label);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    values = speciesData.(fieldName);
    if strcmp(fieldName, 'pH')
        faceData = values(:);
        levelSetOnTriangles = levelSetData(grid.V0T);
        poreTriangles = all(levelSetOnTriangles < 0, 2);
        faceData(~poreTriangles) = NaN;
    else
        faceData = PrepareConcentrationFaceData(values(:), grid.V0T, levelSetData);
    end

    fig = figure('Visible', 'off', 'Position', [100, 100, 800, 700]);
    ax = axes(fig);
    patch(ax, 'Faces', grid.V0T, ...
        'Vertices', [grid.coordV, zeros(grid.numV, 1)], ...
        'FaceVertexCData', faceData, ...
        'FaceColor', 'flat', 'EdgeColor', 'none');
    title(ax, sprintf('%s at t = %.2f s', titleText, timeSeconds), ...
        'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    xlabel(ax, 'X (cm)', 'FontSize', fontSize, 'FontName', fontName);
    ylabel(ax, 'Y (cm)', 'FontSize', fontSize, 'FontName', fontName);
    set(ax, 'FontSize', fontSize, 'FontName', fontName);
    view(ax, 2);
    axis(ax, 'equal', 'tight');
    colorbar(ax);
    print(fig, fullfile(outDir, sprintf('%s_%04d.png', label, stepIndex)), '-dpng', '-r300');
    close(fig);
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
