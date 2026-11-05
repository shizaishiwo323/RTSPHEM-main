function summary = RunPerturbedHexFaceDissolution(geometryDir, resultsRoot)
% Run the fixed perturbed-hex sample geometries with the face-dissolution setup.
if nargin < 1 || isempty(geometryDir)
    error('RTM:PerturbedHex:GeometryDir', 'geometryDir is required.');
end
if nargin < 2 || isempty(resultsRoot)
    resultsRoot = fullfile(fileparts(geometryDir), 'rtm_face_dissolution');
end
files = dir(fullfile(geometryDir, '*_perturbed_hex_600x400um.mat'));
assert(numel(files) == 4, 'RTM:PerturbedHex:GeometryCount', 'Expected four MAT geometries.');
rows = cell(numel(files), 6);
for i = 1:numel(files)
    geometryFile = fullfile(files(i).folder, files(i).name);
    sampleName = erase(string(files(i).name), '_perturbed_hex_600x400um.mat');
    resultDir = fullfile(resultsRoot, char(sampleName));
    metadataFile = fullfile(resultDir, 'run_metadata.json');
    chemistryFile = fullfile(resultDir, 'phreeqc_results', 'seismoelectric_chemistry.csv');
    if exist(metadataFile, 'file') == 2 && exist(chemistryFile, 'file') == 2
        metadata = jsondecode(fileread(metadataFile));
        if isfield(metadata, 'final') && ~metadata.final.timedOut
            chemistry = readtable(chemistryFile);
            rows(i,:) = {sampleName, "SkippedComplete", height(chemistry), ...
                metadata.final.finalPorosity, metadata.final.finalPermeability_mD, string(resultDir)};
            continue;
        end
    end
    cfg = struct( ...
        'layoutType', 'random', 'loadExistingGeometry', true, 'geometryLoadFile', geometryFile, ...
        'resultsDir', resultDir, ...
        'inletVelocity', 1e-4, 'acidConcentration_mol_L', 0.02, ...
        'rateCoefficientTST', 0.00135501355013550, ...
        'phreeqcKineticsCorrectionFactor', 0.00135501355013550 / 8.9e-5, ...
        'numPartitionsMicroscale', 256, 'meshNumPartitionsX', 256, 'meshNumPartitionsY', 171, ...
        'targetDissolutionSlices', 34, 'targetDissolutionFinalPorosity', 0.78, ...
        'phreeqcExportEvery', Inf, 'phreeqcExportSeismoelectric', true, ...
        'saveRealtimePlot', false, 'showDebugFigures', false, ...
        'transportSolverBackend', 'cut_cell_fv', 'phreeqcTransportBasis', 'conserved_totals');
    result = run_phreeqc_se(cfg);
    chemistry = readtable(result.seismoelectricChemistryCsv);
    status = "Success";
    if result.timedOut || height(chemistry) < cfg.targetDissolutionSlices
        status = "Incomplete";
    end
    rows(i,:) = {sampleName, status, height(chemistry), result.finalPorosity, ...
        result.finalPermeability, string(result.resultsDir)};
end
summary = cell2table(rows, 'VariableNames', {'sampleName','status','exportedSlices', ...
    'finalPorosity','finalKxx_mD','resultsDir'});
writetable(summary, fullfile(resultsRoot, 'face_dissolution_summary.csv'));
writetable(summary, fullfile(resultsRoot, 'face_dissolution_summary.xlsx'));
assert(all(summary.status == "Success" | summary.status == "SkippedComplete"), ...
    'RTM:PerturbedHex:Incomplete', 'One or more face-dissolution cases did not finish.');
end
