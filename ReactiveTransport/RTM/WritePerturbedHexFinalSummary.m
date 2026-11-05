function summary = WritePerturbedHexFinalSummary(geometryDir, smokeSummaryCsv)
% Attach experimental provenance and RTM-loaded values to the geometry files.
samples = ["M14_24";"M14_25";"M14_26";"M14_34"];
experimentalPorosity = [0.382;0.374;0.342;0.356];
experimentalGasK_mD = [124.3;107.6;59.7;108.7];
smoke = readtable(smokeSummaryCsv, 'TextType','string');
geometry = readtable(fullfile(geometryDir, 'geometry_validation.csv'), 'TextType','string');
assert(isequal(string(smoke.sampleName), samples) && isequal(string(geometry.sampleName), samples), ...
    'RTM:GeometrySummary:Samples', 'Expected the four M14 samples in order.');
summary = table(samples, experimentalPorosity, smoke.porosity, ...
    smoke.porosity-experimentalPorosity, experimentalGasK_mD, smoke.initialKxx_mD, ...
    geometry.Iv, geometry.minimumGap_cm*1e4, geometry.cellsPerMinRadius, ...
    geometry.leftRightConnected, ...
    'VariableNames', {'sampleName','experimentalPorosity','rtmPorosity', ...
    'rtmPorosityError','experimentalGasK_mD','twoDimensionalKxx_mD','Iv', ...
    'minimumGap_um','cellsPerMinimumRadius','leftRightConnected'});
for i = 1:height(summary)
    file = fullfile(geometryDir, sprintf('%s_perturbed_hex_600x400um.mat', samples(i)));
    s = load(file); %#ok<NASGU>
    designPorosity = s.targetInitialPorosity; %#ok<NASGU>
    rtmPorosity = summary.rtmPorosity(i); %#ok<NASGU>
    rtmKxx_mD = summary.twoDimensionalKxx_mD(i); %#ok<NASGU>
    experimentalPorosityValue = experimentalPorosity(i); %#ok<NASGU>
    experimentalGasK_mDValue = experimentalGasK_mD(i); %#ok<NASGU>
    calibrationScope = "porosity_only_2d_k_not_fitted_to_gas_permeability"; %#ok<NASGU>
    save(file, 'designPorosity','rtmPorosity','rtmKxx_mD', ...
        'experimentalPorosityValue','experimentalGasK_mDValue','calibrationScope','-append');
end
writetable(summary, fullfile(geometryDir, 'final_calibration_summary.csv'));
writetable(summary, fullfile(geometryDir, 'final_calibration_summary.xlsx'));
end
