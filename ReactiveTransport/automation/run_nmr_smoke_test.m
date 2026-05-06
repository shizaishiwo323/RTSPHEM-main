function result = run_nmr_smoke_test(folder_path, timestep_index, calibration_factor)
% run_nmr_smoke_test - Non-interactive NMR smoke test for one RTM output.
%
% Inputs
%   folder_path     RTM result folder containing dxf_pore/dxf_solid.
%   timestep_index      1-based paired DXF index to process. Defaults to 1.
%   calibration_factor  Optional T2 calibration factor reused from step 1.
%
% Outputs
%   result          Struct with COMSOL/T2 status and output paths.

if nargin < 1 || strlength(string(folder_path)) == 0
    error('folder_path is required.');
end
if nargin < 2 || isempty(timestep_index)
    timestep_index = 1;
end
if nargin < 3
    calibration_factor = [];
end

automation_path = fileparts(mfilename('fullpath'));
addpath(automation_path);

folder_path = char(folder_path);
config = AutomationConfig();
config.overwrite_existing = true;
config.export_mph = false;
config.enable_gif = false;

result = struct();
result.folder_path = folder_path;
result.timestep_index = timestep_index;
result.comsol_success = false;
result.inversion_success = false;
result.excel_output = '';
result.inversion_dir = '';
result.total_water = NaN;
result.raw_spectrum_sum = NaN;
result.calibration_factor = NaN;

if ~exist(folder_path, 'dir')
    error('RTM folder does not exist: %s', folder_path);
end

fprintf('NMR smoke test folder: %s\n', folder_path);

params = parse_folder_name(folder_path);
fprintf('Parsed params: Da=%.4g, Pe=%.4g, L=%.4g, X=%.6g, Y=%.6g, layout=%s\n', ...
    params.Da, params.Pe, params.L, params.lengthXAxis, params.lengthYAxis, params.layoutType);

[pore_files, solid_files] = get_dxf_files(folder_path);
if isempty(pore_files)
    error('No paired DXF files found in: %s', folder_path);
end
if timestep_index < 1 || timestep_index > numel(pore_files)
    error('timestep_index %d is out of range 1..%d.', timestep_index, numel(pore_files));
end

comsol_output_dir = fullfile(folder_path, 'comsol_results');
inversion_output_dir = fullfile(folder_path, 'inversion_results');
if ~exist(comsol_output_dir, 'dir')
    mkdir(comsol_output_dir);
end
if ~exist(inversion_output_dir, 'dir')
    mkdir(inversion_output_dir);
end

pore_dxf = fullfile(folder_path, 'dxf_pore', pore_files(timestep_index).name);
solid_dxf = fullfile(folder_path, 'dxf_solid', solid_files(timestep_index).name);
timestep = extract_timestep(pore_files(timestep_index).name);

excel_filename = sprintf('T2_Da%.4f_Pe%.4f_X%.4f_Y%.4f_t%s.xlsx', ...
    params.Da, params.Pe, params.lengthXAxis, params.lengthYAxis, timestep);
excel_output = fullfile(comsol_output_dir, excel_filename);

result.excel_output = excel_output;
result.inversion_dir = inversion_output_dir;

fprintf('Processing DXF pair %d/%d: %s <-> %s\n', ...
    timestep_index, numel(pore_files), pore_files(timestep_index).name, solid_files(timestep_index).name);
fprintf('COMSOL output: %s\n', excel_output);

result.comsol_success = run_comsol_processing( ...
    config.mph_file, pore_dxf, solid_dxf, ...
    params.lengthXAxis, params.lengthYAxis, excel_output, config);

if ~result.comsol_success
    fprintf('NMR smoke test stopped: COMSOL processing failed.\n');
    return;
end

first_porosity = NaN;
global_evolution_file = fullfile(folder_path, 'global_evolution.xlsx');
if exist(global_evolution_file, 'file')
    try
        global_data = readtable(global_evolution_file);
        porosity_col_idx = find(contains(lower(global_data.Properties.VariableNames), 'porosity'), 1);
        if ~isempty(porosity_col_idx) && height(global_data) >= 1
            first_porosity = global_data{1, porosity_col_idx};
            fprintf('First RTM porosity for calibration reference: %.6f\n', first_porosity);
        end
    catch ME
        fprintf('Could not read global_evolution.xlsx for calibration reference: %s\n', ME.message);
    end
end

[result.inversion_success, result.total_water, result.raw_spectrum_sum, calibration_factor] = ...
    run_python_inversion(excel_output, inversion_output_dir, config, calibration_factor);

if result.inversion_success && timestep_index == 1 && ...
        isfinite(first_porosity) && isfinite(result.raw_spectrum_sum) && result.raw_spectrum_sum > 0
    calibration_factor = first_porosity / result.raw_spectrum_sum;
    fprintf('Re-running T2 inversion with RTM porosity calibration factor: %.6e\n', calibration_factor);
    [result.inversion_success, result.total_water, result.raw_spectrum_sum, calibration_factor] = ...
        run_python_inversion(excel_output, inversion_output_dir, config, calibration_factor);
end
result.calibration_factor = calibration_factor;

fprintf('NMR smoke test summary: COMSOL=%d, inversion=%d, total_water=%.6g\n', ...
    result.comsol_success, result.inversion_success, result.total_water);
end
