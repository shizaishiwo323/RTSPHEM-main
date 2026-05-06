%% 处理单个Dissolution结果文件夹
% ============================================================================
% 功能: 处理指定的单个dissolution结果文件夹
%       适用于测试或单独处理特定数据
%
% 用法:
%   修改下方的 target_folder 变量指定要处理的文件夹
%
% 作者: 自动化框架
% 日期: 2025-11-29
% ============================================================================

clear all;
close all;
clc;

%% ========== 添加自动化脚本路径 ==========
automation_path = fileparts(mfilename('fullpath'));
addpath(automation_path);

%% ========== 配置 ==========
config = AutomationConfig();

% ===== 修改这里指定要处理的文件夹 =====
% 支持两种格式:
%   1. 绝对路径: 'C:\...\outputs\rtm_runs\rtm_...'
%   2. 相对路径(文件夹名): 'rtm_...' 或旧版 'dissolution_results-Da_...'
% 留空 "" 时自动选择 config.data_root 中按名称排序的最后一个RTM结果文件夹。
target_folder = "";

% ===== 采样参数设置 =====
max_samples = 10;  % 最多处理的DXF对数量 (设为inf表示处理全部)
% =========================================

% 判断是绝对路径还是相对路径
if strlength(string(target_folder)) == 0
    folders = scan_dissolution_folders(config.data_root);
    if isempty(folders)
        error('data_root 中未找到RTM结果文件夹: %s', config.data_root);
    end
    folder_name = folders(end).name;
    folder_path = fullfile(config.data_root, folder_name);
elseif exist(target_folder, 'dir')
    % 绝对路径，直接使用
    folder_path = target_folder;
    % 注意：fileparts 会把小数点后的部分当作扩展名，需要合并 name 和 ext
    [~, name_part, ext_part] = fileparts(target_folder);
    folder_name = [name_part, ext_part];  % 合并避免截断
else
    % 相对路径，拼接数据根目录
    folder_path = fullfile(config.data_root, target_folder);
    folder_name = target_folder;
end

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║  单文件夹处理模式                                          ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% ========== 检查文件夹 ==========
if ~exist(folder_path, 'dir')
    error('文件夹不存在: %s', folder_path);
end

fprintf('处理文件夹: %s\n', folder_path);
fprintf('文件夹名称: %s\n\n', folder_name);

%% ========== 解析参数 ==========
params = parse_folder_name(folder_path, config.metadata_filename);

fprintf('解析参数:\n');
fprintf('  参数来源:          %s\n', describe_parameter_source(params));
fprintf('  Damkohler数 (Da): %.4f\n', params.Da);
fprintf('  Peclet数 (Pe):    %.4f\n', params.Pe);
fprintf('  特征长度 (L):     %.4f\n', params.L);
fprintf('  X轴长度:          %.4f cm\n', params.lengthXAxis);
fprintf('  Y轴长度:          %.4f cm\n', params.lengthYAxis);
fprintf('  布局类型:         %s\n\n', params.layoutType);

validate_nmr_parameters(params, folder_path, config);

%% ========== 获取DXF文件 ==========
[pore_files, solid_files] = get_dxf_files(folder_path);

if isempty(pore_files)
    error('未找到DXF文件');
end

fprintf('找到 %d 对 DXF 文件:\n', length(pore_files));
for i = 1:min(5, length(pore_files))
    fprintf('  [%d] %s <-> %s\n', i, pore_files(i).name, solid_files(i).name);
end
if length(pore_files) > 5
    fprintf('  ... 还有 %d 对文件\n', length(pore_files) - 5);
end
fprintf('\n');

%% ========== 计算采样策略 ==========
total_dxf_pairs = length(pore_files);

if isinf(max_samples) || total_dxf_pairs <= max_samples
    % 数量不超过最大值，全部处理
    sample_indices = 1:total_dxf_pairs;
else
    % 使用 linspace 实现全局均匀分布采样
    % 从第1个到最后一个均匀选取 max_samples 个点
    sample_indices = unique(round(linspace(1, total_dxf_pairs, max_samples)));
end

fprintf('采样策略:\n');
fprintf('  总DXF对数: %d\n', total_dxf_pairs);
fprintf('  最大采样数: %d\n', max_samples);
fprintf('  实际处理数: %d\n', length(sample_indices));
fprintf('  采样索引: %s\n\n', mat2str(sample_indices));

%% ========== 创建输出目录 ==========
comsol_output_dir = fullfile(folder_path, 'comsol_results');
inversion_output_dir = fullfile(folder_path, 'inversion_results');

if ~exist(comsol_output_dir, 'dir')
    mkdir(comsol_output_dir);
    fprintf('创建COMSOL结果目录: %s\n', comsol_output_dir);
end

if ~exist(inversion_output_dir, 'dir')
    mkdir(inversion_output_dir);
    fprintf('创建反演结果目录: %s\n', inversion_output_dir);
end
fprintf('\n');

%% ========== 读取 global_evolution.xlsx 获取孔隙率数据用于校准 ==========
global_evolution_file = fullfile(folder_path, 'global_evolution.xlsx');
global_porosity_data = [];
global_time_data = [];
first_timestep_porosity = NaN;
calibration_factor = [];  % 初始为空，将在第一次反演后计算

if exist(global_evolution_file, 'file')
    try
        global_data = readtable(global_evolution_file);
        col_names = global_data.Properties.VariableNames;
        
        % 查找时间列
        time_col_idx = find(contains(lower(col_names), 'time'), 1);
        if isempty(time_col_idx)
            time_col_idx = find(contains(lower(col_names), 'step'), 1);
        end
        if isempty(time_col_idx)
            time_col_idx = 1;
        end
        
        % 查找孔隙率列
        porosity_col_idx = find(contains(lower(col_names), 'porosity'), 1);
        if isempty(porosity_col_idx)
            porosity_col_idx = 3;
        end
        
        global_time_data = global_data{:, time_col_idx};
        global_porosity_data = global_data{:, porosity_col_idx};
        first_timestep_porosity = global_porosity_data(1);
        
        fprintf('✓ 读取 global_evolution.xlsx 成功\n');
        fprintf('  第一时间步孔隙率: %.6f\n', first_timestep_porosity);
        fprintf('  总时间步数: %d\n\n', length(global_porosity_data));
    catch ME
        fprintf('⚠ 读取 global_evolution.xlsx 失败: %s\n\n', ME.message);
    end
else
    fprintf('⚠ 未找到 global_evolution.xlsx，将使用默认校准因子\n\n');
end

%% ========== 选择要处理的时间步 ==========
fprintf('选择处理模式:\n');
fprintf('  [1] 使用采样策略 (处理 %d 个时间步)\n', length(sample_indices));
fprintf('  [2] 处理所有时间步 (%d 个)\n', total_dxf_pairs);
fprintf('  [3] 仅处理第一个时间步 (测试)\n');
fprintf('  [4] 指定时间步范围\n');
fprintf('  [5] 仅可视化已有反演结果 (跳过COMSOL计算)\n');
fprintf('\n');

mode = input('请选择 (1-5, 默认1): ');
if isempty(mode)
    mode = 1;
end

% 标记是否跳过计算，仅可视化
skip_computation = false;

switch mode
    case 1
        % 使用采样策略
        process_indices = sample_indices;
    case 2
        % 处理所有
        process_indices = 1:total_dxf_pairs;
    case 3
        % 仅第一个
        process_indices = 1;
    case 4
        start_idx = input('起始时间步 (1): ');
        if isempty(start_idx), start_idx = 1; end
        end_idx = input(sprintf('结束时间步 (%d): ', total_dxf_pairs));
        if isempty(end_idx), end_idx = total_dxf_pairs; end
        process_indices = start_idx:end_idx;
    case 5
        % 仅可视化已有反演结果
        skip_computation = true;
        process_indices = [];  % 不需要处理索引
        fprintf('\n[模式5] 跳过COMSOL计算，直接读取已有反演结果...\n');
    otherwise
        process_indices = sample_indices;
end

%% ========== 模式5: 读取已有反演结果 ==========
if skip_computation
    fprintf('\n========== 读取已有反演结果 ==========\n');
    
    % 检查反演结果目录
    if ~exist(inversion_output_dir, 'dir')
        error('反演结果目录不存在: %s', inversion_output_dir);
    end
    
    % 查找所有 MAT 文件
    mat_files = dir(fullfile(inversion_output_dir, '*_T2.mat'));
    
    if isempty(mat_files)
        error('未找到反演结果文件 (*_T2.mat)');
    end
    
    fprintf('找到 %d 个反演结果文件\n', length(mat_files));
    
    % 初始化反演结果存储
    inversion_results = struct();
    inversion_results.timesteps = [];
    inversion_results.water_contents = [];
    inversion_results.timestep_strings = {};
    
    % 遍历读取每个 MAT 文件
    for k = 1:length(mat_files)
        mat_path = fullfile(inversion_output_dir, mat_files(k).name);
        
        try
            % 从文件名提取时间步
            % 文件名格式: T2_Da0.0185_Pe0.0500_X0.0600_Y0.0390_t0001_T2.mat
            tokens = regexp(mat_files(k).name, '_t(\d{4})_', 'tokens');
            if ~isempty(tokens)
                timestep_str = tokens{1}{1};
                timestep_num = str2double(timestep_str);
            else
                fprintf('  跳过 (无法解析时间步): %s\n', mat_files(k).name);
                continue;
            end
            
            % 加载 MAT 文件
            data = load(mat_path);
            
            if isfield(data, 'total_water')
                inversion_results.timesteps(end+1) = timestep_num;
                inversion_results.water_contents(end+1) = data.total_water;
                inversion_results.timestep_strings{end+1} = timestep_str;
                fprintf('  ✓ 读取 t%s: 含水率 = %.6f\n', timestep_str, data.total_water);
            else
                fprintf('  ⚠ 跳过 (无 total_water 字段): %s\n', mat_files(k).name);
            end
        catch ME
            fprintf('  × 读取失败: %s (%s)\n', mat_files(k).name, ME.message);
        end
    end
    
    % 按时间步排序
    if ~isempty(inversion_results.timesteps)
        [inversion_results.timesteps, sort_idx] = sort(inversion_results.timesteps);
        inversion_results.water_contents = inversion_results.water_contents(sort_idx);
        inversion_results.timestep_strings = inversion_results.timestep_strings(sort_idx);
    end
    
    fprintf('\n成功读取 %d 个反演结果\n', length(inversion_results.timesteps));
    
    % 更新 process_indices 用于后续可视化
    process_indices = inversion_results.timesteps;
    
    % 跳过处理循环，直接进入可视化
    success_count = length(inversion_results.timesteps);
    fail_count = 0;
    
else
    % 正常处理模式
    fprintf('\n将处理 %d 个时间步\n', length(process_indices));
    if length(process_indices) <= 20
        fprintf('索引: %s\n', mat2str(process_indices));
    end
    fprintf('按Enter继续...\n');
    pause;

    %% ========== 初始化反演结果存储 ==========
    inversion_results = struct();
    inversion_results.timesteps = [];
    inversion_results.water_contents = [];
    inversion_results.timestep_strings = {};

    %% ========== 处理循环 ==========
    success_count = 0;
    fail_count = 0;
    total_to_process = length(process_indices);

    for idx = 1:total_to_process
        j = process_indices(idx);  % 实际的文件索引
        pore_dxf = fullfile(folder_path, 'dxf_pore', pore_files(j).name);
        solid_dxf = fullfile(folder_path, 'dxf_solid', solid_files(j).name);
        
        timestep = extract_timestep(pore_files(j).name);
        
        fprintf('\n════════════════════════════════════════════════════════════\n');
        fprintf('处理时间步 %s [%d/%d] (原索引 %d/%d)\n', timestep, idx, total_to_process, j, total_dxf_pairs);
        fprintf('════════════════════════════════════════════════════════════\n');
        
        % 生成输出文件名
        excel_filename = sprintf('T2_Da%.4f_Pe%.4f_X%.4f_Y%.4f_t%s.xlsx', ...
            params.Da, params.Pe, params.lengthXAxis, params.lengthYAxis, timestep);
        excel_output = fullfile(comsol_output_dir, excel_filename);
        
        fprintf('孔隙DXF: %s\n', pore_files(j).name);
        fprintf('固体DXF: %s\n', solid_files(j).name);
        fprintf('输出文件: %s\n\n', excel_filename);
        
        % COMSOL处理
        fprintf('[COMSOL] 处理几何并求解...\n');
        
        comsol_success = run_comsol_processing(...
            config.mph_file, ...
            pore_dxf, ...
            solid_dxf, ...
            params.lengthXAxis, ...
            params.lengthYAxis, ...
            excel_output, ...
            config);
        
        if ~comsol_success
            fprintf('× COMSOL处理失败\n');
            fail_count = fail_count + 1;
            continue;
        end
        
        fprintf('✓ COMSOL处理完成\n\n');
        
        % 检查Excel文件
        if ~exist(excel_output, 'file')
            fprintf('× Excel文件未生成\n');
            fail_count = fail_count + 1;
            continue;
        end
        
        % MATLAB T2反演
        fprintf('[MATLAB] 运行T2反演...\n');
        
        [inversion_success, inv_water_content, raw_spectrum_sum, calibration_factor] = run_python_inversion(...
            excel_output, ...
            inversion_output_dir, ...
            config, ...
            calibration_factor);
        
        % 如果是第一个时间步且校准因子未设置，计算校准因子
        if idx == 1 && ~isnan(first_timestep_porosity) && ~isnan(raw_spectrum_sum) && raw_spectrum_sum > 0
            calibration_factor = first_timestep_porosity / raw_spectrum_sum;
            fprintf('[校准] 基于第一时间步孔隙率 %.6f 计算校准因子: %.6e\n', ...
                first_timestep_porosity, calibration_factor);
            
            % 重新运行第一个时间步的反演
            fprintf('[重新反演] 使用校准因子重新计算...\n');
            [inversion_success, inv_water_content, ~, ~] = run_python_inversion(...
                excel_output, ...
                inversion_output_dir, ...
                config, ...
                calibration_factor);
        end
        
        if inversion_success
            fprintf('✓ T2反演完成\n');
            success_count = success_count + 1;
            
            % 存储反演结果
            inversion_results.timesteps(end+1) = j;
            inversion_results.water_contents(end+1) = inv_water_content;
            inversion_results.timestep_strings{end+1} = timestep;
        else
            fprintf('× T2反演失败\n');
            fail_count = fail_count + 1;
        end
    end
end  % end if skip_computation

%% ========== 生成对比图 ==========
if ~isempty(inversion_results.timesteps) && ~isempty(global_porosity_data)
    fprintf('\n[对比图] 生成反演含水率与原始孔隙率对比图...\n');
    try
        generate_comparison_plot(...
            inversion_results, ...
            global_time_data, ...
            global_porosity_data, ...
            process_indices, ...
            inversion_output_dir, ...
            folder_name, ...
            config.show_nmr_porosity);
        fprintf('✓ 对比图生成完成\n');
    catch ME
        fprintf('× 对比图生成失败: %s\n', ME.message);
    end
end

%% ========== 生成GIF动画 ==========
if config.enable_gif
    fprintf('\n[GIF] 生成动画...\n');
    try
        generate_gif(folder_path, config, process_indices);
    catch ME
        fprintf('× GIF生成失败: %s\n', ME.message);
    end
end

%% ========== 总结 ==========
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf('处理完成总结:\n');
fprintf('  成功: %d\n', success_count);
fprintf('  失败: %d\n', fail_count);
fprintf('  COMSOL结果: %s\n', comsol_output_dir);
fprintf('  反演结果:   %s\n', inversion_output_dir);
fprintf('════════════════════════════════════════════════════════════\n');

function validate_nmr_parameters(params, folder_path, config)
    if config.require_metadata_json && ~strcmp(params.parameterSource, 'run_metadata.json')
        error('当前配置要求从JSON读取参数，但未能读取 %s: %s', ...
            config.metadata_filename, folder_path);
    end

    missing = {};
    if isnan(params.Da); missing{end+1} = 'Da'; end
    if isnan(params.Pe); missing{end+1} = 'Pe'; end
    if isnan(params.L); missing{end+1} = 'L/characteristicLength_cm'; end
    if isnan(params.lengthXAxis) || params.lengthXAxis <= 0
        missing{end+1} = 'lengthXAxis_cm';
    end
    if isnan(params.lengthYAxis) || params.lengthYAxis <= 0
        missing{end+1} = 'lengthYAxis_cm';
    end

    if ~isempty(missing)
        error('NMR模拟缺少关键几何/输运参数 (%s)。请检查 %s 或旧版文件夹名。', ...
            strjoin(missing, ', '), fullfile(folder_path, config.metadata_filename));
    end

    if ~strcmp(params.parameterSource, 'run_metadata.json')
        warning('当前样本使用旧版文件夹名解析参数。建议生成 %s 后再运行NMR模拟: %s', ...
            config.metadata_filename, folder_path);
    end
end

function text = describe_parameter_source(params)
    if strcmp(params.parameterSource, 'run_metadata.json')
        text = sprintf('JSON (%s)', params.metadataFile);
    elseif strcmp(params.parameterSource, 'folder_name')
        text = '旧版文件夹名';
    else
        text = params.parameterSource;
    end
end
