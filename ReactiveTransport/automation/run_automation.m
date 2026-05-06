
%% 主自动化脚本 - NMR T2 反演批处理
% ============================================================================
% 功能: 自动处理ALLDissolutionResults中所有子文件夹的DXF文件
%       1. 调用COMSOL处理DXF几何并导出结果表
%       2. 调用Python脚本进行NMR T2反演生成结果图
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

%% ========== 配置区域 ==========
config = AutomationConfig();

% 显示配置信息
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║  NMR T2 反演批处理自动化系统                               ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

fprintf('配置信息:\n');
fprintf('  数据源目录: %s\n', config.data_root);
fprintf('  元数据文件: %s\n', config.metadata_filename);
fprintf('  强制JSON参数: %s\n', mat2str(config.require_metadata_json));
fprintf('  COMSOL模型: %s\n', config.mph_file);
fprintf('  Python解释器: %s\n', config.python_exe);
fprintf('  每文件夹最大采样数: %d\n', config.max_samples_per_folder);
fprintf('\n');

%% ========== 扫描所有子文件夹 ==========
fprintf('[1/4] 扫描数据目录...\n');

folders = scan_dissolution_folders(config.data_root);

if isempty(folders)
    fprintf('  × 未找到任何dissolution结果文件夹\n');
    return;
end

fprintf('  ✓ 找到 %d 个dissolution结果文件夹\n\n', length(folders));

%% ========== 显示所有文件夹信息 ==========
fprintf('[2/4] 解析文件夹参数...\n');
folder_params = cell(length(folders), 1);

for i = 1:length(folders)
    folder_name = folders(i).name;
    folder_path = fullfile(config.data_root, folder_name);
    params = parse_folder_name(folder_path, config.metadata_filename);
    [params.nmrReady, params.validationMessage] = validate_nmr_parameters(params, folder_path, config);
    folder_params{i} = params;
    
    
    fprintf('  [%d] %s\n', i, folder_name);
    fprintf('      Da=%.4f, Pe=%.4f, L=%.4f, X=%.4f, Y=%.4f, Type=%s\n', ...
        params.Da, params.Pe, params.L, params.lengthXAxis, params.lengthYAxis, params.layoutType);
    fprintf('      参数来源: %s\n', describe_parameter_source(params));
    if ~params.nmrReady
        fprintf('      ⚠ %s\n', params.validationMessage);
    elseif ~strcmp(params.parameterSource, 'run_metadata.json')
        fprintf('      ⚠ 使用旧版文件夹名解析参数，建议补齐 %s\n', config.metadata_filename);
    end
end
fprintf('\n');

%% ========== 用户确认 ==========
fprintf('是否开始批处理? (按Enter继续, Ctrl+C取消)\n');
pause;

%% ========== 批处理主循环 ==========
fprintf('[3/4] 开始批处理...\n\n');

total_folders = length(folders);
success_count = 0;
fail_count = 0;
skip_count = 0;

% 创建日志文件
log_file = fullfile(config.data_root, 'batch_logs', ...
    sprintf('batch_log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
if ~exist(fullfile(config.data_root, 'batch_logs'), 'dir')
    mkdir(fullfile(config.data_root, 'batch_logs'));
end
log_fid = fopen(log_file, 'w');
fprintf(log_fid, '批处理日志 - %s\n', datestr(now));
fprintf(log_fid, '================================================\n\n');

for i = 1:total_folders
    folder_path = fullfile(config.data_root, folders(i).name);
    params = folder_params{i};
    
    fprintf('════════════════════════════════════════════════════════════\n');
    fprintf('处理 [%d/%d]: %s\n', i, total_folders, folders(i).name);
    fprintf('════════════════════════════════════════════════════════════\n');
    
    try
        if ~params.nmrReady
            fprintf('  ⚠ 参数不完整，跳过: %s\n', params.validationMessage);
            skip_count = skip_count + 1;
            fprintf(log_fid, '[SKIP] %s - 参数不完整: %s\n', folders(i).name, params.validationMessage);
            continue;
        end

        % 获取DXF文件列表
        [pore_files, solid_files] = get_dxf_files(folder_path);
        
        if isempty(pore_files)
            fprintf('  ⚠ 未找到DXF文件,跳过\n');
            skip_count = skip_count + 1;
            fprintf(log_fid, '[SKIP] %s - 无DXF文件\n', folders(i).name);
            continue;
        end
        
        fprintf('  找到 %d 对 DXF 文件\n', length(pore_files));
        
        % ===== 动态采样: 最多处理 max_samples 个，全局均匀分布 =====
        max_samples = config.max_samples_per_folder;  % 从配置读取
        total_dxf_pairs = length(pore_files);
        
        if total_dxf_pairs <= max_samples
            % 数量不超过最大值，全部处理
            sample_indices = 1:total_dxf_pairs;
        else
            % 使用 linspace 实现全局均匀分布采样
            % 从第1个到最后一个均匀选取 max_samples 个点
            sample_indices = unique(round(linspace(1, total_dxf_pairs, max_samples)));
        end
        
        fprintf('  采样策略: 全局均匀分布, 将处理 %d 个 (索引: %s)\n', ...
            length(sample_indices), mat2str(sample_indices));
        
        % 创建结果输出目录
        comsol_output_dir = fullfile(folder_path, 'comsol_results');
        inversion_output_dir = fullfile(folder_path, 'inversion_results');
        
        if ~exist(comsol_output_dir, 'dir')
            mkdir(comsol_output_dir);
        end
        if ~exist(inversion_output_dir, 'dir')
            mkdir(inversion_output_dir);
        end
        
        % ===== 读取 global_evolution.xlsx 获取孔隙率数据用于校准 =====
        global_evolution_file = fullfile(folder_path, 'global_evolution.xlsx');
        global_porosity_data = [];
        global_time_data = [];
        first_timestep_porosity = NaN;
        calibration_factor = [];  % 初始为空，将在第一次反演后计算
        
        if exist(global_evolution_file, 'file')
            try
                global_data = readtable(global_evolution_file);
                % 假设列名包含 'Time' 和 'Porosity' (或类似名称)
                col_names = global_data.Properties.VariableNames;
                
                % 查找时间列
                time_col_idx = find(contains(lower(col_names), 'time'), 1);
                if isempty(time_col_idx)
                    time_col_idx = find(contains(lower(col_names), 'step'), 1);
                end
                if isempty(time_col_idx)
                    time_col_idx = 1;  % 默认第一列
                end
                
                % 查找孔隙率列
                porosity_col_idx = find(contains(lower(col_names), 'porosity'), 1);
                if isempty(porosity_col_idx)
                    porosity_col_idx = 3;  % 默认第三列 (根据截图)
                end
                
                global_time_data = global_data{:, time_col_idx};
                global_porosity_data = global_data{:, porosity_col_idx};
                first_timestep_porosity = global_porosity_data(1);
                
                fprintf('  ✓ 读取 global_evolution.xlsx 成功\n');
                fprintf('    第一时间步孔隙率: %.6f\n', first_timestep_porosity);
            catch ME
                fprintf('  ⚠ 读取 global_evolution.xlsx 失败: %s\n', ME.message);
            end
        else
            fprintf('  ⚠ 未找到 global_evolution.xlsx\n');
        end
        
        % 存储所有时间步的反演结果用于对比图
        inversion_results = struct();
        inversion_results.timesteps = [];
        inversion_results.water_contents = [];
        inversion_results.timestep_strings = {};
        
        % 处理采样的DXF文件
        for idx = 1:length(sample_indices)
            j = sample_indices(idx);  % 实际的文件索引
            pore_dxf = fullfile(folder_path, 'dxf_pore', pore_files(j).name);
            solid_dxf = fullfile(folder_path, 'dxf_solid', solid_files(j).name);
            
            % 提取时间步编号
            timestep = extract_timestep(pore_files(j).name);
            
            fprintf('\n  [%d/%d] 时间步 %s (原索引 %d/%d)\n', ...
                idx, length(sample_indices), timestep, j, total_dxf_pairs);
            fprintf('    孔隙: %s\n', pore_files(j).name);
            fprintf('    固体: %s\n', solid_files(j).name);
            
            % 生成输出文件名 (使用短文件名避免Windows路径长度限制)
            excel_filename_short = sprintf('T2_t%s.xlsx', timestep);
            excel_filename_long = sprintf('T2_Da%.4f_Pe%.4f_X%.4f_Y%.4f_t%s.xlsx', ...
                params.Da, params.Pe, params.lengthXAxis, params.lengthYAxis, timestep);
            
            % 检查是否已存在结果（兼容长短两种文件名）
            excel_output_short = fullfile(comsol_output_dir, excel_filename_short);
            excel_output_long = fullfile(comsol_output_dir, excel_filename_long);
            
            if exist(excel_output_long, 'file')
                % 优先使用已存在的长文件名
                excel_output = excel_output_long;
                excel_filename = excel_filename_long;
            else
                % 使用短文件名
                excel_output = excel_output_short;
                excel_filename = excel_filename_short;
            end
            
            % 步骤1: 调用COMSOL处理
            fprintf('    [COMSOL] 处理几何并求解...\n');
            
            comsol_success = run_comsol_processing(...
                config.mph_file, ...
                pore_dxf, ...
                solid_dxf, ...
                params.lengthXAxis, ...
                params.lengthYAxis, ...
                excel_output, ...
                config);
            
            if ~comsol_success
                fprintf('    × COMSOL处理失败\n');
                fprintf(log_fid, '[FAIL] %s/t%s - COMSOL处理失败\n', folders(i).name, timestep);
                continue;
            end
            
            fprintf('    ✓ COMSOL处理完成: %s\n', excel_filename);
            
            % 检查Excel文件是否生成
            if ~exist(excel_output, 'file')
                fprintf('    × Excel文件未生成\n');
                fprintf(log_fid, '[FAIL] %s/t%s - Excel未生成\n', folders(i).name, timestep);
                continue;
            end
            
            % 步骤2: 调用Python进行T2反演
            fprintf('    [Python] 运行T2反演...\n');
            
            [python_success, inv_water_content, raw_spectrum_sum, calibration_factor] = run_python_inversion(...
                excel_output, ...
                inversion_output_dir, ...
                config, ...
                calibration_factor);  % 传递校准因子
            
            % 如果是第一个时间步且校准因子未设置，计算校准因子
            if idx == 1 && ~isnan(first_timestep_porosity) && ~isnan(raw_spectrum_sum) && raw_spectrum_sum > 0
                % calibration_factor = porosity / raw_spectrum_sum
                calibration_factor = first_timestep_porosity / raw_spectrum_sum;
                fprintf('    [校准] 基于第一时间步孔隙率 %.6f 计算校准因子: %.6e\n', ...
                    first_timestep_porosity, calibration_factor);
                
                % 重新运行第一个时间步的反演，使用正确的校准因子
                fprintf('    [重新反演] 使用校准因子重新计算第一时间步...\n');
                [python_success, inv_water_content, ~, ~] = run_python_inversion(...
                    excel_output, ...
                    inversion_output_dir, ...
                    config, ...
                    calibration_factor);
            end
            
            if python_success
                fprintf('    ✓ T2反演完成\n');
                fprintf(log_fid, '[OK] %s/t%s\n', folders(i).name, timestep);
                
                % 存储反演结果
                inversion_results.timesteps(end+1) = j;
                inversion_results.water_contents(end+1) = inv_water_content;
                inversion_results.timestep_strings{end+1} = timestep;
            else
                fprintf('    × T2反演失败\n');
                fprintf(log_fid, '[FAIL] %s/t%s - Python反演失败\n', folders(i).name, timestep);
            end
        end
        
        % ===== 生成反演含水率与原始孔隙率对比图 =====
        if ~isempty(inversion_results.timesteps) && ~isempty(global_porosity_data)
            fprintf('\n  [对比图] 生成反演含水率与原始孔隙率对比图...\n');
            try
                generate_comparison_plot(...
                    inversion_results, ...
                    global_time_data, ...
                    global_porosity_data, ...
                    sample_indices, ...
                    inversion_output_dir, ...
                    folders(i).name, ...
                    config.show_nmr_porosity);
                fprintf('  ✓ 对比图生成完成\n');
            catch ME
                fprintf('  × 对比图生成失败: %s\n', ME.message);
            end
        end
        
        % ===== 生成GIF动画 =====
        if config.enable_gif
            try
                generate_gif(folder_path, config, sample_indices);
            catch ME
                fprintf('  × GIF生成失败: %s\n', ME.message);
            end
        end
        
        success_count = success_count + 1;
        
    catch ME
        fprintf('  × 处理失败: %s\n', ME.message);
        fail_count = fail_count + 1;
        fprintf(log_fid, '[ERROR] %s - %s\n', folders(i).name, ME.message);
    end
    
    fprintf('\n');
end

fclose(log_fid);

%% ========== 总结报告 ==========
fprintf('[4/4] 批处理完成\n');
fprintf('════════════════════════════════════════════════════════════\n');
fprintf('批处理总结:\n');
fprintf('  总文件夹数: %d\n', total_folders);
fprintf('  成功: %d\n', success_count);
fprintf('  失败: %d\n', fail_count);
fprintf('  跳过: %d\n', skip_count);
fprintf('  日志文件: %s\n', log_file);
fprintf('════════════════════════════════════════════════════════════\n');

% 生成总结报告
generate_summary_report(config.data_root, folders, folder_params);

fprintf('\n✓ 自动化批处理完成!\n');

function [is_valid, message] = validate_nmr_parameters(params, folder_path, config)
    is_valid = true;
    message = '';

    if config.require_metadata_json && ~strcmp(params.parameterSource, 'run_metadata.json')
        is_valid = false;
        message = sprintf('当前配置要求从JSON读取参数，但未能读取 %s: %s', ...
            config.metadata_filename, folder_path);
        return;
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
        is_valid = false;
        message = sprintf('NMR模拟缺少关键几何/输运参数 (%s)，请检查 %s。', ...
            strjoin(missing, ', '), fullfile(folder_path, config.metadata_filename));
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
