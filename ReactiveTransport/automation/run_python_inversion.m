function [success, total_water, raw_spectrum_sum, calibration_factor] = run_python_inversion(excel_file, output_dir, config, calibration_factor)
    % run_python_inversion - 使用T2_process标准工具包执行T2反演
    %
    % 输入:
    %   excel_file         - Excel输入文件路径
    %   output_dir         - 输出目录
    %   config             - 配置对象
    %   calibration_factor - (可选) 校准因子，如果未提供则由反演脚本使用默认值
    %
    % 输出:
    %   success            - 处理是否成功
    %   total_water        - 反演得到的总含水率
    %   raw_spectrum_sum   - 原始谱的积分值（用于校准）
    %   calibration_factor - 使用的校准因子

    success = false;
    total_water = NaN;
    raw_spectrum_sum = NaN;

    if nargin < 4 || isempty(calibration_factor)
        calibration_factor = [];
    end

    if ~config.enable_inversion
        fprintf('      ⚠ T2反演已禁用\n');
        success = true;
        return;
    end

    if ~exist(excel_file, 'file')
        fprintf('      × Excel文件不存在: %s\n', excel_file);
        return;
    end

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    [~, excel_filename, ~] = fileparts(excel_file);
    excel_filename = char(excel_filename);
    output_dir_char = char(output_dir);
    output_png = fullfile(output_dir_char, [excel_filename, '_T2.png']);
    output_mat = fullfile(output_dir_char, [excel_filename, '_T2.mat']);

    png_exists = ~isempty(dir(output_png)) || check_file_exists_long(output_png);
    mat_exists = ~isempty(dir(output_mat)) || check_file_exists_long(output_mat);

    if png_exists && mat_exists && ~config.overwrite_existing
        fprintf('      ⚠ 反演结果已存在,跳过T2反演\n');
        try
            loaded = load_cached_inversion(output_mat, excel_filename);
            total_water = loaded.total_water;
            raw_spectrum_sum = loaded.raw_spectrum_sum;
            if isfield(loaded, 'calibration_factor')
                calibration_factor = loaded.calibration_factor;
            end
            fprintf('      ✓ 从缓存加载: total_water=%.4f\n', total_water);
        catch ME
            fprintf('      ⚠ 无法加载缓存数据: %s\n', ME.message);
        end
        success = true;
        return;
    end

    % 数值模拟场景统一使用固定平滑/正则化因子10000，不再使用L-curve搜索。
    fprintf('      [T2_process] 执行固定正则化T2反演 (regularization=10000)...\n');
    [success, total_water, raw_spectrum_sum, calibration_factor] = run_t2_process_package_inversion( ...
        excel_file, output_dir, config, calibration_factor);
end

function loaded = load_cached_inversion(output_mat, excel_filename)
    % load_cached_inversion - 读取旧/新反演缓存，兼容Windows长路径

    if ispc && length(output_mat) > 200
        temp_mat = fullfile(tempdir, [excel_filename, '_T2.mat']);
        try
            unc_mat = ['\\?\' output_mat];
            copyfile(unc_mat, temp_mat);
            loaded = load(temp_mat, 'total_water', 'raw_spectrum_sum', 'calibration_factor');
            delete(temp_mat);
        catch
            loaded = load(output_mat, 'total_water', 'raw_spectrum_sum', 'calibration_factor');
        end
    else
        loaded = load(output_mat, 'total_water', 'raw_spectrum_sum', 'calibration_factor');
    end
end

function [success, total_water, raw_spectrum_sum, calibration_factor] = run_t2_process_package_inversion(excel_file, output_dir, config, calibration_factor)
    % run_t2_process_package_inversion - 调用ReactiveTransport/T2_process工具包反演
    %
    % 该桥接函数保留自动化框架既有返回值和*_T2.mat缓存字段，同时将
    % 数值反演和可视化交给T2_process/nmr_t2包完成。

    success = false;
    total_water = NaN;
    raw_spectrum_sum = NaN;

    automation_dir = fileparts(mfilename('fullpath'));
    default_script = fullfile(automation_dir, 'run_t2_process_inversion.py');
    script_file = default_script;

    try
        if isprop(config, 'inversion_script') && exist(config.inversion_script, 'file')
            script_file = config.inversion_script;
        end
    catch
        script_file = default_script;
    end

    if ~exist(script_file, 'file')
        fprintf('        x T2_process反演脚本不存在: %s\n', script_file);
        return;
    end

    python_exe = 'python';
    try
        if isprop(config, 'python_exe') && ~isempty(config.python_exe)
            python_exe = config.python_exe;
        end
    catch
        python_exe = 'python';
    end

    cmd = sprintf('%s %s --input-excel %s --output-dir %s --regularization %.17g --time-to-ms-scale %.17g', ...
        quote_cmd_arg(python_exe), ...
        quote_cmd_arg(script_file), ...
        quote_cmd_arg(excel_file), ...
        quote_cmd_arg(output_dir), ...
        10000, ...
        1000.0);

    if ~isempty(calibration_factor) && isnumeric(calibration_factor) && isfinite(calibration_factor)
        cmd = sprintf('%s --calibration-factor %.17g', cmd, calibration_factor);
    end

    [status, cmdout] = system(cmd);
    if ~isempty(strtrim(cmdout))
        fprintf('%s\n', strtrim(cmdout));
    end

    result = parse_result_json(cmdout);
    if isempty(result)
        fprintf('        x 无法解析T2_process反演返回结果\n');
        return;
    end

    if isfield(result, 'total_water') && ~isempty(result.total_water)
        total_water = result.total_water;
    end
    if isfield(result, 'raw_spectrum_sum') && ~isempty(result.raw_spectrum_sum)
        raw_spectrum_sum = result.raw_spectrum_sum;
    end
    if isfield(result, 'calibration_factor') && ~isempty(result.calibration_factor)
        calibration_factor = result.calibration_factor;
    end

    if status == 0 && isfield(result, 'success') && logical(result.success)
        success = true;
    else
        if isfield(result, 'error')
            fprintf('        x T2_process反演失败: %s\n', result.error);
        else
            fprintf('        x T2_process反演失败，退出码: %d\n', status);
        end
    end
end

function quoted = quote_cmd_arg(value)
    % quote_cmd_arg - Windows/POSIX system命令参数加引号
    value = char(value);
    value = strrep(value, '"', '\"');
    quoted = ['"' value '"'];
end

function result = parse_result_json(cmdout)
    % parse_result_json - 从Python输出中提取RESULT_JSON行
    result = [];
    marker = 'RESULT_JSON=';
    lines = regexp(cmdout, '\r?\n', 'split');
    json_text = '';

    for i = 1:length(lines)
        line = strtrim(lines{i});
        if startsWith(line, marker)
            json_text = extractAfter(line, strlength(marker));
        end
    end

    if strlength(json_text) == 0
        return;
    end

    try
        result = jsondecode(char(json_text));
    catch ME
        fprintf('        x JSON解析失败: %s\n', ME.message);
        result = [];
    end
end

function exists = check_file_exists_long(filepath)
    % check_file_exists_long - 检查文件是否存在，兼容Windows长路径

    exists = false;
    filepath = char(filepath);

    if ~ispc
        exists = exist(filepath, 'file') == 2;
        return;
    end

    [status, ~] = system(sprintf('if exist "%s" (echo 1) else (echo 0)', filepath));
    if status == 0
        [status, result] = system(sprintf('dir /b "%s" 2>nul', filepath));
        exists = (status == 0 && ~isempty(strtrim(result)));
    end
end
