function success = generate_gif(folder_path, config, valid_steps)
% GENERATE_GIF 调用Python脚本生成GIF动画
%
% 输入参数:
%   folder_path  - dissolution结果文件夹路径
%   config       - AutomationConfig配置对象
%   valid_steps  - (可选) 有效的时间步编号数组，用于同步timestep和inversion动画
%
% 输出参数:
%   success - 是否成功生成GIF
%
% 功能:
%   1. 从主文件夹的timestep图片生成 animation_timestep.gif
%   2. 从inversion_results子文件夹生成 animation_inversion.gif
%
% 作者: 自动化框架
% 日期: 2025-12-02

    success = false;
    
    % 检查文件夹是否存在
    if ~exist(folder_path, 'dir')
        fprintf('  × GIF生成失败: 文件夹不存在 %s\n', folder_path);
        return;
    end
    
    % 获取配置参数
    python_exe = config.python_exe;
    gif_speed = config.gif_speed;
    gif_format = config.gif_format;
    
    % 构建Python脚本路径
    automation_path = fileparts(mfilename('fullpath'));
    gif_script = fullfile(automation_path, 'generate_gif_helper.py');
    
    % 检查Python脚本是否存在
    if ~exist(gif_script, 'file')
        fprintf('  × GIF生成失败: Python脚本不存在 %s\n', gif_script);
        return;
    end
    
    % 构建命令参数
    % 处理长路径：使用短路径或UNC路径前缀
    folder_arg = folder_path;
    if ispc && length(folder_path) > 200
        % 尝试获取短路径名
        try
            [~, short_path] = system(['for %I in ("' folder_path '") do @echo %~sI']);
            short_path = strtrim(short_path);
            if ~isempty(short_path) && exist(short_path, 'dir')
                folder_arg = short_path;
                fprintf('  [长路径] 使用短路径: %s\n', short_path);
            else
                % 使用UNC前缀（Python脚本会处理）
                folder_arg = folder_path;
            end
        catch
            folder_arg = folder_path;
        end
    end
    
    if nargin >= 3 && ~isempty(valid_steps)
        % 将valid_steps转换为逗号分隔的字符串
        steps_str = strjoin(arrayfun(@num2str, valid_steps, 'UniformOutput', false), ',');
        cmd = sprintf('"%s" "%s" "%s" --speed %.2f --format %s --valid-steps "%s"', ...
            python_exe, gif_script, folder_arg, gif_speed, gif_format, steps_str);
    else
        cmd = sprintf('"%s" "%s" "%s" --speed %.2f --format %s', ...
            python_exe, gif_script, folder_arg, gif_speed, gif_format);
    end
    
    fprintf('  [GIF] 生成动画...\n');
    fprintf('  命令: %s\n', cmd);
    
    % 执行Python脚本
    [status, output] = system(cmd);
    
    % 显示Python输出
    if ~isempty(output)
        fprintf('  Python输出:\n');
        lines = strsplit(output, newline);
        for i = 1:length(lines)
            if ~isempty(strtrim(lines{i}))
                fprintf('    %s\n', lines{i});
            end
        end
    end
    
    if status == 0
        fprintf('  ✓ GIF动画生成完成\n');
        
        % 检查生成的文件
        timestep_gif = fullfile(folder_path, sprintf('animation_timestep.%s', gif_format));
        inversion_gif = fullfile(folder_path, sprintf('animation_inversion.%s', gif_format));
        
        if exist(timestep_gif, 'file')
            fprintf('    - %s\n', sprintf('animation_timestep.%s', gif_format));
        end
        if exist(inversion_gif, 'file')
            fprintf('    - %s\n', sprintf('animation_inversion.%s', gif_format));
        end
        
        success = true;
    else
        fprintf('  × GIF生成失败 (exit code: %d)\n', status);
    end
end
