classdef AutomationConfig
    % AutomationConfig - 自动化配置类
    % 集中管理所有路径和参数配置
    
    properties
        % 根目录
        project_root
        
        % 数据目录
        data_root
        metadata_filename
        require_metadata_json
        
        % COMSOL配置
        comsol_path
        mph_file
        
        % Python配置
        python_exe
        inversion_script
        
        % 处理选项
        enable_comsol      % 是否启用COMSOL处理
        enable_inversion   % 是否启用Python反演
        overwrite_existing % 是否覆盖已存在的结果
        
        % 采样配置
        max_samples_per_folder  % 每个文件夹最多处理的DXF对数量
        
        % 几何缩放配置
        scale_factor  % 几何缩放系数 (用于缩放导入的DXF几何)
        
        % 导出配置
        export_mph  % 是否导出mph文件
        
        % GIF动画配置
        enable_gif      % 是否启用GIF生成
        gif_speed       % GIF播放速度 (1.0为正常速度)
        gif_format      % 输出格式: 'gif' 或 'mp4'
        
        % 对比图配置
        show_nmr_porosity  % 是否在对比图中显示NMR反演孔隙率结果
        
        % 日志配置
        enable_logging
        log_level  % 'debug', 'info', 'warning', 'error'
    end
    
    methods
        function obj = AutomationConfig()
            % 构造函数 - 初始化默认配置
            
            % 项目根目录：根据本文件位置自动定位到 ReactiveTransport
            automation_dir = fileparts(mfilename('fullpath'));
            obj.project_root = fileparts(automation_dir);
            
            % 数据根目录 (包含所有RTM结果文件夹)
            obj.data_root = fullfile(fileparts(obj.project_root), 'outputs', 'rtm_tests');
            obj.metadata_filename = 'run_metadata.json';
            obj.require_metadata_json = true;  % true时禁止旧版文件夹名参数解析
            % COMSOL配置
            obj.comsol_path = 'C:\Program Files\COMSOL\COMSOL63\Multiphysics';
            obj.mph_file = fullfile(obj.project_root, 'NMR', 'CT-simulation.mph');
            
            % Python配置
            obj.python_exe = 'C:\ProgramData\anaconda3\python.exe';
            % 使用T2_process标准工具包的反演桥接脚本
            obj.inversion_script = fullfile(obj.project_root, 'automation', 'run_t2_process_inversion.py');
            
            % 处理选项
            obj.enable_comsol = true;
            obj.enable_inversion = true;
            obj.overwrite_existing = false;
            
            % 采样配置
            obj.max_samples_per_folder = 100;  % 每个文件夹最多处理100个DXF对
            
            % 几何缩放配置
            obj.scale_factor = 10000;  % 几何缩放系数，默认值为10000
            
            % 导出配置
            obj.export_mph = false;  % 是否导出mph文件，默认不导出
            
            % GIF动画配置
            obj.enable_gif = true;    % 是否启用GIF生成
            obj.gif_speed = 0.5;      % GIF播放速度 (0.5表示较慢)
            obj.gif_format = 'gif';   % 输出格式: 'gif' 或 'mp4'
            
            % 对比图配置
            obj.show_nmr_porosity = true;  % 是否在对比图中显示NMR反演孔隙率结果
            
            % 日志配置
            obj.enable_logging = true;
            obj.log_level = 'info';
        end
        
        function validate(obj)
            % 验证配置有效性
            
            fprintf('验证配置...\n');
            
            % 检查数据目录
            if ~exist(obj.data_root, 'dir')
                error('数据目录不存在: %s', obj.data_root);
            end
            fprintf('  ✓ 数据目录存在\n');
            
            % 检查COMSOL
            if obj.enable_comsol
                if ~exist(obj.comsol_path, 'dir')
                    error('COMSOL安装目录不存在: %s', obj.comsol_path);
                end
                if ~exist(obj.mph_file, 'file')
                    error('COMSOL模型文件不存在: %s', obj.mph_file);
                end
                fprintf('  ✓ COMSOL配置有效\n');
            end
            
            % 检查Python
            if obj.enable_inversion
                if ~exist(obj.python_exe, 'file')
                    error('Python解释器不存在: %s', obj.python_exe);
                end
                if ~exist(obj.inversion_script, 'file')
                    error('反演脚本不存在: %s', obj.inversion_script);
                end
                fprintf('  ✓ Python配置有效\n');
            end
            
            fprintf('配置验证完成!\n\n');
        end
        
        function disp_config(obj)
            % 显示当前配置
            
            fprintf('当前配置:\n');
            fprintf('─────────────────────────────────────────\n');
            fprintf('  项目根目录: %s\n', obj.project_root);
            fprintf('  数据目录:   %s\n', obj.data_root);
            fprintf('  元数据文件: %s\n', obj.metadata_filename);
            fprintf('  强制JSON:   %s\n', mat2str(obj.require_metadata_json));
            fprintf('  COMSOL路径: %s\n', obj.comsol_path);
            fprintf('  模型文件:   %s\n', obj.mph_file);
            fprintf('  Python:     %s\n', obj.python_exe);
            fprintf('  反演脚本:   %s\n', obj.inversion_script);
            fprintf('─────────────────────────────────────────\n');
            fprintf('  启用COMSOL: %s\n', mat2str(obj.enable_comsol));
            fprintf('  启用反演:   %s\n', mat2str(obj.enable_inversion));
            fprintf('  覆盖结果:   %s\n', mat2str(obj.overwrite_existing));
            fprintf('  最大采样数: %d\n', obj.max_samples_per_folder);
            fprintf('  缩放系数:   %d\n', obj.scale_factor);
            fprintf('  导出mph:   %s\n', mat2str(obj.export_mph));
            fprintf('  启用GIF:    %s\n', mat2str(obj.enable_gif));
            fprintf('  GIF速度:    %.2f\n', obj.gif_speed);
            fprintf('  GIF格式:    %s\n', obj.gif_format);
            fprintf('  显示NMR:    %s\n', mat2str(obj.show_nmr_porosity));
            fprintf('─────────────────────────────────────────\n\n');
        end
    end
end
