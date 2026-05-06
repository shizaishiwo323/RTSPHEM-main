%% COMSOL初始化和连接设置
% 此脚本用于初始化COMSOL服务器并建立MATLAB与COMSOL的连接

function model = comsol_init()
    % COMSOL安装路径
    comsol_path = 'C:\Program Files\COMSOL\COMSOL63\Multiphysics';
    
    % 添加COMSOL MATLAB接口到路径
    mli_path = fullfile(comsol_path, 'mli');
    if ~exist(mli_path, 'dir')
        error('COMSOL MATLAB接口路径不存在: %s\n请检查COMSOL安装路径是否正确', mli_path);
    end
    addpath(mli_path);
    fprintf('✓ 已添加COMSOL MATLAB接口路径\n');
    
    % 导入COMSOL类
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    fprintf('正在连接到COMSOL服务器...\n');
    
    % 尝试多种连接方式
    connected = false;
    
    % 方法1: 连接到本地已运行的服务器(默认端口2036)
    try
        fprintf('  [方法1] 尝试连接到localhost:2036...\n');
        mphstart('localhost', 2036);
        model = ModelUtil.create('Model');
        fprintf('✓ 已连接到本地COMSOL服务器(端口2036)\n');
        connected = true;
    catch ME1
        fprintf('  ✗ 方法1失败: %s\n', ME1.message);
    end
    
    % 方法2: 尝试其他常用端口
    if ~connected
        ports = [2037, 2038, 3036];
        for port = ports
            try
                fprintf('  [方法2] 尝试连接到localhost:%d...\n', port);
                mphstart('localhost', port);
                model = ModelUtil.create('Model');
                fprintf('✓ 已连接到本地COMSOL服务器(端口%d)\n', port);
                connected = true;
                break;
            catch
                fprintf('  ✗ 端口%d连接失败\n', port);
            end
        end
    end
    
    % 方法3: 不指定服务器直接创建(使用内置模式)
    if ~connected
        try
            fprintf('  [方法3] 尝试直接创建模型(内置模式)...\n');
            model = ModelUtil.create('Model');
            fprintf('✓ 模型已创建(内置模式)\n');
            connected = true;
        catch ME3
            fprintf('  ✗ 方法3失败: %s\n', ME3.message);
        end
    end
    
    % 方法4: 启动新的服务器实例
    if ~connected
        try
            fprintf('  [方法4] 启动新的COMSOL服务器实例...\n');
            mphstart();
            pause(5); % 等待服务器启动
            model = ModelUtil.create('Model');
            fprintf('✓ 新服务器已启动并创建模型\n');
            connected = true;
        catch ME4
            fprintf('  ✗ 方法4失败: %s\n', ME4.message);
        end
    end
    
    % 如果所有方法都失败
    if ~connected
        fprintf('\n=== 故障排除步骤 ===\n');
        fprintf('所有自动连接方法均失败。请按以下步骤操作:\n\n');
        fprintf('步骤1: 确认COMSOL服务器是否在运行\n');
        fprintf('   - 查看任务管理器中是否有 comsolmphserver.exe 进程\n\n');
        
        fprintf('步骤2: 手动启动COMSOL服务器\n');
        fprintf('   方式A: 双击运行 start_comsol_server.m\n');
        fprintf('   方式B: 在CMD中运行:\n');
        fprintf('   "%s\\bin\\win64\\comsolmphserver.exe"\n\n', comsol_path);
        
        fprintf('步骤3: 检查端口占用\n');
        fprintf('   在CMD中运行: netstat -ano | findstr "2036"\n');
        fprintf('   如果端口被占用,杀死进程或使用其他端口\n\n');
        
        fprintf('步骤4: 尝试使用COMSOL Multiphysics桌面版\n');
        fprintf('   如果服务器模式有问题,可以使用桌面版:\n');
        fprintf('   - 打开COMSOL GUI\n');
        fprintf('   - 服务器 -> 连接到服务器 -> 本地主机\n\n');
        
        error('无法连接到COMSOL服务器。请按上述步骤排查。');
    end
    
    fprintf('\n=== COMSOL初始化成功 ===\n');
end
