function report = check_rtsphem_environment()
%% check_rtsphem_environment - 检查本机能否运行 RTM + PHREEQC 耦合模拟
%
% 用法（在 MATLAB 中，任意当前目录）：
%   report = check_rtsphem_environment();
%
% 检查内容：
%   1) MATLAB 版本（参考环境 R2025b）
%   2) 工具箱（Parallel Computing Toolbox 为可选）
%   3) IPhreeqcCOM COM 组件能否创建（actxserver('IPhreeqcCOM.Object')）
%   4) 能否加载仓库自带 PHREEQC 数据库并完成一次方解石平衡计算
%   5) HyPHM / RTM 路径能否找到 PNM_beauty3.m
% 输出：
%   report  struct 数组，字段 name / ok / detail；同时在命令行打印 PASS/FAIL。
% 依赖：
%   仅 MATLAB 本体；不修改任何文件，不安装任何东西。

rtmDir = fileparts(mfilename('fullpath'));
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));

report = struct('name', {}, 'ok', {}, 'detail', {});

%% 1) MATLAB 版本
rel = version('-release');
report(end+1) = item('MATLAB 版本', true, ...
    sprintf('R%s（参考环境 R2025b%s）', rel, ternary(strcmp(rel, '2025b'), '', '，版本不同可能导致数值有微小差异')));

%% 2) 工具箱
v = ver;
names = {v.Name};
optional = {'Parallel Computing Toolbox'};
for i = 1:numel(optional)
    has = any(strcmp(names, optional{i}));
    report(end+1) = item(['工具箱(可选): ' optional{i}], true, ...
        ternary(has, '已安装', '未安装，parfor 会退化为串行，结果不受影响')); %#ok<AGROW>
end

%% 3) Windows + IPhreeqcCOM
if ~ispc
    report(end+1) = item('操作系统', false, 'IPhreeqcCOM 只能在 Windows 上使用');
    printReport(report);
    return;
end
engine = [];
try
    engine = actxserver('IPhreeqcCOM.Object');
    report(end+1) = item('IPhreeqcCOM COM 组件', true, ...
        sprintf('可创建（MATLAB 位数: %s）', computer('arch')));
catch ME
    report(end+1) = item('IPhreeqcCOM COM 组件', false, ...
        ['actxserver 失败: ' ME.message ...
        '。请安装 IPhreeqcCOM-3.8.6-17100-x64.msi（见仓库 README.md）']);
end

%% 4) 数据库 + 方解石平衡试算
if ~isempty(engine)
    dbPath = fullfile(rtmDir, 'phreeqc', 'database', 'backup_legacy_20260626', 'phreeqc.dat');
    try
        if exist(dbPath, 'file') ~= 2
            error('找不到数据库文件: %s', dbPath);
        end
        if engine.LoadDatabase(dbPath) ~= 0
            error('LoadDatabase 失败: %s', engine.GetErrorString);
        end
        input = sprintf([ ...
            'SOLUTION 1\n  temp 25\n  pH 3 charge\n  Cl 1 mmol/kgw\n' ...
            'EQUILIBRIUM_PHASES 1\n  Calcite 0 10\n' ...
            'SELECTED_OUTPUT 1\n  -reset false\n  -pH true\n  -totals Ca\nEND\n']);
        if engine.RunString(input) ~= 0
            error('RunString 失败: %s', engine.GetErrorString);
        end
        out = engine.GetSelectedOutputArray;
        report(end+1) = item('PHREEQC 方解石平衡试算', true, ...
            sprintf('pH = %.3f, Ca = %.3e mol/kgw', out{end, 1}, out{end, 2}));
    catch ME
        report(end+1) = item('PHREEQC 方解石平衡试算', false, ME.message);
    end
    try
        delete(engine);
    catch
    end
end

%% 5) 代码路径
report(end+1) = item('RTM 主程序 PNM_beauty3.m', exist('PNM_beauty3', 'file') == 2, rtmDir);

printReport(report);
end

function s = item(name, ok, detail)
s = struct('name', name, 'ok', logical(ok), 'detail', detail);
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function printReport(report)
fprintf('\n===== RTSPHEM 环境检查 =====\n');
for i = 1:numel(report)
    fprintf('[%s] %s: %s\n', ternary(report(i).ok, 'PASS', 'FAIL'), report(i).name, report(i).detail);
end
if all([report.ok])
    fprintf('全部通过，可以运行 run_single_circle_phreeqc。\n');
else
    fprintf('存在 FAIL 项，请按仓库 README.md 的“环境配置”一节处理。\n');
end
fprintf('============================\n\n');
end
