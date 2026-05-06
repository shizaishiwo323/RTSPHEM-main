%% 批次实验脚本 - 自动循环不同参数组合
% 通过调用 PNM_batch 函数自动运行多组参数的溶解模拟
% 生成总结表包含 PBTime（渗透率达到初始值100倍的时间）
%
% 参数说明：
%   L (cm)       - 特征长度，即 circleSpacing 或 targetAvgSpacing [cm]
%   u_in (cm/s)  - 入口流速 [cm/s]
%   c_in         - 入口浓度 [mol/cm³]
%   Geometry     - 几何布局 ('hex' | 'square' | 'random')
%
% 时间步长控制逻辑（根据 Pe 和 Da 自动设置 Time_stepmax ）：
%   第一步：根据 Pe 确定基础值
%     Pe >= 10        -> base = 1   s
%     1   <= Pe < 10  -> base = 30  s
%     0.01<= Pe < 1   -> base = 90  s
%     Pe < 0.01       -> base = 300 s
%   第二步：根据 Da 确定倍率
%     Da >= 0.1       -> multiplier = 1
%     0.01<= Da < 0.1 -> multiplier = 5
%     Da < 0.01       -> multiplier = 10
%   最终：Time_stepmax = base * multiplier
%
% 输出表头：
%   Pe, Da, Pe/Da, Character, Time_stepmax, L(cm), u_in(cm/s), 
%   c_in(mol/cm³), D(cm²/s), molarVolume(cm3/mol), k_m(mol/cm²/s), 
%   k_m(mol/dm²/s), Geometry, PBTimeStep, PBTime
%
% 断点续做功能：
%   设置 resumeFromDir 为已有的批次结果文件夹路径，将跳过已完成的实验
%   留空 '' 则创建新的结果文件夹，全部重新运行
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
batchTimer = tic;

batchScriptDir = fileparts(mfilename('fullpath'));
rtmDir = fileparts(batchScriptDir);
reactiveRoot = fileparts(rtmDir);
projectRoot = fileparts(reactiveRoot);
addpath(batchScriptDir);
addpath(rtmDir);

batchOutputRoot = fullfile(projectRoot, 'outputs', 'rtm_batches');
if ~exist(batchOutputRoot, 'dir')
    mkdir(batchOutputRoot);
end

%% ===================== 断点续做设置 =====================
% 如果要从已有的批次结果继续运行，设置为已有文件夹路径
% 例如: resumeFromDir = fullfile(batchOutputRoot, 'batch_20260502_120000');
% 留空 '' 则创建新的结果文件夹，全部重新运行

resumeFromDir = '';
% resumeFromDir = fullfile(batchOutputRoot, 'batch_20260502_120000');

%% ===================== 输出策略 =====================
% exportEvery=1 表示每个时间步都导出 DXF，便于后续 COMSOL/NMR 全量处理。
% 如果只想先快速看趋势，可临时设为 5 或 10，并减少图片输出。
batchOutputOptions = struct();
batchOutputOptions.exportEvery = 1;
batchOutputOptions.exportDXF = true;
batchOutputOptions.saveMainPlot = true;
batchOutputOptions.saveIndividualPlots = false;
batchOutputOptions.saveInterfaceMask = true;
batchOutputOptions.saveRealtimePlot = false;
batchOutputOptions.saveFigureFiles = false;
batchOutputOptions.writeExcel = true;
batchOutputOptions.saveFinalPlot = true;
% true 时，每次导出 pore/solid DXF 后立即同步运行 COMSOL NMR + T2反演。
% NMR路径、COMSOL模型、Python解释器、覆盖策略等在 ReactiveTransport/automation/AutomationConfig.m 中设置。
batchOutputOptions.enableNMRSimulation = false;
%% ===================== 批次实验参数设置 =====================
% 是否使用手动定义的实验方案（与你的表格对应）
useManualSchemes = true;

if useManualSchemes
    % ========== 手动定义精确的实验方案（根据你的表格） ==========
    % 格式: {L_cm, u_in, c_in, 'layoutType'} 或 {L_cm, u_in, c_in, 'layoutType', Time_stepmax}
    % 如果指定第5个参数 Time_stepmax，则使用手动值；否则根据 Pe/Da 自动计算

    manualSchemes = {
    %     % === 方案6: Pe=0.05/0.10/0.20, hex ===
    %     % L=0.0005/0.001/0.002 cm, u=0.001 cm/s, c=1e-4 mol/cm³, hex
        {0.001, 0.001, 1e-4, 'random'};    % Pe=0.05, Da=0.01845
        {0.001,  0.01, 1e-4, 'random'};    % Pe=0.10, Da=0.0369
        {0.001,  0.1, 1e-4, 'random'};    % Pe=0.20, Da=0.0738
        {0.001, 0.001, 1e-4, 'hex'};    % Pe=0.05, Da=0.01845
        {0.001,  0.01, 1e-4, 'hex'};    % Pe=0.10, Da=0.0369
        {0.001,  0.1, 1e-4, 'hex'};    % Pe=0.20, Da=0.0738
      
        
    };

    
    paramCombinations = cell(size(manualSchemes, 1), 1);
    for i = 1:size(manualSchemes, 1)
        scheme = manualSchemes{i};
        paramCombinations{i} = struct(...
            'L_cm', scheme{1}, ...
            'u_in', scheme{2}, ...
            'c_in', scheme{3}, ...
            'layoutType', scheme{4} ...
        );
        % 如果指定了第5个参数（Time_stepmax），则标记为手动设置
        if length(scheme) >= 5 && ~isempty(scheme{5})
            paramCombinations{i}.Time_stepmax_manual = scheme{5};
        end
    end
else
    % ========== 使用笛卡尔积自动生成组合 ==========
    % 可根据需要调整以下列表（Time_stepmax 将由 Pe 自动计算）
    L_cm_list = [0.0005, 0.002];
    u_in_list = [0.001, 0.1];
    c_in_list = [1e-4];
    layoutType_list = {'hex', 'random'};
    
    paramCombinations = {};
    comboIndex = 0;
    for iL = 1:length(L_cm_list)
        for iU = 1:length(u_in_list)
            for iC = 1:length(c_in_list)
                for iG = 1:length(layoutType_list)
                    comboIndex = comboIndex + 1;
                    paramCombinations{comboIndex} = struct(...
                        'L_cm', L_cm_list(iL), ...
                        'u_in', u_in_list(iU), ...
                        'c_in', c_in_list(iC), ...
                        'layoutType', layoutType_list{iG} ...
                    );
                end
            end
        end
    end
end

%% ===================== 固定物理参数（用于 Pe/Da 及 Time_stepmax 计算） =====================
D_fixed = 1e-5;           % 扩散系数 [cm²/s]
molarVolume_fixed = 36.9; % 摩尔体积 [cm³/mol]
k_m_dm2_fixed = 1e-4;     % 反应速率系数 [mol/dm²/s]
k_m_cm2_fixed = k_m_dm2_fixed * 0.01;  % 转换为 [mol/cm²/s]

%% ===================== 根据 Pe 和 Da 自动设置 Time_stepmax =====================
% 逻辑说明：
%   1. 根据 Pe 确定基础值：
%      Pe >= 10       -> base = 1 s
%      1 <= Pe < 10   -> base = 30 s
%      0.01 <= Pe < 1 -> base = 90 s
%      Pe < 0.01      -> base = 300 s
%   2. 根据 Da 确定倍率：
%      Da >= 0.1      -> multiplier = 1
%      0.01 <= Da < 0.1 -> multiplier = 5
%      Da < 0.01      -> multiplier = 10
%   3. Time_stepmax = base * multiplier

for i = 1:numel(paramCombinations)
    p = paramCombinations{i};
    
    % 如果已手动指定 Time_stepmax，直接使用手动值
    if isfield(p, 'Time_stepmax_manual')
        p.Time_stepmax = p.Time_stepmax_manual;
        p = rmfield(p, 'Time_stepmax_manual');  % 移除临时字段
        paramCombinations{i} = p;
        continue;
    end
    
    % 否则根据 Pe 和 Da 自动计算
    Pe_i = p.u_in * p.L_cm / D_fixed;
    Da_i = p.c_in * molarVolume_fixed * k_m_dm2_fixed * 1000 * p.L_cm / D_fixed;
    
    % Pe-based 基础值
    if Pe_i >= 10
        base_step = 1;
    elseif Pe_i >= 1
        base_step = 5;
    elseif Pe_i >= 0.01
        base_step = 90;
    else
        base_step = 300;
    end
    
    % Da-based 倍率
    if Da_i >= 0.1
        da_mult = 1;
    elseif Da_i >= 0.01
        da_mult = 5;
    else
        da_mult = 10;
    end
    
    p.Time_stepmax = base_step * da_mult;
    paramCombinations{i} = p;
end

outputOptionNames = fieldnames(batchOutputOptions);
for i = 1:numel(paramCombinations)
    p = paramCombinations{i};
    for iOpt = 1:numel(outputOptionNames)
        optName = outputOptionNames{iOpt};
        if ~isfield(p, optName)
            p.(optName) = batchOutputOptions.(optName);
        end
    end
    paramCombinations{i} = p;
end

numExperiments = length(paramCombinations);
fprintf('========================================\n');
fprintf('批次实验：共 %d 组参数组合\n', numExperiments);
fprintf('========================================\n');
%% ===================== 创建或复用批次结果文件夹 =====================
if ~isempty(resumeFromDir) && exist(resumeFromDir, 'dir')
    % 断点续做模式：使用已有文件夹
    batchResultsDir = resumeFromDir;
    isResumeMode = true;
    fprintf('>>> 断点续做模式：使用已有文件夹 <<<\n');
    fprintf('批次结果目录: %s\n', batchResultsDir);
    
    % 扫描已完成的实验（通过检查 experiment_completed.log 文件）
    completedConfigs = scanCompletedExperiments(batchResultsDir);
    fprintf('已检测到 %d 个已完成的实验\n', length(completedConfigs));
else
    % 全新开始模式：创建新文件夹
    batchResultsDir = fullfile(batchOutputRoot, sprintf('batch_%s', datestr(now, 'yyyymmdd_HHMMSS')));
    if ~exist(batchResultsDir, 'dir')
        mkdir(batchResultsDir);
    end
    isResumeMode = false;
    completedConfigs = {};
    fprintf('>>> 全新开始模式 <<<\n');
    fprintf('批次结果保存目录: %s\n', batchResultsDir);
end

% PB（穿透）判断阈值
permeabilityRatioThreshold = 100;  % 渗透率达到初始值的100倍

%% ===================== 初始化总结表 =====================
% 预分配结果表（增加 Status 列记录成功/失败）
varNames = {'ExpIdx', 'Pe', 'Da', 'Pe_Da', 'Character', 'Time_stepmax', 'L_cm', ...
    'u_in_cm_s', 'c_in_mol_cm3', 'D_cm2_s', 'molarVolume_cm3_mol', ...
    'k_m_mol_cm2_s', 'k_m_mol_dm2_s', 'Geometry', 'PBTimeStep', 'PBTime', ...
    'InitialPerm_mD', 'FinalPerm_mD', 'Status'};
varTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', ...
    'double', 'double', 'cell', 'double', 'double', ...
    'double', 'double', 'cell'};
summaryTable = table('Size', [0, length(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);

%% ===================== 初始化错误日志 =====================
errorLogFile = fullfile(batchResultsDir, 'error_log.txt');
fid = fopen(errorLogFile, 'w');
fprintf(fid, '========================================\n');
fprintf(fid, '批次实验错误日志\n');
fprintf(fid, '开始时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '========================================\n\n');
fclose(fid);

% 统计变量
successCount = 0;
failCount = 0;
skippedCount = 0;  % 新增：跳过的实验数
failedExperiments = [];  % 记录失败的实验索引
skippedExperiments = []; % 新增：记录跳过的实验索引

%% ===================== 打印实验方案预览 =====================
fprintf('\n=== 实验方案预览 ===\n');
fprintf('%-5s %-8s %-10s %-10s %-10s %-10s %-10s %-8s\n', ...
    'Exp', 'Layout', 'L(cm)', 'u_in', 'c_in', 'Pe', 'Da', 'Status');
fprintf('%s\n', repmat('-', 1, 85));

for expIdx = 1:numExperiments
    params = paramCombinations{expIdx};
    Pe = params.u_in * params.L_cm / D_fixed;
    Da = params.c_in * molarVolume_fixed * k_m_dm2_fixed * 1000 * params.L_cm / D_fixed;
    
    % 检查是否已完成
    configKey = generateConfigKey(params, D_fixed, molarVolume_fixed, k_m_dm2_fixed);
    if isResumeMode && any(strcmp(completedConfigs, configKey))
        statusStr = '[已完成]';
    else
        statusStr = '[待运行]';
    end
    
    fprintf('%-5d %-8s %-10.4f %-10.4f %-10.2e %-10.4f %-10.6f %s\n', ...
        expIdx, params.layoutType, params.L_cm, params.u_in, params.c_in, Pe, Da, statusStr);
end
fprintf('%s\n\n', repmat('-', 1, 85));

%% ===================== 循环运行每组实验 =====================
for expIdx = 1:numExperiments
    fprintf('\n');
    fprintf('########################################\n');
    fprintf('## 实验 %d / %d\n', expIdx, numExperiments);
    fprintf('########################################\n');
    
    % 获取当前参数组
    params = paramCombinations{expIdx};
    
    % 计算 Pe, Da
    Pe = params.u_in * params.L_cm / D_fixed;
    Da = params.c_in * molarVolume_fixed * k_m_dm2_fixed * 1000 * params.L_cm / D_fixed;
    PeDaRatio = Pe / Da;
    
    % === 检查是否已完成（断点续做） ===
    configKey = generateConfigKey(params, D_fixed, molarVolume_fixed, k_m_dm2_fixed);
    if isResumeMode && any(strcmp(completedConfigs, configKey))
        fprintf('>>> 实验已完成，跳过 <<<\n');
        fprintf('配置: Pe=%.4f, Da=%.6f, Geometry=%s\n', Pe, Da, params.layoutType);
        skippedCount = skippedCount + 1;
        skippedExperiments = [skippedExperiments, expIdx];
        continue;  % 跳过此实验
    end
    
    fprintf('参数: Time_stepmax=%d, L=%.4f cm, u_in=%.4f cm/s, c_in=%.2e mol/cm³, Geometry=%s\n', ...
        params.Time_stepmax, params.L_cm, params.u_in, params.c_in, params.layoutType);
    fprintf('计算: Pe=%.4f, Da=%.6f, Pe/Da=%.4f\n', Pe, Da, PeDaRatio);
    
    % 调用单次模拟函数
    try
        result = PNM_batch(params, batchResultsDir, expIdx, permeabilityRatioThreshold);
        
        % 构建成功行
        newRow = {expIdx, Pe, Da, PeDaRatio, params.L_cm, params.Time_stepmax, params.L_cm, ...
            params.u_in, params.c_in, D_fixed, molarVolume_fixed, ...
            k_m_cm2_fixed, k_m_dm2_fixed, ...
            {params.layoutType}, result.PBTimeStep, result.PBTime, ...
            result.initialPermeability, result.finalPermeability, {'Success'}};
        
        summaryTable = [summaryTable; newRow];
        successCount = successCount + 1;
        
        % === 生成子实验完成日志 ===
        writeExperimentCompletedLog(result.resultsDir, params, Pe, Da, result, D_fixed, molarVolume_fixed, k_m_dm2_fixed);
        
        fprintf('## 实验 %d 完成: PBTimeStep=%d, PBTime=%.2f s\n', expIdx, result.PBTimeStep, result.PBTime);
        
    catch ME
        failCount = failCount + 1;
        failedExperiments = [failedExperiments, expIdx];
        
        % 控制台输出错误信息
        fprintf('\n');
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        fprintf('!! 实验 %d 执行失败，已跳过\n', expIdx);
        fprintf('!! 错误信息: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('!! 错误位置: %s (行 %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        
        % 写入错误日志文件
        fid = fopen(errorLogFile, 'a');
        fprintf(fid, '----------------------------------------\n');
        fprintf(fid, '实验 %d 失败\n', expIdx);
        fprintf(fid, '时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, '参数:\n');
        fprintf(fid, '  Time_stepmax = %d\n', params.Time_stepmax);
        fprintf(fid, '  L_cm = %.6f\n', params.L_cm);
        fprintf(fid, '  u_in = %.6f cm/s\n', params.u_in);
        fprintf(fid, '  c_in = %.2e mol/cm³\n', params.c_in);
        fprintf(fid, '  Geometry = %s\n', params.layoutType);
        fprintf(fid, '  Pe = %.6f, Da = %.6f\n', Pe, Da);
        fprintf(fid, '错误信息: %s\n', ME.message);
        fprintf(fid, '错误标识: %s\n', ME.identifier);
        if ~isempty(ME.stack)
            fprintf(fid, '调用栈:\n');
            for stackIdx = 1:min(5, length(ME.stack))  % 最多记录5层调用栈
                fprintf(fid, '  [%d] %s (行 %d)\n', stackIdx, ...
                    ME.stack(stackIdx).name, ME.stack(stackIdx).line);
            end
        end
        fprintf(fid, '\n');
        fclose(fid);
        
        % 记录失败行（PBTimeStep 和 PBTime 为 NaN）
        newRow = {expIdx, Pe, Da, PeDaRatio, params.L_cm, params.Time_stepmax, params.L_cm, ...
            params.u_in, params.c_in, D_fixed, molarVolume_fixed, ...
            k_m_cm2_fixed, k_m_dm2_fixed, ...
            {params.layoutType}, NaN, NaN, NaN, NaN, {'Failed'}};
        
        summaryTable = [summaryTable; newRow];
        
        % 继续下一个实验
        fprintf('>> 继续执行下一个实验...\n');
    end
    
    % 每次实验后保存中间结果（防止意外中断丢失数据）
    writetable(summaryTable, fullfile(batchResultsDir, 'batch_summary_partial.xlsx'));
    writetable(summaryTable, fullfile(batchResultsDir, 'batch_summary_partial.csv'));
end

%% ===================== 保存最终总结表 =====================
summaryXlsxPath = fullfile(batchResultsDir, 'batch_summary_final.xlsx');
summaryCsvPath = fullfile(batchResultsDir, 'batch_summary_final.csv');

writetable(summaryTable, summaryXlsxPath);
writetable(summaryTable, summaryCsvPath);

%% ===================== 生成汇总图表 =====================
try
    % 创建汇总可视化
    summaryFig = figure('Position', [100, 100, 1400, 1000]);
    
    % 提取成功实验的数据
    successIdx = strcmp(summaryTable.Status, 'Success');
    if any(successIdx)
        successData = summaryTable(successIdx, :);
        
        % 子图1: PBTime vs Pe/Da
        subplot(2,2,1);
        scatter(successData.Pe_Da, successData.PBTime, 100, 'filled');
        xlabel('Pe/Da'); ylabel('PBTime [s]');
        title('Breakthrough Time vs Pe/Da');
        grid on; set(gca, 'XScale', 'log');
        
        % 子图2: PBTime vs Pe (按 Geometry 分组)
        subplot(2,2,2);
        hold on;
        geoTypes = unique(successData.Geometry);
        colors = {'b', 'r', 'g', 'm'};
        for iGeo = 1:length(geoTypes)
            idx = strcmp(successData.Geometry, geoTypes{iGeo});
            scatter(successData.Pe(idx), successData.PBTime(idx), 100, colors{iGeo}, 'filled', 'DisplayName', geoTypes{iGeo});
        end
        xlabel('Pe'); ylabel('PBTime [s]');
        title('Breakthrough Time vs Pe by Geometry');
        legend('Location', 'best'); grid on; set(gca, 'XScale', 'log');
        
        % 子图3: Initial vs Final Permeability
        subplot(2,2,3);
        validPerm = ~isnan(successData.InitialPerm_mD) & ~isnan(successData.FinalPerm_mD);
        if any(validPerm)
            scatter(successData.InitialPerm_mD(validPerm), successData.FinalPerm_mD(validPerm), 100, 'filled');
            hold on;
            xRange = xlim;
            plot(xRange, xRange * 100, 'r--', 'LineWidth', 2, 'DisplayName', '100x threshold');
            xlabel('Initial Permeability [mD]'); ylabel('Final Permeability [mD]');
            title('Permeability Evolution');
            legend('Location', 'best'); grid on;
            set(gca, 'XScale', 'log', 'YScale', 'log');
        end
        
        % 子图4: 实验结果统计
        subplot(2,2,4);
        bar([successCount, skippedCount, failCount]);
        xticklabels({'Success', 'Skipped', 'Failed'});
        ylabel('Number of Experiments');
        title(sprintf('Batch Results Summary (Total: %d)', numExperiments));
        text(1, successCount + 0.5, num2str(successCount), 'HorizontalAlignment', 'center');
        text(2, skippedCount + 0.5, num2str(skippedCount), 'HorizontalAlignment', 'center');
        text(3, failCount + 0.5, num2str(failCount), 'HorizontalAlignment', 'center');
        grid on;
    end
    
    saveas(summaryFig, fullfile(batchResultsDir, 'batch_summary_plots.png'));
    saveas(summaryFig, fullfile(batchResultsDir, 'batch_summary_plots.fig'));
    close(summaryFig);
catch ME
    warning(ME.identifier, '汇总图表生成失败: %s', ME.message);
end

%% ===================== 写入最终错误日志摘要 =====================
fid = fopen(errorLogFile, 'a');
fprintf(fid, '\n========================================\n');
fprintf(fid, '批次实验完成摘要\n');
fprintf(fid, '结束时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '========================================\n');
fprintf(fid, '总实验数: %d\n', numExperiments);
fprintf(fid, '成功: %d\n', successCount);
fprintf(fid, '跳过(已完成): %d\n', skippedCount);
fprintf(fid, '失败: %d\n', failCount);
if ~isempty(failedExperiments)
    fprintf(fid, '\n失败的实验编号: ');
    fprintf(fid, '%d ', failedExperiments);
    fprintf(fid, '\n\n');
    fprintf(fid, '失败实验详情:\n');
    for i = 1:length(failedExperiments)
        idx = failedExperiments(i);
        p = paramCombinations{idx};
        fprintf(fid, '  实验 %d: Time_stepmax=%d, L=%.4f, u_in=%.4f, c_in=%.2e, Geometry=%s\n', ...
            idx, p.Time_stepmax, p.L_cm, p.u_in, p.c_in, p.layoutType);
    end
end
fprintf(fid, '========================================\n');
fclose(fid);

%% ===================== 最终控制台输出 =====================
fprintf('\n');
fprintf('========================================\n');
fprintf('批次实验完成！\n');
fprintf('总耗时: %.2f 分钟\n', toc(batchTimer)/60);
fprintf('----------------------------------------\n');
fprintf('成功: %d / %d\n', successCount, numExperiments);
fprintf('跳过(已完成): %d / %d\n', skippedCount, numExperiments);
fprintf('失败: %d / %d\n', failCount, numExperiments);
if ~isempty(failedExperiments)
    fprintf('失败的实验编号: ');
    fprintf('%d ', failedExperiments);
    fprintf('\n');
    fprintf('详细错误信息请查看: %s\n', errorLogFile);
end
fprintf('----------------------------------------\n');
fprintf('结果保存至: %s\n', batchResultsDir);
fprintf('总结表: %s\n', summaryXlsxPath);
fprintf('========================================\n');

% 显示总结表
disp(' ');
disp('=== 批次实验总结表 ===');
disp(summaryTable);

%% ===================== 保存工作区变量 =====================
save(fullfile(batchResultsDir, 'batch_workspace.mat'), ...
    'summaryTable', 'paramCombinations', 'batchResultsDir', ...
    'successCount', 'failCount', 'skippedCount', 'failedExperiments', 'skippedExperiments', ...
    'D_fixed', 'molarVolume_fixed', 'k_m_dm2_fixed', 'k_m_cm2_fixed', ...
    'permeabilityRatioThreshold', 'isResumeMode');

fprintf('\n工作区变量已保存至: %s\n', fullfile(batchResultsDir, 'batch_workspace.mat'));

%% ===================== 辅助函数定义 =====================

function configKey = generateConfigKey(params, D_fixed, molarVolume_fixed, k_m_dm2_fixed)
% 生成唯一的实验配置键（用于识别已完成的实验）
% 基于 Pe, Da 和布局类型生成唯一标识
    Pe = params.u_in * params.L_cm / D_fixed;
    Da = params.c_in * molarVolume_fixed * k_m_dm2_fixed * 1000 * params.L_cm / D_fixed;
    configKey = sprintf('Pe_%.6f_Da_%.6f_%s', Pe, Da, params.layoutType);
end

function completedConfigs = scanCompletedExperiments(batchResultsDir)
% 扫描批次结果文件夹，识别已完成的实验
% 通过检查每个子文件夹中的 experiment_completed.log 文件
    completedConfigs = {};
    
    % 获取所有子文件夹
    items = dir(batchResultsDir);
    subDirs = items([items.isdir]);
    
    for i = 1:length(subDirs)
        subDirName = subDirs(i).name;
        if strcmp(subDirName, '.') || strcmp(subDirName, '..')
            continue;
        end
        
        % 检查是否存在完成日志
        logFile = fullfile(batchResultsDir, subDirName, 'experiment_completed.log');
        if exist(logFile, 'file')
            % 读取配置键
            try
                fid = fopen(logFile, 'r');
                content = fread(fid, '*char')';
                fclose(fid);
                
                % 提取配置键（在 CONFIG_KEY: 行）
                keyMatch = regexp(content, 'CONFIG_KEY:\s*(\S+)', 'tokens');
                if ~isempty(keyMatch)
                    completedConfigs{end+1} = keyMatch{1}{1}; %#ok<AGROW>
                end
            catch
                % 忽略读取错误
            end
        end
    end
end

function writeExperimentCompletedLog(resultsDir, params, Pe, Da, result, D_fixed, molarVolume_fixed, k_m_dm2_fixed)
% 为已完成的子实验写入完成日志
    logFile = fullfile(resultsDir, 'experiment_completed.log');
    configKey = sprintf('Pe_%.6f_Da_%.6f_%s', Pe, Da, params.layoutType);
    
    fid = fopen(logFile, 'w');
    fprintf(fid, '========================================\n');
    fprintf(fid, '实验完成日志\n');
    fprintf(fid, '========================================\n');
    fprintf(fid, '完成时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '\n');
    fprintf(fid, '--- 配置标识 ---\n');
    fprintf(fid, 'CONFIG_KEY: %s\n', configKey);
    fprintf(fid, '\n');
    fprintf(fid, '--- 输入参数 ---\n');
    fprintf(fid, 'Time_stepmax: %d s\n', params.Time_stepmax);
    fprintf(fid, 'L_cm: %.6f cm\n', params.L_cm);
    fprintf(fid, 'u_in: %.6f cm/s\n', params.u_in);
    fprintf(fid, 'c_in: %.2e mol/cm³\n', params.c_in);
    fprintf(fid, 'layoutType: %s\n', params.layoutType);
    fprintf(fid, '\n');
    fprintf(fid, '--- 计算参数 ---\n');
    fprintf(fid, 'Pe: %.6f\n', Pe);
    fprintf(fid, 'Da: %.6f\n', Da);
    fprintf(fid, 'Pe/Da: %.6f\n', Pe/Da);
    fprintf(fid, 'D: %.2e cm²/s\n', D_fixed);
    fprintf(fid, 'molarVolume: %.2f cm³/mol\n', molarVolume_fixed);
    fprintf(fid, 'k_m: %.2e mol/dm²/s\n', k_m_dm2_fixed);
    fprintf(fid, '\n');
    fprintf(fid, '--- 结果摘要 ---\n');
    fprintf(fid, 'PBTimeStep: %d\n', result.PBTimeStep);
    fprintf(fid, 'PBTime: %.4f s\n', result.PBTime);
    fprintf(fid, 'Initial Permeability: %.4f mD\n', result.initialPermeability);
    fprintf(fid, 'Final Permeability: %.4f mD\n', result.finalPermeability);
    fprintf(fid, 'Results Directory: %s\n', resultsDir);
    fprintf(fid, '========================================\n');
    fclose(fid);
end
