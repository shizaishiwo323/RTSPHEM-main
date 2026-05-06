%% RunBatchExperiments_simple - 批量运行 RTM + 可选同步 NMR 的精简入口
%
% 用法：
%   1. 只修改本文件顶部的“用户可调参数”。
%   2. 在 MATLAB 中运行本文件。
%   3. geometryCases 定义几何尺寸，peList 定义每个几何要跑的 Pe regime。
%   4. 如果 batchOptions.enableNMRSimulation=true，NMR 细节在
%      ReactiveTransport/automation/AutomationConfig.m 中设置。
%
% 输出：
%   outputs/rtm_batches/batch_YYYYMMDD_HHMMSS/
%     exp_001/, exp_002/, ...
%     batch_summary_simple.xlsx
%     batch_summary_simple.csv
%     batch_error_log.txt
 
clear; clc;
batchTimer = tic;

%% ===================== 路径与批次目录 =====================
batchScriptDir = fileparts(mfilename('fullpath'));
rtmDir = fileparts(batchScriptDir);
reactiveRoot = fileparts(rtmDir);
projectRoot = fileparts(reactiveRoot);
addpath(batchScriptDir);
addpath(rtmDir);

% 批量输出根目录。
batchOutputRoot = fullfile(projectRoot, 'outputs', 'rtm_batches');

% 批次名。留空 '' 时自动使用 batch_时间戳。
batchName = '';

% 如需固定写入某个已有/指定目录，取消下一行注释。
% batchResultsDir = fullfile(batchOutputRoot, 'my_batch_test');

%% ===================== 全局物理参数 =====================
% 这些参数用于实际 RTM 计算，也用于脚本预估 Pe、Da 和自动 Time_stepmax。
physics = struct();
physics.diffusionCoefficient = 1e-5;  % [cm^2/s]
physics.molarVolume = 36.9;           % [cm^3/mol]
physics.rateCoefficientTST = 1e-4;    % [mol/dm^2/s]

%% ===================== 第一版数据集设计 =====================
% 目标：为后续 NMR-agent 本地训练准备一版小而有判别力的 RTM 样本。
%
% 全局固定参数：
%   L_cm      = 0.001 cm，固定孔喉/特征长度，避免 L 与几何尺寸同时变化。
%   Da        = 0.0369，固定反应强度；c_in 由 Da、L 和物理参数反推。
%   Pe_list   = [0.1, 1, 10]，只改变入口速度 u_in 来覆盖扩散/过渡/对流 regime。
%
% 每个轮次变化参数：
%   shapeFamily        - 数据标签用形状：square | rectangle | random。
%   layoutType         - PNM_beauty3 求解器实际布局：square | random。
%   targetLengthX/Y_um - 目标模拟域尺寸，范围覆盖 200-1000 um。
%   Pe_target          - 当前轮次目标 Pe；u_in = Pe_target * D / L。
%
% 说明：
%   PNM_beauty3 没有单独的 'rectangle' layoutType。rectangle 使用规则 square
%   排布加非 1:1 的 targetAspectRatio 表示，并在 summary 中保留 shapeFamily。
fixedDesign = struct();
fixedDesign.L_cm = 0.001;
fixedDesign.Da = 0.0369;
fixedDesign.c_in = concentrationForDa(fixedDesign.Da, fixedDesign.L_cm, physics);

peList = [0.1, 1, 10];

% 格式：
%   {'geometryCase', 'shapeFamily', 'layoutType', targetLengthX_um, targetLengthY_um}
geometryCases = {
    'square_200',         'square',    'square',  200,  200;
    'square_400',         'square',    'square',  400,  400;
    'square_600',         'square',    'square',  600,  600;
    'square_800',         'square',    'square',  800,  800;
    'square_1000',        'square',    'square', 1000, 1000;
    'rectangle_600x400',  'rectangle', 'square',  600,  400;
    'rectangle_1000x600', 'rectangle', 'square', 1000,  600;
    'rectangle_400x800',  'rectangle', 'square',  400,  800;
    'random_600x400',     'random',    'random',  600,  400;
    'random_1000x600',    'random',    'random', 1000,  600;
    'random_400x800',     'random',    'random',  400,  800;
};

%% ===================== 全局几何参数 =====================
% 所有实验共用的几何控制。目标 X/Y 尺寸由 geometryCases 逐轮覆盖。
geometry = struct();
geometry.loadExistingGeometry = false;   % 批量生成数据时通常为 false
geometry.geometryLoadFile = '';
geometry.geometrySaveFile = '';
geometry.useExternalGeometry = false;
geometry.tifPath = "";

%% ===================== 时间步与终止条件 =====================
timeControl = struct();
timeControl.initialMacroscaleTimeStepSize = 0.10;  % [s]
timeControl.endTime = [];                          % [] 表示由 PNM_beauty3 自动估算
timeControl.permeabilityRatioThreshold = 100;      % k/k0 达到该倍数后完成当前步并停止

%% ===================== 网格与导出精度 =====================
meshExport = struct();
meshExport.numPartitionsMicroscale = 2 * 64;
meshExport.dxfMicronsPerPixel = 4;  % 派生 DXF 规则网格分辨率，约 4 um/像素
meshExport.minDxfResolution = 80;   % 小尺寸样本也保留最低导出分辨率

%% ===================== 输出与 NMR 控制 =====================
batchOptions = struct();

% 每隔多少个 RTM 时间步导出一次结构。同步 NMR 时建议保持 1。
batchOptions.exportEvery = 1;

% 同步 NMR 需要 DXF，打开 enableNMRSimulation 时应保持 true。
batchOptions.exportDXF = true;

% 图像和表格输出。
batchOptions.saveMainPlot = true;
batchOptions.saveIndividualPlots = true;
batchOptions.saveInterfaceMask = true;
batchOptions.saveRealtimePlot = false;
batchOptions.saveFigureFiles = false;
batchOptions.writeExcel = true;
batchOptions.saveFinalPlot = true;

% false：只批量跑 RTM；true：每次导出 DXF 后同步跑 COMSOL NMR + T2 反演。
batchOptions.enableNMRSimulation = true;

%% ===================== 由参数表构建实验列表 =====================
if ~exist(batchOutputRoot, 'dir')
    mkdir(batchOutputRoot);
end

if exist('batchResultsDir', 'var') && strlength(string(batchResultsDir)) > 0
    batchResultsDir = char(batchResultsDir);
elseif strlength(string(batchName)) > 0
    batchResultsDir = fullfile(batchOutputRoot, char(batchName));
else
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    batchResultsDir = fullfile(batchOutputRoot, sprintf('batch_%s', timestamp));
end
if ~exist(batchResultsDir, 'dir')
    mkdir(batchResultsDir);
end

paramList = {};
expCounter = 0;
for iCase = 1:size(geometryCases, 1)
    row = geometryCases(iCase, :);
    targetX_um = row{4};
    targetY_um = row{5};
    targetX_cm = micronToCm(targetX_um);
    targetY_cm = micronToCm(targetY_um);

    for iPe = 1:numel(peList)
        expCounter = expCounter + 1;
        Pe_target = peList(iPe);

        p = struct();
        p.geometryCase = row{1};
        p.shapeFamily = row{2};
        p.layoutType = row{3};
        p.targetLengthXAxis_um = targetX_um;
        p.targetLengthYAxis_um = targetY_um;
        p.targetLengthXAxis = targetX_cm;
        p.targetLengthYAxis = targetY_cm;
        p.targetAspectRatio = targetX_cm / targetY_cm;
        p.Pe_target = Pe_target;
        p.Da_target = fixedDesign.Da;
        p.L_cm = fixedDesign.L_cm;
        p.u_in = Pe_target * physics.diffusionCoefficient / fixedDesign.L_cm;
        p.c_in = fixedDesign.c_in;
        p.Time_stepmax = estimateTimeStepmax(p, physics);
        [p.dxfResolutionX, p.dxfResolutionY] = estimateDxfResolution( ...
            targetX_um, targetY_um, meshExport);

        p = mergeStructs(p, physics);
        p = mergeStructs(p, geometry);
        p = mergeStructs(p, timeControl);
        p = mergeStructs(p, meshExport);
        p = mergeStructs(p, batchOptions);
        paramList{expCounter, 1} = p; %#ok<SAGROW>
    end
end

%% ===================== 打印预览 =====================
fprintf('========================================\n');
fprintf('精简批量 RTM/NMR 运行\n');
fprintf('  批次目录: %s\n', batchResultsDir);
fprintf('  实验数量: %d\n', numel(paramList));
fprintf('  sync NMR: %s\n', mat2str(batchOptions.enableNMRSimulation));
fprintf('  固定 L: %.4g cm | 固定 Da: %.4g | 固定 c_in: %.4g mol/cm^3\n', ...
    fixedDesign.L_cm, fixedDesign.Da, fixedDesign.c_in);
fprintf('========================================\n\n');

fprintf('%-5s %-18s %-10s %-8s %-9s %-9s %-8s %-10s %-10s %-8s\n', ...
    'Exp', 'Case', 'Shape', 'Layout', 'X(um)', 'Y(um)', 'Pe', 'u_in', 'Da', 'dtmax');
fprintf('%s\n', repmat('-', 1, 116));
for i = 1:numel(paramList)
    p = paramList{i};
    [Pe, Da] = calcPeDa(p, physics);
    fprintf('%-5d %-18s %-10s %-8s %-9.0f %-9.0f %-8.4g %-10.4g %-10.4g %-8.4g\n', ...
        i, p.geometryCase, p.shapeFamily, p.layoutType, p.targetLengthXAxis_um, ...
        p.targetLengthYAxis_um, Pe, p.u_in, Da, p.Time_stepmax);
end
fprintf('%s\n\n', repmat('-', 1, 116));

designTable = buildDesignTable(paramList, physics);
writetable(designTable, fullfile(batchResultsDir, 'batch_design_table.xlsx'));
writetable(designTable, fullfile(batchResultsDir, 'batch_design_table.csv'));

%% ===================== 执行批量实验 =====================
summary = table();
errorLogFile = fullfile(batchResultsDir, 'batch_error_log.txt');
fid = fopen(errorLogFile, 'w');
fprintf(fid, 'Batch error log - %s\n\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fclose(fid);

for expIdx = 1:numel(paramList)
    p = paramList{expIdx};
    [Pe, Da] = calcPeDa(p, physics);
    PeDa = Pe / Da;

    fprintf('\n########################################\n');
    fprintf('实验 %d/%d | %s | %s | %gx%g um | L=%.4g cm\n', ...
        expIdx, numel(paramList), p.geometryCase, p.shapeFamily, ...
        p.targetLengthXAxis_um, p.targetLengthYAxis_um, p.L_cm);
    fprintf('u=%.4g cm/s, c=%.2e mol/cm^3\n', p.u_in, p.c_in);
    fprintf('Pe=%.6g, Da=%.6g, Time_stepmax=%.6g s\n', Pe, Da, p.Time_stepmax);
    fprintf('########################################\n');

    try
        result = PNM_batch(p, batchResultsDir, expIdx, timeControl.permeabilityRatioThreshold);
        status = "Success";
        message = "";
        resultsDir = string(result.resultsDir);
        finalPorosity = result.finalPorosity;
        initialPerm = result.initialPermeability;
        finalPerm = result.finalPermeability;
        pbStep = result.PBTimeStep;
        pbTime = result.PBTime;
        [actualX_um, actualY_um] = getActualSizeMicrons(result);
    catch ME
        status = "Failed";
        message = string(ME.message);
        resultsDir = "";
        finalPorosity = NaN;
        initialPerm = NaN;
        finalPerm = NaN;
        pbStep = NaN;
        pbTime = NaN;
        actualX_um = NaN;
        actualY_um = NaN;

        fid = fopen(errorLogFile, 'a');
        if fid ~= -1
            fprintf(fid, '[exp_%03d] %s\n', expIdx, ME.message);
            fclose(fid);
        end
        warning('BatchSimple:ExperimentFailed', '实验 %d 失败: %s', expIdx, ME.message);
    end

    newRow = table(expIdx, string(p.geometryCase), string(p.shapeFamily), string(p.layoutType), ...
        p.targetLengthXAxis_um, p.targetLengthYAxis_um, actualX_um, actualY_um, ...
        p.targetAspectRatio, p.dxfResolutionX, p.dxfResolutionY, p.L_cm, ...
        p.Pe_target, p.Da_target, p.u_in, p.c_in, Pe, Da, PeDa, ...
        p.Time_stepmax, finalPorosity, initialPerm, finalPerm, pbStep, pbTime, ...
        status, message, resultsDir, ...
        'VariableNames', {'ExpIdx', 'GeometryCase', 'ShapeFamily', 'SolverLayout', ...
        'TargetX_um', 'TargetY_um', 'ActualX_um', 'ActualY_um', ...
        'TargetAspectRatio', 'DxfResolutionX', 'DxfResolutionY', 'L_cm', ...
        'Pe_target', 'Da_target', 'u_in_cm_s', 'c_in_mol_cm3', ...
        'Pe', 'Da', 'Pe_Da', 'Time_stepmax_s', 'FinalPorosity', ...
        'InitialPerm_mD', 'FinalPerm_mD', 'PBTimeStep', 'PBTime_s', ...
        'Status', 'Message', 'ResultsDir'});
    summary = [summary; newRow]; %#ok<AGROW>

    writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple_partial.xlsx'));
    writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple_partial.csv'));
end

%% ===================== 保存最终总结 =====================
writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple.xlsx'));
writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple.csv'));
save(fullfile(batchResultsDir, 'batch_workspace_simple.mat'), ...
    'paramList', 'summary', 'designTable', 'batchResultsDir', 'physics', ...
    'fixedDesign', 'peList', 'geometryCases', 'geometry', ...
    'timeControl', 'meshExport', 'batchOptions');

fprintf('\n========================================\n');
fprintf('批量运行完成，用时 %.2f 秒\n', toc(batchTimer));
fprintf('结果目录: %s\n', batchResultsDir);
fprintf('总结表: %s\n', fullfile(batchResultsDir, 'batch_summary_simple.xlsx'));
fprintf('========================================\n');

%% ===================== 本脚本辅助函数 =====================
function p = mergeStructs(p, extra)
names = fieldnames(extra);
for i = 1:numel(names)
    name = names{i};
    p.(name) = extra.(name);
end
end

function [Pe, Da] = calcPeDa(p, physics)
Pe = p.u_in * p.L_cm / physics.diffusionCoefficient;
Da = p.c_in * physics.molarVolume * physics.rateCoefficientTST * 1000 * p.L_cm / physics.diffusionCoefficient;
end

function cm = micronToCm(micronValue)
cm = micronValue * 1e-4;
end

function c_in = concentrationForDa(targetDa, L_cm, physics)
c_in = targetDa * physics.diffusionCoefficient / ...
    (physics.molarVolume * physics.rateCoefficientTST * 1000 * L_cm);
end

function [resolutionX, resolutionY] = estimateDxfResolution(targetX_um, targetY_um, meshExport)
resolutionX = max(meshExport.minDxfResolution, round(targetX_um / meshExport.dxfMicronsPerPixel));
resolutionY = max(meshExport.minDxfResolution, round(targetY_um / meshExport.dxfMicronsPerPixel));
end

function designTable = buildDesignTable(paramList, physics)
designTable = table();
for i = 1:numel(paramList)
    p = paramList{i};
    [Pe, Da] = calcPeDa(p, physics);
    row = table(i, string(p.geometryCase), string(p.shapeFamily), string(p.layoutType), ...
        p.targetLengthXAxis_um, p.targetLengthYAxis_um, p.targetAspectRatio, ...
        p.dxfResolutionX, p.dxfResolutionY, p.L_cm, p.Pe_target, p.Da_target, ...
        p.u_in, p.c_in, Pe, Da, p.Time_stepmax, ...
        'VariableNames', {'ExpIdx', 'GeometryCase', 'ShapeFamily', 'SolverLayout', ...
        'TargetX_um', 'TargetY_um', 'TargetAspectRatio', 'DxfResolutionX', ...
        'DxfResolutionY', 'L_cm', 'Pe_target', 'Da_target', 'u_in_cm_s', ...
        'c_in_mol_cm3', 'Pe', 'Da', 'Time_stepmax_s'});
    designTable = [designTable; row]; %#ok<AGROW>
end
end

function [actualX_um, actualY_um] = getActualSizeMicrons(result)
actualX_um = NaN;
actualY_um = NaN;
if isfield(result, 'metadata') && isfield(result.metadata, 'parameters')
    params = result.metadata.parameters;
    if isfield(params, 'lengthXAxis_cm')
        actualX_um = params.lengthXAxis_cm * 1e4;
    end
    if isfield(params, 'lengthYAxis_cm')
        actualY_um = params.lengthYAxis_cm * 1e4;
    end
end
end

function dt = estimateTimeStepmax(p, physics)
[Pe, Da] = calcPeDa(p, physics);

if Pe >= 10
    baseStep = 1;
elseif Pe >= 1
    baseStep = 5;
elseif Pe >= 0.01
    baseStep = 90;
else
    baseStep = 300;
end

if Da >= 0.1
    daMultiplier = 1;
elseif Da >= 0.01
    daMultiplier = 5;
else
    daMultiplier = 10;
end

dt = baseStep * daMultiplier;
end
