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

% 运行模式：
%   'new'    - 新建一个 batch_时间戳目录，从 exp_001 开始。
%   'resume' - 复用已有批次目录，自动跳过目录中已经存在的 exp_NNN。
batchRunMode = 'resume';

% 续跑时填写已有批次目录。例如：
% resumeBatchResultsDir = fullfile(batchOutputRoot, 'batch_20260508_074936');
resumeBatchResultsDir = 'C:\Users\imgw\Documents\Codex\RTSPHEM-main\outputs\rtm_batches\batch_20260508_074936';

% 如需固定写入某个新目录，可设置 batchName，或取消下一行注释。
% batchResultsDir = fullfile(batchOutputRoot, 'my_batch_test');

%% ===================== 全局物理参数 =====================
% 这些参数提供默认值；每个 Da 组会覆盖 D_H+ 和 k_TST。
physics = struct();
physics.diffusionCoefficient = 1e-4;  % [cm^2/s], per-case D_H+ will override this.
physics.molarVolume = 36.9;           % [cm^3/mol]
physics.rateCoefficientTST = 1e-4;    % Solver convention; treated consistently in calcPeDa.

%% ===================== 600x400 um Pe-Da 批次设计 =====================
% 目标：固定模拟域约 600 um x 400 um，系统扫描 Pe 与 Da。
%
% 全局固定参数：
%   targetLengthX/Y = 600/400 um，固定几何尺寸。
%   L_cm            = 孔喉/平均孔隙间距，作为 Pe/Da/Re 使用的统一特征长度。
%   c_in            = 1 / molarVolume，使本脚本的 Da 表与现有求解器定义一致：
%                     Da = c_in * molarVolume * k_TST * 1000 * L / D
%                        = 1000 * k_TST * L / D
%   Pe_list         = 0.01 到 50，对数式覆盖扩散主导到强对流。
%
% 每个轮次变化参数：
%   geometryCase      - hex、random_1、random_2 三个几何因子。
%   layoutType        - PNM_beauty3 求解器布局：hex | random。
%   Pe_target         - 当前目标 Pe；u_in = Pe_target * D_H+ / L_cm。
%   Da_target         - 当前目标 Da；由 D_H+ 和 k_TST 共同指定。
%
% 说明：
%   Da=0.01-1 为低 Da 数值敏感性扩展，k_TST 会低于常规文献范围。
%   每个 Pe-Da 组都跑 1 个 hex 和 2 个多粒径 random 几何，共 15*9*3=405 组。
%   D_H+ 统一取 3e-5 cm^2/s。若某个 Pe 在当前 L_cm 下导致 u>0.3 cm/s，
%   该组合会写入 skipped 设计表并跳过，不启动模拟。
fixedDesign = struct();
fixedDesign.targetLengthX_um = 600;
fixedDesign.targetLengthY_um = 400;
fixedDesign.c_in = 1 / physics.molarVolume;
fixedDesign.maxInletVelocity = 0.3;  % [cm/s], Stokes-flow safeguard requested by user.

peList = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 50];

% 格式：
%   {'geometryCase', 'shapeFamily', 'layoutType', targetLengthX_um, targetLengthY_um,
%    geometryFactor, randomSeed, particleRadius_um, throatSpacing_um, minThroatRandom_um,
%    randomDensityFactor, useRandomParticleRadii, randomRadiusMin_um, randomRadiusMax_um, targetInitialPorosity}
geometryCases = {
    'hex_600x400_R50_throat20',          'hex',    'hex',    600, 400, 1, NaN,        50.0, 20.0, NaN, 1.0, false, NaN,  NaN,  NaN;
    'random_600x400_poly_R25_50_throat20','random', 'random', 600, 400, 1, 2026050701, 37.5, 20.0, 3.0, 1.0, true,  25.0, 50.0, 0.40;
    'random_600x400_poly_R30_60_throat20','random', 'random', 600, 400, 2, 2026050702, 45.0, 20.0, 3.0, 1.0, true,  30.0, 60.0, 0.38;
};

% 格式：
%   {'daCase', Da_target, D_H+ [cm^2/s], 'category'}
daCases = {
    'Da_0p01', 0.01, 3.0e-5, 'low_Da_sensitivity';
    'Da_0p03', 0.03, 3.0e-5, 'low_Da_sensitivity';
    'Da_0p1',  0.1,  3.0e-5, 'low_Da_sensitivity';
    'Da_0p3',  0.3,  3.0e-5, 'low_Da_sensitivity';
    'Da_1',    1,    3.0e-5, 'low_Da_boundary';
    'Da_2p1',  2.1,  3.0e-5, 'physical_lower_bound';
    'Da_3',    3,    3.0e-5, 'physical_core';
    'Da_5',    5,    3.0e-5, 'physical_core';
    'Da_10',   10,   3.0e-5, 'physical_core';
    'Da_20',   20,   3.0e-5, 'physical_core';
    'Da_30',   30,   3.0e-5, 'physical_core';
    'Da_50',   50,   3.0e-5, 'physical_core';
    'Da_75',   75,   3.0e-5, 'high_Da_extension';
    'Da_100',  100,  3.0e-5, 'high_Da_extension';
    'Da_150',  150,  3.0e-5, 'high_Da_extension';
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
timeControl.endTime = [];                          % adaptive_porosity 下仅作为内部安全时长自动生成
timeControl.timeStepperType = 'adaptive_porosity'; % 根据每步孔隙率增量动态调整下一步
timeControl.maxTotalTimeSteps = 1000;              % 仅作异常保护/预分配，不作为目标步数
timeControl.porosityStepTarget = 0.01;             % 目标：每个 RTM 步孔隙率约增加 1%
timeControl.porosityStepTolerance = 0.0;           % 低于目标就继续增长；上限由 porosityStepUpperFactor 控制
timeControl.porosityStepUpperFactor = 2.0;         % 允许最多约 2%/step，超过后才回调
timeControl.adaptiveGrowthFactor = 2.0;            % 未达到目标时按 0.1,0.2,0.4,0.8... 指数增长
timeControl.adaptiveShrinkSafety = 0.85;           % 超过目标时，按比例回缩并留安全余量
timeControl.adaptiveMinTimeStep = 1e-5;            % 高 Da 快反应时允许更小的起始时间步
timeControl.adaptiveMaxTimeStepCap = 60;           % [s] 只是增长上限，真正步长由孔隙率反馈决定
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
batchOptions.enableNMRSimulation = false;

% true：使用 NMR-agent 图像替代模型预测归一化弛豫曲线，再做 T2 反演。
% 与 enableNMRSimulation 互斥，二者最多只能打开一个。
batchOptions.enableNMRSurrogate = true;
batchOptions.nmrSurrogateModelPath = 'C:\Users\imgw\Documents\Codex\NMR-agent\runs\IMGW_256_300_20260507-130311_3a583275\latest_model.pt';
batchOptions.nmrSurrogateRoot = 'C:\Users\imgw\Documents\Codex\NMR-agent';
batchOptions.nmrSurrogatePythonExe = 'C:\Users\imgw\Documents\Codex\NMR-agent\.venv\Scripts\python.exe';
batchOptions.nmrSurrogateResolution = 256;
batchOptions.nmrSurrogateDevice = 'auto';

%% ===================== 由参数表构建实验列表 =====================
if ~exist(batchOutputRoot, 'dir')
    mkdir(batchOutputRoot);
end

isResumeMode = strcmpi(strtrim(string(batchRunMode)), "resume");
existingExperimentIndices = [];
lastExistingExperimentIdx = 0;

if isResumeMode
    if strlength(strtrim(string(resumeBatchResultsDir))) == 0
        error('BatchSimple:MissingResumeDir', ...
            'batchRunMode=''resume'' 时必须填写 resumeBatchResultsDir。');
    end
    batchResultsDir = char(resumeBatchResultsDir);
    if ~exist(batchResultsDir, 'dir')
        error('BatchSimple:ResumeDirNotFound', ...
            '续跑目录不存在: %s', batchResultsDir);
    end
elseif exist('batchResultsDir', 'var') && strlength(string(batchResultsDir)) > 0
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
if isResumeMode
    existingExperimentIndices = scanExistingExperimentIndices(batchResultsDir);
    if ~isempty(existingExperimentIndices)
        lastExistingExperimentIdx = max(existingExperimentIndices);
    end
end

paramList = {};
skippedParamList = {};
expCounter = 0;
skippedCounter = 0;
for iCase = 1:size(geometryCases, 1)
    row = geometryCases(iCase, :);
    targetX_um = row{4};
    targetY_um = row{5};
    targetX_cm = micronToCm(targetX_um);
    targetY_cm = micronToCm(targetY_um);
    geometryFactor = row{6};
    randomSeed = row{7};
    particleRadius_um = row{8};
    throatSpacing_um = row{9};
    minThroatRandom_um = row{10};
    randomDensityFactor = row{11};
    useRandomParticleRadii = logical(row{12});
    randomRadiusMin_um = row{13};
    randomRadiusMax_um = row{14};
    targetInitialPorosity = row{15};
    particleRadius_cm = micronToCm(particleRadius_um);
    throatSpacing_cm = micronToCm(throatSpacing_um);
    minThroatRandom_cm = micronToCm(minThroatRandom_um);
    randomRadiusMin_cm = micronToCm(randomRadiusMin_um);
    randomRadiusMax_cm = micronToCm(randomRadiusMax_um);
    characteristicLength_cm = throatSpacing_cm;
    estimatedInitialPorosity = estimateInitialPorosity( ...
        row{3}, targetX_um, targetY_um, particleRadius_um, throatSpacing_um, ...
        randomDensityFactor, useRandomParticleRadii, randomRadiusMin_um, ...
        randomRadiusMax_um, targetInitialPorosity);
    if estimatedInitialPorosity >= 0.50
        error('BatchSimple:InitialPorosityTooHigh', ...
            '%s estimated initial porosity is %.3f; adjust particle radius/spacing below 0.50.', ...
            row{1}, estimatedInitialPorosity);
    end

    for iDa = 1:size(daCases, 1)
        daRow = daCases(iDa, :);
        daCase = daRow{1};
        Da_target = daRow{2};
        diffusionCoefficient = daRow{3};
        daCategory = daRow{4};
        rateCoefficientTST = rateCoefficientForDa( ...
            Da_target, characteristicLength_cm, diffusionCoefficient, physics, fixedDesign.c_in);

        for iPe = 1:numel(peList)
            Pe_target = peList(iPe);

            p = struct();
            p.geometryCase = row{1};
            p.shapeFamily = row{2};
            p.layoutType = row{3};
            p.geometryFactor = geometryFactor;
            p.randomSeed = randomSeed;
            p.daCase = daCase;
            p.daCategory = daCategory;
            p.targetLengthXAxis_um = targetX_um;
            p.targetLengthYAxis_um = targetY_um;
            p.targetLengthXAxis = targetX_cm;
            p.targetLengthYAxis = targetY_cm;
            p.targetAspectRatio = targetX_cm / targetY_cm;
            p.Pe_target = Pe_target;
            p.Da_target = Da_target;
            p.L_cm = characteristicLength_cm;
            p.particleRadius_um = particleRadius_um;
            p.particleRadius_cm = particleRadius_cm;
            p.circleRadius = particleRadius_cm;
            p.circleSpacing = throatSpacing_cm;
            p.targetAvgSpacing = throatSpacing_cm;
            p.minThroatRandom = minThroatRandom_cm;
            p.randomDensityFactor = randomDensityFactor;
            p.useRandomParticleRadii = useRandomParticleRadii;
            p.randomParticleRadiusMin = randomRadiusMin_cm;
            p.randomParticleRadiusMax = randomRadiusMax_cm;
            p.randomParticleRadiusMin_um = randomRadiusMin_um;
            p.randomParticleRadiusMax_um = randomRadiusMax_um;
            p.targetInitialPorosity = targetInitialPorosity;
            p.estimatedInitialPorosity = estimatedInitialPorosity;
            p.diffusionCoefficient = diffusionCoefficient;
            p.rateCoefficientTST = rateCoefficientTST;
            p.u_in = Pe_target * diffusionCoefficient / p.L_cm;
            p.c_in = fixedDesign.c_in;
            if p.u_in > fixedDesign.maxInletVelocity
                skippedCounter = skippedCounter + 1;
                p.skipReason = sprintf('u_in %.4g cm/s exceeds %.4g cm/s for Pe=%g, D=%g, L=%g', ...
                    p.u_in, fixedDesign.maxInletVelocity, Pe_target, diffusionCoefficient, p.L_cm);
                [p.dxfResolutionX, p.dxfResolutionY] = estimateDxfResolution( ...
                    targetX_um, targetY_um, meshExport);
                p.Time_stepmax = NaN;
                p.initialMacroscaleTimeStepSize = NaN;
                p.adaptiveMaxTimeStep = NaN;
                p.porosityStepTarget = NaN;
                p.endTime = NaN;
                p.estimatedTotalTimeSteps = NaN;
                skippedParamList{skippedCounter, 1} = p; %#ok<SAGROW>
                continue;
            end
            timeStepSettings = estimateTimeStepSettings(p, physics, timeControl);
            p.Time_stepmax = timeStepSettings.maximalStep;
            [p.dxfResolutionX, p.dxfResolutionY] = estimateDxfResolution( ...
                targetX_um, targetY_um, meshExport);

            p = mergeStructs(p, physics);
            p.diffusionCoefficient = diffusionCoefficient;
            p.rateCoefficientTST = rateCoefficientTST;
            p = mergeStructs(p, geometry);
            p = mergeStructs(p, timeControl);
            p.initialMacroscaleTimeStepSize = timeStepSettings.initialStep;
            p.adaptiveMaxTimeStep = timeStepSettings.adaptiveMaxStep;
            p.porosityStepTarget = timeStepSettings.porosityStepTarget;
            if isempty(p.endTime)
                [p.endTime, p.estimatedTotalTimeSteps] = estimateEndTimeForMaxSteps( ...
                    p.initialMacroscaleTimeStepSize, p.Time_stepmax, p.maxTotalTimeSteps);
            else
                p.estimatedTotalTimeSteps = estimateNumTimeSteps( ...
                    p.initialMacroscaleTimeStepSize, p.Time_stepmax, p.endTime);
                if p.estimatedTotalTimeSteps > p.maxTotalTimeSteps
                    error('BatchSimple:TimeStepLimitExceeded', ...
                        'Configured endTime=%g creates %d time steps, exceeding maxTotalTimeSteps=%d.', ...
                        p.endTime, p.estimatedTotalTimeSteps, p.maxTotalTimeSteps);
                end
            end
            p = mergeStructs(p, meshExport);
            p = mergeStructs(p, batchOptions);
            expCounter = expCounter + 1;
            paramList{expCounter, 1} = p; %#ok<SAGROW>
        end
    end
end

%% ===================== 打印预览 =====================
fprintf('========================================\n');
fprintf('精简批量 RTM/NMR 运行\n');
fprintf('  运行模式: %s\n', char(batchRunMode));
fprintf('  批次目录: %s\n', batchResultsDir);
fprintf('  实验数量: %d\n', numel(paramList));
fprintf('  跳过数量: %d (velocity limit)\n', numel(skippedParamList));
if isResumeMode
    fprintf('  已存在 exp_NNN 目录: %d 个', numel(existingExperimentIndices));
    if ~isempty(existingExperimentIndices)
        fprintf('，最大编号 exp_%03d，续跑将跳过这些编号', lastExistingExperimentIdx);
    end
    fprintf('\n');
end
fprintf('  sync NMR: %s\n', mat2str(batchOptions.enableNMRSimulation));
fprintf('  surrogate NMR: %s\n', mat2str(batchOptions.enableNMRSurrogate));
fprintf('  固定几何: %.0fx%.0f um | L=throat/avg spacing per geometry | c_in: %.4g mol/cm^3 | u_max: %.4g cm/s\n', ...
    fixedDesign.targetLengthX_um, fixedDesign.targetLengthY_um, ...
    fixedDesign.c_in, fixedDesign.maxInletVelocity);
fprintf('  time-step mode: %s | target dPorosity/step %.3f-%.3f | guard steps <= %d\n', ...
    timeControl.timeStepperType, timeControl.porosityStepTarget, ...
    timeControl.porosityStepTarget * timeControl.porosityStepUpperFactor, ...
    timeControl.maxTotalTimeSteps);
fprintf('========================================\n\n');

fprintf('%-5s %-20s %-8s %-8s %-8s %-9s %-10s %-8s %-9s %-9s %-8s %-10s %-8s %-10s %-8s\n', ...
    'Exp', 'Case', 'R(um)', 'L(cm)', 'DaCase', 'Pe', 'phi0_est', 'D_H+', 'k_TST', 'u_in', 'Da', 'dtmax', 'dt0', 'nSteps', 'GeomFac');
fprintf('%s\n', repmat('-', 1, 184));
for i = 1:numel(paramList)
    p = paramList{i};
    [Pe, Da] = calcPeDa(p, physics);
    fprintf('%-5d %-20s %-8.1f %-8.4g %-8s %-9.4g %-10.3f %-8.2e %-9.2e %-9.4g %-8.4g %-10.4g %-8.3g %-10d %-8g\n', ...
        i, p.geometryCase, p.particleRadius_um, p.L_cm, p.daCase, Pe, ...
        p.estimatedInitialPorosity, p.diffusionCoefficient, p.rateCoefficientTST, ...
        p.u_in, Da, p.Time_stepmax, p.initialMacroscaleTimeStepSize, ...
        p.estimatedTotalTimeSteps, p.geometryFactor);
end
fprintf('%s\n\n', repmat('-', 1, 184));

designTable = buildDesignTable(paramList, physics);
writetable(designTable, fullfile(batchResultsDir, 'batch_design_table.xlsx'));
writetable(designTable, fullfile(batchResultsDir, 'batch_design_table.csv'));
skippedDesignTable = buildSkippedDesignTable(skippedParamList, physics);
if ~isempty(skippedDesignTable)
    writetable(skippedDesignTable, fullfile(batchResultsDir, 'batch_skipped_design_table.xlsx'));
    writetable(skippedDesignTable, fullfile(batchResultsDir, 'batch_skipped_design_table.csv'));
end

%% ===================== 执行批量实验 =====================
if isResumeMode
    summary = readExistingSummaryTable(batchResultsDir);
else
    summary = table();
end
errorLogFile = fullfile(batchResultsDir, 'batch_error_log.txt');
if isResumeMode
    fid = fopen(errorLogFile, 'a');
    fprintf(fid, '\nResume batch error log - %s\n\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
else
    fid = fopen(errorLogFile, 'w');
    fprintf(fid, 'Batch error log - %s\n\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
end
fclose(fid);

for expIdx = 1:numel(paramList)
    if isResumeMode && ismember(expIdx, existingExperimentIndices)
        fprintf('跳过已存在实验 exp_%03d，保留已有结果。\n', expIdx);
        continue;
    end

    p = paramList{expIdx};
    [Pe, Da] = calcPeDa(p, physics);
    PeDa = Pe / Da;

    fprintf('\n########################################\n');
    fprintf('实验 %d/%d | %s | %s | %gx%g um | L=%.4g cm\n', ...
        expIdx, numel(paramList), p.geometryCase, p.shapeFamily, ...
        p.targetLengthXAxis_um, p.targetLengthYAxis_um, p.L_cm);
    if p.useRandomParticleRadii
        fprintf('R=%.1f-%.1f um, throat=%.1f um, estimated initial porosity=%.3f\n', ...
            p.randomParticleRadiusMin_um, p.randomParticleRadiusMax_um, ...
            p.circleSpacing * 1e4, p.estimatedInitialPorosity);
    else
        fprintf('R=%.1f um, throat=%.1f um, estimated initial porosity=%.3f\n', ...
            p.particleRadius_um, p.circleSpacing * 1e4, p.estimatedInitialPorosity);
    end
    fprintf('D_H+=%.4g cm^2/s, k_TST=%.4g, u=%.4g cm/s, c=%.2e mol/cm^3\n', ...
        p.diffusionCoefficient, p.rateCoefficientTST, p.u_in, p.c_in);
    fprintf('Pe=%.6g, Da=%.6g (%s), Time_stepmax=%.6g s, endTime=%.6g s, nSteps<=%d\n', ...
        Pe, Da, p.daCategory, p.Time_stepmax, p.endTime, p.estimatedTotalTimeSteps);
    fprintf('########################################\n');

    try
        if isfield(p, 'randomSeed') && ~isnan(p.randomSeed)
            rng(p.randomSeed);
        end
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
        p.geometryFactor, p.randomSeed, string(p.daCase), string(p.daCategory), ...
        p.particleRadius_um, p.particleRadius_cm, p.circleSpacing, p.targetAvgSpacing, ...
        p.minThroatRandom, p.randomDensityFactor, p.useRandomParticleRadii, ...
        p.randomParticleRadiusMin_um, p.randomParticleRadiusMax_um, ...
        p.targetInitialPorosity, p.estimatedInitialPorosity, ...
        p.targetLengthXAxis_um, p.targetLengthYAxis_um, actualX_um, actualY_um, ...
        p.targetAspectRatio, p.dxfResolutionX, p.dxfResolutionY, p.L_cm, ...
        p.Pe_target, p.Da_target, p.diffusionCoefficient, p.rateCoefficientTST, ...
        p.u_in, p.c_in, Pe, Da, PeDa, ...
        p.Time_stepmax, p.initialMacroscaleTimeStepSize, p.adaptiveMaxTimeStep, ...
        p.porosityStepTarget, p.endTime, p.estimatedTotalTimeSteps, ...
        finalPorosity, initialPerm, finalPerm, pbStep, pbTime, ...
        status, message, resultsDir, ...
        'VariableNames', {'ExpIdx', 'GeometryCase', 'ShapeFamily', 'SolverLayout', ...
        'GeometryFactor', 'RandomSeed', 'DaCase', 'DaCategory', ...
        'ParticleRadius_um', 'ParticleRadius_cm', 'CircleSpacing_cm', 'TargetAvgSpacing_cm', ...
        'MinThroatRandom_cm', 'RandomDensityFactor', 'UseRandomParticleRadii', ...
        'RandomParticleRadiusMin_um', 'RandomParticleRadiusMax_um', ...
        'TargetInitialPorosity', 'EstimatedInitialPorosity', ...
        'TargetX_um', 'TargetY_um', 'ActualX_um', 'ActualY_um', ...
        'TargetAspectRatio', 'DxfResolutionX', 'DxfResolutionY', 'L_cm', ...
        'Pe_target', 'Da_target', 'D_H_cm2_s', 'k_TST', 'u_in_cm_s', 'c_in_mol_cm3', ...
        'Pe', 'Da', 'Pe_Da', 'Time_stepmax_s', 'InitialTimeStep_s', ...
        'AdaptiveMaxTimeStep_s', 'PorosityStepTarget', ...
        'EndTime_s', 'EstimatedTotalTimeSteps', 'FinalPorosity', ...
        'InitialPerm_mD', 'FinalPerm_mD', 'PBTimeStep', 'PBTime_s', ...
        'Status', 'Message', 'ResultsDir'});
    summary = upsertSummaryRow(summary, newRow, expIdx);

    writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple_partial.xlsx'));
    writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple_partial.csv'));
end

%% ===================== 保存最终总结 =====================
writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple.xlsx'));
writetable(summary, fullfile(batchResultsDir, 'batch_summary_simple.csv'));
save(fullfile(batchResultsDir, 'batch_workspace_simple.mat'), ...
    'paramList', 'skippedParamList', 'summary', 'designTable', 'skippedDesignTable', ...
    'batchResultsDir', 'batchRunMode', 'isResumeMode', ...
    'existingExperimentIndices', 'lastExistingExperimentIdx', 'physics', ...
    'fixedDesign', 'peList', 'daCases', 'geometryCases', 'geometry', ...
    'timeControl', 'meshExport', 'batchOptions');

fprintf('\n========================================\n');
fprintf('批量运行完成，用时 %.2f 秒\n', toc(batchTimer));
fprintf('结果目录: %s\n', batchResultsDir);
fprintf('总结表: %s\n', fullfile(batchResultsDir, 'batch_summary_simple.xlsx'));
fprintf('========================================\n');

%% ===================== 本脚本辅助函数 =====================
function expIndices = scanExistingExperimentIndices(batchResultsDir)
% 扫描已有 exp_NNN 子目录。只用于跳过已有结果，不判断完成质量。
items = dir(fullfile(batchResultsDir, 'exp_*'));
expIndices = [];
for i = 1:numel(items)
    if ~items(i).isdir
        continue;
    end
    token = regexp(items(i).name, '^exp_(\d+)$', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    expIndices(end + 1) = str2double(token{1}); %#ok<AGROW>
end
expIndices = unique(sort(expIndices));
end

function summary = readExistingSummaryTable(batchResultsDir)
summary = table();
summaryCandidates = {
    fullfile(batchResultsDir, 'batch_summary_simple_partial.xlsx')
    fullfile(batchResultsDir, 'batch_summary_simple.xlsx')
    fullfile(batchResultsDir, 'batch_summary_simple_partial.csv')
    fullfile(batchResultsDir, 'batch_summary_simple.csv')
};
for i = 1:numel(summaryCandidates)
    summaryPath = summaryCandidates{i};
    if ~exist(summaryPath, 'file')
        continue;
    end
    try
        summary = readtable(summaryPath, 'TextType', 'string');
        fprintf('已载入已有总结表: %s (%d 行)\n', summaryPath, height(summary));
        return;
    catch ME
        warning('BatchSimple:ReadExistingSummaryFailed', ...
            '读取已有总结表失败，将只记录本次新增实验: %s\n%s', summaryPath, ME.message);
        summary = table();
        return;
    end
end
fprintf('未找到已有 batch_summary_simple，总结表将只包含本次新增实验。\n');
end

function summary = upsertSummaryRow(summary, newRow, expIdx)
if ~isempty(summary) && ismember('ExpIdx', summary.Properties.VariableNames)
    duplicateMask = summary.ExpIdx == expIdx;
    if any(duplicateMask)
        fprintf('更新已有总结记录 exp_%03d，覆盖旧路径/旧状态。\n', expIdx);
        summary(duplicateMask, :) = [];
    end
end
summary = [summary; newRow];
if ismember('ExpIdx', summary.Properties.VariableNames)
    summary = sortrows(summary, 'ExpIdx');
end
end

function p = mergeStructs(p, extra)
names = fieldnames(extra);
for i = 1:numel(names)
    name = names{i};
    p.(name) = extra.(name);
end
end

function [Pe, Da] = calcPeDa(p, physics)
D = getStructValue(p, 'diffusionCoefficient', physics.diffusionCoefficient);
kTST = getStructValue(p, 'rateCoefficientTST', physics.rateCoefficientTST);
molarVolume = getStructValue(p, 'molarVolume', physics.molarVolume);
Pe = p.u_in * p.L_cm / D;
Da = p.c_in * molarVolume * kTST * 1000 * p.L_cm / D;
end

function cm = micronToCm(micronValue)
cm = micronValue * 1e-4;
end

function kTST = rateCoefficientForDa(targetDa, L_cm, diffusionCoefficient, physics, c_in)
kTST = targetDa * diffusionCoefficient / ...
    (c_in * physics.molarVolume * 1000 * L_cm);
end

function value = getStructValue(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function [resolutionX, resolutionY] = estimateDxfResolution(targetX_um, targetY_um, meshExport)
resolutionX = max(meshExport.minDxfResolution, round(targetX_um / meshExport.dxfMicronsPerPixel));
resolutionY = max(meshExport.minDxfResolution, round(targetY_um / meshExport.dxfMicronsPerPixel));
end

function phi0 = estimateInitialPorosity(layoutType, targetX_um, targetY_um, radius_um, throat_um, ...
    randomDensityFactor, useRandomParticleRadii, randomRadiusMin_um, randomRadiusMax_um, targetInitialPorosity)
marginLeft_um = 0.5 * radius_um;

switch char(layoutType)
    case 'hex'
        stepX_um = 2 * radius_um + throat_um;
        stepY_um = sqrt(3) * radius_um + throat_um;
        actualY_um = max(1, round(targetY_um / stepY_um)) * stepY_um;
        actualX_um = marginLeft_um + max(1, round((targetX_um - marginLeft_um) / stepX_um)) * stepX_um;
        xStart_um = marginLeft_um + radius_um;
        xCenters1 = xStart_um : stepX_um : actualX_um;
        xCenters2 = (xStart_um + stepX_um / 2) : stepX_um : actualX_um;
        yCenters = 0 : stepY_um : actualY_um;
        nCircles = 0;
        for i = 1:numel(yCenters)
            if mod(i, 2) == 1
                nCircles = nCircles + numel(xCenters1);
            else
                nCircles = nCircles + numel(xCenters2);
            end
        end
        domainArea_um2 = actualX_um * actualY_um;

    case 'random'
        if useRandomParticleRadii && ~isnan(targetInitialPorosity)
            phi0 = targetInitialPorosity;
            return;
        end
        randomArea_um2 = max(eps, (targetX_um - marginLeft_um) * targetY_um);
        targetAvgDist_um = 2 * radius_um + throat_um;
        if useRandomParticleRadii
            meanRadius_um = 0.5 * (randomRadiusMin_um + randomRadiusMax_um);
            targetAvgDist_um = 2 * meanRadius_um + throat_um;
            meanSolidArea_um2 = pi * mean([randomRadiusMin_um, randomRadiusMax_um].^2);
            nCircles = floor(randomArea_um2 / (randomDensityFactor * targetAvgDist_um^2));
            solidArea_um2 = nCircles * meanSolidArea_um2;
            phi0 = max(0, min(1, 1 - solidArea_um2 / randomArea_um2));
            return;
        end
        nCircles = floor(randomArea_um2 / (randomDensityFactor * targetAvgDist_um^2));
        domainArea_um2 = randomArea_um2;

    otherwise
        error('BatchSimple:UnknownLayoutForPorosity', 'Unknown layout type: %s', layoutType);
end

solidArea_um2 = nCircles * pi * radius_um^2;
phi0 = max(0, min(1, 1 - solidArea_um2 / domainArea_um2));
end

function designTable = buildDesignTable(paramList, physics)
designTable = table();
for i = 1:numel(paramList)
    p = paramList{i};
    [Pe, Da] = calcPeDa(p, physics);
    row = table(i, string(p.geometryCase), string(p.shapeFamily), string(p.layoutType), ...
        p.geometryFactor, p.randomSeed, string(p.daCase), string(p.daCategory), ...
        p.particleRadius_um, p.particleRadius_cm, p.circleSpacing, p.targetAvgSpacing, ...
        p.minThroatRandom, p.randomDensityFactor, p.useRandomParticleRadii, ...
        p.randomParticleRadiusMin_um, p.randomParticleRadiusMax_um, ...
        p.targetInitialPorosity, p.estimatedInitialPorosity, ...
        p.targetLengthXAxis_um, p.targetLengthYAxis_um, p.targetAspectRatio, ...
        p.dxfResolutionX, p.dxfResolutionY, p.L_cm, p.Pe_target, p.Da_target, ...
        p.diffusionCoefficient, p.rateCoefficientTST, p.u_in, p.c_in, Pe, Da, ...
        p.Time_stepmax, p.initialMacroscaleTimeStepSize, p.adaptiveMaxTimeStep, ...
        p.porosityStepTarget, p.endTime, p.estimatedTotalTimeSteps, ...
        'VariableNames', {'ExpIdx', 'GeometryCase', 'ShapeFamily', 'SolverLayout', ...
        'GeometryFactor', 'RandomSeed', 'DaCase', 'DaCategory', ...
        'ParticleRadius_um', 'ParticleRadius_cm', 'CircleSpacing_cm', 'TargetAvgSpacing_cm', ...
        'MinThroatRandom_cm', 'RandomDensityFactor', 'UseRandomParticleRadii', ...
        'RandomParticleRadiusMin_um', 'RandomParticleRadiusMax_um', ...
        'TargetInitialPorosity', 'EstimatedInitialPorosity', ...
        'TargetX_um', 'TargetY_um', 'TargetAspectRatio', 'DxfResolutionX', ...
        'DxfResolutionY', 'L_cm', 'Pe_target', 'Da_target', 'D_H_cm2_s', 'k_TST', 'u_in_cm_s', ...
        'c_in_mol_cm3', 'Pe', 'Da', 'Time_stepmax_s', 'InitialTimeStep_s', ...
        'AdaptiveMaxTimeStep_s', 'PorosityStepTarget', ...
        'EndTime_s', 'EstimatedTotalTimeSteps'});
    designTable = [designTable; row]; %#ok<AGROW>
end
end

function skippedTable = buildSkippedDesignTable(skippedParamList, physics)
skippedTable = table();
for i = 1:numel(skippedParamList)
    p = skippedParamList{i};
    [Pe, Da] = calcPeDa(p, physics);
    row = table(i, string(p.geometryCase), string(p.shapeFamily), string(p.layoutType), ...
        p.geometryFactor, p.randomSeed, string(p.daCase), string(p.daCategory), ...
        p.particleRadius_um, p.particleRadius_cm, p.circleSpacing, p.targetAvgSpacing, ...
        p.minThroatRandom, p.randomDensityFactor, p.useRandomParticleRadii, ...
        p.randomParticleRadiusMin_um, p.randomParticleRadiusMax_um, ...
        p.targetInitialPorosity, p.estimatedInitialPorosity, ...
        p.targetLengthXAxis_um, p.targetLengthYAxis_um, p.targetAspectRatio, ...
        p.dxfResolutionX, p.dxfResolutionY, p.L_cm, p.Pe_target, p.Da_target, ...
        p.diffusionCoefficient, p.rateCoefficientTST, p.u_in, p.c_in, Pe, Da, ...
        p.Time_stepmax, p.initialMacroscaleTimeStepSize, p.adaptiveMaxTimeStep, ...
        p.porosityStepTarget, p.endTime, p.estimatedTotalTimeSteps, string(p.skipReason), ...
        'VariableNames', {'SkippedIdx', 'GeometryCase', 'ShapeFamily', 'SolverLayout', ...
        'GeometryFactor', 'RandomSeed', 'DaCase', 'DaCategory', ...
        'ParticleRadius_um', 'ParticleRadius_cm', 'CircleSpacing_cm', 'TargetAvgSpacing_cm', ...
        'MinThroatRandom_cm', 'RandomDensityFactor', 'UseRandomParticleRadii', ...
        'RandomParticleRadiusMin_um', 'RandomParticleRadiusMax_um', ...
        'TargetInitialPorosity', 'EstimatedInitialPorosity', ...
        'TargetX_um', 'TargetY_um', 'TargetAspectRatio', 'DxfResolutionX', ...
        'DxfResolutionY', 'L_cm', 'Pe_target', 'Da_target', 'D_H_cm2_s', 'k_TST', 'u_in_cm_s', ...
        'c_in_mol_cm3', 'Pe', 'Da', 'Time_stepmax_s', 'InitialTimeStep_s', ...
        'AdaptiveMaxTimeStep_s', 'PorosityStepTarget', ...
        'EndTime_s', 'EstimatedTotalTimeSteps', ...
        'SkipReason'});
    skippedTable = [skippedTable; row]; %#ok<AGROW>
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

function [endTime, numTimeSteps] = estimateEndTimeForMaxSteps(initialStep, maximalStep, maxTotalTimeSteps)
if maxTotalTimeSteps < 2
    error('BatchSimple:InvalidTimeStepCap', 'maxTotalTimeSteps must be at least 2.');
end

timeSteps = [0, initialStep];
while numel(timeSteps) < maxTotalTimeSteps
    nextTime = min(timeSteps(end) * 2, timeSteps(end) + maximalStep);
    if nextTime <= timeSteps(end)
        nextTime = timeSteps(end) + maximalStep;
    end
    timeSteps(end + 1) = nextTime; %#ok<AGROW>
end

endTime = timeSteps(end);
numTimeSteps = numel(timeSteps);
end

function numTimeSteps = estimateNumTimeSteps(initialStep, maximalStep, endTime)
if endTime <= 0
    numTimeSteps = 1;
    return;
end

timeSteps = [0, initialStep];
while timeSteps(end) < endTime
    nextTime = min(timeSteps(end) * 2, timeSteps(end) + maximalStep);
    if nextTime <= timeSteps(end)
        nextTime = timeSteps(end) + maximalStep;
    end
    timeSteps(end + 1) = nextTime; %#ok<AGROW>
end
numTimeSteps = numel(timeSteps);
end

function settings = estimateTimeStepSettings(p, physics, timeControl)
[~, Da] = calcPeDa(p, physics);

defaultInitialStep = getStructValue(timeControl, 'initialMacroscaleTimeStepSize', 0.10);
minStep = getStructValue(timeControl, 'adaptiveMinTimeStep', 1e-5);
maxStepCap = getStructValue(timeControl, 'adaptiveMaxTimeStepCap', 60);
porosityStepTarget = getStructValue(timeControl, 'porosityStepTarget', 0.01);

% 这里只估计自适应步长的上下界，不估计最终溶解时间。
% Da 越大，第一步越小；最大步长不再跟第一步绑定，避免低 Pe/高 Da 被 dtmax 卡死。
referenceDa = 0.01;
reactionScale = min(1, referenceDa / max(Da, eps));
initialStep = min(defaultInitialStep, max(minStep, defaultInitialStep * reactionScale));

adaptiveMaxStep = max(maxStepCap, initialStep);
adaptiveMaxStep = max(adaptiveMaxStep, initialStep);

settings = struct();
settings.maximalStep = adaptiveMaxStep;
settings.initialStep = initialStep;
settings.adaptiveMaxStep = adaptiveMaxStep;
settings.porosityStepTarget = porosityStepTarget;
end
