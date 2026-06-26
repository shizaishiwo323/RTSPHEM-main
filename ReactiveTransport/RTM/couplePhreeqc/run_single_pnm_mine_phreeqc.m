%% run_single_pnm - 单组参数运行 RTM + 可选同步 NMR 的精简入口
%
% 用法：
%   1. 只修改本文件顶部的“用户可调参数”。
%   2. 在 MATLAB 中运行本文件。
%   3. NMR 路径统一在 ReactiveTransport/RTM/NMRSimulationConfig.m 中选择：
%      none | comsol | surrogate | png_pixel_cpu | png_pixel_gpu | png_mesh。
%      COMSOL 仍复用 AutomationConfig.m 的类结构，但常用覆盖参数集中写在
%      NMRSimulationConfig.m，避免入口脚本里分散配置。
%
% 输出：
%   outputs/rtm_runs/<runName>/
%     run_metadata.json
%     global_evolution.xlsx
%     dxf_pore/, dxf_solid/
%     comsol_results/        仅 enableNMRSimulation=true 时生成
%     inversion_results/     仅 enableNMRSimulation=true 时生成
%     nmr_sync_log.csv       仅 enableNMRSimulation=true 时生成
%     png_nmr_results/      仅 enablePNGSimulation=true 时生成
%     png_nmr_sync_log.csv  仅 enablePNGSimulation=true 时生成
%     surrogate_results/     仅 enableNMRSurrogate=true 时生成
%     surrogate_inversion_results/ 仅 enableNMRSurrogate=true 时生成
%     nmr_surrogate_sync_log.csv   仅 enableNMRSurrogate=true 时生成

clear; clc;

%% ===================== 路径与运行命名 =====================
scriptDir = fileparts(mfilename('fullpath'));
rtmDir = fileparts(scriptDir);
projectRoot = fileparts(fileparts(rtmDir));
addpath(rtmDir);
addpath(scriptDir);

cfg = struct();

% 输出根目录。每次运行会在这里生成一个结果文件夹。
cfg.outputRoot = fullfile(projectRoot, 'outputs', 'rtm_runs');

% 运行名。留空 '' 时，PNM_beauty3 会自动生成 rtm_时间戳_layout 名称。
cfg.runName = '';

% 如需固定输出到某个绝对路径，取消下一行注释并填写；设置后会覆盖 outputRoot/runName。
% cfg.resultsDir = fullfile(projectRoot, 'outputs', 'rtm_runs', 'Pe10_phreeqc_limited');

%% ===================== 几何参数 =====================
% 布局类型：
%   'random' 随机颗粒；characteristicLength 表示目标平均孔喉
%   'hex'    六角排布；characteristicLength 表示最小孔喉
%   'square' 方形排布；characteristicLength 表示最小孔喉
cfg.layoutType = 'random';

% 与 Study process/PNM_beauty3.m 对齐的随机颗粒几何参数。
cfg.sizeScale = 1;
cfg.circleRadius = cfg.sizeScale * 0.003;
cfg.circleSpacing = cfg.sizeScale * 0.0005;
cfg.circleSpacingXLeft = 0.5;

% 特征长度/孔喉尺度 [cm]
cfg.characteristicLength = 0.001;
cfg.targetAvgSpacing = cfg.characteristicLength;
cfg.minThroatRandom = cfg.targetAvgSpacing / 3.0;

% 目标 Y 方向宽度 [cm]，以及目标长宽比 X/Y。
cfg.targetLengthYAxis = 0.04;
cfg.targetAspectRatio = 0.6 / 0.4;

% random 几何复用控制：
%   false 每次生成新随机几何
%   true  尝试读取 geometryLoadFile；为空时使用结果目录中的 random_geometry_config.mat
cfg.loadExistingGeometry = true;
cfg.geometrySaveFile = "C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\ReactiveTransport\Analysis\孔隙耦合参数\original_data\T2lm-single\dissolution_results-Da_0.0369_Pe_0.1000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random\random_geometry_config.mat";
cfg.geometryLoadFile = cfg.geometrySaveFile;

% 外部 TIF 几何。通常保持 false。
cfg.useExternalGeometry = false;
cfg.tifPath = "";

% 外部 DXF 几何。domain 多段线和 calcite 实体/边缘需要放在不同图层。
% 这里按 1500 DXF 单位 = 0.15 cm 换算到模拟长度。
cfg.useExternalDxfGeometry = false;
cfg.externalGeometryType = 'dxf';
cfg.externalDxfPath = "C:\Users\imgw\Documents\Codex\论文复现\pore-scale-simulation-reproduction\organized_gpu_ac3d_reproduction_20260527\data\dissolution_results-Da_40.4424_Pe_4.1640_L_0.1200_square\validation-domin.dxf";
cfg.externalDxfDomainLayerNames = {'domin', 'DOMAIN'};
cfg.externalDxfSolidLayerNames = {'calcite'};
cfg.externalDxfReferenceLength = 1500;
cfg.externalDxfReferenceLengthCm = 0.15;
% DXF 导入方向：as_is, rotate90_cw, rotate90_ccw, rotate180, flip_x, flip_y。
cfg.externalDxfImportDirection = 'rotate90_cw';

%% ===================== 物理参数 =====================
% 入口流速 [cm/s]。这里设为 Pe0.1 对照算例的流速：
% Pe = u * L / D = 0.001 * 0.001 / 1e-5 = 0.1。
cfg.inletVelocity = 0.001;

% 流体方向：left_to_right 或 bottom_to_top。
cfg.flowDirection = 'left_to_right';

% 入口 H+ 浓度 [mol/cm^3]
cfg.initialHydrogenConcentration = 1e-4;

% 化学反应求解器：PHREEQC 单矿物 calcite / HCl-NaCl 体系。
% 旧 TST 求解器仍可通过 cfg.reactionModel = 'tst' 使用。
cfg.reactionModel = 'phreeqc';
% 默认跑两组：
%   phreeqc_database_calcite : phreeqc_rates.dat 中真实 Calcite carbonate chemistry
%   external_tst_phreeqc     : 输运后的 H+ 场定义一级 TST 反应量，再交给 PHREEQC 做物种平衡
% 如只想跑单组，可改为 {'external_tst_phreeqc'} 或 {'phreeqc_database_calcite'}。
cfg.phreeqcRunGroups = {'phreeqc_database_calcite', 'external_tst_phreeqc'};
cfg.phreeqcDatabasePath = ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc_rates.dat');
cfg.phreeqcTemperatureC = 25;
cfg.phreeqcKineticsCorrectionFactor = 1;
cfg.phreeqcMaxSpecificSurfaceArea = 10;
cfg.phreeqcCalciteKineticsParameterConvention = 'phreeqc_rates_cm2_per_mol';
cfg.phreeqcCalciteSurfaceAreaExponent = 1;
cfg.phreeqcBadStepMax = 5000;
cfg.phreeqcKineticsTolerance = 1e-8;
cfg.phreeqcMinHForPHMolL = 1e-7;
cfg.phreeqcMinActiveWaterVolumeFraction = 0;
cfg.phreeqcMinActiveWaterVolumeCm3 = 0;
cfg.phreeqcReactNeutralInterfaceCells = false;
% PHREEQC 中每个 solution 用 1 kg 参考水量；回写 PNM 时再按真实单元
% water kg 缩放 dk_Calcite。
cfg.phreeqcSolutionWaterKg = 1;
cfg.phreeqcWriteSolutionWaterLine = false;
cfg.phreeqcKineticsReservoirMoles = 1;
cfg.initialCalciumConcentration = 0;
cfg.initialCarbonConcentration = 0;
cfg.initialSodiumConcentration = 0;
cfg.initialChlorideConcentration = 0;
cfg.inletCalciumConcentration = cfg.initialCalciumConcentration;
cfg.inletCarbonConcentration = cfg.initialCarbonConcentration;
cfg.inletSodiumConcentration = cfg.initialSodiumConcentration;
cfg.inletChlorideConcentration = cfg.initialHydrogenConcentration;
% PHREEQC summary 每步记录；空间分布 CSV/图按此频率导出，避免完整 128x128
% 案例生成过大的逐步组分文件。若需要每步组分分布，改回 1。
cfg.phreeqcExportEvery = 10;

% 扩散系数 [cm^2/s]
cfg.diffusionCoefficient = 1e-5;

% 碳酸钙摩尔体积 [cm^3/mol]
cfg.molarVolume = 36.9;

% 反应速率常数 [mol/dm^2/s]，与 PNM_beauty3.m 中的 TST 公式单位保持一致。
cfg.rateCoefficientTST = 1e-4;

%% ===================== 时间步与终止条件 =====================
% 初始宏观时间步 [s]
cfg.initialMacroscaleTimeStepSize = 0.10;

% 最大时间步 [s]。如果留空 []，PNM_beauty3 会根据 Pe 自动估算。
cfg.maximalStep = [];

% 结束时间 [s]。如果留空 []，PNM_beauty3 会根据流速等参数自动估算。
cfg.endTime = [];
% cfg.endTime = 10;

% 目标溶解过程切片数。设置为 100 时，会自动按溶解进度调整时间步，
% 预计从初始结构到接近完全溶解约导出 100 个 RTM/DXF 过程切片。
% 留空 [] 时沿用 timeStepperType/exportEvery 的原始行为。
cfg.targetDissolutionSlices = 100;

% 当渗透率达到初始值的多少倍时，完成当前步后停止并导出最终结构。
cfg.permeabilityRatioThreshold = 10000000;

%% ===================== 网格与 DXF 导出精度 =====================
% 微尺度网格分区数。越大越精细，也越慢。
cfg.numPartitionsMicroscale = 2 * 64;

% RTM/HyPHM 三角网格密度控制，优先级：
%   1) meshTargetElementSizeCm 非空时，按目标网格尺寸 [cm] 自动算 X/Y 分区数；
%   2) meshNumPartitionsX/Y 非空时，直接指定 X/Y 方向分区数；
%   3) 否则沿用 numPartitionsMicroscale。
cfg.meshTargetElementSizeCm =  [];
cfg.meshNumPartitionsX = [];
cfg.meshNumPartitionsY = [];

% DXF/掩膜导出的规则网格分辨率。越大 DXF 越细，也越慢。
cfg.dxfResolutionX = 200;
cfg.dxfResolutionY = 100;

%% ===================== 输出控制 =====================
% 每隔多少个 RTM 时间步导出一次结构。同步 NMR 时建议保持 1。
cfg.exportEvery = 1;

% 是否导出 pore/solid DXF。同步 NMR 会自动要求 true。
cfg.exportDXF = true;

% 是否保存主图、单独子图、固液结构图、实时总览、fig 文件、Excel、最终总结图。
cfg.saveMainPlot = true;
cfg.saveIndividualPlots = true;
cfg.saveInterfaceMask = true;
cfg.saveRealtimePlot = false;
cfg.saveFigureFiles = false;
cfg.writeExcel = true;
cfg.saveFinalPlot = true;
cfg.saveMeshDiagnostics = true;

% 是否显示调试图。批量数据生成通常保持 false。
cfg.showDebugFigures = false;

%% ===================== 同步 NMR 模拟 =====================
% 统一由 ReactiveTransport/RTM/NMRSimulationConfig.m 控制。
% 在该配置文件中修改 cfg.nmr_method：
%   none | comsol | surrogate | png_pixel_cpu | png_pixel_gpu | png_mesh
nmrWorkflowConfig = NMRSimulationConfig();
% 本脚本用于 PHREEQC 化学耦合与 Pe0p1 RTM 基线对比，默认不在同一轮
% 同步跑 NMR，避免把 NMR 后处理耗时/失败混入化学求解调试。
% 如需同步 NMR，把此项改为 true 后会继续遵循 NMRSimulationConfig.m。
cfg.syncNMRForThisRun = false;
if ~cfg.syncNMRForThisRun
    nmrWorkflowConfig.nmr_method = 'none';
    nmrWorkflowConfig.enableNMRSimulation = false;
    nmrWorkflowConfig.enableNMRSurrogate = false;
    nmrWorkflowConfig.enablePNGSimulation = false;
end
nmrOptions = ResolveNMRSimulationOptions(nmrWorkflowConfig);
cfg.enableNMRSimulation = nmrOptions.enableNMRSimulation;
cfg.enableNMRSurrogate = nmrOptions.enableNMRSurrogate;
cfg.enablePNGSimulation = nmrOptions.enablePNGSimulation;
cfg.pngNMRMethod = nmrOptions.pngNMRMethod;
cfg.pngNMRConfig = nmrOptions.config;

% 下面这些字段保留在 cfg 顶层，方便旧代码和 run_metadata.json 读取；
% 实际默认值来自 NMRSimulationConfig.m。
cfg.nmrSurrogateModelPath = nmrWorkflowConfig.nmrSurrogateModelPath;
cfg.nmrSurrogateRoot = nmrWorkflowConfig.nmrSurrogateRoot;
cfg.nmrSurrogatePythonExe = nmrWorkflowConfig.nmrSurrogatePythonExe;
cfg.nmrSurrogateDatasetPath = nmrWorkflowConfig.nmrSurrogateDatasetPath;
cfg.nmrSurrogateResolution = nmrWorkflowConfig.nmrSurrogateResolution;
cfg.nmrSurrogateDevice = nmrWorkflowConfig.nmrSurrogateDevice;

% 注意：
%   NMR 方法选择、COMSOL 常用覆盖参数、NMR-agent surrogate 参数、
%   PNG NMR 参数和 T2 反演参数统一在 ReactiveTransport/RTM/NMRSimulationConfig.m 中设置。

%% ===================== 开始运行 =====================
baseCfg = cfg;
runStamp = datestr(now, 'yyyymmdd_HHMMSS');
runGroups = baseCfg.phreeqcRunGroups;
if ischar(runGroups) || isstring(runGroups)
    runGroups = cellstr(runGroups);
end

results = struct();
for iGroup = 1:numel(runGroups)
    groupCfg = ConfigurePhreeqcRunGroup(baseCfg, runGroups{iGroup}, runStamp);
    if isfield(baseCfg, 'resultsDir') && ~isempty(baseCfg.resultsDir)
        groupCfg.resultsDir = fullfile(char(baseCfg.resultsDir), groupCfg.phreeqcRunGroup);
    end

    fprintf('========================================\n');
    fprintf('PHREEQC RTM/NMR 运行组 %d/%d\n', iGroup, numel(runGroups));
    fprintf('  group  = %s\n', groupCfg.phreeqcRunGroup);
    fprintf('  rate law = %s\n', groupCfg.phreeqcRateLaw);
    fprintf('  layout = %s\n', groupCfg.layoutType);
    fprintf('  L      = %.6g cm\n', groupCfg.characteristicLength);
    fprintf('  u_in   = %.6g cm/s\n', groupCfg.inletVelocity);
    fprintf('  D      = %.6g cm^2/s\n', groupCfg.diffusionCoefficient);
    fprintf('  c_in   = %.6g mol/cm^3\n', groupCfg.initialHydrogenConcentration);
    fprintf('  reaction model = %s\n', groupCfg.reactionModel);
    fprintf('  PHREEQC DB = %s\n', groupCfg.phreeqcDatabasePath);
    if strcmpi(groupCfg.phreeqcRateLaw, 'tst_match')
        fprintf('  TST-match k = %.6g mol/dm^2/s\n', groupCfg.phreeqcTstRateCoefficient);
    end
    if ~isempty(groupCfg.targetDissolutionSlices)
        fprintf('  target dissolution slices = %d\n', groupCfg.targetDissolutionSlices);
    end
    if strcmpi(groupCfg.layoutType, 'random') && groupCfg.loadExistingGeometry
        fprintf('  random geometry = %s\n', groupCfg.geometryLoadFile);
    end
    fprintf('  sync NMR = %s\n', mat2str(groupCfg.enableNMRSimulation));
    fprintf('  surrogate NMR = %s\n', mat2str(groupCfg.enableNMRSurrogate));
    fprintf('  NMR method = %s\n', nmrOptions.nmr_method);
    fprintf('  PNG NMR = %s (%s)\n', mat2str(groupCfg.enablePNGSimulation), groupCfg.pngNMRMethod);
    fprintf('========================================\n\n');

    result = PNM_beauty3(groupCfg);
    results.(matlab.lang.makeValidName(groupCfg.phreeqcRunGroup)) = result;

    fprintf('\n运行完成 [%s]:\n', groupCfg.phreeqcRunGroup);
    fprintf('  结果目录: %s\n', result.resultsDir);
    fprintf('  最终孔隙率: %.6f\n', result.finalPorosity);
    fprintf('  最终渗透率: %.6f mD\n', result.finalPermeability);

    try
        comparison = ComparePhreeqcRunToPe0p1(result.resultsDir);
        results.(matlab.lang.makeValidName(groupCfg.phreeqcRunGroup)).comparisonToPe0p1 = comparison;
        fprintf('  Pe0p1 对比目录: %s\n', fullfile(result.resultsDir, 'comparison_to_Pe0p1'));
    catch compareErr
        warning('RTSPHEM:Phreeqc:ComparisonFailed', ...
            'Pe0p1 comparison failed for %s: %s', groupCfg.phreeqcRunGroup, compareErr.message);
    end
end

fprintf('\n全部 PHREEQC 运行组完成:\n');
for iGroup = 1:numel(runGroups)
    fieldName = matlab.lang.makeValidName(ConfigurePhreeqcRunGroup(struct(), runGroups{iGroup}).phreeqcRunGroup);
    if isfield(results, fieldName)
        fprintf('  %s -> %s\n', fieldName, results.(fieldName).resultsDir);
    end
end
