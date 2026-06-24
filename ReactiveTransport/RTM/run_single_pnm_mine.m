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
rtmDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(rtmDir));
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));

cfg = struct();

% 输出根目录。每次运行会在这里生成一个结果文件夹。
cfg.outputRoot = fullfile(projectRoot, 'outputs', 'rtm_runs');

% 运行名。留空 '' 时，PNM_beauty3 会自动生成 rtm_时间戳_layout 名称。
cfg.runName = '';

% 如需固定输出到某个绝对路径，取消下一行注释并填写；设置后会覆盖 outputRoot/runName。
% cfg.resultsDir = fullfile(projectRoot, 'outputs', 'rtm_runs', 'my_test_run');

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
% 入口流速 [cm/s]
cfg.inletVelocity = 0.001;

% 流体方向：left_to_right 或 bottom_to_top。
cfg.flowDirection = 'left_to_right';

% 入口 H+ 浓度 [mol/cm^3]
cfg.initialHydrogenConcentration = 1e-4;

% 化学反应求解器：
%   'phreeqc' 使用 phreeqc-m.dat 中的 Calcite RATES + 每单元 KINETICS；
%   'tst'     使用原有简化 TST 速率模型。
cfg.reactionModel = 'phreeqc';
% PHREEQC 运行组：
%   phreeqc_database_calcite : 真实 PHREEQC carbonate chemistry；
%   phreeqc_tst_match        : PHREEQC 框架中复刻旧 TST 速率律。
cfg.phreeqcRunGroup = 'phreeqc_database_calcite';
cfg.phreeqcDatabasePath = ResolvePhreeqcDatabasePath(rtmDir, 'phreeqc-m.dat');
cfg.phreeqcTemperatureC = 25;
cfg.phreeqcKineticsCorrectionFactor = 1;
cfg.phreeqcMaxSpecificSurfaceArea = 10;
cfg.phreeqcBadStepMax = 5000;
cfg.phreeqcMinHForPHMolL = 1e-7;
cfg.phreeqcMinActiveWaterVolumeFraction = 1e-2;
cfg.phreeqcMinActiveWaterVolumeCm3 = 0;
cfg.phreeqcReactNeutralInterfaceCells = false;
cfg.initialCalciumConcentration = 0;
cfg.initialCarbonConcentration = 0;
cfg.initialSodiumConcentration = 0;
cfg.initialChlorideConcentration = 0;
cfg.inletCalciumConcentration = cfg.initialCalciumConcentration;
cfg.inletCarbonConcentration = cfg.initialCarbonConcentration;
cfg.inletSodiumConcentration = cfg.initialSodiumConcentration;
cfg.inletChlorideConcentration = cfg.initialHydrogenConcentration;
cfg.phreeqcExportEvery = 1;

% 扩散系数 [cm^2/s]
cfg.diffusionCoefficient = 1e-5;

% 碳酸钙摩尔体积 [cm^3/mol]
cfg.molarVolume = 36.9;

% 反应速率常数 [mol/dm^2/s]，与 PNM_beauty3.m 中的 TST 公式单位保持一致。
cfg.rateCoefficientTST = 1e-4;
if strcmpi(cfg.reactionModel, 'phreeqc')
    cfg = ConfigurePhreeqcRunGroup(cfg, cfg.phreeqcRunGroup);
end

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
fprintf('========================================\n');
fprintf('单组 RTM/NMR 运行\n');
fprintf('  layout = %s\n', cfg.layoutType);
fprintf('  L      = %.6g cm\n', cfg.characteristicLength);
fprintf('  u_in   = %.6g cm/s\n', cfg.inletVelocity);
fprintf('  D      = %.6g cm^2/s\n', cfg.diffusionCoefficient);
fprintf('  c_in   = %.6g mol/cm^3\n', cfg.initialHydrogenConcentration);
fprintf('  reaction model = %s\n', cfg.reactionModel);
if strcmpi(cfg.reactionModel, 'phreeqc') && isfield(cfg, 'phreeqcRunGroup')
    fprintf('  PHREEQC group = %s\n', cfg.phreeqcRunGroup);
    fprintf('  PHREEQC rate law = %s\n', cfg.phreeqcRateLaw);
    fprintf('  PHREEQC DB = %s\n', cfg.phreeqcDatabasePath);
end
if ~isempty(cfg.targetDissolutionSlices)
    fprintf('  target dissolution slices = %d\n', cfg.targetDissolutionSlices);
end
if strcmpi(cfg.layoutType, 'random') && cfg.loadExistingGeometry
    fprintf('  random geometry = %s\n', cfg.geometryLoadFile);
end
fprintf('  sync NMR = %s\n', mat2str(cfg.enableNMRSimulation));
fprintf('  surrogate NMR = %s\n', mat2str(cfg.enableNMRSurrogate));
fprintf('  NMR method = %s\n', nmrOptions.nmr_method);
fprintf('  PNG NMR = %s (%s)\n', mat2str(cfg.enablePNGSimulation), cfg.pngNMRMethod);
fprintf('========================================\n\n');

result = PNM_beauty3(cfg);

fprintf('\n运行完成:\n');
fprintf('  结果目录: %s\n', result.resultsDir);
fprintf('  最终孔隙率: %.6f\n', result.finalPorosity);
fprintf('  最终渗透率: %.6f mD\n', result.finalPermeability);
