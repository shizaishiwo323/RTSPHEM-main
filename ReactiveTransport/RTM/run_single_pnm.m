%% run_single_pnm - 单组参数运行 RTM + 可选同步 NMR 的精简入口
%
% 用法：
%   1. 只修改本文件顶部的“用户可调参数”。
%   2. 在 MATLAB 中运行本文件。
%   3. 如果 enableNMRSimulation=true，NMR 的 COMSOL/Python 路径等细节在
%      ReactiveTransport/automation/AutomationConfig.m 中设置。
%
% 输出：
%   outputs/rtm_runs/<runName>/
%     run_metadata.json
%     global_evolution.xlsx
%     dxf_pore/, dxf_solid/
%     comsol_results/        仅 enableNMRSimulation=true 时生成
%     inversion_results/     仅 enableNMRSimulation=true 时生成
%     nmr_sync_log.csv       仅 enableNMRSimulation=true 时生成

clear; clc;

%% ===================== 路径与运行命名 =====================
rtmDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(rtmDir));
addpath(rtmDir);

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
cfg.layoutType = 'hex';

% 特征长度/孔喉尺度 [cm]
cfg.characteristicLength = 0.001;

% 目标 Y 方向宽度 [cm]，以及目标长宽比 X/Y。
cfg.targetLengthYAxis = 0.04;
cfg.targetAspectRatio = 0.6 / 0.4;

% random 几何复用控制：
%   false 每次生成新随机几何
%   true  尝试读取 geometryLoadFile；为空时使用结果目录中的 random_geometry_config.mat
cfg.loadExistingGeometry = false;
cfg.geometryLoadFile = '';
cfg.geometrySaveFile = '';

% 外部 TIF 几何。通常保持 false。
cfg.useExternalGeometry = false;
cfg.tifPath = "";

%% ===================== 物理参数 =====================
% 入口流速 [cm/s]
cfg.inletVelocity = 0.01;

% 入口 H+ 浓度 [mol/cm^3]
cfg.initialHydrogenConcentration = 1e-4;

% 扩散系数 [cm^2/s]
cfg.diffusionCoefficient = 1e-5;

% 碳酸钙摩尔体积 [cm^3/mol]
cfg.molarVolume = 36.9;

% 反应速率常数 [mol/dm^2/s]
cfg.rateCoefficientTST = 1e-4;

%% ===================== 时间步与终止条件 =====================
% 初始宏观时间步 [s]
cfg.initialMacroscaleTimeStepSize = 0.10;

% 最大时间步 [s]。如果留空 []，PNM_beauty3 会根据 Pe 自动估算。
cfg.maximalStep = [];

% 结束时间 [s]。如果留空 []，PNM_beauty3 会根据流速等参数自动估算。
cfg.endTime = [];

% 当渗透率达到初始值的多少倍时，完成当前步后停止并导出最终结构。
cfg.permeabilityRatioThreshold = 1000;

%% ===================== 网格与 DXF 导出精度 =====================
% 微尺度网格分区数。越大越精细，也越慢。
cfg.numPartitionsMicroscale = 2 * 64;

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

% 是否显示调试图。批量数据生成通常保持 false。
cfg.showDebugFigures = false;

%% ===================== 同步 NMR 模拟 =====================
% false：只跑 RTM，后续再单独跑 NMR。
% true ：每次导出 pore/solid DXF 后立即运行 COMSOL NMR + T2 反演。
cfg.enableNMRSimulation = true;

% 注意：
%   NMR 的 COMSOL 模型、Python解释器、是否覆盖已有结果、是否启用COMSOL/反演等
%   在 ReactiveTransport/automation/AutomationConfig.m 中设置。

%% ===================== 开始运行 =====================
fprintf('========================================\n');
fprintf('单组 RTM/NMR 运行\n');
fprintf('  layout = %s\n', cfg.layoutType);
fprintf('  L      = %.6g cm\n', cfg.characteristicLength);
fprintf('  u_in   = %.6g cm/s\n', cfg.inletVelocity);
fprintf('  c_in   = %.6g mol/cm^3\n', cfg.initialHydrogenConcentration);
fprintf('  sync NMR = %s\n', mat2str(cfg.enableNMRSimulation));
fprintf('========================================\n\n');

result = PNM_beauty3(cfg);

fprintf('\n运行完成:\n');
fprintf('  结果目录: %s\n', result.resultsDir);
fprintf('  最终孔隙率: %.6f\n', result.finalPorosity);
fprintf('  最终渗透率: %.6f mD\n', result.finalPermeability);
