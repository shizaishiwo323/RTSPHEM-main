function result = run_phreeqc_se(overrides)
%% run_phreeqc_se - 导出震电正演所需的反应输运与 PHREEQC 水化学
%
% 用法：
%   result = run_phreeqc_se();
%   result = run_phreeqc_se(struct('inletVelocity', 0.01, ...
%       'acidConcentration_mol_L', 0.02, 'rateCoefficientTST', 1e-3));
% 输出：
%   outputs/rtm_runs/<runName>/
%     run_metadata.json
%     global_evolution.xlsx
%     phreeqc_results/phreeqc_species_step_*.csv
%     phreeqc_results/seismoelectric_chemistry.csv/.xlsx

if nargin < 1 || isempty(overrides)
    overrides = struct();
end
validateattributes(overrides, {'struct'}, {'scalar'});

%% ===================== 路径与运行命名 =====================
rtmDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(rtmDir));
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));

cfg = struct();

% 输出根目录。每次运行会在这里生成一个结果文件夹。
cfg.outputRoot = fullfile(projectRoot, 'outputs', 'rtm_runs');

% 新目录避免覆盖未包含早期过渡状态的旧运行。
cfg.runName = 'M14_24_bottom_left_600x400_literal_wormhole_Pe12_DaII20_200x100_early_transition_complete';

% 如需固定输出到某个绝对路径，取消下一行注释并填写；设置后会覆盖 outputRoot/runName。
% cfg.resultsDir = fullfile(projectRoot, 'outputs', 'rtm_runs', 'my_test_run');

%% ===================== 几何参数 =====================
% 布局类型：
%   'random' 随机颗粒；characteristicLength 表示目标平均孔喉
%   'hex'    六角排布；characteristicLength 表示最小孔喉
%   'square' 方形排布；characteristicLength 表示最小孔喉
cfg.layoutType = 'random';

% M14/24 600 x 400 um literal 方解石圆粒几何。
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
cfg.geometrySaveFile = "C:\Users\imgw\Documents\Codex\震电模拟\Wang_SE_experiment\generated_geometries\M14_24_bottom_left_600x400\M14_24_bottom_left_600x400_literal.mat";
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
% 分支 wormhole：Pe=12、DaII=20、DaI=1.667（L=0.001 cm，D=1e-5 cm^2/s）。
cfg.inletVelocity = 0.12; % [cm/s]

% 流体方向：left_to_right 或 bottom_to_top。
cfg.flowDirection = 'left_to_right';

% 酸浓度以 mol/L 输入，再显式换算为求解器使用的 mol/cm^3。
% 原脚本的 1e-2 mol/cm^3 等于 10 mol/L，超出本模型的稀溶液换算适用范围。
cfg.acidConcentration_mol_L = 0.100;
validateattributes(cfg.acidConcentration_mol_L, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', '<=', 0.5});
cfg.initialHydrogenConcentration = cfg.acidConcentration_mol_L * 1e-3;

% 化学反应求解器：
%   'phreeqc' 使用 PHREEQC Calcite RATES + 每单元 KINETICS；
%   'tst'     使用原有简化 TST 速率模型。
cfg.reactionModel = 'phreeqc';
% PHREEQC 运行组：
%   phreeqc_database_calcite : 真实 PHREEQC carbonate chemistry；
%   phreeqc_tst_match        : PHREEQC 框架中复刻旧 TST 速率律。
cfg.phreeqcRunGroup = 'phreeqc_database_calcite';
cfg.phreeqcDatabasePath = ResolvePhreeqcDatabasePath( ...
    rtmDir, fullfile('backup_legacy_20260626', 'phreeqc.dat'));
cfg.phreeqcTemperatureC = 25;
cfg.phreeqcKineticsCorrectionFactor = 0.0542005420054201 / 8.9e-5;
cfg.phreeqcMaxSpecificSurfaceArea = Inf;
cfg.phreeqcCalciteKineticsParameterConvention = 'legacy_m2_per_kgw';
cfg.phreeqcCalciteSurfaceAreaExponent = 1;
% PHREEQC database 界面推进：当前 database_calcite 模式按旧脚本尺度使用
% PHREEQC 原始 dk_Calcite/dt，即 normalSpeed ≈ rawRate * phreeqcLegacySpeedScaleFactor * 39.63e-3。
cfg.phreeqcInterfaceRateMode = 'per_area_from_dissolved_moles';
cfg.phreeqcLegacySpeedScaleFactor = 1;
cfg.phreeqcLegacyCalciteMolarVolumeForSpeed = 39.63e-3;
% PHREEQC 网格模式：
%   'separate'  使用 PHREEQC 专用反应聚合网格；
%   'level_set' 使用与 level-set / transport 相同的 gridHyPHM 网格。
cfg.phreeqcGridMode = 'level-set';
cfg.useSeparatePhreeqcGrid = strcmpi(cfg.phreeqcGridMode, 'separate');
% PHREEQC 专用反应聚合网格尺寸。仅 phreeqcGridMode='separate' 时生效。
% 20 um = 20e-4 cm = 0.002 cm。
cfg.phreeqcGridHmax = 20e-4;
cfg.exportPhreeqcGridPlot = true;
cfg.phreeqcBadStepMax = 5000;
cfg.phreeqcMinHForPHMolL = 1e-7;
cfg.phreeqcMinActiveWaterVolumeFraction = 1e-2;
cfg.phreeqcMinActiveWaterVolumeCm3 = 0;
cfg.phreeqcReactNeutralInterfaceCells = false;
cfg.initialCalciumConcentration = 0;
cfg.initialCarbonConcentration = 0;
cfg.initialSodiumConcentration = 0;
% 初始孔隙水不预装注入酸；HCl 只通过左侧入口进入。
cfg.initialChlorideConcentration = 0;
cfg.inletCalciumConcentration = cfg.initialCalciumConcentration;
cfg.inletCarbonConcentration = cfg.initialCarbonConcentration;
cfg.inletSodiumConcentration = cfg.initialSodiumConcentration;
cfg.inletChlorideConcentration = ...
    cfg.initialHydrogenConcentration + cfg.initialChlorideConcentration;
% 逐单元 PHREEQC 表很大；保留最初 6 步以解析酸突破前沿，之后仍按目标切片导出。
cfg.phreeqcExportEvery = Inf;
cfg.phreeqcEarlyExportSteps = 6;
cfg.phreeqcExportSeismoelectric = true;

% 扩散系数 [cm^2/s]
cfg.diffusionCoefficient = 1e-5;

% 碳酸钙摩尔体积 [cm^3/mol]
cfg.molarVolume = 36.9;

% DaII=20 对应的方解石速率；PHREEQC 修正因子同步改变实际动力学。
cfg.rateCoefficientTST = 0.0542005420054201; % [mol/dm^2/s]

%% ===================== 时间步与终止条件 =====================
% 高 Pe、快反应 wormhole 的保守初始宏观时间步 [s]
cfg.initialMacroscaleTimeStepSize = 5e-5;
% HyPHM 输运为隐式求解；硬 CFL 会把本算例锁死在约 0.01 s 小步。
% 保留浓度越界限制和自适应孔隙率时间步控制。
cfg.enableHardCflLimit = false;
cfg.hardCflLimitInitialGrid = false;
cfg.advectiveCflSafety = 0.3;
cfg.reactiveCflSafety = 0.2;
% PHREEQC 组分采用守恒隐式切单元有限体积输运，可在高 Pe 下保持非负和质量守恒。
cfg.transportSolverBackend = 'cut_cell_fv';
% 仅输运 Ca/C/Na/Cl 总量；pH、自由 H+ 和碳酸根物种均由 PHREEQC 反算。
cfg.phreeqcTransportBasis = 'conserved_totals';

% 最大时间步 [s]。如果留空 []，PNM_beauty3 会根据 Pe 自动估算。
cfg.maximalStep = [];

% 结束时间 [s]。如果留空 []，PNM_beauty3 会根据流速等参数自动估算。
cfg.endTime = [];
% cfg.endTime = 10;

% 从初始结构到完全溶解导出 100 个 RTM/DXF 过程切片。
% 留空 [] 时沿用 timeStepperType/exportEvery 的原始行为。
cfg.targetDissolutionSlices = 100;
% 目标切片逻辑会为粒体积归零的最终步保留第 100 个导出位。
cfg.targetDissolutionFinalPorosity = 1.0;

% 当渗透率达到初始值的多少倍时，完成当前步后停止并导出最终结构。
cfg.permeabilityRatioThreshold = 10000000;

%% ===================== 网格与 DXF 导出精度 =====================
% M14/24 的默认高分辨率诊断网格；正式运行前必须通过颗粒间隙分辨率检查。
cfg.numPartitionsMicroscale = 200;

% RTM/HyPHM 三角网格密度控制，优先级：
%   1) meshTargetElementSizeCm 非空时，按目标网格尺寸 [cm] 自动算 X/Y 分区数；
%   2) meshNumPartitionsX/Y 非空时，直接指定 X/Y 方向分区数；
%   3) 否则沿用 numPartitionsMicroscale。
cfg.meshTargetElementSizeCm =  [];
cfg.meshNumPartitionsX = 200;
cfg.meshNumPartitionsY = 100;

% DXF/掩膜导出的规则网格分辨率。越大 DXF 越细，也越慢。
cfg.dxfResolutionX = 512;
cfg.dxfResolutionY = 384;

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

%% ===================== 本入口不运行 NMR =====================
cfg.enableNMRSimulation = false;
cfg.enableNMRSurrogate = false;
cfg.enablePNGSimulation = false;
cfg.pngNMRMethod = 'none';
cfg.pngNMRConfig = struct();
nmrOptions = struct('nmr_method', 'none');

%% ===================== 开始运行 =====================
cfg = applyOverrides(cfg, overrides);
if strcmpi(cfg.reactionModel, 'phreeqc')
    cfg = ConfigurePhreeqcRunGroup(cfg, cfg.phreeqcRunGroup);
end

fprintf('========================================\n');
fprintf('单组 RTM/NMR 运行\n');
fprintf('  layout = %s\n', cfg.layoutType);
fprintf('  L      = %.6g cm\n', cfg.characteristicLength);
fprintf('  u_in   = %.6g cm/s\n', cfg.inletVelocity);
fprintf('  D      = %.6g cm^2/s\n', cfg.diffusionCoefficient);
fprintf('  c_in   = %.6g mol/cm^3\n', cfg.initialHydrogenConcentration);
fprintf('         = %.6g mol/L HCl\n', cfg.acidConcentration_mol_L);
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
assertEarlyChemistryExport(result.seismoelectricChemistryCsv, ...
    cfg.phreeqcEarlyExportSteps);

fprintf('\n运行完成:\n');
fprintf('  结果目录: %s\n', result.resultsDir);
fprintf('  最终孔隙率: %.6f\n', result.finalPorosity);
fprintf('  最终渗透率: %.6f mD\n', result.finalPermeability);
fprintf('  震电水化学 CSV: %s\n', result.seismoelectricChemistryCsv);
fprintf('  震电水化学 XLSX: %s\n', result.seismoelectricChemistryXlsx);
end

function assertEarlyChemistryExport(csvPath, requiredSteps)
if requiredSteps <= 0
    return;
end
if exist(csvPath, 'file') ~= 2
    error('RTSPHEM:SE:MissingChemistryCsv', ...
        'PHREEQC seismoelectric chemistry CSV was not created: %s', csvPath);
end
chemistry = readtable(csvPath, 'VariableNamingRule', 'preserve');
if ~ismember('TimeStep', chemistry.Properties.VariableNames)
    error('RTSPHEM:SE:MissingTimeStepColumn', ...
        'Chemistry CSV has no TimeStep column: %s', csvPath);
end
missing = setdiff((1:requiredSteps)', chemistry.TimeStep(:));
if ~isempty(missing)
    error('RTSPHEM:SE:MissingEarlyChemistry', ...
        'Chemistry CSV is missing required early TimeStep values: %s', ...
        strjoin(string(missing), ', '));
end
end

function cfg = applyOverrides(cfg, overrides)
names = fieldnames(overrides);
for i = 1:numel(names)
    cfg.(names{i}) = overrides.(names{i});
end

if isfield(overrides, 'acidConcentration_mol_L')
    validateattributes(cfg.acidConcentration_mol_L, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive', '<=', 0.5});
    if ~isfield(overrides, 'initialHydrogenConcentration')
        cfg.initialHydrogenConcentration = cfg.acidConcentration_mol_L * 1e-3;
    end
    if ~isfield(overrides, 'inletChlorideConcentration')
        cfg.inletChlorideConcentration = ...
            cfg.initialHydrogenConcentration + cfg.initialChlorideConcentration;
    end
elseif isfield(overrides, 'initialHydrogenConcentration')
    cfg.acidConcentration_mol_L = cfg.initialHydrogenConcentration * 1e3;
    if ~isfield(overrides, 'inletChlorideConcentration')
        cfg.inletChlorideConcentration = ...
            cfg.initialHydrogenConcentration + cfg.initialChlorideConcentration;
    end
end

validateattributes(cfg.inletVelocity, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(cfg.rateCoefficientTST, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(cfg.phreeqcKineticsCorrectionFactor, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
if ~isempty(cfg.targetDissolutionSlices)
    validateattributes(cfg.targetDissolutionSlices, {'numeric'}, ...
        {'scalar', 'integer', '>=', 2});
end
end
