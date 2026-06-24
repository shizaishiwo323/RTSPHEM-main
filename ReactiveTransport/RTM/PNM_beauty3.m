function result = PNM_beauty3(config)
% PNM_beauty3 - Reactive-transport pore-network dissolution simulation.
%
% Inputs:
%   config (optional struct) overrides the defaults below. Important fields:
%     outputRoot, resultsDir, runName, layoutType, characteristicLength,
%     inletVelocity, initialHydrogenConcentration, exportEvery,
%     targetDissolutionSlices,
%     saveMainPlot, saveIndividualPlots, exportDXF, saveRealtimePlot,
%     enableNMRSimulation, enableNMRSurrogate.
%
% Outputs:
%   result struct with the run directory, metadata, and final scalar metrics.
%
% Test of transport/geometry coupling against second problem in the benchmark of
%'Simulation of mineral dissolution at the pore scale with evolving
%fluid-solid interfaces: review of approaches and benchmark problem set',
%Molins et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General setting variables
if nargin < 1 || isempty(config)
    config = struct();
end

rtmDir = fileparts(mfilename('fullpath'));
reactiveRoot = fileparts(rtmDir);
projectRoot = fileparts(reactiveRoot);
setupPNMPaths(reactiveRoot, rtmDir);

global Solver
Solver = 'StandardDirect';
global EPS;
EPS = eps;
simulationTimer = tic;

%% 新增：外部几何输入设置
useExternalGeometry = logical(cfgget(config, 'useExternalGeometry', false)); % 使用外部几何
externalGeometryType = char(cfgget(config, 'externalGeometryType', 'tif'));
useExternalDxfGeometry = logical(cfgget(config, 'useExternalDxfGeometry', false)) || ...
    (useExternalGeometry && strcmpi(externalGeometryType, 'dxf'));
useExternalGeometry = useExternalGeometry || useExternalDxfGeometry;
useExternalTifGeometry = useExternalGeometry && ~useExternalDxfGeometry;
% tifPath = '/Users/wangbin/Desktop/pore copy-1-2.tif'; % TIF 路径
tifPath = cfgget(config, 'tifPath', "");
externalTifSolidValue = cfgget(config, 'externalTifSolidValue', []);
externalTifLiquidValue = cfgget(config, 'externalTifLiquidValue', []);
externalTifPixelSizeMicron = cfgget(config, 'externalTifPixelSizeMicron', []);
externalTifPixelSizeCm = cfgget(config, 'externalTifPixelSizeCm', []);
externalTifSmoothingSigmaPixels = cfgget(config, 'externalTifSmoothingSigmaPixels', 1.5);
externalTifCropTopRows = cfgget(config, 'externalTifCropTopRows', 0);
externalTifCropBottomRows = cfgget(config, 'externalTifCropBottomRows', 0);
externalTifCropLeftCols = cfgget(config, 'externalTifCropLeftCols', 0);
externalTifCropRightCols = cfgget(config, 'externalTifCropRightCols', 0);
externalDxfPath = cfgget(config, 'externalDxfPath', cfgget(config, 'dxfPath', ""));
externalDxfDomainLayerNames = cfgget(config, 'externalDxfDomainLayerNames', {'domin', 'domain', 'DOMAIN'});
externalDxfSolidLayerNames = cfgget(config, 'externalDxfSolidLayerNames', {'calcite'});
externalDxfReferenceLength = cfgget(config, 'externalDxfReferenceLength', 1500);
externalDxfReferenceLengthCm = cfgget(config, 'externalDxfReferenceLengthCm', 0.15);
externalDxfScaleCmPerUnit = cfgget(config, 'externalDxfScaleCmPerUnit', []);
externalDxfOutsideDomainIsSolid = cfgget(config, 'externalDxfOutsideDomainIsSolid', true);
externalDxfImportDirection = cfgget(config, 'externalDxfImportDirection', ...
    cfgget(config, 'dxfImportDirection', 'as_is'));
numPartitionsMicroscale = cfgget(config, 'numPartitionsMicroscale', 2 * 64); % Number of partitions
meshNumPartitionsX = cfgget(config, 'meshNumPartitionsX', []);
meshNumPartitionsY = cfgget(config, 'meshNumPartitionsY', []);
meshTargetElementSizeCm = cfgget(config, 'meshTargetElementSizeCm', []);

% 时间步长控制参数（maximalStep 将根据 Pe 自动设置）
initialMacroscaleTimeStepSize = cfgget(config, 'initialMacroscaleTimeStepSize', 0.10); % 初始宏观时间步长 [s]
timeStepperType = char(cfgget(config, 'timeStepperType', 'expmax'));
maxTotalTimeSteps = cfgget(config, 'maxTotalTimeSteps', []);
hasExplicitMaxTotalTimeSteps = isstruct(config) && isfield(config, 'maxTotalTimeSteps') && ...
    ~isempty(config.maxTotalTimeSteps);
targetDissolutionSlices = cfgget(config, 'targetDissolutionSlices', []);
targetDissolutionSliceSafetyFactor = cfgget(config, 'targetDissolutionSliceSafetyFactor', 2.0);
targetDissolutionSliceMinExtraSteps = cfgget(config, 'targetDissolutionSliceMinExtraSteps', 20);
porosityStepTarget = cfgget(config, 'porosityStepTarget', 0.01);
porosityStepTolerance = cfgget(config, 'porosityStepTolerance', 0.0);
adaptivePorosityDefaults = ResolveAdaptivePorosityStepSettings(config);
porosityStepUpperFactor = adaptivePorosityDefaults.porosityStepUpperFactor;
adaptiveGrowthFactor = adaptivePorosityDefaults.adaptiveGrowthFactor;
adaptiveShrinkSafety = cfgget(config, 'adaptiveShrinkSafety', 0.85);
targetSliceStepBoundaryFactor = cfgget(config, 'targetSliceStepBoundaryFactor', 0.95);
adaptiveMinTimeStep = cfgget(config, 'adaptiveMinTimeStep', 1e-5);
targetSliceSettings = ResolveTargetDissolutionSlices( ...
    targetDissolutionSlices, NaN, maxTotalTimeSteps, porosityStepTarget, ...
    targetDissolutionSliceSafetyFactor, targetDissolutionSliceMinExtraSteps);
if targetSliceSettings.enabled
    timeStepperType = 'adaptive_porosity';
    maxTotalTimeSteps = targetSliceSettings.maxTotalTimeSteps;
    porosityStepTarget = targetSliceSettings.porosityStepTarget;
end
enableHardCflLimit = logical(cfgget(config, 'enableHardCflLimit', false));
hardCflLimitInitialGrid = logical(cfgget(config, 'hardCflLimitInitialGrid', false));
advectiveCflSafety = cfgget(config, 'advectiveCflSafety', 0.3);
reactiveCflSafety = cfgget(config, 'reactiveCflSafety', 0.2);
hardCflMinTimeStep = cfgget(config, 'hardCflMinTimeStep', adaptiveMinTimeStep);
enableConcentrationCflLimit = logical(cfgget(config, 'enableConcentrationCflLimit', true));
concentrationOvershootThreshold = cfgget(config, 'concentrationOvershootThreshold', 1.05);
concentrationHardOvershootThreshold = cfgget(config, 'concentrationHardOvershootThreshold', 10.0);
concentrationShrinkFactor = cfgget(config, 'concentrationShrinkFactor', 0.5);
hardCflMaxInitialSteps = cfgget(config, 'hardCflMaxInitialSteps', 50000);

% ----------------------------------
% 可调参数
% ----------------------------------
SizeMicron = cfgget(config, 'sizeScale', 1);  % cm 的缩放系数
circleRadius  = cfgget(config, 'circleRadius', SizeMicron * 0.003);   % 颗粒半径 [cm]
circleSpacing = cfgget(config, 'circleSpacing', SizeMicron * 0.0005);   % 最小孔喉间距（用于 square/hex，或 random 的下限） [cm]
% circleSpacing = SizeMicron * 0.0005; 
circleSpacingX_left = cfgget(config, 'circleSpacingXLeft', 0.5);   % X 方向左边界留白倍数

layoutType = cfgget(config, 'layoutType', 'random');   % 'square' | 'hex' | 'random'

% === 新增：目标尺寸与长宽比设置 ===
targetLengthYAxis = cfgget(config, 'targetLengthYAxis', SizeMicron * 0.04);  % 目标宽度 [cm]（Y方向）
targetAspectRatio = cfgget(config, 'targetAspectRatio', 0.6 / 0.4);          % 目标长宽比 X/Y

% === 新增：random 模式的孔喉控制参数 ===
% 1) targetAvgSpacing 控制期望的平均孔喉（用于特征长度、优化目标）
% 2) minThroatRandom 用作随机排布时的最小孔喉约束（通常取 targetAvgSpacing/3 左右）
targetAvgSpacing = cfgget(config, 'targetAvgSpacing', cfgget(config, 'characteristicLength', SizeMicron * 0.001));      % 目标平均孔喉 [cm]
minThroatRandom = cfgget(config, 'minThroatRandom', targetAvgSpacing / 3.0);   % random 模式的最小孔喉 [cm]
useRandomParticleRadii = logical(cfgget(config, 'useRandomParticleRadii', false));
randomParticleRadiusMin = cfgget(config, 'randomParticleRadiusMin', circleRadius);
randomParticleRadiusMax = cfgget(config, 'randomParticleRadiusMax', circleRadius);
targetInitialPorosity = cfgget(config, 'targetInitialPorosity', []);
if ~strcmp(layoutType, 'random')
    circleSpacing = cfgget(config, 'circleSpacing', cfgget(config, 'characteristicLength', circleSpacing));
end

% === 随机几何加载/保存选项 ===
loadExistingGeometry = cfgget(config, 'loadExistingGeometry', true);  % true=加载已有几何, false=生成新随机几何
geometrySaveFile = cfgget(config, 'geometrySaveFile', ""); % 随机几何配置保存路径（仅 random 模式有效）
geometryLoadFile = cfgget(config, 'geometryLoadFile', geometrySaveFile); % 随机几何配置加载路径
loadedRandomGeometry = [];
randomGeometryPreloadAttempted = false;
if ~useExternalGeometry && strcmp(layoutType, 'random') && loadExistingGeometry && ...
        ~isempty(strtrim(char(geometryLoadFile))) && exist(char(geometryLoadFile), 'file') == 2
    randomGeometryPreloadAttempted = true;
    geometryDefaults = struct( ...
        'circleRadius', circleRadius, ...
        'circleSpacing', circleSpacing, ...
        'targetAvgSpacing', targetAvgSpacing, ...
        'minThroatRandom', minThroatRandom, ...
        'targetLengthYAxis', targetLengthYAxis, ...
        'targetAspectRatio', targetAspectRatio, ...
        'useRandomParticleRadii', useRandomParticleRadii, ...
        'randomParticleRadiusMin', randomParticleRadiusMin, ...
        'randomParticleRadiusMax', randomParticleRadiusMax, ...
        'targetInitialPorosity', targetInitialPorosity);
    try
        [loadedRandomGeometry, loadedGeometryParams] = LoadRandomGeometryConfig(geometryLoadFile, geometryDefaults);
        circleRadius = loadedGeometryParams.circleRadius;
        circleSpacing = loadedGeometryParams.circleSpacing;
        targetAvgSpacing = loadedGeometryParams.targetAvgSpacing;
        minThroatRandom = loadedGeometryParams.minThroatRandom;
        targetLengthYAxis = loadedGeometryParams.targetLengthYAxis;
        targetAspectRatio = loadedGeometryParams.targetAspectRatio;
        useRandomParticleRadii = loadedGeometryParams.useRandomParticleRadii;
        randomParticleRadiusMin = loadedGeometryParams.randomParticleRadiusMin;
        randomParticleRadiusMax = loadedGeometryParams.randomParticleRadiusMax;
        if isfield(loadedGeometryParams, 'targetInitialPorosity')
            targetInitialPorosity = loadedGeometryParams.targetInitialPorosity;
        end
        fprintf('=== 已预加载随机几何配置，几何尺寸和颗粒参数以配置文件为准 ===\n');
        fprintf('  文件: %s\n', geometryLoadFile);
        fprintf('  lengthXAxis = %.6f cm, lengthYAxis = %.6f cm\n', ...
            loadedGeometryParams.lengthXAxis, loadedGeometryParams.lengthYAxis);
        fprintf('  circleRadius = %.6f cm\n', circleRadius);
        fprintf('============================================================\n');
    catch ME
        warning('MATLAB:GeometryLoad', '加载几何配置失败: %s', ME.message);
    end
end

% === DXF 导出分辨率控制 ===
dxfResolutionX = cfgget(config, 'dxfResolutionX', 200);   % DXF 导出 X 方向网格点数
dxfResolutionY = cfgget(config, 'dxfResolutionY', 100);   % DXF 导出 Y 方向网格点数
% 注意：分辨率越高，DXF 文件越大，导出时间越长

% === 输出控制：提高效率时可降低导出频率或关闭重型图形 ===
exportEvery = max(1, cfgget(config, 'exportEvery', 1));
if targetSliceSettings.enabled && exportEvery ~= 1
    warning('MATLAB:PNMTargetDissolutionSlicesExportEvery', ...
        'targetDissolutionSlices controls exported process slices, so exportEvery was set from %d to 1.', ...
        exportEvery);
    exportEvery = 1;
end
saveMainPlot = cfgget(config, 'saveMainPlot', true);
saveIndividualPlots = cfgget(config, 'saveIndividualPlots', true);
saveInterfaceMask = cfgget(config, 'saveInterfaceMask', true);
exportDXF = cfgget(config, 'exportDXF', true);
saveRealtimePlot = cfgget(config, 'saveRealtimePlot', true);
saveFigureFiles = cfgget(config, 'saveFigureFiles', true);
writeExcel = cfgget(config, 'writeExcel', true);
saveFinalPlot = cfgget(config, 'saveFinalPlot', true);
permeabilityRatioThreshold = cfgget(config, 'permeabilityRatioThreshold', 100);
showDebugFigures = cfgget(config, 'showDebugFigures', false);
saveMeshDiagnostics = logical(cfgget(config, 'saveMeshDiagnostics', true));
maxExperimentWallSeconds = cfgget(config, 'maxExperimentWallSeconds', Inf);
if isempty(maxExperimentWallSeconds) || ~isnumeric(maxExperimentWallSeconds)
    maxExperimentWallSeconds = Inf;
else
    maxExperimentWallSeconds = maxExperimentWallSeconds(1);
    if isnan(maxExperimentWallSeconds) || maxExperimentWallSeconds <= 0
        maxExperimentWallSeconds = Inf;
    end
end

% === 同步 NMR 模拟 ===
% NMRSimulationConfig.m 统一控制 COMSOL、NMR-agent surrogate 和 PNG NMR。
nmrWorkflowConfig = loadNMRWorkflowConfig(rtmDir, config);
nmrOptions = ResolveNMRSimulationOptions(nmrWorkflowConfig);
enableNMRSimulation = nmrOptions.enableNMRSimulation;
enableNMRSurrogate = nmrOptions.enableNMRSurrogate;
enablePNGSimulation = nmrOptions.enablePNGSimulation;
pngNMRConfig = nmrOptions.config;
pngNMRConfig.method = char(cfgget(nmrOptions, 'pngNMRMethod', cfgget(pngNMRConfig, 'method', 'png_mesh')));
if enableNMRSimulation && ~exportDXF
    warning('MATLAB:PNMNMRSync', 'enableNMRSimulation=true 需要导出 pore/solid DXF，已自动启用 exportDXF。');
    exportDXF = true;
end
if enablePNGSimulation && ~saveInterfaceMask
    warning('MATLAB:PNMPNGNMRSync', 'enablePNGSimulation=true 需要 interface PNG，已自动启用 saveInterfaceMask。');
    saveInterfaceMask = true;
end
nmrConfig = [];
nmrComsolOutputDir = '';
nmrInversionOutputDir = '';
nmrSyncLogFile = '';
nmrCalibrationFactor = [];
nmrSurrogateConfig = [];
nmrSurrogateOutputDir = '';
nmrSurrogateInversionOutputDir = '';
nmrSurrogateMaskDir = '';
nmrSurrogateSyncLogFile = '';
nmrSurrogateCalibrationFactor = [];
pngNMROutputDir = '';
pngNMRSyncLogFile = '';
pngNMRRuntimeConfigFile = '';
 
% 计算比例供参考
minSpacingRatio = circleSpacing / circleRadius;      % square/hex 的最小孔喉相对半径
avgSpacingRatio = targetAvgSpacing / circleRadius;   % random 目标平均相对半径
minRandomRatio  = minThroatRandom / circleRadius;    % random 最小孔喉相对半径

% === 根据布局类型自动计算实际 lengthXAxis 和 lengthYAxis ===
% 保证上下左右边界条件满足
if ~useExternalGeometry
    marginLeft = circleRadius * circleSpacingX_left;  % 左边入口留白
    circleRadii = [];
    
    switch layoutType
        case 'square'
            % 方形网格：步长 = 2R + spacing
            stepY = 2 * circleRadius + circleSpacing;
            stepX = stepY;  % 方形网格 X/Y 步长相同
            
            % Y方向：上下边界都是半圆，所以从 -R 开始到 lengthY + R
            % 实际颗粒中心 Y 范围：0, stepY, 2*stepY, ...
            % 为使上下边界恰好有半圆，lengthY 应为 N_y * stepY
            N_y = max(1, round((targetLengthYAxis) / stepY));
            lengthYAxis = N_y * stepY;
            
            % X方向：左边留白，右边有半圆
            % X方向颗粒中心范围：marginLeft, marginLeft+stepX, ..., 到 lengthX
            targetLengthXAxis = targetLengthYAxis * targetAspectRatio;
            N_x = max(1, round((targetLengthXAxis - marginLeft) / stepX));
            lengthXAxis = marginLeft + N_x * stepX;
            
        case 'hex'
            % 六角密堆
            stepX = 2 * circleRadius + circleSpacing;
            stepY = sqrt(3) * circleRadius + circleSpacing;
            
            % Y方向：上下边界有半圆
            N_y = max(1, round((targetLengthYAxis) / stepY));
            lengthYAxis = N_y * stepY;
            
            % X方向：左边留白，右边有半圆
            targetLengthXAxis = targetLengthYAxis * targetAspectRatio;
            N_x = max(1, round((targetLengthXAxis - marginLeft) / stepX));
            lengthXAxis = marginLeft + N_x * stepX;
            
        case 'random'
            % random 模式保持原有逻辑，使用目标尺寸
            lengthYAxis = targetLengthYAxis;
            lengthXAxis = targetLengthYAxis * targetAspectRatio;

        case 'single_circle'
            % Single centered grain benchmark geometry, e.g. Molins et al.
            % parts I/II: 0.1 cm x 0.05 cm channel with one calcite disk.
            lengthYAxis = targetLengthYAxis;
            lengthXAxis = targetLengthYAxis * targetAspectRatio;
            
        otherwise
            error('Unknown layoutType: %s', layoutType);
    end
    
    fprintf('=== 尺寸自动调整 ===\n');
    fprintf('目标宽度 (Y): %.6f cm -> 实际: %.6f cm\n', targetLengthYAxis, lengthYAxis);
    fprintf('目标长度 (X): %.6f cm -> 实际: %.6f cm\n', targetLengthYAxis * targetAspectRatio, lengthXAxis);
    fprintf('实际长宽比: %.4f (目标: %.4f)\n', lengthXAxis/lengthYAxis, targetAspectRatio);
    fprintf('====================\n');
end

% 输出几何参数供参考
fprintf('=== 几何参数设置 ===\n');
fprintf('颗粒半径: %.4f cm (%.1f μm)\n', circleRadius, circleRadius*1e4);
fprintf('最小孔喉(square/hex): %.4f cm (%.1f μm, %.2f×R)\n', circleSpacing, circleSpacing*1e4, minSpacingRatio);
fprintf('random目标平均孔喉   : %.4f cm (%.1f μm, %.2f×R)\n', targetAvgSpacing, targetAvgSpacing*1e4, avgSpacingRatio);
fprintf('random最小孔喉约束   : %.4f cm (%.1f μm, %.2f×R)\n', minThroatRandom, minThroatRandom*1e4, minRandomRatio);
fprintf('布局类型: %s\n', layoutType);
fprintf('=====================\n');
% ----------------------------------


tic; % Preprocessing

% Physical parameters
dimension = 2;
spaceScaleFactor = 1; % spaceScaleFactor [length of Y] = 1 [cm]

% Parameters in dm and s
pixelSizeMicron = cfgget(config, 'pixelSizeMicron', 1); % 每个像素对应的实际长度（微米，例：1 μm/px）
if ~isempty(externalTifPixelSizeMicron)
    pixelSizeMicron = externalTifPixelSizeMicron;
end
pixelSizeCm = pixelSizeMicron * 1e-4; % 将微米转换为厘米（1 μm = 1e-4 cm）
if ~isempty(externalTifPixelSizeCm)
    pixelSizeCm = externalTifPixelSizeCm;
    pixelSizeMicron = pixelSizeCm * 1e4;
end
externalGeometryXCenters = [];
externalGeometryYCenters = [];
externalDxfGeometry = [];

diffusionCoefficient = cfgget(config, 'diffusionCoefficient', 1e-5); % [ cm^2 s^(-1) ]
molarVolume = cfgget(config, 'molarVolume', 36.9); % [ cm^3 mol^(-1) ]
rateCoefficientTST = cfgget(config, 'rateCoefficientTST', 10^(-4)); % [ mol dm^(-2) s^(-1) ]
inletVelocity = cfgget(config, 'inletVelocity', 0.01); % [ cm s^(-1) ]
flowDirection = cfgget(config, 'flowDirection', 'left_to_right');
initialHydrogenConcentration = cfgget(config, 'initialHydrogenConcentration', 1e-4); % [ mol cm^(-3) ]
reactionModel = lower(char(cfgget(config, 'reactionModel', 'tst')));
usePhreeqc = strcmpi(reactionModel, 'phreeqc');
initialCalciumConcentration = cfgget(config, 'initialCalciumConcentration', 0); % [ mol cm^(-3) ]
initialCarbonConcentration = cfgget(config, 'initialCarbonConcentration', 0); % [ mol cm^(-3) ]
initialSodiumConcentration = cfgget(config, 'initialSodiumConcentration', 0); % [ mol cm^(-3) ]
initialChlorideConcentration = cfgget(config, 'initialChlorideConcentration', 0); % [ mol cm^(-3) ]
inletCalciumConcentration = cfgget(config, 'inletCalciumConcentration', initialCalciumConcentration); % [ mol cm^(-3) ]
inletCarbonConcentration = cfgget(config, 'inletCarbonConcentration', initialCarbonConcentration); % [ mol cm^(-3) ]
inletSodiumConcentration = cfgget(config, 'inletSodiumConcentration', initialSodiumConcentration); % [ mol cm^(-3) ]
inletChlorideConcentration = cfgget(config, 'inletChlorideConcentration', initialHydrogenConcentration + initialChlorideConcentration); % [ mol cm^(-3) ]
phreeqcTemperatureC = cfgget(config, 'phreeqcTemperatureC', 25);
phreeqcDatabasePath = cfgget(config, 'phreeqcDatabasePath', defaultPhreeqcDatabasePath());
phreeqcExportEvery = cfgget(config, 'phreeqcExportEvery', exportEvery);
phreeqcRunGroup = char(cfgget(config, 'phreeqcRunGroup', 'phreeqc_database_calcite'));
if usePhreeqc
    phreeqcRateLaw = normalizePhreeqcRateLaw( ...
        cfgget(config, 'phreeqcRateLaw', defaultPhreeqcRateLaw(phreeqcRunGroup)));
else
    phreeqcRateLaw = normalizePhreeqcRateLaw( ...
        cfgget(config, 'phreeqcRateLaw', 'database_calcite'));
end
phreeqcTransportFirstSplit = usePhreeqc && strcmpi(phreeqcRateLaw, 'tst_match');
phreeqcTstRateCoefficient = cfgget(config, 'phreeqcTstRateCoefficient', rateCoefficientTST);
phreeqcKineticsCorrectionFactor = cfgget(config, 'phreeqcKineticsCorrectionFactor', 1);
phreeqcMaxSpecificSurfaceArea = cfgget(config, 'phreeqcMaxSpecificSurfaceArea', 10);
if usePhreeqc && strcmpi(phreeqcRateLaw, 'tst_match') && ...
        ~cfgHasNonEmptyField(config, 'phreeqcMaxSpecificSurfaceArea')
    phreeqcMaxSpecificSurfaceArea = Inf;
end
phreeqcBadStepMax = cfgget(config, 'phreeqcBadStepMax', 5000);
phreeqcMinHForPHMolL = cfgget(config, 'phreeqcMinHForPHMolL', 1e-7);
phreeqcMinActiveWaterVolumeFraction = cfgget(config, 'phreeqcMinActiveWaterVolumeFraction', 0);
phreeqcMinActiveWaterVolumeCm3 = cfgget(config, 'phreeqcMinActiveWaterVolumeCm3', 0);
phreeqcReactionWaterVolumeFloorFraction = cfgget(config, 'phreeqcReactionWaterVolumeFloorFraction', 0);
phreeqcReactionWaterVolumeFloorCm3 = cfgget(config, 'phreeqcReactionWaterVolumeFloorCm3', 0);
phreeqcReactionProjectionMinCells = cfgget(config, 'phreeqcReactionProjectionMinCells', 6);
phreeqcReactionProjectionMaxCells = cfgget(config, 'phreeqcReactionProjectionMaxCells', 18);
phreeqcReactionProjectionMaxRings = cfgget(config, 'phreeqcReactionProjectionMaxRings', 3);
phreeqcReactionProjectionFlowBias = cfgget(config, 'phreeqcReactionProjectionFlowBias', 0);
phreeqcReactionProjectionMode = lower(char(cfgget(config, ...
    'phreeqcReactionProjectionMode', 'interface_cell')));
phreeqcReactionWaterVolumeMode = lower(char(cfgget(config, ...
    'phreeqcReactionWaterVolumeMode', 'cut_cell_agglomerated')));
phreeqcReactNeutralInterfaceCells = logical(cfgget(config, 'phreeqcReactNeutralInterfaceCells', false));
phreeqcSolutionWaterKg = cfgget(config, 'phreeqcSolutionWaterKg', 1);
phreeqcWriteSolutionWaterLine = logical(cfgget(config, 'phreeqcWriteSolutionWaterLine', false));
phreeqcKineticsReservoirMoles = cfgget(config, 'phreeqcKineticsReservoirMoles', 1);
phreeqcKineticsTolerance = cfgget(config, 'phreeqcKineticsTolerance', 1e-8);
phreeqcTransportMaxFactor = cfgget(config, 'phreeqcTransportMaxFactor', Inf);
phreeqcTransportSpeciesNames = {'H', 'Ca_total', 'C_total', 'Na_total', 'Cl_total'};
phreeqcTransportUpperBounds = [ ...
    max(initialHydrogenConcentration, 0) * phreeqcTransportMaxFactor; ...
    Inf; ...
    Inf; ...
    max([initialSodiumConcentration, inletSodiumConcentration, 0]) * phreeqcTransportMaxFactor; ...
    max([initialChlorideConcentration, inletChlorideConcentration, 0]) * phreeqcTransportMaxFactor];

dissolutionReactionRate = @(cHydrogen) ...
    (cHydrogen * 1000 * rateCoefficientTST);


% ========== 新增：计算雷诺数 ==========
kinematicViscosityWater = 0.01;  % [cm²/s]
rhoWater = 1;                                  % 密度 [g/cm^3]（cgs）
mu = rhoWater * kinematicViscosityWater;   
thickness = cfgget(config, 'thicknessCm', cfgget(config, 'partIIIHeightCm', 1));

if strcmp(layoutType, 'random')
    defaultCharacteristicLength = targetAvgSpacing; % [cm]
else
    defaultCharacteristicLength = circleSpacing; % [cm]
end
characteristicLength = cfgget(config, 'characteristicLength', defaultCharacteristicLength); % [cm]

reynoldsNumber = inletVelocity * characteristicLength / kinematicViscosityWater;

fprintf('================================================\n');
fprintf('Rate coefficient (TST)          : %.4e mol dm^(-2) s^(-1)\n', rateCoefficientTST);
fprintf('Reaction model                 : %s\n', reactionModel);
if usePhreeqc
    fprintf('PHREEQC database               : %s\n', phreeqcDatabasePath);
    fprintf('PHREEQC run group              : %s\n', phreeqcRunGroup);
    fprintf('PHREEQC rate law               : %s\n', phreeqcRateLaw);
    if strcmpi(phreeqcRateLaw, 'tst_match')
        fprintf('PHREEQC TST-match coefficient  : %.4e mol dm^(-2) s^(-1)\n', ...
            phreeqcTstRateCoefficient);
    end
    fprintf('Na/Cl initial concentration    : %.4e / %.4e mol cm^(-3)\n', ...
        initialSodiumConcentration, initialChlorideConcentration);
end
fprintf('Feature velocity (inlet)       : %.4f cm/s\n', inletVelocity);
fprintf('Characteristic length (grain)  : %.4f cm\n', characteristicLength);
fprintf('Kinematic viscosity (water)    : %.4e cm²/s\n', kinematicViscosityWater);
fprintf('>>> Reynolds number (Re)       : %.6f\n', reynoldsNumber);
if reynoldsNumber < 1
    fprintf('Flow regime: Creeping flow (Stokes flow valid)\n');
else
    fprintf('Flow regime: Inertial effects may appear\n');
end

fprintf('dissolutionReactionRate: %s\n', dissolutionReactionRate(initialHydrogenConcentration));
fprintf('================================================\n');    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation and discretization parameters (comment for convergence tests)
% 读取外部几何并构造初始水平集（signed distance，固体为正，孔隙为负）
if useExternalDxfGeometry
    dxfOptions = struct();
    dxfOptions.domainLayerNames = externalDxfDomainLayerNames;
    dxfOptions.solidLayerNames = externalDxfSolidLayerNames;
    dxfOptions.referenceLengthDxf = externalDxfReferenceLength;
    dxfOptions.referenceLengthCm = externalDxfReferenceLengthCm;
    dxfOptions.scaleCmPerDxfUnit = externalDxfScaleCmPerUnit;
    dxfOptions.resolutionX = dxfResolutionX;
    dxfOptions.resolutionY = dxfResolutionY;
    dxfOptions.outsideDomainIsSolid = externalDxfOutsideDomainIsSolid;
    dxfOptions.smoothingSigmaPixels = cfgget(config, 'externalDxfSmoothingSigmaPixels', 1.5);
    dxfOptions.importDirection = externalDxfImportDirection;

    externalDxfGeometry = ReadLayeredDxfGeometry(char(externalDxfPath), dxfOptions);
    phi_cm = externalDxfGeometry.phiCm;
    BW_solid = externalDxfGeometry.solidMask;
    [ny_img, nx_img] = size(BW_solid);
    lengthXAxis = externalDxfGeometry.lengthXAxis; % [cm]
    lengthYAxis = externalDxfGeometry.lengthYAxis; % [cm]
    externalGeometryXCenters = externalDxfGeometry.xCentersCm;
    externalGeometryYCenters = externalDxfGeometry.yCentersCm;
    pixelSizeCm = mean([lengthXAxis / nx_img, lengthYAxis / ny_img]);

    fprintf('DXF geometry: %s\n', char(externalDxfPath));
    fprintf('DXF domain layer: %s, solid layers: %s\n', ...
        externalDxfGeometry.domainLayer, strjoin(cellstr(externalDxfGeometry.solidLayers), ', '));
    fprintf('DXF scale: %.8g cm/unit\n', externalDxfGeometry.scaleCmPerDxfUnit);
    fprintf('DXF import direction: %s\n', externalDxfGeometry.importDirection);
    fprintf('DXF raster size: %d x %d\n', nx_img, ny_img);
    fprintf('Length X Axis: %.6f cm, Length Y Axis: %.6f cm\n', lengthXAxis, lengthYAxis);
elseif useExternalTifGeometry
    Iraw = imread(tifPath);
    if size(Iraw, 3) == 3
        Iraw = rgb2gray(Iraw);
    elseif size(Iraw, 3) == 4
        Iraw = rgb2gray(Iraw(:, :, 1:3));
    elseif ndims(Iraw) > 2
        Iraw = Iraw(:, :, 1);
    end
    externalTifCropTopRows = max(0, round(externalTifCropTopRows));
    externalTifCropBottomRows = max(0, round(externalTifCropBottomRows));
    externalTifCropLeftCols = max(0, round(externalTifCropLeftCols));
    externalTifCropRightCols = max(0, round(externalTifCropRightCols));
    rawRows = size(Iraw, 1);
    rawCols = size(Iraw, 2);
    rowStart = 1 + externalTifCropTopRows;
    rowEnd = rawRows - externalTifCropBottomRows;
    colStart = 1 + externalTifCropLeftCols;
    colEnd = rawCols - externalTifCropRightCols;
    if rowStart > rowEnd || colStart > colEnd
        error('External TIF crop removes the whole image: top=%d, bottom=%d, left=%d, right=%d.', ...
            externalTifCropTopRows, externalTifCropBottomRows, ...
            externalTifCropLeftCols, externalTifCropRightCols);
    end
    if externalTifCropTopRows > 0 || externalTifCropBottomRows > 0 || ...
            externalTifCropLeftCols > 0 || externalTifCropRightCols > 0
        fprintf('Cropping external TIF: top=%d rows, bottom=%d rows, left=%d cols, right=%d cols\n', ...
            externalTifCropTopRows, externalTifCropBottomRows, ...
            externalTifCropLeftCols, externalTifCropRightCols);
        Iraw = Iraw(rowStart:rowEnd, colStart:colEnd);
    end
    % 将图像转为二值：BW_solid=true 表示固体，false 表示孔隙/液体。
    if ~isempty(externalTifSolidValue) || ~isempty(externalTifLiquidValue)
        if isempty(externalTifSolidValue) || isempty(externalTifLiquidValue)
            error('externalTifSolidValue and externalTifLiquidValue must be set together.');
        end

        BW_solid = Iraw == externalTifSolidValue;
        BW_liquid = Iraw == externalTifLiquidValue;
        unknownPixels = ~(BW_solid | BW_liquid);
        if any(unknownPixels(:))
            uniqueValues = unique(Iraw(:));
            error('External TIF contains values outside solid=%g and liquid=%g. Unique values: %s', ...
                externalTifSolidValue, externalTifLiquidValue, mat2str(double(uniqueValues(:)')));
        end
    else
        % 兼容旧数据：非零像素作为固体，0 作为孔隙。
        BW_solid = Iraw > 0;
    end
    % 计算有符号距离场（像素为单位）
    % phi > 0: solid; phi < 0: pore; phi = 0: 界面
    phi_pixels = bwdist(~BW_solid) - bwdist(BW_solid);
    sigma = externalTifSmoothingSigmaPixels; % 平滑半径，通常取 1 到 2 个像素
    phi_pixels = imgaussfilt(phi_pixels, sigma);
    % 将 y 方向从"图像上->下"转换为"几何下->上"
    phi_pixels = flipud(phi_pixels);
    % 换算单位：像素 -> 厘米
    phi_cm = double(phi_pixels) * pixelSizeCm;
    % 计算物理尺寸
    [ny_img, nx_img] = size(BW_solid);
    lengthXAxis = nx_img * pixelSizeCm; % [cm]
    lengthYAxis = ny_img * pixelSizeCm; % [cm]
    externalGeometryXCenters = (0.5:1:(nx_img-0.5)) * pixelSizeCm; % [cm]
    externalGeometryYCenters = (0.5:1:(ny_img-0.5)) * pixelSizeCm; % [cm] bottom-up after flip
    fprintf('Image size (pixels): %d x %d\n', nx_img, ny_img);
    fprintf('TIF phase values: solid=%s, liquid=%s\n', ...
        mat2str(externalTifSolidValue), mat2str(externalTifLiquidValue));
    fprintf('TIF pixel size: %.6g micron/px\n', pixelSizeMicron);
    fprintf('Length X Axis: %.4f cm, Length Y Axis: %.4f cm\n', lengthXAxis, lengthYAxis);
else
    % 如果不用外部几何，lengthXAxis 和 lengthYAxis 已在前面自动计算
    fprintf('Length X Axis: %.6f cm, Length Y Axis: %.6f cm\n', lengthXAxis, lengthYAxis);
end


numberOfSlices = 1;
disp(['numSlices = ', num2str(numberOfSlices)]);
pecletNumber = (inletVelocity * characteristicLength) / diffusionCoefficient;
fprintf(['Peclet number: ', num2str(pecletNumber), '\n']);

DamkohlerNumber =  initialHydrogenConcentration * molarVolume * rateCoefficientTST * 1000 * characteristicLength/ diffusionCoefficient; 
fprintf(['Damkohler number: ', num2str(DamkohlerNumber), '\n']);

% 根据 Pe 与 Da 自动设置最大时间步长（同步 RunBatchExperiments 逻辑）
% ---------- Time_stepmax 自动策略 ----------
% Pe-based 基础值（Pe 越小，步长上限越大）
if pecletNumber >= 10
    base_step = 5;
elseif pecletNumber >= 1
    base_step = 30;
elseif pecletNumber >= 0.1
    base_step = 120;
elseif pecletNumber >= 0.01
    base_step = 300;
else
    base_step = 900;
end

% Da-based 倍率（Da 越小，步长上限倍率越大）
if DamkohlerNumber >= 10
    da_mult = 1;
elseif DamkohlerNumber >= 1
    da_mult = 2;
elseif DamkohlerNumber >= 0.1
    da_mult = 5;
elseif DamkohlerNumber >= 0.01
    da_mult = 12;
else
    da_mult = 25;
end

% 联合放大策略：当 Pe 和 Da 同时较小时，进一步增加时间步上限
lowPeDa_step_boost = 1.0;
if pecletNumber < 0.01 && DamkohlerNumber < 0.01
    lowPeDa_step_boost = 90.0;
elseif pecletNumber < 0.01 && DamkohlerNumber < 0.1
    lowPeDa_step_boost = 60.0;
elseif pecletNumber < 0.1 && DamkohlerNumber < 0.01
    lowPeDa_step_boost = 52.0;
end

maximalStepAuto = min(4320000, max(1, round(base_step * da_mult * lowPeDa_step_boost)));
maximalStep = cfgget(config, 'maximalStep', maximalStepAuto);
adaptiveMaxTimeStep = cfgget(config, 'adaptiveMaxTimeStep', maximalStep);
initialMacroscaleTimeStepSize = min(max(initialMacroscaleTimeStepSize, adaptiveMinTimeStep), maximalStep);
fprintf('maximalStep = %g s (base: %g, da_mult: %g, boost: %g)\n', ...
    maximalStep, base_step, da_mult, lowPeDa_step_boost);

%% 创建保存结果的文件夹
outputRoot = cfgget(config, 'outputRoot', fullfile(projectRoot, 'outputs', 'rtm_runs'));
runName = cfgget(config, 'runName', '');
if isempty(runName)
    runName = sprintf('rtm_%s_%s', datestr(now, 'yyyymmdd_HHMMSS_FFF'), layoutType);
end
resultsDir = cfgget(config, 'resultsDir', fullfile(outputRoot, runName));
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% 如果是random模式且不使用外部几何，更新几何配置文件路径到结果文件夹
if ~useExternalGeometry && strcmp(layoutType, 'random')
    % 说明：
    % geometryLoadFile 保持用户配置（用于“加载已有几何”）
    % geometrySaveFile 指向本次结果目录（用于“保存新生成几何”）
    geometrySaveFile = fullfile(resultsDir, 'random_geometry_config.mat');
    if loadExistingGeometry
        fprintf('随机几何配置加载路径: %s\n', geometryLoadFile);
    end
    fprintf('随机几何配置保存路径: %s\n', geometrySaveFile);
end


numMicroscaleTimeSteps = 1;
% 根据 Pe 与 Da 自动设置总模拟时长（同步 RunBatchExperiments 逻辑）
% ---------- endTime 自动策略 ----------
% Pe-based 基础总时长（Pe 越小，总时长越大）
if pecletNumber >= 10
    base_endtime = 2700 * 10;
elseif pecletNumber >= 1
    base_endtime = 2700 * 20;
elseif pecletNumber >= 0.1
    base_endtime = 2700 * 40;
elseif pecletNumber >= 0.01
    base_endtime = 2700 * 100;
else
    base_endtime = 2700 * 300;
end

% Da-based 倍率（Da 越小，溶解越慢，总时长倍率越大）
if DamkohlerNumber >= 10
    da_end_mult = 0.8;
elseif DamkohlerNumber >= 1
    da_end_mult = 1;
elseif DamkohlerNumber >= 0.1
    da_end_mult = 2;
elseif DamkohlerNumber >= 0.01
    da_end_mult = 5;
else
    da_end_mult = 10;
end

% 联合放大策略：当 Pe 和 Da 同时较小时，进一步增加总模拟时长
lowPeDa_end_boost = 1.0;
if pecletNumber < 0.01 && DamkohlerNumber < 0.01
    lowPeDa_end_boost = 120.0;
elseif pecletNumber < 0.01 && DamkohlerNumber < 0.1
    lowPeDa_end_boost = 80.0;
elseif pecletNumber < 0.1 && DamkohlerNumber < 0.01
    lowPeDa_end_boost = 64.0;
end

endTimeAuto = min(2700 * 1000000, max(2700 * 10, round(base_endtime * da_end_mult * lowPeDa_end_boost))); % [s]
endTime = cfgget(config, 'endTime', endTimeAuto);
fprintf('endTime = %g s (base: %g, da_mult: %g, boost: %g)\n', ...
    endTime, base_endtime, da_end_mult, lowPeDa_end_boost);

%% Computation of time steps in simulation
switch (timeStepperType)
    case 'linear'
        if (mod(endTime, initialMacroscaleTimeStepSize) < EPS)
            timeSteps = 0:initialMacroscaleTimeStepSize:endTime;
        else
            timeSteps = [0:initialMacroscaleTimeStepSize:endTime, endTime];
        end
        timeStepSizeFactor = 1;
    case 'exp'
        timeStepSizeFactor = 2;
        numTimeSteps = floor(log(endTime / initialMacroscaleTimeStepSize) ...
            /log(timeStepSizeFactor));
        timeSteps = [0; timeStepSizeFactor.^(0:numTimeSteps)'] ...
            /initialMacroscaleTimeStepSize;
        if (timeSteps(end) < endTime)
            timeSteps(end+1) = endTime;
        end
    case 'expmax'
        
        timeStepSizeFactor = 2;
        timeSteps = [0, initialMacroscaleTimeStepSize];
        while timeSteps(end) < endTime
            timeSteps = [timeSteps, min(timeSteps(end) * timeStepSizeFactor, timeSteps(end) + maximalStep)];
        end
        timeSteps = [timeSteps(1:(end -1)), endTime];
        size(timeSteps); %#ok<SIZINT>
    case 'adaptive_porosity'
        if isempty(maxTotalTimeSteps)
            maxTotalTimeSteps = 120;
        end
        maxTotalTimeSteps = max(2, round(maxTotalTimeSteps));
        adaptiveMaxTimeStep = max(adaptiveMinTimeStep, adaptiveMaxTimeStep);
        initialMacroscaleTimeStepSize = min(max(initialMacroscaleTimeStepSize, adaptiveMinTimeStep), adaptiveMaxTimeStep);

        timeSteps = zeros(1, maxTotalTimeSteps + 1);
        timeSteps(1) = 0;
        nextStepSize = initialMacroscaleTimeStepSize;
        for kStep = 1:maxTotalTimeSteps
            timeSteps(kStep + 1) = min(timeSteps(kStep) + nextStepSize, endTime);
            if timeSteps(kStep + 1) >= endTime - EPS
                timeSteps = timeSteps(1:(kStep + 1));
                break;
            end
            nextStepSize = min(adaptiveMaxTimeStep, nextStepSize * adaptiveGrowthFactor);
        end
        endTime = timeSteps(end);
        fprintf(['adaptive_porosity: target dPorosity=[%.4g, %.4g], ', ...
            'dt0=%g s, dt range=[%g, %g] s, safety steps=%d\n'], ...
            porosityStepTarget, porosityStepTarget * porosityStepUpperFactor, initialMacroscaleTimeStepSize, ...
            adaptiveMinTimeStep, adaptiveMaxTimeStep, numel(timeSteps)-1);
    otherwise
        error('Time stepper type not implemented.');
end

%% Initialization of level set method variables

% 自适应网格分块。可用 meshTargetElementSizeCm 或 meshNumPartitionsX/Y 显式控制。
[nxParts, nyParts, meshPartitionMode] = ResolveMeshPartitions(lengthXAxis, lengthYAxis, ...
    numPartitionsMicroscale, meshNumPartitionsX, meshNumPartitionsY, meshTargetElementSizeCm);
nxNodes = nxParts + 1;
nyNodes = nyParts + 1;
fprintf('Mesh partitions: nx=%d, ny=%d (%s)\n', nxParts, nyParts, meshPartitionMode);

microscaleGrid = FoldedCartesianGrid(dimension, ...
    [0, lengthXAxis, 0, lengthYAxis], ...
    [nxParts, nyParts]);

coord = microscaleGrid.coordinates;
coordCell = mat2cell(coord, ones(1, microscaleGrid.nodes), dimension);
meshStatsForCfl = ComputeMeshDiagnostics(coord, microscaleGrid.triangles);
meshMinLengthForCfl = meshStatsForCfl.edge_length_min;
inletCflMaxStep = advectiveCflSafety * meshMinLengthForCfl / max(abs(inletVelocity), eps);
hardCflInitialMaxStep = max(hardCflMinTimeStep, inletCflMaxStep);
if hardCflLimitInitialGrid
    oldNumTimeSteps = numel(timeSteps) - 1;
    projectedNumTimeSteps = countLimitedTimeGridSteps(timeSteps, hardCflInitialMaxStep);
    if projectedNumTimeSteps > hardCflMaxInitialSteps
        warning('MATLAB:PNMHardCFLInitialGrid', ...
            ['Hard CFL initial grid was skipped because it would create %d steps ', ...
            '(limit hardCflMaxInitialSteps=%d). Use a shorter endTime, larger safety/min dt, ', ...
            'or raise hardCflMaxInitialSteps if this is intentional.'], ...
            projectedNumTimeSteps, hardCflMaxInitialSteps);
    else
        timeSteps = limitTimeGridByMaxStep(timeSteps, hardCflInitialMaxStep);
        maximalStep = min(maximalStep, hardCflInitialMaxStep);
        adaptiveMaxTimeStep = min(adaptiveMaxTimeStep, hardCflInitialMaxStep);
        initialMacroscaleTimeStepSize = min(initialMacroscaleTimeStepSize, hardCflInitialMaxStep);
        fprintf(['Hard CFL initial cap: dt <= %.6g s ', ...
            '(h_min=%.6g cm, u_ref=%.6g cm/s, C_adv=%.3g, steps %d -> %d)\n'], ...
            hardCflInitialMaxStep, meshMinLengthForCfl, abs(inletVelocity), ...
            advectiveCflSafety, oldNumTimeSteps, numel(timeSteps)-1);
    end
else
    fprintf(['Hard CFL initial grid cap disabled: candidate dt <= %.6g s ', ...
        '(enable with hardCflLimitInitialGrid=true)\n'], hardCflInitialMaxStep);
end

numTimeSlices = numel(timeSteps);
levelSetEvolutionTime = NaN(numTimeSlices, 1);
cellProblemTime = NaN(numTimeSlices, 1);

fprintf('Total number of time steps: %d\n', numTimeSlices);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 初始水平集：如果使用外部几何，则用 TIF 的有符号距离；否则用原圆形示例
if useExternalGeometry
    % 用像素中心构造插值器
    if isempty(externalGeometryXCenters) || isempty(externalGeometryYCenters)
        xCenters = (0.5:1:(nx_img-0.5)) * pixelSizeCm; % [cm]
        yCenters = (0.5:1:(ny_img-0.5)) * pixelSizeCm; % [cm] bottom-up after flip
    else
        xCenters = externalGeometryXCenters;
        yCenters = externalGeometryYCenters;
    end
    % 注意 griddedInterpolant 需要 V(i,j) 对应 x(i), y(j)，因此需要转置
    Fphi = griddedInterpolant({xCenters, yCenters}, phi_cm', 'linear', 'nearest');
    phi_at_nodes = Fphi(coord(:,1), coord(:,2));
    initialLevelSetDataCells = cell(numberOfSlices, 1);
    initialLevelSetDataCells{1} = phi_at_nodes;
else

    % 左边入口留白
    marginLeft = circleRadius * circleSpacingX_left;

    % 区域包围盒
    minX = min(coord(:,1));  maxX = max(coord(:,1));
    minY = min(coord(:,2));  maxY = max(coord(:,2));

    % ----------------------------------
    % 生成圆心（根据 layoutType）
    % 上下边界和右边界都放置截断半圆
    % 左边界保留入口间距
    % ----------------------------------
    switch layoutType
        case 'square'
            % 方形规则网格
            step = 2 * circleRadius + circleSpacing;
            
            % X方向：圆心最左边位置应使圆片左边界距离计算域左边界 >= marginLeft
            % 即圆心 x >= marginLeft + circleRadius
            xStart = marginLeft + circleRadius;
            xCenters = xStart : step : maxX;
            
            % Y方向：从 0 开始（下边界有半圆），到 maxY（上边界有半圆）
            yCenters = 0 : step : maxY;
            
            [Xc, Yc] = meshgrid(xCenters, yCenters);
            circleCenters = [Xc(:), Yc(:)];

        case 'hex'
            % 六角密堆（等边三角形）
            stepX = 2 * circleRadius + circleSpacing;
            stepY = sqrt(3) * circleRadius + circleSpacing;

            % X方向：圆心最左边位置应使圆片左边界距离计算域左边界 >= marginLeft
            % 即圆心 x >= marginLeft + circleRadius
            xStart = marginLeft + circleRadius;
            
            % 奇数行（第1、3、5...行）从 xStart 开始
            xCenters1 = xStart : stepX : maxX;
            % 偶数行（第2、4、6...行）有 stepX/2 的偏移，但也必须满足 x >= xStart
            xCenters2_raw = (xStart + stepX/2) : stepX : maxX;
            % 确保偶数行也从 xStart 以后开始（过滤掉太靠左的点）
            xCenters2 = xCenters2_raw(xCenters2_raw >= xStart);

            % Y方向：从 0 开始（下边界有半圆），到 maxY（上边界有半圆）
            yCenters = 0 : stepY : maxY;

            circleCenters = [];
            for i = 1:length(yCenters)
                if mod(i,2)==1
                    Xrow = xCenters1;
                else
                    Xrow = xCenters2;
                end
                Yrow = yCenters(i) * ones(size(Xrow));
                circleCenters = [circleCenters; [Xrow(:), Yrow(:)]];
            end

        case 'random'
            % 随机分布，精确控制平均孔喉间距（迭代优化算法）
            % random 模式：左边留白，其他三边允许圆心靠近边界（截断效果）
            randomArea = (maxX - marginLeft) * maxY;
            targetAvgDist = 2 * circleRadius + targetAvgSpacing; % 目标平均圆心间距
            
            % === 检查是否加载已有几何 ===
            geometryLoaded = false;
            if loadExistingGeometry && ~isempty(loadedRandomGeometry)
                fprintf('=== 使用已加载随机几何配置 ===\n');
                loadedData = loadedRandomGeometry;
                circleCenters = loadedData.circleCenters;
                if isfield(loadedData, 'circleRadii')
                    circleRadii = loadedData.circleRadii(:);
                else
                    circleRadii = repmat(circleRadius, size(circleCenters, 1), 1);
                end
                geometryLoaded = true;

                if isfield(loadedData, 'finalAvgSpacing')
                    savedAvgSpacing = loadedData.finalAvgSpacing;
                elseif isfield(loadedData, 'targetAvgSpacing')
                    savedAvgSpacing = loadedData.targetAvgSpacing;
                else
                    savedAvgSpacing = NaN;
                end
                spacingRadius = mean(circleRadii);
                actualAvgSpacing = calculateAverageSpacing(circleCenters, spacingRadius);
                fprintf('✓ 成功加载几何配置\n');
                fprintf('  颗粒数: %d\n', size(circleCenters,1));
                fprintf('  几何尺寸: %.6f × %.6f cm\n', lengthXAxis, lengthYAxis);
                fprintf('  保存的平均孔喉: %.4f cm\n', savedAvgSpacing);
                fprintf('  实际平均孔喉: %.4f cm (%.2f×R_mean)\n', ...
                    actualAvgSpacing, actualAvgSpacing / spacingRadius);
                fprintf('  文件: %s\n', geometryLoadFile);
            elseif loadExistingGeometry && randomGeometryPreloadAttempted
                fprintf('未使用已有随机几何，将生成新的随机几何。\n');
            elseif loadExistingGeometry
                warning('MATLAB:GeometryLoad', '未找到几何配置文件，将生成新的随机几何: %s', geometryLoadFile);
            end
            
            % === 如果未加载，则生成新的随机几何 ===
            if ~geometryLoaded
                if useRandomParticleRadii
                    fprintf('=== 多粒径随机分布 ===\n');
                    fprintf('颗粒半径范围: %.4f-%.4f cm (%.1f-%.1f μm)\n', ...
                        randomParticleRadiusMin, randomParticleRadiusMax, ...
                        randomParticleRadiusMin * 1e4, randomParticleRadiusMax * 1e4);
                    if isempty(targetInitialPorosity)
                        targetInitialPorosity = 0.35;
                    end
                    targetSolidArea = randomArea * (1 - targetInitialPorosity);
                    maxTrials = 80000;
                    circleCenters = [];
                    circleRadii = [];
                    solidArea = 0;
                    trial = 0;
                    while trial < maxTrials && solidArea < targetSolidArea
                        trial = trial + 1;
                        fillFraction = solidArea / max(targetSolidArea, eps);
                        if fillFraction < 0.75
                            radiusNow = randomParticleRadiusMin + ...
                                rand * (randomParticleRadiusMax - randomParticleRadiusMin);
                        else
                            radiusNow = randomParticleRadiusMin + ...
                                rand^2 * (randomParticleRadiusMax - randomParticleRadiusMin);
                        end
                        xMin = marginLeft + radiusNow;
                        if xMin >= maxX
                            continue;
                        end
                        cx = xMin + rand * (maxX - xMin);
                        cy = rand * maxY;
                        candidate = [cx, cy];
                        if isempty(circleCenters)
                            accept = true;
                        else
                            distances = sqrt(sum((circleCenters - candidate).^2, 2));
                            accept = all(distances >= (circleRadii + radiusNow + minThroatRandom));
                        end
                        if accept
                            circleCenters = [circleCenters; candidate]; %#ok<AGROW>
                            circleRadii = [circleRadii; radiusNow]; %#ok<AGROW>
                            solidArea = solidArea + pi * radiusNow^2;
                        end
                    end
                    finalPorosity = 1 - solidArea / randomArea;
                    fprintf('颗粒数: %d\n', size(circleCenters, 1));
                    fprintf('估算孔隙率: %.3f [目标: %.3f]\n', finalPorosity, targetInitialPorosity);
                    if finalPorosity >= 0.50
                        warning('MATLAB:RandomPorosityHigh', ...
                            '多粒径随机几何估算孔隙率 %.3f >= 0.50；可降低 targetInitialPorosity 或减小 minThroatRandom。', ...
                            finalPorosity);
                    end
                    finalAvgSpacing = calculateAverageSpacing(circleCenters, mean(circleRadii));
                    try
                        save(geometrySaveFile, 'circleCenters', 'circleRadii', 'circleRadius', ...
                             'targetAvgSpacing', 'circleSpacing', 'lengthXAxis', 'lengthYAxis', ...
                             'finalAvgSpacing', 'finalPorosity', 'randomParticleRadiusMin', ...
                             'randomParticleRadiusMax', 'targetInitialPorosity', '-v7.3');
                        fprintf('✓ 几何配置已保存至: %s\n', geometrySaveFile);
                    catch ME
                        warning('MATLAB:GeometrySave', '保存几何配置失败: %s', ME.message);
                    end
                else
                fprintf('=== 迭代优化随机分布 ===\n');
                fprintf('目标平均孔喉: %.4f cm (%.2f×R)\n', targetAvgSpacing, avgSpacingRatio);
            
            % 迭代调整参数
            tolerance = targetAvgSpacing * 0.05; % 允许5%误差
            maxIterations = 8;
            bestResult = [];
            bestError = inf;
            
            % 初始参数估算
            densityFactor = cfgget(config, 'randomDensityFactor', 1.8); % 密度调整因子
            
            for iteration = 1:maxIterations
                % 根据密度因子估算目标颗粒数
                targetN = floor(randomArea / (densityFactor * targetAvgDist^2));
                
                % 动态调整最小距离约束（保证基本不重叠）
                % 对于random，最小表面间距由 minThroatRandom 控制
                adaptiveMinDist = 2 * circleRadius + minThroatRandom; % 中心距 = 2R + 最小孔喉
                
                circleCenters = [];
                maxTrials = targetN * 150; % 增加尝试次数
                
                fprintf('迭代 %d: 目标颗粒数=%d, 密度因子=%.2f\n', iteration, targetN, densityFactor);
                
                % 生成随机分布（左边留白，上下右边可截断）
                % X方向：圆心最左边位置应使圆片左边界距离计算域左边界 >= marginLeft
                % 即圆心 x >= marginLeft + circleRadius
                xStartRandom = marginLeft + circleRadius;
                for trial = 1:maxTrials
                    % X: 从 xStartRandom 到 maxX（右边界可截断）
                    cx = xStartRandom + rand * (maxX - xStartRandom);
                    % Y: 从 0 到 maxY（上下边界可截断）
                    cy = rand * maxY;
                    candidate = [cx, cy];
                    
                    if isempty(circleCenters)
                        circleCenters = candidate;
                    else
                        distances = sqrt(sum((circleCenters - candidate).^2, 2));
                        
                        % 只需满足自适应最小距离
                        if all(distances >= adaptiveMinDist)
                            circleCenters = [circleCenters; candidate];
                            
                            % 达到目标数量时停止
                            if size(circleCenters,1) >= targetN
                                break;
                            end
                        end
                    end
                end
                
                % 评估当前结果
                if size(circleCenters,1) > 1
                    actualAvgSpacing = calculateAverageSpacing(circleCenters, circleRadius);
                    currentError = abs(actualAvgSpacing - targetAvgSpacing);
                    
                    fprintf('  实际颗粒数: %d, 平均孔喉: %.4f cm, 误差: %.4f cm\n', ...
                            size(circleCenters,1), actualAvgSpacing, currentError);
                    
                    % 记录最佳结果
                    if currentError < bestError
                        bestError = currentError;
                        bestResult = circleCenters;
                    end
                    
                    % 检查是否达到目标精度
                    if currentError < tolerance
                        fprintf('✓ 达到目标精度！\n');
                        break;
                    end
                    
                    % 根据误差调整密度因子
                    if actualAvgSpacing > targetAvgSpacing
                        % 平均间距过大，需要增加颗粒密度
                        densityFactor = densityFactor * 0.85;
                    else
                        % 平均间距过小，需要减少颗粒密度
                        densityFactor = densityFactor * 1.15;
                    end
                    
                else
                    fprintf('  生成失败，调整参数...\n');
                    densityFactor = densityFactor * 1.5; % 大幅减少密度
                end
            end
            
                % 使用最佳结果
                if ~isempty(bestResult)
                    circleCenters = bestResult;
                    finalAvgSpacing = calculateAverageSpacing(circleCenters, circleRadius);
                    finalPorosity = 1 - size(circleCenters,1) * pi * circleRadius^2 / randomArea;
                    
                    fprintf('=== 最终结果 ===\n');
                    fprintf('颗粒数: %d\n', size(circleCenters,1));
                    fprintf('平均孔喉: %.4f cm (%.2f×R) [目标: %.4f cm]\n', ...
                            finalAvgSpacing, finalAvgSpacing/circleRadius, targetAvgSpacing);
                    fprintf('相对误差: %.2f%%\n', abs(finalAvgSpacing-targetAvgSpacing)/targetAvgSpacing*100);
                    fprintf('估算孔隙率: %.3f\n', finalPorosity);
                    
                    % 分析间距分布
                    spacings = [];
                    n = size(circleCenters,1);
                    for i = 1:n
                        for j = i+1:n
                            centerDist = norm(circleCenters(i,:) - circleCenters(j,:));
                            surfaceSpacing = centerDist - 2 * circleRadius;
                            if surfaceSpacing > 0
                                spacings = [spacings; surfaceSpacing];
                            end
                        end
                    end
                    
                    if ~isempty(spacings)
                        fprintf('间距范围: [%.4f, %.4f] cm\n', min(spacings), max(spacings));
                        fprintf('间距标准差: %.4f cm\n', std(spacings));
                    end
                    
                    % === 保存几何配置到文件 ===
                    try
                        save(geometrySaveFile, 'circleCenters', 'circleRadius', ...
                             'targetAvgSpacing', 'circleSpacing', 'lengthXAxis', 'lengthYAxis', ...
                             'finalAvgSpacing', 'finalPorosity', '-v7.3');
                        fprintf('✓ 几何配置已保存至: %s\n', geometrySaveFile);
                    catch ME
                        warning('MATLAB:GeometrySave', '保存几何配置失败: %s', ME.message);
                    end
                    
                else
                    warning('所有迭代均失败，回退到规则网格');
                    step = targetAvgDist;
                    xCenters = marginLeft : step : maxX;
                    yCenters = 0 : step : maxY;
                    [Xc, Yc] = meshgrid(xCenters, yCenters);
                    circleCenters = [Xc(:), Yc(:)];
                end
                end
            end  % end if ~geometryLoaded

        case 'single_circle'
            defaultCenter = [0.5 * lengthXAxis, 0.5 * lengthYAxis];
            circleCenters = cfgget(config, 'singleCircleCenter', defaultCenter);
            circleCenters = double(circleCenters);
            if numel(circleCenters) ~= 2
                error('singleCircleCenter must be a two-element [x, y] coordinate in cm.');
            end
            circleCenters = reshape(circleCenters, 1, 2);
            circleRadii = circleRadius;

        otherwise
            error('Unknown layoutType. Use square | hex | random | single_circle.');
    end

    if isempty(circleRadii)
        circleRadii = repmat(circleRadius, size(circleCenters, 1), 1);
    end

    % ----------------------------------
    % 构造水平集 φ(x)
    % Solid inside: φ = R - dist (>0 inside)
    % union 用 max
    % ----------------------------------
    phi_at_nodes = -inf(size(coord,1), 1);

    for k = 1:size(circleCenters, 1)
        center = circleCenters(k,:);
        radiusK = circleRadii(k);
        dist = sqrt( (coord(:,1) - center(1)).^2 + (coord(:,2) - center(2)).^2 );
        phi_k = radiusK - dist;
        phi_at_nodes = max(phi_at_nodes, phi_k);
    end

    initialLevelSetDataCells = cell(numberOfSlices, 1);
    initialLevelSetDataCells{1} = phi_at_nodes;
end



% 预览初始界面（保持原有结构，但自适应网格维度）
xlin = linspace(0, lengthXAxis, nxNodes);
ylin = linspace(0, lengthYAxis, nyNodes);
[a, b] = meshgrid(xlin, ylin);
if showDebugFigures
    figure;
    contour(a, b, reshape(initialLevelSetDataCells{1}, nxNodes, nyNodes)', [0, 1]);
    axis equal
    title('Initial Level Set (Interface)'); xlabel('X [cm]'); ylabel('Y [cm]');
end

interfaceNormalVelocity = @(cHydrogen) ...
    molarVolume * dissolutionReactionRate(cHydrogen);

currentTime = 0;

currentLevelSetDataCells = cell(numberOfSlices, 1);
[currentLevelSetDataCells{:}] = deal(initialLevelSetDataCells{1});
oldLevelSetDataCells = currentLevelSetDataCells;

levelSetData = NaN(numel(initialLevelSetDataCells{1}), numTimeSlices);
levelSetData(:, 1) = currentLevelSetDataCells{1};
levelSet = cell(numberOfSlices, 1);
[levelSet{:}] = deal(levelSetData);
clear levelSetData;

% 保存初始水平集数据用于绘制参考轮廓线
initialLevelSetForContour = initialLevelSetDataCells{1};

% 固定的“初始界面轮廓线”绘制网格（用于可视化参考线，必须始终存在）
% 说明：避免在每个时间步用 scatteredInterpolant 重新插值（可能产生 NaN 导致轮廓缺失）
X0_initContour = a; % meshgrid(xlin, ylin)
Y0_initContour = b;
phi0_initContour = reshape(initialLevelSetForContour, nxNodes, nyNodes)';

% 使用矢量形式的坐标计算初始界面 0 等值线，避免 contourc 对矩阵坐标的兼容性问题
% 预先提取初始界面的 0 等值线段，后续直接重用，避免随溶解过程消失
C0_initContour = contourc(xlin, ylin, phi0_initContour, [0 0]);
initInterfaceSegments = {};
kC = 1;
while kC < size(C0_initContour, 2)
    numPts = C0_initContour(2, kC);
    initInterfaceSegments{end+1} = C0_initContour(:, kC + 1 : kC + numPts); %#ok<SAGROW>
    kC = kC + numPts + 1;
end

meshDiagnosticsDir = fullfile(resultsDir, 'mesh_diagnostics');
if saveMeshDiagnostics
    try
        WriteMeshDiagnostics(meshDiagnosticsDir, coord, microscaleGrid.triangles, ...
            initialLevelSetDataCells{1}, initInterfaceSegments, lengthXAxis, lengthYAxis);
        fprintf('Mesh diagnostics saved to: %s\n', meshDiagnosticsDir);
    catch ME
        warning('MATLAB:MeshDiagnostics', 'Failed to write mesh diagnostics: %s', ME.message);
    end
end

% assemble HyPHM grid and label outer edges
gridHyPHM = Grid(coord, microscaleGrid.triangles);

% 修复边界条件判断 - 确保索引在有效范围内
numV = size(gridHyPHM.coordV, 1);
for i = 1:gridHyPHM.numE
    vIndices = gridHyPHM.V0E(i, :);
    if any(vIndices < 1) || any(vIndices > numV)
        error('Invalid vertex indices in edge %d: [%d, %d]', i, vIndices(1), vIndices(2));
    end
    if (all(gridHyPHM.coordV(vIndices, 1) < eps))
        gridHyPHM.idE(i) = 4;
    end
    if (all(gridHyPHM.coordV(vIndices, 1) > lengthXAxis - eps))
        gridHyPHM.idE(i) = 2;
    end
    if (all(gridHyPHM.coordV(vIndices, 2) < eps))
        gridHyPHM.idE(i) = 1;
    end
    if (all(gridHyPHM.coordV(vIndices, 2) > lengthYAxis - eps))
        gridHyPHM.idE(i) = 3;
    end
end

macroCoordCell = mat2cell(gridHyPHM.baryT, ones(gridHyPHM.numT, 1), 2);
flowConfig = GetFlowBoundaryConfig(flowDirection, lengthXAxis, lengthYAxis, nxParts, nyParts);
fprintf('Flow direction: %s\n', flowConfig.direction);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of initial effective parameters

diffusionTensors = cell(numberOfSlices, 1);
[diffusionTensors{:}] = deal(NaN(4, numTimeSlices));
porosities = cell(numberOfSlices, 1);
[porosities{1:end}] = deal(NaN(numTimeSlices, 1));
clear numTimeSlices;

surfaceArea = cell(numberOfSlices, 1);

for i = 1:numberOfSlices
    [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, triangleSurfaces] ...
        = assembleCellProblem(microscaleGrid, levelSet{i}(:, 1));
    surfaceArea{i}(1) = sum(triangleSurfaces);
    
    % 计算初始孔隙体积 (PV)
    if i == 1
        InitialPoreVolume = sum(triangleVolumes); % 假设 triangleVolumes 是孔隙空间的单元体积
        fprintf('Initial Pore Volume (Area): %.6f cm^2\n', InitialPoreVolume);
    end
end

initialPorosityForTargetSlices = InitialPoreVolume / max(lengthXAxis * lengthYAxis, eps);
targetSliceSettings = ResolveTargetDissolutionSlices( ...
    targetDissolutionSlices, initialPorosityForTargetSlices, maxTotalTimeSteps, ...
    porosityStepTarget, targetDissolutionSliceSafetyFactor, targetDissolutionSliceMinExtraSteps);
if targetSliceSettings.enabled
    porosityStepTarget = targetSliceSettings.porosityStepTarget;
    maxTotalTimeSteps = targetSliceSettings.maxTotalTimeSteps;
    fprintf(['targetDissolutionSlices=%d: initial porosity=%.4f, ', ...
        'target dPorosity per exported RTM step=%.4g, safety steps=%d\n'], ...
        targetSliceSettings.targetSliceCount, initialPorosityForTargetSlices, ...
        porosityStepTarget, maxTotalTimeSteps);
end
allowAdaptiveTimeGridExtension = targetSliceSettings.enabled && ~hasExplicitMaxTotalTimeSteps;
targetExportEnabled = targetSliceSettings.enabled;
targetExportInitialPorosity = initialPorosityForTargetSlices;
targetExportLastIndex = 0;
targetExportCount = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of HyPHM variables (cont.)

flowStepper = Stepper(0:1);
transportStepper = Stepper(timeSteps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stokes
%Calculate Pressure distribution and Stokes velocity
disp(' ');
disp('Initializing Stokes problem...');

phead = Variable(gridHyPHM, flowStepper, 'pressure head', 'P1');
StokesVelocity = Variable(gridHyPHM, flowStepper, 'Stokes Velocity', 'P2P2');
phead.setdata(0, @(t, x) 0.0);
StokesVelocity.setdata(0, @(t, x) 0.0);

StokesL = StokesLEVEL(gridHyPHM, flowStepper, 'Stokes problem');
StokesL.L.setdata(levelSet{1}(:, 1))
StokesL.id2D = num2cell(flowConfig.stokesDirichletIds);
StokesL.uD.setdata(@(t, x) double(flowConfig.inletPredicate(x)) .* flowConfig.velocityVector(inletVelocity));
StokesL.F.setdata(@(t, x) 0);
StokesL.U = StokesVelocity;
StokesL.P = phead;
try
    StokesL.N.setdata(mu);
catch
    % 如果 setdata 需要 step-index 形式，尝试使用 step 0
    StokesL.N.setdata(0, mu);
end

flowStepper.next;
StokesL.computeLevel('s');
flowStepper.prev;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flow = Variable(gridHyPHM, flowStepper, 'Flow', 'RT0');
helper = StokesL.U.getdata(1);
flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));
porosityFunc = @(t, x) porosityHelperFun(t, x, porosities, numberOfSlices, 1);

% Boundary IDs: 1 = down, 2 = right, 3 = up, 4 = left

disp(' ');
disp('Initializing hydrogen transport...');

hydrogenConcentration = Variable(gridHyPHM, transportStepper, ...
    'H^+_pH', 'P0');
% hydrogenConcentration.setdata(0, @(t, x) initialHydrogenConcentration);
hydrogenConcentration.setdata(0, @(t, x) 0);
hydrogenTransport = TransportLEVEL(gridHyPHM, transportStepper, 'H^+ Transport');
hydrogenTransport.id2N = {1, 2, 3};
hydrogenTransport.id2F = {4};
hydrogenTransport.U = hydrogenConcentration;
hydrogenTransport.D.setdata(diffusionCoefficient*eye(2));
hydrogenTransport.gF.setdata( ...
    @(t, x) -initialHydrogenConcentration*inletVelocity*(x(1) < EPS));
hydrogenTransport.A.setdata(0, @(t, x) 1);
hydrogenTransport.C.setdata(0, flow.getdata(1));
hydrogenTransport.isUpwind = 'exp';
fprintf('Transport upwinding mode: %s\n', hydrogenTransport.isUpwind);

if usePhreeqc
    disp('Initializing PHREEQC transport components...');
    calciumTransport = createPhreeqcTransport(gridHyPHM, transportStepper, flow, ...
        'Ca_total', initialCalciumConcentration, inletCalciumConcentration, ...
        diffusionCoefficient, inletVelocity, EPS, 'exp');
    carbonTransport = createPhreeqcTransport(gridHyPHM, transportStepper, flow, ...
        'C_total', initialCarbonConcentration, inletCarbonConcentration, ...
        diffusionCoefficient, inletVelocity, EPS, 'exp');
    sodiumTransport = createPhreeqcTransport(gridHyPHM, transportStepper, flow, ...
        'Na_total', initialSodiumConcentration, inletSodiumConcentration, ...
        diffusionCoefficient, inletVelocity, EPS, 'exp');
    chlorideTransport = createPhreeqcTransport(gridHyPHM, transportStepper, flow, ...
        'Cl_total', initialChlorideConcentration, inletChlorideConcentration, ...
        diffusionCoefficient, inletVelocity, EPS, 'exp');
end

hydrogenDataOld = hydrogenConcentration.getdata(0);

preprocessingTime = toc; % Preprocessing

disp(' ');
disp(['Preprocessing done in ', num2str(preprocessingTime), ' seconds.']);

hydrogenTransport.L.setdata(0, [levelSet{1}(:, 1)]);
if usePhreeqc
    calciumTransport.L.setdata(0, [levelSet{1}(:, 1)]);
    carbonTransport.L.setdata(0, [levelSet{1}(:, 1)]);
    sodiumTransport.L.setdata(0, [levelSet{1}(:, 1)]);
    chlorideTransport.L.setdata(0, [levelSet{1}(:, 1)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping Vertices --> triangles

VertexTriMatrix = cell(1, gridHyPHM.numV);
for i = 1:gridHyPHM.numT
    for j = 1:3
        VertexTriMatrix{1, gridHyPHM.V0T(i, j)} = [VertexTriMatrix{1, gridHyPHM.V0T(i, j)}, i];
    end
end
Rate = nan(1, numel(timeSteps));
Volume = nan(1, numel(timeSteps));
% --- 新增：初始化孔隙率与渗透率存储数组（插入点：Volume = nan(...) 之后） ---
Porosity = nan(1, numel(timeSteps));
Permeability = nan(1, numel(timeSteps));
PV_injected = nan(1, numel(timeSteps)); % 新增：PV 数组
OutletHConc = nan(1, numel(timeSteps)); % 新增：出口处平均 H+ 浓度
Tortuosity  = nan(1, numel(timeSteps)); % 迂曲度 τ = <|v|>/<vx>（基于孔隙速度场）
numTortuositySegments = 5; % 新增：沿 X 方向等分的分段数量
TortuositySegments = nan(numTortuositySegments, numel(timeSteps)); % 新增：每段迂曲度
k0 = NaN; % 初始渗透率（用于计算 k/k0）

% ----- 日志文件 & 边界缓存（只做一次） -----
logfile = fullfile(resultsDir, 'global_evolution_log.csv');
xlsxFile = fullfile(resultsDir, 'global_evolution.xlsx'); % 新增：Excel 文件路径
tortuositySegCsvFile = fullfile(resultsDir, 'tortuosity_segments_log.csv'); % 新增：分段迂曲度 CSV
tortuositySegXlsxFile = fullfile(resultsDir, 'tortuosity_segments.xlsx');   % 新增：分段迂曲度 Excel
stabilityDiagnosticsCsvFile = fullfile(resultsDir, 'stability_diagnostics_log.csv');
if writeExcel && exist(xlsxFile, 'file'), delete(xlsxFile); end % 删除旧文件
if writeExcel && exist(tortuositySegXlsxFile, 'file'), delete(tortuositySegXlsxFile); end % 删除旧文件

if ~exist(logfile, 'file')
    fid = fopen(logfile, 'w');
    fprintf(fid, ['timestep,time_s,porosity,permeability_mD,k_k0,avg_dissolution_rate,', ...
        'interface_tst_rate_mol_cm2_s,flux_apparent_rate_mol_cm2_s,', ...
        'surface_area_cm2,grain_volume_cm3,injected_pv,outlet_H_conc,tortuosity\n']);
    fclose(fid);
end

if ~exist(tortuositySegCsvFile, 'file')
    fid = fopen(tortuositySegCsvFile, 'w');
    fprintf(fid, 'timestep,time_s,tortuosity_all,tortuosity_seg1,tortuosity_seg2,tortuosity_seg3,tortuosity_seg4,tortuosity_seg5\n');
    fclose(fid);
end

if ~exist(stabilityDiagnosticsCsvFile, 'file')
    fid = fopen(stabilityDiagnosticsCsvFile, 'w');
    fprintf(fid, ['timestep,time_s,dt_s,c_min,c_max,total_mass_mol,mass_balance_residual_mol,', ...
        'mass_balance_relative,local_u_max_cm_s,local_pe_max,advective_cfl,reaction_cfl,', ...
        'normal_speed_max_cm_s,pore_components,largest_pore_component_fraction,porosity,permeability_mD,k_k0,', ...
        'transport_upwind,diagnostic_flag\n']);
    fclose(fid);
end

% 新增：预计算分段边界（沿主流方向平均 5 段）
tortuositySegEdges = flowConfig.segmentEdges;

% 边界索引缓存（避免在每步重复计算）
edgesRightCache = flowConfig.outletEdgePredicate(gridHyPHM.baryE);
leftNodesCache  = flowConfig.inletNodePredicate(gridHyPHM.coordV);
rightNodesCache = flowConfig.outletNodePredicate(gridHyPHM.coordV);
% ----------------------------------------------



solidDXFDir = fullfile(resultsDir, 'dxf_solid');
if ~exist(solidDXFDir, 'dir'); mkdir(solidDXFDir); end
poreDXFDir = fullfile(resultsDir, 'dxf_pore');
if ~exist(poreDXFDir, 'dir'); mkdir(poreDXFDir); end
domainDXFPath = fullfile(resultsDir, 'simulation_domain.dxf');
ExportDomainRectangleToDXF(lengthXAxis, lengthYAxis, domainDXFPath, 'DOMAIN');
imageDir = fullfile(resultsDir, 'interface_images');
if ~exist(imageDir, 'dir');    mkdir(imageDir); end
phreeqcOutputDir = '';
phreeqcWorkDir = '';
phreeqcSummaryCsvFile = '';
phreeqcSpeciesPlotDir = '';
if usePhreeqc
    phreeqcOutputDir = fullfile(resultsDir, 'phreeqc_results');
    phreeqcWorkDir = fullfile(phreeqcOutputDir, 'inputs');
    if ~exist(phreeqcOutputDir, 'dir'); mkdir(phreeqcOutputDir); end
    if ~exist(phreeqcWorkDir, 'dir'); mkdir(phreeqcWorkDir); end
    phreeqcSummaryCsvFile = fullfile(phreeqcOutputDir, 'phreeqc_summary_log.csv');
    initializePhreeqcSummaryLog(phreeqcSummaryCsvFile);
    phreeqcOptions = struct( ...
        'databasePath', char(phreeqcDatabasePath), ...
        'workDir', phreeqcWorkDir, ...
        'temperatureC', phreeqcTemperatureC, ...
        'mineralName', 'Calcite', ...
        'mineralFormula', 'CaCO3', ...
        'rateLaw', phreeqcRateLaw, ...
        'phreeqcRunGroup', phreeqcRunGroup, ...
        'rateCoefficientTST', phreeqcTstRateCoefficient, ...
        'phreeqcTstRateCoefficient', phreeqcTstRateCoefficient, ...
        'kineticsCorrectionFactor', phreeqcKineticsCorrectionFactor, ...
        'maxSpecificSurfaceArea', phreeqcMaxSpecificSurfaceArea, ...
        'badStepMax', phreeqcBadStepMax, ...
        'minHForPHMolL', phreeqcMinHForPHMolL, ...
        'minActiveWaterVolumeFraction', phreeqcMinActiveWaterVolumeFraction, ...
        'minActiveWaterVolumeCm3', phreeqcMinActiveWaterVolumeCm3, ...
        'reactionWaterVolumeFloorFraction', phreeqcReactionWaterVolumeFloorFraction, ...
        'reactionWaterVolumeFloorCm3', phreeqcReactionWaterVolumeFloorCm3, ...
        'reactionProjectionMinCells', phreeqcReactionProjectionMinCells, ...
        'reactionProjectionMaxCells', phreeqcReactionProjectionMaxCells, ...
        'reactionProjectionMaxRings', phreeqcReactionProjectionMaxRings, ...
        'reactionProjectionFlowBias', phreeqcReactionProjectionFlowBias, ...
        'reactionProjectionMode', phreeqcReactionProjectionMode, ...
        'reactionWaterVolumeMode', phreeqcReactionWaterVolumeMode, ...
        'reactNeutralInterfaceCells', phreeqcReactNeutralInterfaceCells, ...
        'solutionWaterKg', phreeqcSolutionWaterKg, ...
        'writeSolutionWaterLine', phreeqcWriteSolutionWaterLine, ...
        'kineticsReservoirMoles', phreeqcKineticsReservoirMoles, ...
        'kineticsTolerance', phreeqcKineticsTolerance, ...
        'timeStepSize', 1, ...
        'timeStepIndex', 0);
    phreeqcSpeciesData = initializePhreeqcSpeciesDataFromTransports( ...
        hydrogenTransport, calciumTransport, carbonTransport, sodiumTransport, chlorideTransport, 0);
    phreeqcCalciteRateData = zeros(gridHyPHM.numT, 1);
    if strcmpi(phreeqcRateLaw, 'tst_match')
        initialPhreeqcGeometryState = computePhreeqcGeometryState( ...
            gridHyPHM, levelSet{1}(:, 1), molarVolume, thickness);
        initialPhreeqcState = collectPhreeqcState(hydrogenTransport, calciumTransport, ...
            carbonTransport, sodiumTransport, chlorideTransport, 0, initialPhreeqcGeometryState);
        initialPhreeqcState.prescribed_calcite_dissolved_moles = ...
            zeros(size(initialPhreeqcGeometryState.water_volume_cm3(:)));
        initialPhreeqcState.reaction_water_volume_cm3 = ...
            initialPhreeqcGeometryState.water_volume_cm3(:);
        initialEquilibriumOptions = phreeqcOptions;
        initialEquilibriumOptions.timeStepIndex = -1;
        initialEquilibriumData = RunPhreeqcCalciteBatch( ...
            initialPhreeqcState, initialEquilibriumOptions);
        phreeqcSpeciesData = updatePhreeqcStateFromResult( ...
            phreeqcSpeciesData, initialEquilibriumData);
        initialPhreeqcStepSize = transportStepper.timepts(2) - transportStepper.timepts(1);
        initialRateDriver = getPhreeqcTstRateDriver(phreeqcSpeciesData);
        initialTstStep = ComputeTstMatchStepDissolution( ...
            initialRateDriver, initialPhreeqcGeometryState.interface_area_cm2, ...
            initialPhreeqcGeometryState.calcite_moles, initialPhreeqcStepSize, ...
            gridHyPHM.V0T, VertexTriMatrix, gridHyPHM.numV, ...
            phreeqcTstRateCoefficient, molarVolume);
        phreeqcCalciteRateData = initialTstStep.ratePerArea_mol_cm2_s;
    end
    phreeqcSpeciesData.calciteRatePerArea_mol_cm2_s = phreeqcCalciteRateData;
    phreeqcSpeciesData.calciteDissolvedMoles = zeros(gridHyPHM.numT, 1);
end

if enableNMRSimulation
    [nmrConfig, nmrComsolOutputDir, nmrInversionOutputDir, nmrSyncLogFile] = ...
        initializeNMRSync(reactiveRoot, resultsDir, pngNMRConfig);
end
if enableNMRSurrogate
    [nmrSurrogateConfig, nmrSurrogateOutputDir, nmrSurrogateInversionOutputDir, ...
        nmrSurrogateMaskDir, nmrSurrogateSyncLogFile] = initializeNMRSurrogateSync( ...
        reactiveRoot, resultsDir, pngNMRConfig);
end
if enablePNGSimulation
    [pngNMRConfig, pngNMROutputDir, pngNMRSyncLogFile, pngNMRRuntimeConfigFile] = ...
        initializePNGNMRSync(rtmDir, reactiveRoot, resultsDir, pngNMRConfig);
end

% 创建单独子图保存文件夹
individualPlotsDir = fullfile(resultsDir, 'individual_plots');
if ~exist(individualPlotsDir, 'dir'); mkdir(individualPlotsDir); end
if usePhreeqc
    phreeqcSpeciesPlotDir = fullfile(individualPlotsDir, 'phreeqc_species');
    if ~exist(phreeqcSpeciesPlotDir, 'dir'); mkdir(phreeqcSpeciesPlotDir); end
end
subplot1Dir = fullfile(individualPlotsDir, 'concentration');
subplot2Dir = fullfile(individualPlotsDir, 'interface');
subplot3Dir = fullfile(individualPlotsDir, 'velocity');
if ~exist(subplot1Dir, 'dir'); mkdir(subplot1Dir); end
if ~exist(subplot2Dir, 'dir'); mkdir(subplot2Dir); end
if ~exist(subplot3Dir, 'dir'); mkdir(subplot3Dir); end

metadata = struct();
metadata.schema_version = "rtm_run_metadata_v1";
metadata.run_id = string(runName);
metadata.created_at = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
metadata.results_dir = string(resultsDir);
metadata.layoutType = string(layoutType);
metadata.useExternalGeometry = useExternalGeometry;
metadata.externalGeometry = struct( ...
    'type', string(externalGeometryType), ...
    'useExternalDxfGeometry', useExternalDxfGeometry, ...
    'dxfPath', string(externalDxfPath), ...
    'dxfDomainLayerNames', string(externalDxfDomainLayerNames), ...
    'dxfSolidLayerNames', string(externalDxfSolidLayerNames), ...
    'dxfReferenceLength', externalDxfReferenceLength, ...
    'dxfReferenceLengthCm', externalDxfReferenceLengthCm, ...
    'dxfImportDirection', string(externalDxfImportDirection));
metadata.externalGeometry.tifPath = string(tifPath);
metadata.externalGeometry.tifSolidValue = externalTifSolidValue;
metadata.externalGeometry.tifLiquidValue = externalTifLiquidValue;
metadata.externalGeometry.tifPixelSizeMicron = pixelSizeMicron;
metadata.externalGeometry.tifSmoothingSigmaPixels = externalTifSmoothingSigmaPixels;
metadata.externalGeometry.tifCropTopRows = externalTifCropTopRows;
metadata.externalGeometry.tifCropBottomRows = externalTifCropBottomRows;
metadata.externalGeometry.tifCropLeftCols = externalTifCropLeftCols;
metadata.externalGeometry.tifCropRightCols = externalTifCropRightCols;
if useExternalDxfGeometry && ~isempty(externalDxfGeometry)
    metadata.externalGeometry.dxfScaleCmPerUnit = externalDxfGeometry.scaleCmPerDxfUnit;
    metadata.externalGeometry.dxfDomainLayer = externalDxfGeometry.domainLayer;
    metadata.externalGeometry.dxfSolidLayers = externalDxfGeometry.solidLayers;
    metadata.externalGeometry.dxfDomainBoundsRaw = externalDxfGeometry.domainBoundsRaw;
    metadata.externalGeometry.dxfResolutionX = externalDxfGeometry.resolutionX;
    metadata.externalGeometry.dxfResolutionY = externalDxfGeometry.resolutionY;
    metadata.externalGeometry.dxfImportDirectionApplied = externalDxfGeometry.importDirection;
end
metadata.parameters = struct( ...
    'Da', DamkohlerNumber, ...
    'Pe', pecletNumber, ...
    'Re', reynoldsNumber, ...
    'characteristicLength_cm', characteristicLength, ...
    'lengthXAxis_cm', lengthXAxis, ...
    'lengthYAxis_cm', lengthYAxis, ...
    'circleRadius_cm', circleRadius, ...
    'circleSpacing_cm', circleSpacing, ...
    'targetAvgSpacing_cm', targetAvgSpacing, ...
    'minThroatRandom_cm', minThroatRandom, ...
    'useRandomParticleRadii', useRandomParticleRadii, ...
    'randomParticleRadiusMin_cm', randomParticleRadiusMin, ...
    'randomParticleRadiusMax_cm', randomParticleRadiusMax, ...
    'targetInitialPorosity', targetInitialPorosity, ...
    'inletVelocity_cm_s', inletVelocity, ...
    'flowDirection', string(flowConfig.direction), ...
    'initialHydrogenConcentration_mol_cm3', initialHydrogenConcentration, ...
    'diffusionCoefficient_cm2_s', diffusionCoefficient, ...
    'molarVolume_cm3_mol', molarVolume, ...
    'rateCoefficientTST_mol_dm2_s', rateCoefficientTST, ...
    'reactionModel', string(reactionModel), ...
    'phreeqcRunGroup', string(phreeqcRunGroup), ...
    'phreeqcRateLaw', string(phreeqcRateLaw), ...
    'phreeqcTstRateCoefficient_mol_dm2_s', phreeqcTstRateCoefficient, ...
    'phreeqcDatabasePath', string(phreeqcDatabasePath), ...
    'phreeqcTemperatureC', phreeqcTemperatureC, ...
    'phreeqcKineticsCorrectionFactor', phreeqcKineticsCorrectionFactor, ...
    'phreeqcMaxSpecificSurfaceArea_m2_kgw', phreeqcMaxSpecificSurfaceArea, ...
    'phreeqcBadStepMax', phreeqcBadStepMax, ...
    'phreeqcMinHForPHMolL', phreeqcMinHForPHMolL, ...
    'phreeqcMinActiveWaterVolumeFraction', phreeqcMinActiveWaterVolumeFraction, ...
    'phreeqcMinActiveWaterVolumeCm3', phreeqcMinActiveWaterVolumeCm3, ...
    'phreeqcReactionWaterVolumeFloorFraction', phreeqcReactionWaterVolumeFloorFraction, ...
    'phreeqcReactionWaterVolumeFloorCm3', phreeqcReactionWaterVolumeFloorCm3, ...
    'phreeqcReactionProjectionMinCells', phreeqcReactionProjectionMinCells, ...
    'phreeqcReactionProjectionMaxCells', phreeqcReactionProjectionMaxCells, ...
    'phreeqcReactionProjectionMaxRings', phreeqcReactionProjectionMaxRings, ...
    'phreeqcReactionProjectionFlowBias', phreeqcReactionProjectionFlowBias, ...
    'phreeqcReactionProjectionMode', string(phreeqcReactionProjectionMode), ...
    'phreeqcReactionWaterVolumeMode', string(phreeqcReactionWaterVolumeMode), ...
    'phreeqcReactNeutralInterfaceCells', phreeqcReactNeutralInterfaceCells, ...
    'phreeqcSolutionWaterKg', phreeqcSolutionWaterKg, ...
    'phreeqcWriteSolutionWaterLine', phreeqcWriteSolutionWaterLine, ...
    'phreeqcKineticsReservoirMoles', phreeqcKineticsReservoirMoles, ...
    'phreeqcKineticsTolerance', phreeqcKineticsTolerance, ...
    'phreeqcTransportMaxFactor', phreeqcTransportMaxFactor, ...
    'initialCalciumConcentration_mol_cm3', initialCalciumConcentration, ...
    'initialCarbonConcentration_mol_cm3', initialCarbonConcentration, ...
    'initialSodiumConcentration_mol_cm3', initialSodiumConcentration, ...
    'initialChlorideConcentration_mol_cm3', initialChlorideConcentration, ...
    'inletCalciumConcentration_mol_cm3', inletCalciumConcentration, ...
    'inletCarbonConcentration_mol_cm3', inletCarbonConcentration, ...
    'inletSodiumConcentration_mol_cm3', inletSodiumConcentration, ...
    'inletChlorideConcentration_mol_cm3', inletChlorideConcentration, ...
    'timeStepperType', string(timeStepperType), ...
    'initialMacroscaleTimeStepSize_s', initialMacroscaleTimeStepSize, ...
    'maximalStep_s', maximalStep, ...
    'adaptiveMaxTimeStep_s', adaptiveMaxTimeStep, ...
    'targetDissolutionSlices', targetDissolutionSlices, ...
    'targetDissolutionSliceSafetyFactor', targetDissolutionSliceSafetyFactor, ...
    'targetDissolutionSliceMinExtraSteps', targetDissolutionSliceMinExtraSteps, ...
    'allowAdaptiveTimeGridExtension', allowAdaptiveTimeGridExtension, ...
    'initialPorosityForTargetSlices', initialPorosityForTargetSlices, ...
    'porosityStepTarget', porosityStepTarget, ...
    'porosityStepTolerance', porosityStepTolerance, ...
    'porosityStepUpperFactor', porosityStepUpperFactor, ...
    'adaptiveGrowthFactor', adaptiveGrowthFactor, ...
    'adaptiveShrinkSafety', adaptiveShrinkSafety, ...
    'targetSliceStepBoundaryFactor', targetSliceStepBoundaryFactor, ...
    'transportUpwindMode', string(hydrogenTransport.isUpwind), ...
    'enableHardCflLimit', enableHardCflLimit, ...
    'hardCflLimitInitialGrid', hardCflLimitInitialGrid, ...
    'enableConcentrationCflLimit', enableConcentrationCflLimit, ...
    'concentrationOvershootThreshold', concentrationOvershootThreshold, ...
    'concentrationHardOvershootThreshold', concentrationHardOvershootThreshold, ...
    'concentrationShrinkFactor', concentrationShrinkFactor, ...
    'advectiveCflSafety', advectiveCflSafety, ...
    'reactiveCflSafety', reactiveCflSafety, ...
    'hardCflMinTimeStep_s', hardCflMinTimeStep, ...
    'hardCflMaxInitialSteps', hardCflMaxInitialSteps, ...
    'hardCflInitialCandidateStep_s', hardCflInitialMaxStep, ...
    'meshMinLengthForCfl_cm', meshMinLengthForCfl, ...
    'endTime_s', endTime, ...
    'maxTotalTimeSteps', maxTotalTimeSteps, ...
    'maxExperimentWallSeconds', maxExperimentWallSeconds, ...
    'numPartitionsMicroscale', numPartitionsMicroscale, ...
    'meshNumPartitionsX', nxParts, ...
    'meshNumPartitionsY', nyParts, ...
    'meshPartitionMode', string(meshPartitionMode), ...
    'meshTargetElementSizeCm', meshTargetElementSizeCm, ...
    'saveMeshDiagnostics', saveMeshDiagnostics);
metadata.outputs = struct( ...
    'global_evolution', string(xlsxFile), ...
    'global_evolution_csv', string(logfile), ...
    'tortuosity_segments', string(tortuositySegXlsxFile), ...
    'tortuosity_segments_csv', string(tortuositySegCsvFile), ...
    'stability_diagnostics_csv', string(stabilityDiagnosticsCsvFile), ...
    'dxf_pore_dir', string(poreDXFDir), ...
    'dxf_solid_dir', string(solidDXFDir), ...
    'mesh_diagnostics_dir', string(meshDiagnosticsDir), ...
    'phreeqc_results_dir', string(phreeqcOutputDir), ...
    'phreeqc_inputs_dir', string(phreeqcWorkDir), ...
    'phreeqc_summary_log_csv', string(phreeqcSummaryCsvFile), ...
    'phreeqc_species_plot_dir', string(phreeqcSpeciesPlotDir));
metadata.nmr = struct( ...
    'enabled', enableNMRSimulation, ...
    'method', string(nmrOptions.nmr_method), ...
    'trigger', "after_exported_pore_solid_dxf_pair", ...
    'config_source', string(fullfile(reactiveRoot, 'automation', 'AutomationConfig.m')), ...
    'workflow_config_source', string(fullfile(rtmDir, 'NMRSimulationConfig.m')), ...
    'comsol_mph_file', string(cfgget(pngNMRConfig, 'mph_file', '')), ...
    'inversion_script', string(cfgget(pngNMRConfig, 'inversion_script', '')), ...
    'comsol_results_dir', string(nmrComsolOutputDir), ...
    'inversion_results_dir', string(nmrInversionOutputDir), ...
    'sync_log', string(nmrSyncLogFile), ...
    'surrogate_enabled', enableNMRSurrogate, ...
    'surrogate_model_path', string(cfgget(pngNMRConfig, 'nmrSurrogateModelPath', '')), ...
    'surrogate_results_dir', string(nmrSurrogateOutputDir), ...
    'surrogate_inversion_results_dir', string(nmrSurrogateInversionOutputDir), ...
    'surrogate_sync_log', string(nmrSurrogateSyncLogFile), ...
    'png_enabled', enablePNGSimulation, ...
    'png_method', string(cfgget(pngNMRConfig, 'method', 'png_mesh')), ...
    'png_results_dir', string(pngNMROutputDir), ...
    'png_sync_log', string(pngNMRSyncLogFile), ...
    'png_runtime_config', string(pngNMRRuntimeConfigFile), ...
    'png_config_source', string(fullfile(rtmDir, 'NMRSimulationConfig.m')));
writeJsonFile(fullfile(resultsDir, 'run_metadata.json'), metadata);

%% Time iteration
useAdaptivePorositySteps = strcmpi(timeStepperType, 'adaptive_porosity');
adaptivePreviousPorosity = NaN;
timedOut = false;
timeoutMessage = "";
while transportStepper.next
    timeIterationStep = transportStepper.curstep;
    macroscaleTimeStepSize = transportStepper.curtau;
    disp('----------------------------------------');
    fprintf('  Total Runtime: %.2f seconds\n', toc(simulationTimer));
    disp(['Time step ', num2str(timeIterationStep)]);
    disp(['  Current time = ', ...
        num2str(currentTime)]);
    disp(['  Current time step size = ', ...
        num2str(macroscaleTimeStepSize)]);

    if (currentTime + macroscaleTimeStepSize > endTime - EPS)
        macroscaleTimeStepSize = endTime - currentTime;
        if (abs(macroscaleTimeStepSize) < EPS)
            break;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   inner corrector
    if timeIterationStep > 0
        if usePhreeqc
            speciesCells = {hydrogenTransport; calciumTransport; carbonTransport; sodiumTransport; chlorideTransport};
        else
            speciesCells = {hydrogenTransport};
        end
        [out, Corrector] = Continuation(speciesCells, timeIterationStep, levelSet{1}(:, timeIterationStep), 20);
        hydrogenTransport = out{1};
        if usePhreeqc
            calciumTransport = out{2};
            carbonTransport = out{3};
            sodiumTransport = out{4};
            chlorideTransport = out{5};
        end
    else
        if usePhreeqc
            Corrector = zeros(gridHyPHM.numT, 5);
        else
            Corrector = zeros(gridHyPHM.numT, 1);
        end
    end

    stopAfterThisStep = false;
    hydrogenData = hydrogenTransport.U.getdata(timeIterationStep-1);

    if ~phreeqcTransportFirstSplit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level set evolution step

    disp('Evolution of level set ...');
    tic;

    phreeqcLevelSetStep = [];

    % Convert reaction-driving data from T --> V.
    if usePhreeqc
        if strcmpi(phreeqcRateLaw, 'tst_match')
            phreeqcLevelSetGeometryState = computePhreeqcGeometryState( ...
                gridHyPHM, levelSet{1}(:, timeIterationStep), molarVolume, thickness);
            phreeqcLevelSetStep = ComputeTstMatchStepDissolution( ...
                getPhreeqcTstRateDriver(phreeqcSpeciesData), ...
                phreeqcLevelSetGeometryState.interface_area_cm2, ...
                phreeqcLevelSetGeometryState.calcite_moles, macroscaleTimeStepSize, ...
                gridHyPHM.V0T, VertexTriMatrix, gridHyPHM.numV, ...
                phreeqcTstRateCoefficient, molarVolume);
            normalSpeed = phreeqcLevelSetStep.normalSpeed_cm_s;
        else
            phreeqcRateV = averageTriangleDataToVertices(phreeqcCalciteRateData, VertexTriMatrix, gridHyPHM.numV);
            normalSpeed = molarVolume * phreeqcRateV;
        end
    else
        normalSpeed = computeLegacyTstNormalSpeedFromTriangleHydrogen( ...
            hydrogenData, VertexTriMatrix, gridHyPHM.numV, ...
            molarVolume, rateCoefficientTST);
    end
    normalSpeedMax = max(abs(normalSpeed));

    % 可视化（保持结构，自动调整维度）
    % figure('Visible', 'off');
    % imagesc(reshape(normalSpeed, nxNodes, nyNodes)');
    % title('Normal speed'); axis equal tight; colorbar;

    CFL = 1 / 4 * 1 / 10 * 1 / numPartitionsMicroscale / max(normalSpeedMax, eps);

    oldMicroscaleTime = currentTime;
    levelSetTargetTime = currentTime + macroscaleTimeStepSize;
    for j = 1:ceil(macroscaleTimeStepSize/CFL)
        microscaleTimeStepSize = min(levelSetTargetTime - oldMicroscaleTime, CFL);
        newMicroscaleTime = oldMicroscaleTime + microscaleTimeStepSize;
        % [] argument is unused in method (needed for implicit methods)
        currentLevelSetDataCells{1} = levelSetEquationTimeStep( ...
            newMicroscaleTime, ...
            oldMicroscaleTime, oldLevelSetDataCells{1}, microscaleGrid, normalSpeed, 1);
        oldMicroscaleTime = newMicroscaleTime;
        oldLevelSetDataCells{1} = currentLevelSetDataCells{1};
    end

    levelSet{1}(:, timeIterationStep + 1) = currentLevelSetDataCells{1};
    currentTime = currentTime + macroscaleTimeStepSize;

    levelSetEvolutionTime(timeIterationStep) = toc;
    disp(['    ... done in ', ...
        num2str(levelSetEvolutionTime(timeIterationStep)), ' seconds.']);

    hydrogenTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
    if usePhreeqc
        calciumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
        carbonTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
        sodiumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
        chlorideTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculation of effective parameters

    disp('Calculation of effective parameters ...');
    tic;

    for i = 1:numberOfSlices
        [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
            triangleSurfaces] = assembleCellProblem(microscaleGrid, ...
            levelSet{i}(:, timeIterationStep + 1));
        SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
        [diffusion, porosities] = computeDiffusionTensor(microscaleGrid, SOL, ...
            triangleVolumes);
        surfaceArea{i}(timeIterationStep + 1) = sum(triangleSurfaces);
    end
    Volume(timeIterationStep+1) = (lengthYAxis * lengthXAxis - porosities);
    cellProblemTime(timeIterationStep) = toc;
    disp(['    ... done in ', num2str(cellProblemTime(timeIterationStep)), ...
        ' seconds.']);

    % === 检查 Grain Volume 是否已溶解完毕 ===
    grainVolumeTolerance = 1e-10;  % 体积阈值，低于此值视为完全溶解
    if Volume(timeIterationStep+1) <= grainVolumeTolerance
        fprintf('\n========================================\n');
        fprintf('>>> Grain Volume reached zero (%.2e cm^3)\n', Volume(timeIterationStep+1));
        fprintf('>>> Will complete this step to export final dissolution images...\n');
        fprintf('========================================\n\n');
        % 设置标志，完成当前步骤后再退出
        stopAfterThisStep = true;
    else
        stopAfterThisStep = false;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Stokes
    %Calculate Pressure distribution and Stokes velocity
    disp(' ');
    disp('Initializing Stokes problem...');

    StokesL.L.setdata(levelSet{1}(:, timeIterationStep + 1))

    flowStepper.next;
    StokesL.computeLevel('s');

    helper = StokesL.U.getdata(1);
    flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));
    flowStepper.prev;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --- 新增：基于边界流量与压降估算等效渗透率（cm^2） ---
    edgesRight = flowConfig.outletEdgePredicate(gridHyPHM.baryE);
    edgeFlux = flow.getdata(1);                     % per-edge flux or normal velocity
    Q = sum(edgeFlux .* gridHyPHM.areaE .* double(edgesRight));  % [cm^3/s]

    pData = phead.getdata(1);
    leftNodes = flowConfig.inletNodePredicate(gridHyPHM.coordV);
    rightNodes = flowConfig.outletNodePredicate(gridHyPHM.coordV);
    dp = mean(pData(leftNodes)) - mean(pData(rightNodes));

    A = flowConfig.crossSectionLength * thickness;  % [cm^2]
    L = flowConfig.flowLength;                      % [cm]

    if abs(dp) < EPS
        k_eff_cm2 = NaN;
    else
        k_eff_cm2 = mu * Q * L / (A * dp);          % 渗透率（cm^2）
    end

    % --- 单位转换：cm^2 → mD ---
    % 正确值：1 Darcy = 9.869233e-9 cm², 1 mD = 9.869233e-12 cm²
    darcy_to_cm2 = 9.869233e-9;          % 1 Darcy = 9.869233e-9 cm^2
    md_to_cm2    = darcy_to_cm2 * 1e-3;  % 1 mD = 9.869233e-12 cm^2
    k_eff_mD = k_eff_cm2 / md_to_cm2;

    % k_eff_mD = 100000;  % 占位符，保持代码结构完整,待替换为真实渗透率计算
    % 存储（单位：mD）
    Permeability(timeIterationStep+1) = k_eff_mD;
    
    % 存储初始渗透率 k0（仅在第一步时设置）
    if isnan(k0) && ~isnan(k_eff_mD)
        k0 = k_eff_mD;
    end
    
    % --- 迂曲度计算：τ = Σ|v| / Σv_flow，仅在孔隙空间节点上求和 ---
    % 方法见：基于流场的迂曲度计算.md （Duda et al. 2011 / Koponen et al.）
    velocityDataTau = StokesL.U.getdata(1);  % P2P2 节点速度，[numV+numE, 2]
    nP2 = gridHyPHM.numV + gridHyPHM.numE;
    if size(velocityDataTau, 2) == 2 && size(velocityDataTau, 1) == nP2
        vx_p2 = velocityDataTau(:, 1);
        vy_p2 = velocityDataTau(:, 2);
    elseif size(velocityDataTau, 2) == 1 && length(velocityDataTau) == 2 * nP2
        vx_p2 = velocityDataTau(1:nP2);
        vy_p2 = velocityDataTau(nP2+1:end);
    else
        vx_p2 = [];
        vy_p2 = [];
    end

    if ~isempty(vx_p2)
        % 当前水平集（P1 顶点值）
        lsNow = levelSet{1}(:, timeIterationStep + 1);  % length = numV
        % 边中点水平集：取两端顶点均值
        lsEdgeMid = 0.5 * (lsNow(gridHyPHM.V0E(:,1)) + lsNow(gridHyPHM.V0E(:,2)));
        ls_p2 = [lsNow; lsEdgeMid];  % length = nP2
        % 孔隙节点掩膜（levelSet < 0 为孔隙）
        poreMask_p2 = ls_p2 < 0;
        sum_speed_pore = sum(sqrt(vx_p2(poreMask_p2).^2 + vy_p2(poreMask_p2).^2));
        if flowConfig.axisIndex == 1
            v_axis_p2 = vx_p2;
        else
            v_axis_p2 = vy_p2;
        end
        sum_vaxis_pore = sum(v_axis_p2(poreMask_p2));
        if sum_vaxis_pore > EPS
            tau_now = sum_speed_pore / sum_vaxis_pore;
        else
            tau_now = NaN;
        end

        % --- 新增：5 段分区迂曲度（沿主流方向等分）---
        x_p2 = flowConfig.axisCoordinate([gridHyPHM.coordV; gridHyPHM.baryE]);
        for iSeg = 1:numTortuositySegments
            xL = tortuositySegEdges(iSeg);
            xR = tortuositySegEdges(iSeg + 1);
            if iSeg < numTortuositySegments
                segMask = (x_p2 >= xL) & (x_p2 < xR);
            else
                segMask = (x_p2 >= xL) & (x_p2 <= xR);
            end
            poreSegMask = poreMask_p2 & segMask;

            sum_speed_seg = sum(sqrt(vx_p2(poreSegMask).^2 + vy_p2(poreSegMask).^2));
            sum_vaxis_seg = sum(v_axis_p2(poreSegMask));
            if sum_vaxis_seg > EPS
                TortuositySegments(iSeg, timeIterationStep+1) = sum_speed_seg / sum_vaxis_seg;
            else
                TortuositySegments(iSeg, timeIterationStep+1) = NaN;
            end
        end
    else
        tau_now = NaN;
        TortuositySegments(:, timeIterationStep+1) = NaN;
    end
    Tortuosity(timeIterationStep+1) = tau_now;
    % -------------------------------------------------------------------------
    end

    %% Macroscopic transport step

    if phreeqcTransportFirstSplit
        hydrogenTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep));
        calciumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep));
        carbonTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep));
        sodiumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep));
        chlorideTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep));
    end

    surfaceAreaFunc = @(x) areaHelperFun(-1, x, surfaceArea, numberOfSlices, ...
        timeIterationStep);

    Levels = hydrogenTransport.L.getdata(timeIterationStep);

    % concentrations P_0(T)--> P_0(E)
    hydrogenDataE = zeros(gridHyPHM.numE, 1);
    for i = 1:gridHyPHM.numE
        if ((Levels(gridHyPHM.V0E(i, 1)) > -eps) & (Levels(gridHyPHM.V0E(i, 2)) > -eps))
            helper = gridHyPHM.T0E(i, :);
            helper = helper(helper > 0); % 去掉无效索引
            hydrogenDataE(i, 1) = max(hydrogenData(helper));
        end
    end

    hydrogenTransportRhsData = zeros(gridHyPHM.numE, 1);
    hydrogenTransport.A.setdata(timeIterationStep, @(t, x) 1);
    hydrogenTransport.gF.setdata(timeIterationStep, hydrogenTransport.gF.getdata(0)+hydrogenTransportRhsData);
    hydrogenTransport.C.setdata(timeIterationStep, flow.getdata(1));

    SorceScale = zeros(gridHyPHM.numT, 1);
    for kT = 1:gridHyPHM.numT
        L = 0;
        for i = 1:3
            if Levels(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 1)) > -eps & Levels(gridHyPHM.V0E(gridHyPHM.E0T(kT, i), 2)) > -eps
                L = L + gridHyPHM.areaE(gridHyPHM.E0T(kT, i)) / gridHyPHM.areaT(kT);
            end
        end
        SorceScale(kT) = L;
    end

    if usePhreeqc
        setPhreeqcTransportStepData(calciumTransport, timeIterationStep, flow.getdata(1));
        setPhreeqcTransportStepData(carbonTransport, timeIterationStep, flow.getdata(1));
        setPhreeqcTransportStepData(sodiumTransport, timeIterationStep, flow.getdata(1));
        setPhreeqcTransportStepData(chlorideTransport, timeIterationStep, flow.getdata(1));
        speciesCells = {hydrogenTransport; calciumTransport; carbonTransport; sodiumTransport; chlorideTransport};
        for iSpeciesSolve = 1:numel(speciesCells)
            [nonlinearFunc, nonlinearJacFunc] = zeroPhreeqcNonlinearFunctions(1, gridHyPHM.numT);
            newtonIteration(speciesCells(iSpeciesSolve), nonlinearFunc, nonlinearJacFunc, 2);
        end
        applyTransportCorrector(speciesCells, Corrector, timeIterationStep);
        speciesCells = limitPhreeqcTransportBounds(speciesCells, phreeqcTransportUpperBounds, ...
            phreeqcTransportSpeciesNames, timeIterationStep);
        hydrogenTransport = speciesCells{1};
        calciumTransport = speciesCells{2};
        carbonTransport = speciesCells{3};
        sodiumTransport = speciesCells{4};
        chlorideTransport = speciesCells{5};

        phreeqcOptions.timeStepSize = macroscaleTimeStepSize;
        phreeqcOptions.timeStepIndex = timeIterationStep;
        phreeqcReactionLevelSetData = levelSet{1}(:, timeIterationStep + 1);
        if phreeqcTransportFirstSplit
            phreeqcReactionLevelSetData = levelSet{1}(:, timeIterationStep);
        end
        phreeqcGeometryState = computePhreeqcGeometryState( ...
            gridHyPHM, phreeqcReactionLevelSetData, molarVolume, thickness);
        phreeqcBaseState = collectPhreeqcState(hydrogenTransport, calciumTransport, carbonTransport, ...
            sodiumTransport, chlorideTransport, timeIterationStep, phreeqcGeometryState);
        phreeqcState = phreeqcBaseState;
        prescribedInterfaceMoles = [];
        projectedMoles = [];
        agglomerateWeightMatrix = [];
        phreeqcReactionInputState = [];
        if strcmpi(phreeqcRateLaw, 'tst_match')
            phreeqcEquilibriumState = phreeqcState;
            phreeqcEquilibriumState.prescribed_calcite_dissolved_moles = ...
                zeros(size(phreeqcGeometryState.water_volume_cm3(:)));
            phreeqcEquilibriumState.reaction_water_volume_cm3 = ...
                phreeqcGeometryState.water_volume_cm3(:);
            phreeqcEquilibriumOptions = phreeqcOptions;
            phreeqcEquilibriumOptions.timeStepIndex = 100000 + timeIterationStep;
            phreeqcEquilibriumData = RunPhreeqcCalciteBatch( ...
                phreeqcEquilibriumState, phreeqcEquilibriumOptions);
            phreeqcState = updatePhreeqcStateFromResult(phreeqcState, phreeqcEquilibriumData);
            phreeqcBaseState = phreeqcState;

            if phreeqcTransportFirstSplit
                phreeqcLevelSetStep = ComputeTstMatchStepDissolution( ...
                    getPhreeqcTstRateDriver(phreeqcState), ...
                    phreeqcGeometryState.interface_area_cm2, ...
                    phreeqcGeometryState.calcite_moles, macroscaleTimeStepSize, ...
                    gridHyPHM.V0T, VertexTriMatrix, gridHyPHM.numV, ...
                    phreeqcTstRateCoefficient, molarVolume);
            end
            if isempty(phreeqcLevelSetStep)
                error('RTSPHEM:Phreeqc:MissingTstStep', ...
                    ['PHREEQC TST-match requires the same dissolution step ', ...
                    'for level-set motion and prescribed CaCO3 reaction.']);
            end
            phreeqcCalciteRateData = phreeqcLevelSetStep.ratePerArea_mol_cm2_s;
            prescribedInterfaceMoles = phreeqcLevelSetStep.prescribedMoles;
            switch phreeqcReactionProjectionMode
                case {'interface_cell', 'cut_cell', 'same_cell'}
                    projectedMoles = prescribedInterfaceMoles;
                    switch phreeqcReactionWaterVolumeMode
                        case {'actual_cell', 'cell', 'none'}
                            projectedReactionWaterCm3 = phreeqcGeometryState.water_volume_cm3(:);
                        case {'cut_cell_agglomerated', 'agglomerated', 'local_agglomerated'}
                            [projectedReactionWaterCm3, agglomerateWeightMatrix] = ...
                                computeAgglomeratedReactionWaterVolumes( ...
                                prescribedInterfaceMoles, gridHyPHM, VertexTriMatrix, ...
                                phreeqcReactionLevelSetData, ...
                                phreeqcGeometryState.water_volume_cm3, ...
                                phreeqcReactionProjectionMinCells, ...
                                phreeqcReactionProjectionMaxCells, ...
                                phreeqcReactionProjectionMaxRings);
                        otherwise
                            error('RTSPHEM:Phreeqc:UnknownReactionWaterVolumeMode', ...
                                'Unknown PHREEQC reaction water volume mode: %s.', ...
                                phreeqcReactionWaterVolumeMode);
                    end
                case {'neighbor_pore', 'pore_neighbor', 'projected_pore'}
                    triangleVelocityData = stokesVelocityToTriangleVelocity(gridHyPHM, StokesL.U.getdata(1));
                    [projectedMoles, projectedReactionWaterCm3] = projectInterfaceReactionToPoreCells( ...
                        prescribedInterfaceMoles, gridHyPHM, VertexTriMatrix, ...
                        phreeqcReactionLevelSetData, phreeqcGeometryState.water_volume_cm3, ...
                        phreeqcReactionProjectionMinCells, phreeqcReactionProjectionMaxCells, ...
                        phreeqcReactionProjectionMaxRings, triangleVelocityData, ...
                        phreeqcReactionProjectionFlowBias);
                otherwise
                    error('RTSPHEM:Phreeqc:UnknownReactionProjectionMode', ...
                        'Unknown PHREEQC reaction projection mode: %s.', ...
                        phreeqcReactionProjectionMode);
            end
            positiveWaterVolumes = projectedReactionWaterCm3(projectedReactionWaterCm3 > 0);
            reactionWaterFloorCm3 = phreeqcReactionWaterVolumeFloorCm3;
            if phreeqcReactionWaterVolumeFloorFraction > 0 && ~isempty(positiveWaterVolumes)
                reactionWaterFloorCm3 = max(reactionWaterFloorCm3, ...
                    max(positiveWaterVolumes) * phreeqcReactionWaterVolumeFloorFraction);
            end
            phreeqcState.prescribed_calcite_dissolved_moles = projectedMoles;
            phreeqcState.reaction_water_volume_cm3 = max(projectedReactionWaterCm3(:), ...
                reactionWaterFloorCm3);
            if ~isempty(agglomerateWeightMatrix) && nnz(agglomerateWeightMatrix) > 0
                phreeqcReactionInputState = BuildPhreeqcAgglomeratedState( ...
                    phreeqcState, phreeqcGeometryState.water_volume_cm3, ...
                    agglomerateWeightMatrix);
                phreeqcState = phreeqcReactionInputState;
            end
        end
        phreeqcSpeciesData = RunPhreeqcCalciteBatch( ...
            phreeqcState, phreeqcOptions);
        if strcmpi(phreeqcRateLaw, 'tst_match') && ...
                ~isempty(agglomerateWeightMatrix) && nnz(agglomerateWeightMatrix) > 0
            phreeqcSpeciesData = applyAgglomeratedPhreeqcUpdate( ...
                phreeqcSpeciesData, phreeqcBaseState, phreeqcReactionInputState, ...
                phreeqcGeometryState.water_volume_cm3, agglomerateWeightMatrix);
        end
        phreeqcSpeciesData = addPhreeqcWaterPhaseMassDiagnostics( ...
            phreeqcSpeciesData, phreeqcBaseState, phreeqcGeometryState.water_volume_cm3);
        applyPhreeqcResultToTransports(phreeqcSpeciesData, hydrogenTransport, calciumTransport, ...
            carbonTransport, sodiumTransport, chlorideTransport, timeIterationStep);
        if ~strcmpi(phreeqcRateLaw, 'tst_match')
            phreeqcCalciteRateData = ComputePhreeqcInterfaceRatePerArea( ...
                phreeqcSpeciesData, phreeqcGeometryState.interface_area_cm2, ...
                macroscaleTimeStepSize, phreeqcOptions, phreeqcState.h_mol_cm3);
        end
        phreeqcSpeciesData.calciteRatePerArea_mol_cm2_s = phreeqcCalciteRateData;
        if strcmpi(phreeqcRateLaw, 'tst_match')
            phreeqcSpeciesData.interfacePrescribedCalciteMoles = prescribedInterfaceMoles;
            phreeqcSpeciesData.projectedPrescribedCalciteMoles = projectedMoles;
            phreeqcSpeciesData.agglomerateMaxRowOverlap = computeAgglomerateMaxRowOverlap(agglomerateWeightMatrix);
            phreeqcSpeciesData.agglomerateMeanRowOverlap = computeAgglomerateMeanRowOverlap(agglomerateWeightMatrix);
            phreeqcSpeciesData.agglomerateMaxVolumeRatio = computeAgglomerateMaxVolumeRatio( ...
                phreeqcState.reaction_water_volume_cm3, phreeqcGeometryState.water_volume_cm3, ...
                projectedMoles);
        end
        phreeqcLogTime = currentTime;
        if phreeqcTransportFirstSplit
            phreeqcLogTime = currentTime + macroscaleTimeStepSize;
        end
        appendPhreeqcSummaryLog(phreeqcSummaryCsvFile, timeIterationStep, phreeqcLogTime, phreeqcSpeciesData);
    else
        speciesCells = {hydrogenTransport};
        nonlinearFunc = cell(1, 1);
        nonlinearFunc{1} = @(x, y) -SorceScale ...
            .* (x * rateCoefficientTST * 1000);

        nonlinearJacFunc = cell(1, 1);
        nonlinearJacFunc{1, 1} = @(x, y) -SorceScale ...
            .* (rateCoefficientTST * 1000);

        newtonIteration(speciesCells, nonlinearFunc, nonlinearJacFunc, 2);

        hydrogenTransport.U.setdata(timeIterationStep-1, hydrogenTransport.U.getdata(timeIterationStep - 1)-Corrector(:, 1));
    end

    if phreeqcTransportFirstSplit
        disp('Evolution of level set ...');
        tic;

        if isempty(phreeqcLevelSetStep)
            error('RTSPHEM:Phreeqc:MissingTransportFirstTstStep', ...
                'PHREEQC transport-first TST split did not compute a dissolution step.');
        end
        normalSpeed = phreeqcLevelSetStep.normalSpeed_cm_s;
        normalSpeedMax = max(abs(normalSpeed));
        CFL = 1 / 4 * 1 / 10 * 1 / numPartitionsMicroscale / max(normalSpeedMax, eps);

        oldMicroscaleTime = currentTime;
        levelSetTargetTime = currentTime + macroscaleTimeStepSize;
        for j = 1:ceil(macroscaleTimeStepSize/CFL)
            microscaleTimeStepSize = min(levelSetTargetTime - oldMicroscaleTime, CFL);
            newMicroscaleTime = oldMicroscaleTime + microscaleTimeStepSize;
            currentLevelSetDataCells{1} = levelSetEquationTimeStep( ...
                newMicroscaleTime, oldMicroscaleTime, oldLevelSetDataCells{1}, ...
                microscaleGrid, normalSpeed, 1);
            oldMicroscaleTime = newMicroscaleTime;
            oldLevelSetDataCells{1} = currentLevelSetDataCells{1};
        end

        levelSet{1}(:, timeIterationStep + 1) = currentLevelSetDataCells{1};
        currentTime = currentTime + macroscaleTimeStepSize;

        levelSetEvolutionTime(timeIterationStep) = toc;
        disp(['    ... done in ', ...
            num2str(levelSetEvolutionTime(timeIterationStep)), ' seconds.']);

        hydrogenTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
        calciumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
        carbonTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
        sodiumTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));
        chlorideTransport.L.setdata(timeIterationStep, levelSet{1}(:, timeIterationStep + 1));

        disp('Calculation of effective parameters ...');
        tic;

        for i = 1:numberOfSlices
            [cellProblemSystemMatrix, rhs, isDoF, triangleVolumes, ...
                triangleSurfaces] = assembleCellProblem(microscaleGrid, ...
                levelSet{i}(:, timeIterationStep + 1));
            SOL = solveSystemFE(microscaleGrid, cellProblemSystemMatrix, rhs, isDoF);
            [diffusion, porosities] = computeDiffusionTensor(microscaleGrid, SOL, ...
                triangleVolumes);
            surfaceArea{i}(timeIterationStep + 1) = sum(triangleSurfaces);
        end
        Volume(timeIterationStep+1) = (lengthYAxis * lengthXAxis - porosities);
        cellProblemTime(timeIterationStep) = toc;
        disp(['    ... done in ', num2str(cellProblemTime(timeIterationStep)), ...
            ' seconds.']);

        grainVolumeTolerance = 1e-10;
        if Volume(timeIterationStep+1) <= grainVolumeTolerance
            fprintf('\n========================================\n');
            fprintf('>>> Grain Volume reached zero (%.2e cm^3)\n', Volume(timeIterationStep+1));
            fprintf('>>> Will complete this step to export final dissolution images...\n');
            fprintf('========================================\n\n');
            stopAfterThisStep = true;
        else
            stopAfterThisStep = false;
        end

        disp(' ');
        disp('Initializing Stokes problem...');

        StokesL.L.setdata(levelSet{1}(:, timeIterationStep + 1))

        flowStepper.next;
        StokesL.computeLevel('s');

        helper = StokesL.U.getdata(1);
        flow.setdata(helper(gridHyPHM.numV + 1:end, 1).*gridHyPHM.nuE(:, 1)+helper(gridHyPHM.numV + 1:end, 2).*gridHyPHM.nuE(:, 2));
        flowStepper.prev;

        edgesRight = flowConfig.outletEdgePredicate(gridHyPHM.baryE);
        edgeFlux = flow.getdata(1);
        Q = sum(edgeFlux .* gridHyPHM.areaE .* double(edgesRight));

        pData = phead.getdata(1);
        leftNodes = flowConfig.inletNodePredicate(gridHyPHM.coordV);
        rightNodes = flowConfig.outletNodePredicate(gridHyPHM.coordV);
        dp = mean(pData(leftNodes)) - mean(pData(rightNodes));

        A = flowConfig.crossSectionLength * thickness;
        L = flowConfig.flowLength;

        if abs(dp) < EPS
            k_eff_cm2 = NaN;
        else
            k_eff_cm2 = mu * Q * L / (A * dp);
        end

        darcy_to_cm2 = 9.869233e-9;
        md_to_cm2 = darcy_to_cm2 * 1e-3;
        k_eff_mD = k_eff_cm2 / md_to_cm2;
        Permeability(timeIterationStep+1) = k_eff_mD;

        if isnan(k0) && ~isnan(k_eff_mD)
            k0 = k_eff_mD;
        end

        velocityDataTau = StokesL.U.getdata(1);
        nP2 = gridHyPHM.numV + gridHyPHM.numE;
        if size(velocityDataTau, 2) == 2 && size(velocityDataTau, 1) == nP2
            vx_p2 = velocityDataTau(:, 1);
            vy_p2 = velocityDataTau(:, 2);
        elseif size(velocityDataTau, 2) == 1 && length(velocityDataTau) == 2 * nP2
            vx_p2 = velocityDataTau(1:nP2);
            vy_p2 = velocityDataTau(nP2+1:end);
        else
            vx_p2 = [];
            vy_p2 = [];
        end

        if ~isempty(vx_p2)
            lsNow = levelSet{1}(:, timeIterationStep + 1);
            lsEdgeMid = 0.5 * (lsNow(gridHyPHM.V0E(:,1)) + lsNow(gridHyPHM.V0E(:,2)));
            ls_p2 = [lsNow; lsEdgeMid];
            poreMask_p2 = ls_p2 < 0;
            sum_speed_pore = sum(sqrt(vx_p2(poreMask_p2).^2 + vy_p2(poreMask_p2).^2));
            if flowConfig.axisIndex == 1
                v_axis_p2 = vx_p2;
            else
                v_axis_p2 = vy_p2;
            end
            sum_vaxis_pore = sum(v_axis_p2(poreMask_p2));
            if sum_vaxis_pore > EPS
                tau_now = sum_speed_pore / sum_vaxis_pore;
            else
                tau_now = NaN;
            end

            x_p2 = flowConfig.axisCoordinate([gridHyPHM.coordV; gridHyPHM.baryE]);
            for iSeg = 1:numTortuositySegments
                xL = tortuositySegEdges(iSeg);
                xR = tortuositySegEdges(iSeg + 1);
                if iSeg < numTortuositySegments
                    segMask = (x_p2 >= xL) & (x_p2 < xR);
                else
                    segMask = (x_p2 >= xL) & (x_p2 <= xR);
                end
                poreSegMask = poreMask_p2 & segMask;

                sum_speed_seg = sum(sqrt(vx_p2(poreSegMask).^2 + vy_p2(poreSegMask).^2));
                sum_vaxis_seg = sum(v_axis_p2(poreSegMask));
                if sum_vaxis_seg > EPS
                    TortuositySegments(iSeg, timeIterationStep+1) = sum_speed_seg / sum_vaxis_seg;
                else
                    TortuositySegments(iSeg, timeIterationStep+1) = NaN;
                end
            end
        else
            tau_now = NaN;
            TortuositySegments(:, timeIterationStep+1) = NaN;
        end
        Tortuosity(timeIterationStep+1) = tau_now;
    end

    % --- 计算出口处平均 H+ 浓度：必须在输运求解之后读取当前步浓度 ---
    hydrogenDataCurrent = hydrogenTransport.U.getdata(timeIterationStep);
    outletTriangles = flowConfig.outletTrianglePredicate(gridHyPHM.baryT);
    OutletHConc(timeIterationStep+1) = ComputeOutletConcentration(hydrogenDataCurrent, outletTriangles);
    fluxApparentRateNow = computeHydrogenFluxApparentRate( ...
        hydrogenTransport, timeIterationStep, initialHydrogenConcentration, ...
        inletVelocity, flowConfig, gridHyPHM, surfaceArea{1}(timeIterationStep));
    interfaceTstRateNow = NaN;

    if usePhreeqc
        if strcmpi(phreeqcRateLaw, 'tst_match')
            interfaceAreaNow = phreeqcGeometryState.interface_area_cm2(:);
            activeRateArea = interfaceAreaNow > 0 & isfinite(phreeqcCalciteRateData(:));
            Rate(timeIterationStep+1) = sum(phreeqcCalciteRateData(activeRateArea) ...
                .* interfaceAreaNow(activeRateArea), 'omitnan') ...
                ./ max(sum(interfaceAreaNow(activeRateArea), 'omitnan'), eps);
            interfaceTstRateNow = Rate(timeIterationStep+1);
        else
            Rate(timeIterationStep+1) = mean(phreeqcCalciteRateData, 'omitnan');
        end
    else
        Rate(timeIterationStep+1) = fluxApparentRateNow;
        legacyRateGeometryState = computePhreeqcGeometryState( ...
            gridHyPHM, levelSet{1}(:, timeIterationStep), molarVolume, thickness);
        legacyInterfaceRateData = ComputeLegacyEquivalentTstTriangleRate( ...
            hydrogenData, legacyRateGeometryState.interface_area_cm2, ...
            gridHyPHM.V0T, VertexTriMatrix, gridHyPHM.numV, rateCoefficientTST);
        legacyInterfaceArea = legacyRateGeometryState.interface_area_cm2(:);
        activeLegacyRateArea = legacyInterfaceArea > 0 & isfinite(legacyInterfaceRateData(:));
        interfaceTstRateNow = sum(legacyInterfaceRateData(activeLegacyRateArea) ...
            .* legacyInterfaceArea(activeLegacyRateArea), 'omitnan') ...
            ./ max(sum(legacyInterfaceArea(activeLegacyRateArea), 'omitnan'), eps);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---------------- 每步输出并记录全局量 ----------------
    t_now = currentTime;                    % 当前时间（s）
    por_now = Porosity(timeIterationStep+1);
    perm_now = Permeability(timeIterationStep+1); % 已经是 mD
    rate_now = Rate(timeIterationStep+1);

    % surfaceArea 存储为 cell，确保索引存在
    if numel(surfaceArea) >= 1 && numel(surfaceArea{1}) >= (timeIterationStep+1)
        sa_now = surfaceArea{1}(timeIterationStep+1);
    else
        sa_now = NaN;
    end

    vol_now = Volume(timeIterationStep+1);

    % 计算 PV
    inletFlux = inletVelocity * flowConfig.crossSectionLength * thickness; % [cm^3/s]
    cumulativeVolume = inletFlux * t_now; % [cm^3]
    % InitialPoreVolume 是面积 [cm^2]，乘以 thickness [cm] 得到体积
    currentPV = cumulativeVolume / (InitialPoreVolume * thickness);
    PV_injected(timeIterationStep+1) = currentPV;
    doExportStep = false;

    %% 在每个时间步骤输出结果的图（保持原结构）
    tic;
    fig1 = figure('Visible', 'off', 'Position', [100, 100, 1200, 1200]);
    
    % === 可视化参数设置 ===
    fontSize = 20;           % 坐标轴标签和刻度字体大小
    titleFontSize = 24;      % 标题字体大小
    fontName = 'Helvetica';  % 字体类型
    
    % 使用 tiledlayout 确保三个子图大小一致
    t_layout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % 添加总标题显示时间信息
    title(t_layout, sprintf('Simulation at {\\itt} = %.2f s', currentTime), ...
          'FontSize', titleFontSize + 2, 'FontName', fontName, 'FontWeight', 'normal');
    
    % 获取浓度数据（三角形 P0 数据），绘图时保持单元常数并屏蔽固体/界面单元。
    hydrogenPlotTitle = 'Hydrogen concentration {\itc} (mol·cm^{-3})';
    hydrogenPlotShortTitle = 'Hydrogen concentration {\itc}';
    hydrogenConcentrationData = hydrogenTransport.U.getdata(timeIterationStep);
    hydrogenUseAdaptiveColorScale = false;
    if phreeqcTransportFirstSplit && isfield(phreeqcSpeciesData, 'h_activity_mol_cm3') && ...
            ~isempty(phreeqcSpeciesData.h_activity_mol_cm3)
        hydrogenConcentrationData = phreeqcSpeciesData.h_activity_mol_cm3(:);
        hydrogenPlotTitle = 'H+ activity {\ita_H} (mol·cm^{-3})';
        hydrogenPlotShortTitle = 'H+ activity {\ita_H}';
        hydrogenUseAdaptiveColorScale = true;
    end
    concentrationFaceData = PrepareConcentrationFaceData( ...
        hydrogenConcentrationData, gridHyPHM.V0T, levelSet{1}(:, timeIterationStep + 1));
    hydrogenColorLimits = [0, initialHydrogenConcentration];
    hydrogenColorbarLabel = '';
    if hydrogenUseAdaptiveColorScale
        finiteHydrogenFaceData = concentrationFaceData(isfinite(concentrationFaceData) & concentrationFaceData >= 0);
        if ~isempty(finiteHydrogenFaceData)
            adaptiveUpperLimit = max(finiteHydrogenFaceData);
        else
            adaptiveUpperLimit = initialHydrogenConcentration;
        end
        if ~isfinite(adaptiveUpperLimit) || adaptiveUpperLimit <= 0
            adaptiveUpperLimit = max(initialHydrogenConcentration, eps);
        end
        hydrogenColorLimits = [0, adaptiveUpperLimit];
        hydrogenColorbarLabel = sprintf('adaptive [0, %.3g]', adaptiveUpperLimit);
    end
    
    % 绘制浓度分布
    ax1 = nexttile(1);
    patch(ax1, 'Faces', gridHyPHM.V0T, ...
        'Vertices', [gridHyPHM.coordV, zeros(gridHyPHM.numV, 1)], ...
        'FaceVertexCData', concentrationFaceData, ...
        'FaceColor', 'flat', 'EdgeColor', 'none');
    title(hydrogenPlotTitle, ...
          'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    set(gca, 'FontSize', fontSize, 'FontName', fontName, 'XTickLabel', []);
    view(2); axis equal tight;
    caxis(hydrogenColorLimits);
    % colormap(ax1, 'plasma');  % 使用Plasma颜色映射
    % 获取子图位置并设置colorbar与子图y轴对齐
    ax1Pos = ax1.Position;
    cb1 = colorbar; cb1.FontSize = fontSize; cb1.FontName = fontName;
    if ~isempty(hydrogenColorbarLabel)
        cb1.Label.String = hydrogenColorbarLabel;
        cb1.Label.FontSize = max(10, fontSize - 4);
        cb1.Label.FontName = fontName;
    end
    cb1.Position = [cb1.Position(1), ax1Pos(2), cb1.Position(3), ax1Pos(4)];  % colorbar上下对齐子图
    cbXPos = cb1.Position(1);  % 记录colorbar的X位置
    
    % 叠加初始界面轮廓线（始终显示在最上层）
    hold(ax1, 'on');
    % 设置较高的Z值，确保轮廓线覆盖在浓度云图上方
    zBase = max(concentrationFaceData, [], 'omitnan');
    if ~isfinite(zBase) || zBase <= 0
        zBase = 1;
    end
    zOffsetBlack = zBase * 1.08;
    zOffsetWhite = zBase * 1.12;

    % 当前（实时）界面轮廓：白色，覆盖在黑色初始轮廓之上
    phi_nowContour = reshape(levelSet{1}(:, timeIterationStep+1), nxNodes, nyNodes)';
    C_nowContour = contourc(xlin, ylin, phi_nowContour, [0 0]);
    nowInterfaceSegments = {};
    kNow = 1;
    while kNow < size(C_nowContour, 2)
        nNow = C_nowContour(2, kNow);
        nowInterfaceSegments{end+1} = C_nowContour(:, kNow + 1 : kNow + nNow); %#ok<SAGROW>
        kNow = kNow + nNow + 1;
    end

    contourHandles = [];  % 存储轮廓线句柄
    for segIdx = 1:numel(initInterfaceSegments)
        seg = initInterfaceSegments{segIdx};
        % 黑色：初始界面参考轮廓
        h = plot3(ax1, seg(1, :), seg(2, :), ones(size(seg(1, :))) * zOffsetBlack, ...
                  'k--', 'LineWidth', 1.5, 'Clipping', 'off');
        contourHandles = [contourHandles; h]; %#ok<AGROW>
    end
    whiteContourHandles = [];  % 白色实时轮廓句柄
    for segIdx = 1:numel(nowInterfaceSegments)
        seg = nowInterfaceSegments{segIdx};
        hW = plot3(ax1, seg(1, :), seg(2, :), ones(size(seg(1, :))) * zOffsetWhite, ...
                   'w-', 'LineWidth', 2.0, 'Clipping', 'off');
        whiteContourHandles = [whiteContourHandles; hW]; %#ok<AGROW>
    end

    % 将白色实时轮廓移到最上层
    uistack(contourHandles, 'top');
    if ~isempty(whiteContourHandles)
        uistack(whiteContourHandles, 'top');
    end
    hold(ax1, 'off');
    
    % 绘制水平集界面
    nexttile(2);
    [X,Y] = meshgrid(linspace(0, lengthXAxis, dxfResolutionX), linspace(0, lengthYAxis, dxfResolutionY));
    F = scatteredInterpolant(microscaleGrid.coordinates(:,1), ...
                             microscaleGrid.coordinates(:,2), ...
                             levelSet{1}(:, timeIterationStep+1), 'natural');
    Z = F(X,Y);
    solidMask = Z >= 0;
    poreMask  = Z <= 0;
    % ExportBinaryMaskToDXF(X, Y, solidfile = fullfile(solidDXFDir, sprintf('solid_t%04d.dxf', timeIterationStep)), 'SOLID'); % ...existing code...
    % ExportBinaryMaskToDXF(X, Y, poreMask,  fullfile(poreDXFDir,  sprintf('pore_t%04d.dxf',  timeIterationStep)), 'PORE');
    % contourf(X, Y, Z, [0,0], 'LineWidth', 2);
    % hold on;
    % contour(X, Y, Z, [0,0], 'k', 'LineWidth', 2);
    % title(['Mineral Interface at t = ', num2str(currentTime), ' s']);
    % 使用细网格 Z 直接估算孔隙率（孔隙：Z < 0）
    % 使用细网格 Z 直接估算孔隙率（孔隙：Z <= 0）
    porosityVal = sum(Z(:) <= 0) / numel(Z); % 无量纲 [0,1]
    % 存储孔隙率（与 timeSteps 对齐）
    Porosity(timeIterationStep+1) = porosityVal;

    
    % 控制台输出（格式化）
    fprintf('---- Time step %d | t = %.2f s ----\n', timeIterationStep, t_now);
    fprintf(' Porosity = %.4f  | Permeability = %g mD | PV = %.4f\n', porosityVal, perm_now, currentPV);
    fprintf(' Avg dissolution rate = %g  | Surface area = %g cm^2  | Volume = %g cm^3\n', ...
            rate_now, sa_now, vol_now);
    fprintf(' Tortuosity (tau) = %.4f\n', tau_now);
    fprintf(' Segment Tortuosity [1..5] = [%.4f, %.4f, %.4f, %.4f, %.4f]\n', ...
        TortuositySegments(1,timeIterationStep+1), TortuositySegments(2,timeIterationStep+1), ...
        TortuositySegments(3,timeIterationStep+1), TortuositySegments(4,timeIterationStep+1), ...
        TortuositySegments(5,timeIterationStep+1));

    % 更新 por_now 使用刚计算的 porosityVal
    por_now = porosityVal;

    if targetExportEnabled
        [doExportStep, nextTargetExportIndex] = ShouldExportTargetDissolutionSlice( ...
            timeIterationStep, porosityVal, targetExportInitialPorosity, ...
            targetSliceSettings.targetSliceCount, targetExportLastIndex, stopAfterThisStep);
        if doExportStep
            targetExportLastIndex = nextTargetExportIndex;
            targetExportCount = targetExportCount + 1;
        end
    else
        doExportStep = shouldExportStep(timeIterationStep, transportStepper.numsteps, ...
            exportEvery, stopAfterThisStep);
    end

    % 计算 k/k0 渗透率比值
    if ~isnan(k0) && k0 > 0
        k_k0_now = perm_now / k0;
    else
        k_k0_now = NaN;
    end
    outletHConc_now = OutletHConc(timeIterationStep+1);

    stabilityDiagnostics = computeStabilityDiagnostics( ...
        gridHyPHM, levelSet{1}(:, timeIterationStep + 1), hydrogenData, ...
        hydrogenDataCurrent, flow.getdata(1), StokesL.U.getdata(1), ...
        diffusionCoefficient, meshMinLengthForCfl, macroscaleTimeStepSize, ...
        SorceScale, rateCoefficientTST, thickness, initialHydrogenConcentration, ...
        inletVelocity, flowConfig, normalSpeedMax, por_now, perm_now, ...
        k_k0_now, hydrogenTransport.isUpwind);
    nextHardCflStepSize = computeHardCflStepSize(stabilityDiagnostics, ...
        meshMinLengthForCfl, advectiveCflSafety, reactiveCflSafety, hardCflMinTimeStep);

    fprintf([' Stability: c=[%.3e, %.3e], mass_res=%.3e mol (rel %.3e), ', ...
        'u_max=%.3e cm/s, Pe_max=%.3g, CFL_adv=%.3g, CFL_react=%.3g, pore_components=%d, flag=%s\n'], ...
        stabilityDiagnostics.cMin, stabilityDiagnostics.cMax, ...
        stabilityDiagnostics.massBalanceResidual, stabilityDiagnostics.massBalanceRelative, ...
        stabilityDiagnostics.localVelocityMax, stabilityDiagnostics.localPeMax, ...
        stabilityDiagnostics.advectiveCfl, stabilityDiagnostics.reactionCfl, ...
        stabilityDiagnostics.poreComponents, stabilityDiagnostics.flag);

    fidStab = fopen(stabilityDiagnosticsCsvFile, 'a');
    if fidStab ~= -1
        fprintf(fidStab, ['%d,%.6f,%.6f,%.12e,%.12e,%.12e,%.12e,%.12e,', ...
            '%.12e,%.12e,%.12e,%.12e,%.12e,%d,%.12e,%.12e,%.12e,%.12e,"%s","%s"\n'], ...
            timeIterationStep, t_now, macroscaleTimeStepSize, ...
            stabilityDiagnostics.cMin, stabilityDiagnostics.cMax, ...
            stabilityDiagnostics.totalMass, stabilityDiagnostics.massBalanceResidual, ...
            stabilityDiagnostics.massBalanceRelative, stabilityDiagnostics.localVelocityMax, ...
            stabilityDiagnostics.localPeMax, stabilityDiagnostics.advectiveCfl, ...
            stabilityDiagnostics.reactionCfl, stabilityDiagnostics.normalSpeedMax, ...
            stabilityDiagnostics.poreComponents, stabilityDiagnostics.largestPoreComponentFraction, ...
            por_now, perm_now, k_k0_now, ...
            char(hydrogenTransport.isUpwind), stabilityDiagnostics.flag);
        fclose(fidStab);
    else
        warning('Cannot open stability diagnostics log file for writing: %s', stabilityDiagnosticsCsvFile);
    end
    
    % 追加写入 CSV（每步一行）
    fid = fopen(logfile, 'a');
    if fid ~= -1
        fprintf(fid, ['%d,%.6f,%.6f,%.6f,%.6f,%.6e,%.6e,%.6e,', ...
            '%.6e,%.6e,%.6f,%.6e,%.6f\n'], ...
            timeIterationStep, t_now, por_now, perm_now, k_k0_now, ...
            rate_now, interfaceTstRateNow, fluxApparentRateNow, ...
            sa_now, vol_now, currentPV, outletHConc_now, tau_now);
        fclose(fid);
    else
        warning('Cannot open log file for writing: %s', logfile);
    end
    
    % 追加写入 Excel
    if writeExcel
        try
            T_row = table(timeIterationStep, t_now, por_now, perm_now, k_k0_now, ...
                rate_now, interfaceTstRateNow, fluxApparentRateNow, ...
                sa_now, vol_now, currentPV, outletHConc_now, tau_now, ...
                'VariableNames', {'TimeStep', 'Time_s', 'Porosity', 'Permeability_mD', ...
                'k_k0', 'DissolutionRate', 'InterfaceTstRate_mol_cm2_s', ...
                'FluxApparentRate_mol_cm2_s', 'SurfaceArea_cm2', 'GrainVolume_cm3', ...
                'InjectedPV', 'OutletHConc', 'Tortuosity'});
            writetable(T_row, xlsxFile, 'WriteMode', 'append');
        catch ME
            warning('MATLAB:PNMExcelWrite', 'Failed to write to Excel: %s', ME.message);
        end
    end

    % 追加写入分段迂曲度 CSV + Excel（额外表）
    tauSeg = TortuositySegments(:, timeIterationStep+1);
    fidSeg = fopen(tortuositySegCsvFile, 'a');
    if fidSeg ~= -1
        fprintf(fidSeg, '%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
            timeIterationStep, t_now, tau_now, tauSeg(1), tauSeg(2), tauSeg(3), tauSeg(4), tauSeg(5));
        fclose(fidSeg);
    else
        warning('Cannot open segment tortuosity log file for writing: %s', tortuositySegCsvFile);
    end

    if writeExcel
        try
            T_seg = table(timeIterationStep, t_now, tau_now, tauSeg(1), tauSeg(2), tauSeg(3), tauSeg(4), tauSeg(5), ...
                'VariableNames', {'TimeStep', 'Time_s', 'Tortuosity_All', 'Tortuosity_Seg1', 'Tortuosity_Seg2', 'Tortuosity_Seg3', 'Tortuosity_Seg4', 'Tortuosity_Seg5'});
            writetable(T_seg, tortuositySegXlsxFile, 'WriteMode', 'append');
        catch ME
            warning('MATLAB:PNMSegmentExcelWrite', 'Failed to write segment tortuosity to Excel: %s', ME.message);
        end
    end
    % -----------------------------------------------------

    %  绘图并在标题中显示孔隙率百分比
    contourf(X, Y, Z, [0,0], 'LineWidth', 2);
    % colormap(gca, 'plasma');  % 使用Plasma颜色映射
    hold on;
    contour(X, Y, Z, [0,0], 'k', 'LineWidth', 2);
    title(sprintf('Mineral interface — Porosity: %.2f%%', porosityVal*100), ...
          'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    ylh = ylabel('{\itY} (cm)', 'FontSize', fontSize, 'FontName', fontName);
    ylh.Position(1) = ylh.Position(1) - 0.007;  % 向左移动Y轴标签
    set(gca, 'FontSize', fontSize, 'FontName', fontName, 'XTickLabel', []);
    axis equal tight;

    % 绘制流场流线
    ax3 = nexttile(3);
    % 获取当前时间步的速度场数据
    velocityData = StokesL.U.getdata(1);
    
    % 确保速度场数据格式正确
    if size(velocityData, 2) == 1
        numComponents = length(velocityData) / (gridHyPHM.numV + gridHyPHM.numE);
        if numComponents == 2
            vx = velocityData(1:gridHyPHM.numV+gridHyPHM.numE);
            vy = velocityData(gridHyPHM.numV+gridHyPHM.numE+1:end);
        else
            error('Unexpected velocity data format');
        end
    else
        vx = velocityData(:,1);
        vy = velocityData(:,2);
    end
    
    % 获取P2节点的坐标（包括顶点和边中点）
    if isprop(gridHyPHM, 'baryE')
        p2Coords = [gridHyPHM.coordV; gridHyPHM.baryE];
    else
        % 如果没有baryE属性，计算边中点
        baryE = zeros(gridHyPHM.numE, 2);
        for i = 1:gridHyPHM.numE
            v1 = gridHyPHM.coordV(gridHyPHM.V0E(i,1), :);
            v2 = gridHyPHM.coordV(gridHyPHM.V0E(i,2), :);
            baryE(i,:) = (v1 + v2) / 2;
        end
        p2Coords = [gridHyPHM.coordV; baryE];
    end
    
    % 确保坐标和速度分量长度匹配
    if length(vx) ~= size(p2Coords, 1) || length(vy) ~= size(p2Coords, 1)
        error('Coordinate and velocity vector lengths do not match');
    end
    
    % 创建规则网格用于流线图
    canshu = 200; % 调整此参数以改变网格密度
    xGrid = linspace(0, lengthXAxis, canshu);
    yGrid = linspace(0, lengthYAxis, canshu/2);
    [XGrid, YGrid] = meshgrid(xGrid, yGrid);
    
    % 插值速度场到规则网格
    VX = griddata(p2Coords(:,1), p2Coords(:,2), vx, XGrid, YGrid);
    VY = griddata(p2Coords(:,1), p2Coords(:,2), vy, XGrid, YGrid);
    
    % 计算流速大小
    speed = sqrt(VX.^2 + VY.^2);
    speedValid = speed(~isnan(speed));
    low = prctile(speedValid, 1);
    high = prctile(speedValid, 99);
    pad = 0.1 * (high - low + eps);
    speedCLim = [max(0, low - pad), high + pad];
    % speedCLim = [0, 0.8]; 
    % 绘制流线和流速大小
    imagesc(xGrid, yGrid, speed); % 显示流速大小
    set(gca, 'YDir', 'normal');              % y轴向上
    caxis(speedCLim); % 固定颜色条范围
    hold on;
    streamslice(XGrid, YGrid, VX, VY, 3); % 流线密度为3
    axis equal tight;
    colormap(jet);
    % colormap(ax3, 'plasma');  % 使用Plasma颜色映射
    % 获取子图位置并设置colorbar与子图y轴对齐
    ax3Pos = ax3.Position;
    cb3 = colorbar; cb3.FontSize = fontSize; cb3.FontName = fontName;
    cb3OffsetY = 0.02;  % 向上偏移量，可调整
    cb3.Position = [cbXPos, ax3Pos(2) + cb3OffsetY, cb3.Position(3), ax3Pos(4)];  % colorbar上下对齐子图，X位置与第一个一致
    title(['Flow velocity field {\bf\itu} (cm·s^{-1})'], ...
          'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    xlabel('{\itX} (cm)', 'FontSize', fontSize, 'FontName', fontName);
    set(gca, 'FontSize', fontSize, 'FontName', fontName);
    
    % 添加矿物表面轮廓
    velContourHandles = [];
    for segIdx = 1:numel(initInterfaceSegments)
        seg = initInterfaceSegments{segIdx};
        hVel = plot(ax3, seg(1, :), seg(2, :), 'k--', 'LineWidth', 1.5, 'Clipping', 'off');
        velContourHandles = [velContourHandles; hVel]; %#ok<AGROW>
    end
    velWhiteContourHandles = [];
    for segIdx = 1:numel(nowInterfaceSegments)
        seg = nowInterfaceSegments{segIdx};
        hVelW = plot(ax3, seg(1, :), seg(2, :), 'w-', 'LineWidth', 2.0, 'Clipping', 'off');
        velWhiteContourHandles = [velWhiteContourHandles; hVelW]; %#ok<AGROW>
    end
    if ~isempty(velContourHandles)
        uistack(velContourHandles, 'top');
    end
    if ~isempty(velWhiteContourHandles)
        uistack(velWhiteContourHandles, 'top');
    end
    
    % 保存图像（300 dpi）
    if doExportStep && saveMainPlot
        safePrintPng(fig1, fullfile(resultsDir, sprintf('timestep_%04d.png', timeIterationStep)), ...
            sprintf('main_timestep_png (step=%d)', timeIterationStep));
        % 保存可编辑的 .fig 文件（含颜色映射、文字等全部样式，可在 MATLAB 中重新打开编辑）
        if saveFigureFiles
            safeSaveFig(fig1, fullfile(resultsDir, sprintf('timestep_%04d.fig', timeIterationStep)), ...
                sprintf('main_timestep_fig (step=%d)', timeIterationStep));
        end
    end
    safeCloseFigure(fig1);
    
    %% === 单独保存每个子图（高清版本）===
    if doExportStep && saveIndividualPlots
    % 子图1：浓度分布
    figSub1 = figure('Visible', 'off', 'Position', [100, 100, 800, 700]);
    axSub1 = axes(figSub1);
    patch(axSub1, 'Faces', gridHyPHM.V0T, ...
        'Vertices', [gridHyPHM.coordV, zeros(gridHyPHM.numV, 1)], ...
        'FaceVertexCData', concentrationFaceData, ...
        'FaceColor', 'flat', 'EdgeColor', 'none');
    titleHandle1 = title(sprintf('%s at {\\itt} = %.2f s (mol·cm^{-3})', hydrogenPlotShortTitle, currentTime), ...
          'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    titleHandle1.Position(2) = titleHandle1.Position(2) * 1.35;  % 增加标题距离
    titleHandle1.Position(1) = titleHandle1.Position(1) - 0.01;  % 向左移动标题
    xlabel('X (cm)', 'FontSize', fontSize, 'FontName', fontName);
    ylabel('Y (cm)', 'FontSize', fontSize, 'FontName', fontName);
    set(gca, 'FontSize', fontSize, 'FontName', fontName);
    view(2); axis equal tight;
    caxis(hydrogenColorLimits);
    cb1 = colorbar; cb1.FontSize = fontSize; cb1.FontName = fontName;
    if ~isempty(hydrogenColorbarLabel)
        cb1.Label.String = hydrogenColorbarLabel;
        cb1.Label.FontSize = max(10, fontSize - 4);
        cb1.Label.FontName = fontName;
    end
    % colormap(figSub1, 'viridis'); % 使用viridis映射
    % 叠加初始界面轮廓线
    hold on;
    zBase = max(concentrationFaceData, [], 'omitnan');
    if ~isfinite(zBase) || zBase <= 0
        zBase = 1;
    end
    zOffsetBlack = zBase * 1.08;
    zOffsetWhite = zBase * 1.12;

    % 实时演化界面（白色）
    phi_nowContour = reshape(levelSet{1}(:, timeIterationStep+1), nxNodes, nyNodes)';
    C_nowContour = contourc(xlin, ylin, phi_nowContour, [0 0]);
    nowInterfaceSegments = {};
    kNow = 1;
    while kNow < size(C_nowContour, 2)
        nNow = C_nowContour(2, kNow);
        nowInterfaceSegments{end+1} = C_nowContour(:, kNow + 1 : kNow + nNow); %#ok<SAGROW>
        kNow = kNow + nNow + 1;
    end

    for segIdx = 1:numel(initInterfaceSegments)
        seg = initInterfaceSegments{segIdx};
        plot3(seg(1, :), seg(2, :), ones(size(seg(1, :))) * zOffsetBlack, ...
              'k--', 'LineWidth', 1.5, 'Clipping', 'off');
    end
    for segIdx = 1:numel(nowInterfaceSegments)
        seg = nowInterfaceSegments{segIdx};
        plot3(seg(1, :), seg(2, :), ones(size(seg(1, :))) * zOffsetWhite, ...
              'w-', 'LineWidth', 2.0, 'Clipping', 'off');
    end
    hold off;
    safePrintPng(figSub1, fullfile(subplot1Dir, sprintf('concentration_%04d.png', timeIterationStep)), ...
        sprintf('subplot_concentration_png (step=%d)', timeIterationStep));
    if saveFigureFiles
        safeSaveFig(figSub1, fullfile(subplot1Dir, sprintf('concentration_%04d.fig', timeIterationStep)), ...
            sprintf('subplot_concentration_fig (step=%d)', timeIterationStep));
    end
    safeCloseFigure(figSub1);

    % 子图2：矿物界面
    figSub2 = figure('Visible', 'off', 'Position', [100, 100, 800, 700]);
    contourf(X, Y, Z, [0,0], 'LineWidth', 2);
    hold on;
    contour(X, Y, Z, [0,0], 'k', 'LineWidth', 2);
    titleHandle2 = title(sprintf('Mineral interface at {\\itt} = %.2f s — Porosity: %.2f%%', currentTime, porosityVal*100), ...
          'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    titleHandle2.Position(2) = titleHandle2.Position(2) * 1.05;  % 增加标题距离
    xlabel('X (cm)', 'FontSize', fontSize, 'FontName', fontName);
    ylabel('Y (cm)', 'FontSize', fontSize, 'FontName', fontName);
    set(gca, 'FontSize', fontSize, 'FontName', fontName);
    axis equal tight;
    hold off;
    safePrintPng(figSub2, fullfile(subplot2Dir, sprintf('interface_%04d.png', timeIterationStep)), ...
        sprintf('subplot_interface_png (step=%d)', timeIterationStep));
    if saveFigureFiles
        safeSaveFig(figSub2, fullfile(subplot2Dir, sprintf('interface_%04d.fig', timeIterationStep)), ...
            sprintf('subplot_interface_fig (step=%d)', timeIterationStep));
    end
    safeCloseFigure(figSub2);
    
    % 子图3：流场流线
    figSub3 = figure('Visible', 'off', 'Position', [100, 100, 800, 700]);
    imagesc(xGrid, yGrid, speed);
    set(gca, 'YDir', 'normal');
    caxis(speedCLim);
    hold on;
    streamslice(XGrid, YGrid, VX, VY, 3);
    axis equal tight;
    colormap(jet);
    cb3 = colorbar; cb3.FontSize = fontSize; cb3.FontName = fontName;
    titleHandle3 = title(sprintf('Flow velocity field at {\\itt} = %.2f s ({\\bf\\itu}, cm·s^{-1})', currentTime), ...
          'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    titleHandle3.Position(2) = titleHandle3.Position(2) * 1.05;  % 增加标题距离
    xlabel('X (cm)', 'FontSize', fontSize, 'FontName', fontName);
    ylabel('Y (cm)', 'FontSize', fontSize, 'FontName', fontName);
    set(gca, 'FontSize', fontSize, 'FontName', fontName);

    % 当前矿物界面（白色）在速度图中的分段（与浓度图一致）
    phi_nowContour_vel = reshape(levelSet{1}(:, timeIterationStep+1), nxNodes, nyNodes)';
    C_nowContour_vel = contourc(xlin, ylin, phi_nowContour_vel, [0 0]);
    nowInterfaceSegmentsVel = {};
    kNowVel = 1;
    while kNowVel < size(C_nowContour_vel, 2)
        nNowVel = C_nowContour_vel(2, kNowVel);
        nowInterfaceSegmentsVel{end+1} = C_nowContour_vel(:, kNowVel + 1 : kNowVel + nNowVel); %#ok<SAGROW>
        kNowVel = kNowVel + nNowVel + 1;
    end

    % 添加矿物表面轮廓：初始黑色虚线 + 当前白色实线
    velSubContourHandles = [];
    for segIdx = 1:numel(initInterfaceSegments)
        seg = initInterfaceSegments{segIdx};
        hSub = plot(seg(1, :), seg(2, :), 'k--', 'LineWidth', 1.5, 'Clipping', 'off');
        velSubContourHandles = [velSubContourHandles; hSub]; %#ok<AGROW>
    end
    velSubWhiteContourHandles = [];
    for segIdx = 1:numel(nowInterfaceSegmentsVel)
        seg = nowInterfaceSegmentsVel{segIdx};
        hSubW = plot(seg(1, :), seg(2, :), 'w-', 'LineWidth', 2.0, 'Clipping', 'off');
        velSubWhiteContourHandles = [velSubWhiteContourHandles; hSubW]; %#ok<AGROW>
    end
    if ~isempty(velSubContourHandles)
        uistack(velSubContourHandles, 'top');
    end
    if ~isempty(velSubWhiteContourHandles)
        uistack(velSubWhiteContourHandles, 'top');
    end
    hold off;
    safePrintPng(figSub3, fullfile(subplot3Dir, sprintf('velocity_%04d.png', timeIterationStep)), ...
        sprintf('subplot_velocity_png (step=%d)', timeIterationStep));
    if saveFigureFiles
        safeSaveFig(figSub3, fullfile(subplot3Dir, sprintf('velocity_%04d.fig', timeIterationStep)), ...
            sprintf('subplot_velocity_fig (step=%d)', timeIterationStep));
    end
    safeCloseFigure(figSub3);
    end

    if usePhreeqc && (mod(timeIterationStep, phreeqcExportEvery) == 0 || doExportStep)
        WritePhreeqcSpeciesTable(phreeqcOutputDir, timeIterationStep, t_now, gridHyPHM, phreeqcSpeciesData);
        if saveIndividualPlots
            ExportPhreeqcSpeciesPlots(phreeqcSpeciesPlotDir, timeIterationStep, t_now, ...
                gridHyPHM, levelSet{1}(:, timeIterationStep + 1), phreeqcSpeciesData, ...
                struct('fontName', fontName, 'fontSize', fontSize, 'titleFontSize', titleFontSize));
        end
    end


    % 新建图窗口,绘制固液边界
    if doExportStep && (exportDXF || saveInterfaceMask || enableNMRSurrogate || enablePNGSimulation)
    fig1 = figure('Visible', 'off');  % 可改为 'on' 调试查看

    [X, Y] = meshgrid(linspace(0, lengthXAxis, dxfResolutionX), linspace(0, lengthYAxis, dxfResolutionY));
    F = scatteredInterpolant(microscaleGrid.coordinates(:,1), ...
                            microscaleGrid.coordinates(:,2), ...
                            levelSet{1}(:, timeIterationStep+1), 'natural');
    Z = F(X,Y);

    % 生成掩膜
    solidMask = Z >= 0;  % 固体区域
    poreMask  = Z < 0;   % 孔隙区域（液体）

    % 导出DXF
    if exportDXF
        solidDxfPath = fullfile(solidDXFDir, sprintf('solid_t%04d.dxf', timeIterationStep));
        poreDxfPath = fullfile(poreDXFDir, sprintf('pore_t%04d.dxf', timeIterationStep));
        ExportBinaryMaskToDXF(X, Y, solidMask, solidDxfPath, 'SOLID');
        ExportBinaryMaskToDXF(X, Y, poreMask, poreDxfPath, 'PORE');

        if enableNMRSimulation
            nmrCalibrationFactor = runSynchronousNMRStep( ...
                nmrConfig, ...
                poreDxfPath, ...
                solidDxfPath, ...
                lengthXAxis, ...
                lengthYAxis, ...
                timeIterationStep, ...
                t_now, ...
                porosityVal, ...
                nmrComsolOutputDir, ...
                nmrInversionOutputDir, ...
                nmrSyncLogFile, ...
                nmrCalibrationFactor);
        end
    end

    if enableNMRSurrogate
        nmrSurrogateCalibrationFactor = runSynchronousNMRSurrogateStep( ...
            nmrSurrogateConfig, ...
            solidMask, ...
            poreMask, ...
            lengthXAxis, ...
            lengthYAxis, ...
            timeIterationStep, ...
            t_now, ...
            porosityVal, ...
            nmrSurrogateOutputDir, ...
            nmrSurrogateInversionOutputDir, ...
            nmrSurrogateMaskDir, ...
            nmrSurrogateSyncLogFile, ...
            nmrSurrogateCalibrationFactor);
    end

    % 绘制彩色图像：固体=黄色，液体=红色
    maskImage = zeros([size(solidMask), 3]);
    maskImage(:,:,1) = poreMask | solidMask;    % 红通道：两者都有红色成分
    maskImage(:,:,2) = solidMask;               % 绿通道：仅固体有 -> 黄色（红+绿）
    % 蓝通道保持为 0

    image(linspace(0, lengthXAxis, size(maskImage,2)), ...
        linspace(0, lengthYAxis, size(maskImage,1)), ...
        maskImage);

    set(gca, 'YDir', 'normal');
    axis equal tight off;
    box off;

    % 保存图像
    if saveInterfaceMask
        interfacePngPath = fullfile(imageDir, sprintf('timestep_%04d.png', timeIterationStep));
        safeSaveAs(fig1, interfacePngPath, ...
            sprintf('interface_mask_png (step=%d)', timeIterationStep));
    end
    if enablePNGSimulation
        interfacePngPath = fullfile(imageDir, sprintf('timestep_%04d.png', timeIterationStep));
        runSynchronousPNGNMRStep( ...
            pngNMRConfig, ...
            pngNMRRuntimeConfigFile, ...
            resultsDir, ...
            interfacePngPath, ...
            timeIterationStep, ...
            t_now, ...
            porosityVal, ...
            pngNMRSyncLogFile);
    end
    safeCloseFigure(fig1);
    % 输出进度信息
    disp(['    Results saved for timestep ', num2str(timeIterationStep), ...
          ' in ', num2str(toc), ' seconds']);
    end
    
    % %% === 新增：导出骨架为 TIF 二值图像（推荐放在上面那段绘图代码的最后）===
    % tifSkeletonDir = fullfile(resultsDir, 'skeleton_tif');
    % if ~exist(tifSkeletonDir, 'dir')
    %     mkdir(tifSkeletonDir);
    % end
    
    % % 目标分辨率：建议和原始 TIF 完全一致，这样像素级对齐
    % targetNx = nx_img;   % 原始图像的列数（x方向像素数）
    % targetNy = ny_img;   % 原始图像的行数（y方向像素数）
    
    % % 构造与原始 TIF 完全对齐的规则网格（像素中心对齐）
    % xEdges = linspace(0, lengthXAxis, targetNx+1);
    % yEdges = linspace(0, lengthYAxis, targetNy+1);
    % [xCent, yCent] = meshgrid( ...
    %     (xEdges(1:end-1) + xEdges(2:end))/2, ...
    %     (yEdges(1:end-1) + yEdges(2:end))/2 );
    
    % % 插值当前 level set 到这个精确像素网格
    % F = scatteredInterpolant(microscaleGrid.coordinates(:,1), ...
    %                          microscaleGrid.coordinates(:,2), ...
    %                          levelSet{1}(:, timeIterationStep+1), ...
    %                          'linear', 'none');   % 界面附近可能会出现 NaN，用 none 更安全
    
    % Z_fine = F(xCent, yCent);
    
    % % 二值化：固体 = 255（白色），孔隙 = 0（黑色）
    % % 注意：你的原始 TIF 是 255=固体，所以保持一致
    % BW_solid_tif = uint8( Z_fine >= 0 ) * 255;
    
    % % 重要：MATLAB 的 imagesc/imwrite 默认 y 方向是"上→下"，而你的 flipud 已经在前面做了
    % % 这里我们已经用了 bottom-up 的 yCent，所以不需要再 flipud
    % % 但为了 100% 保险，再显式翻转一次与原始输入完全一致
    % BW_solid_tif = flipud(BW_solid_tif);
    
    % % 保存为 8-bit TIF
    % tifFilename = fullfile(tifSkeletonDir, sprintf('skeleton_%04d.tif', timeIterationStep));
    % imwrite(BW_solid_tif, tifFilename, 'tiff', ...
    %         'Compression', 'none', ...    % 无损
    %         'Resolution', [1/pixelSizeCm, 1/pixelSizeCm]);  % 单位：像素/cm → dpi ≈ 12700 for 2.02 μm/px
    
    % disp(['    Skeleton TIF saved: ', tifFilename]);
    % %% === 结束：TIF 导出 ===

    % ...existing code...
    % disp(['    Skeleton TIF saved: ', tifFilename]);
    % %% === 结束：TIF 导出 ===

    %% === 新增：实时绘制并保存全局变化曲线 ===
    if saveRealtimePlot && doExportStep
    % 使用不可见 figure 避免弹窗干扰，每次覆盖保存图片以查看实时进度
    rtFig = figure('Visible', 'off', 'Position', [100, 100, 1000, 1800]);
    
    % 截取当前已计算的时间步数据
    currentRange = 1:(timeIterationStep + 1);
    t_hist = transportStepper.timepts(currentRange);
    
    subplot(7,1,1);
    plot(t_hist, Rate(currentRange), 'b.-');
    xlabel('Time [s]'); ylabel('Avg. Dissolution Rate');
    title(['Dissolution Rate (t = ' num2str(currentTime, '%.1f') ' s)']);
    grid on;

    subplot(7,1,2);
    plot(t_hist, Porosity(currentRange), 'k.-');
    xlabel('Time [s]'); ylabel('Porosity [-]');
    title('Porosity Evolution');
    grid on;

    subplot(7,1,3);
    permData = Permeability(currentRange);
    plot(t_hist, permData, 'm.-');
    xlabel('Time [s]'); ylabel('Permeability [mD]');
    title('Permeability Evolution');
    if any(permData > 0) && all(~isnan(permData))
        try
            set(gca, 'YScale', 'log');
        catch
        end
    end
    grid on;

    subplot(7,1,4);
    plot(t_hist, surfaceArea{1}(currentRange), 'r.-');
    xlabel('Time [s]'); ylabel('Surface Area [cm^2]');
    title('Surface Area Evolution');
    grid on;

    subplot(7,1,5);
    plot(t_hist, Volume(currentRange), 'g.-');
    xlabel('Time [s]'); ylabel('Grain Volume [cm^3]');
    title('Grain Volume Evolution');
    grid on;

    subplot(7,1,6);
    plot(t_hist, PV_injected(currentRange), 'c.-');
    xlabel('Time [s]'); ylabel('Injected PV');
    title('Injected Pore Volumes');
    grid on;

    subplot(7,1,7);
    hold on;
    plot(t_hist, Tortuosity(currentRange), 'k.-', 'LineWidth', 1.5, 'DisplayName', '\tau all');
    plot(t_hist, TortuositySegments(1,currentRange), 'r.-', 'DisplayName', '\tau seg1');
    plot(t_hist, TortuositySegments(2,currentRange), 'g.-', 'DisplayName', '\tau seg2');
    plot(t_hist, TortuositySegments(3,currentRange), 'b.-', 'DisplayName', '\tau seg3');
    plot(t_hist, TortuositySegments(4,currentRange), 'm.-', 'DisplayName', '\tau seg4');
    plot(t_hist, TortuositySegments(5,currentRange), 'c.-', 'DisplayName', '\tau seg5');
    xlabel('Time [s]'); ylabel('Tortuosity \tau [-]');
    title('Tortuosity Evolution (all + 5 segments)');
    legend('Location', 'best');
    hold off;
    grid on;

    % 保存图片（覆盖更新名为 global_evolution_realtime.png 的文件）
    safeSaveAs(rtFig, fullfile(resultsDir, 'global_evolution_realtime.png'), ...
        sprintf('global_evolution_realtime (step=%d)', timeIterationStep));
    safeCloseFigure(rtFig);
    end

    elapsedWallSeconds = toc(simulationTimer);
    if elapsedWallSeconds > maxExperimentWallSeconds
        timedOut = true;
        timeoutMessage = sprintf(['Experiment exceeded maxExperimentWallSeconds=%.0f s ' ...
            'after step %d; partial outputs were kept and the batch will continue.'], ...
            maxExperimentWallSeconds, timeIterationStep);
        fprintf('\n========================================\n');
        fprintf('>>> %s\n', timeoutMessage);
        fprintf('========================================\n\n');
        break;
    end

    if ~stopAfterThisStep
        hasPlannedNextStep = timeIterationStep < transportStepper.numsteps;
        if hasPlannedNextStep
            plannedNextStepSize = transportStepper.timepts(timeIterationStep + 2) - ...
                transportStepper.timepts(timeIterationStep + 1);
        else
            plannedNextStepSize = macroscaleTimeStepSize;
        end
        nextMacroscaleStepSize = plannedNextStepSize;
        didUpdateNextStep = false;
        porosityDeltaForStep = NaN;

        if useAdaptivePorositySteps
        if isnan(adaptivePreviousPorosity)
            porosityDeltaForStep = NaN;
        else
            porosityDeltaForStep = max(0, porosityVal - adaptivePreviousPorosity);
        end

            nextMacroscaleStepSize = adaptivePorosityStepSize( ...
            macroscaleTimeStepSize, porosityDeltaForStep, porosityStepTarget, ...
            porosityStepTolerance, porosityStepUpperFactor, adaptiveGrowthFactor, ...
            adaptiveShrinkSafety, adaptiveMinTimeStep, adaptiveMaxTimeStep);
            didUpdateNextStep = true;
        end

        if targetExportEnabled && useAdaptivePorositySteps
            [nextMacroscaleStepSize, wasTargetBoundaryLimited, nextTargetBoundaryPorosity] = ...
                LimitTargetDissolutionStepSize( ...
                nextMacroscaleStepSize, macroscaleTimeStepSize, porosityDeltaForStep, ...
                porosityVal, targetExportInitialPorosity, targetExportLastIndex, ...
                targetSliceSettings.targetSliceCount, adaptiveMinTimeStep, ...
                adaptiveMaxTimeStep, targetSliceStepBoundaryFactor);
            if wasTargetBoundaryLimited
                didUpdateNextStep = true;
                fprintf(' target_slice_dt: limiting next dt to %g s for porosity boundary %.5f\n', ...
                    nextMacroscaleStepSize, nextTargetBoundaryPorosity);
            end
        end

        if enableHardCflLimit && isfinite(nextHardCflStepSize) && nextHardCflStepSize > 0
            nextMacroscaleStepSize = min(nextMacroscaleStepSize, nextHardCflStepSize);
            didUpdateNextStep = true;
        end

        % Concentration-based time step reduction (inspired by GeoChemFoam setDeltaT.H + deltaEpsMax.H)
        if enableConcentrationCflLimit && isfield(stabilityDiagnostics, 'cMax')
            [nextMacroscaleStepSize, wasConcentrationLimited, overshootRatio, activeConcentrationThreshold] = ...
                ResolveConcentrationLimitedStep( ...
                nextMacroscaleStepSize, stabilityDiagnostics.cMax, ...
                initialHydrogenConcentration, concentrationOvershootThreshold, ...
                concentrationHardOvershootThreshold, concentrationShrinkFactor, ...
                adaptiveMinTimeStep, targetExportEnabled);
            if wasConcentrationLimited
                didUpdateNextStep = true;
                fprintf(' [ConcCFL] Overshoot ratio=%.3f > %.3f, reducing next dt to %g s\n', ...
                    overshootRatio, activeConcentrationThreshold, nextMacroscaleStepSize);
            elseif targetExportEnabled && overshootRatio > concentrationOvershootThreshold
                fprintf(' [ConcCFL] Mild overshoot ratio=%.3f ignored in target-slice mode (hard threshold %.3f)\n', ...
                    overshootRatio, activeConcentrationThreshold);
            end
        end

        remainingTime = endTime - currentTime;
        if remainingTime > EPS
            nextMacroscaleStepSize = min(nextMacroscaleStepSize, remainingTime);
        end

        shouldApplyNextStep = didUpdateNextStep && ...
            (abs(nextMacroscaleStepSize - plannedNextStepSize) > eps(plannedNextStepSize) || ...
            (~hasPlannedNextStep && allowAdaptiveTimeGridExtension && remainingTime > EPS));
        if shouldApplyNextStep
            if hasPlannedNextStep
                transportStepper.setTimeStepSize(timeIterationStep + 1, nextMacroscaleStepSize);
            elseif allowAdaptiveTimeGridExtension && remainingTime > EPS
                transportStepper.appendTimeStepSize(nextMacroscaleStepSize);
                extendTransportVariablesToStepper(hydrogenTransport);
                newNumTimePoints = numel(transportStepper.timepts);
                for iSlice = 1:numberOfSlices
                    levelSet{iSlice}(:, newNumTimePoints) = NaN;
                    surfaceArea{iSlice}(newNumTimePoints) = NaN;
                end
                Rate(newNumTimePoints) = NaN;
                Volume(newNumTimePoints) = NaN;
                Porosity(newNumTimePoints) = NaN;
                Permeability(newNumTimePoints) = NaN;
                PV_injected(newNumTimePoints) = NaN;
                OutletHConc(newNumTimePoints) = NaN;
                Tortuosity(newNumTimePoints) = NaN;
                TortuositySegments(:, newNumTimePoints) = NaN;
                fprintf(' adaptive_porosity: extended internal time grid to %d steps.\n', ...
                    transportStepper.numsteps);
            end
            timeSteps = transportStepper.timepts(:)';
            fprintf(' dt update: dPorosity=%s, next dt=%g s (CFL cap=%g s)\n', ...
                formatAdaptiveDelta(porosityDeltaForStep), nextMacroscaleStepSize, nextHardCflStepSize);
        end
    end
    adaptivePreviousPorosity = porosityVal;
    
    % === 检查是否需要在完成当前步骤后停止 ===
    if stopAfterThisStep
        fprintf('\n========================================\n');
        fprintf('>>> Final dissolution images exported.\n');
        fprintf('>>> Stopping simulation at t = %.2f s (step %d)\n', currentTime, timeIterationStep);
        fprintf('========================================\n\n');
        break;  % 跳出 while 循环
    end
    
end % while

% 获取实际完成的时间步数（处理提前终止的情况）
timeSteps = transportStepper.timepts(:)';
nValidSteps = length(surfaceArea{1});
validTimeSteps = timeSteps(1:nValidSteps);
validRate = Rate(1:nValidSteps);
validVolume = Volume(1:nValidSteps);
validTortuosity = Tortuosity(1:nValidSteps);
validTortuositySegments = TortuositySegments(:, 1:nValidSteps);
validPorosity = Porosity(1:nValidSteps);
validPermeability = Permeability(1:nValidSteps);
validPV = PV_injected(1:nValidSteps);

if showDebugFigures
    figure
    plot(validTimeSteps, validRate);
    xlabel('time [s]')
    ylabel('average dissolution rate')

    figure
    plot(validTimeSteps, surfaceArea{1}(:));
    xlabel('time [s]')
    ylabel('surface Area')

    figure
    plot(validTimeSteps, validVolume);
    xlabel('time [s]')
    ylabel('Grain Volume')
end

if saveFinalPlot
    plotEvolutionSummary(validTimeSteps, validRate, validPorosity, validPermeability, ...
        surfaceArea{1}(:), validVolume, validPV, validTortuosity, ...
        validTortuositySegments, fullfile(resultsDir, 'global_evolution_with_porosity_permeability.png'));
end

initialPermCandidates = validPermeability(~isnan(validPermeability) & validPermeability > 0);
if isempty(initialPermCandidates)
    initialPermeability = NaN;
else
    initialPermeability = initialPermCandidates(1);
end

if isnan(initialPermeability)
    PBTimeStep = NaN;
    PBTime = NaN;
else
    breakthroughIdx = find(validPermeability >= initialPermeability * permeabilityRatioThreshold, 1, 'first');
    if isempty(breakthroughIdx)
        PBTimeStep = NaN;
        PBTime = NaN;
    else
        PBTimeStep = breakthroughIdx - 1;
        PBTime = validTimeSteps(breakthroughIdx);
    end
end

result = struct();
result.resultsDir = resultsDir;
result.runName = runName;
result.metadata = metadata;
result.PBTimeStep = PBTimeStep;
result.PBTime = PBTime;
result.finalTimeStep = nValidSteps - 1;
result.finalTime = validTimeSteps(end);
result.initialPermeability = initialPermeability;
result.finalPermeability = validPermeability(end);
result.finalTortuosity = validTortuosity(end);
result.finalPorosity = validPorosity(end);
result.timedOut = timedOut;
result.timeoutMessage = timeoutMessage;
result.exportedDissolutionSlices = targetExportCount;

metadata.final = struct( ...
    'PBTimeStep', PBTimeStep, ...
    'PBTime_s', PBTime, ...
    'finalTime_s', result.finalTime, ...
    'initialPermeability_mD', initialPermeability, ...
    'finalPermeability_mD', result.finalPermeability, ...
    'finalPorosity', result.finalPorosity, ...
    'finalTortuosity', result.finalTortuosity, ...
    'exportedDissolutionSlices', targetExportCount, ...
    'runtimeWallSeconds', toc(simulationTimer), ...
    'timedOut', timedOut, ...
    'timeoutMessage', timeoutMessage);
writeJsonFile(fullfile(resultsDir, 'run_metadata.json'), metadata);

disp('All results (incl. porosity, permeability & tortuosity) saved in results folder');
end

%% Helper functions
function por = porosityHelperFun(t, x, porosityCells, numberOfSlices, timeStep)
por = 0;
for i = 1:numberOfSlices
    por = por + porosityCells{1}(timeStep);
end
end

function diff = diffusionHelperFun(t, x, diffusionCells, numberOfSlices, timeStep)
diff = zeros(2);
for i = 1:numberOfSlices
    diff = diff + reshape(diffusionCells{1}(:, timeStep), 2, 2);
end
end

function area = areaHelperFun(t, x, surfaceAreaCells, numberOfSlices, timeStep)
area = surfaceAreaCells{findSlice(x, numberOfSlices)}(timeStep);
end

function sliceNumber = findSlice(x, numberOfSlices)
sliceNumber = 1;
end

function out = getInterfacePoints(Levels)
% 注意：保留原函数结构，但这里未被调用
% 2*128 是原代码中 x 向分块数的写法，这里如果需要用请替换为实际 nxParts
nxPartsX = 2*128;
out = find((Levels < 0).*(Levels > -1 / (10 * nxPartsX)));
end

function safePrintPng(figHandle, outputPath, contextTag)
try
    print(figHandle, outputPath, '-dpng', '-r300');
catch ME
    warning('MATLAB:FigureExport', 'PNG导出失败 [%s]: %s | file=%s', contextTag, ME.message, outputPath);
end
end

function safeSaveFig(figHandle, outputPath, contextTag)
try
    savefig(figHandle, outputPath);
catch ME
    warning('MATLAB:FigureExport', 'FIG导出失败 [%s]: %s | file=%s', contextTag, ME.message, outputPath);
end
end

function safeSaveAs(figHandle, outputPath, contextTag)
try
    saveas(figHandle, outputPath);
catch ME
    warning('MATLAB:FigureExport', 'saveas导出失败 [%s]: %s | file=%s', contextTag, ME.message, outputPath);
end
end

function safeCloseFigure(figHandle)
try
    if ~isempty(figHandle) && ishghandle(figHandle)
        close(figHandle);
    end
catch
end
end

function out = getClosePoints(Center, dist, grid)
out = find(sqrt((grid.coordV(:, 1)-grid.coordV(Center, 1)).^2 + (grid.coordV(:, 2) - grid.coordV(Center, 2)).^2) < dist);
end

% 新增：计算随机分布中颗粒间的平均表面间距（基于最近邻）
function avgSpacing = calculateAverageSpacing(circleCenters, circleRadius)
% 计算每个颗粒与其最近邻颗粒之间的平均孔喉间距
if size(circleCenters,1) < 2
    avgSpacing = NaN;
    return;
end

n = size(circleCenters,1);
minSpacings = zeros(n, 1);

% 对每个颗粒，找到其最近邻的孔喉间距
for i = 1:n
    minDist = inf;
    for j = 1:n
        if i ~= j
            centerDist = norm(circleCenters(i,:) - circleCenters(j,:));
            surfaceSpacing = centerDist - 2 * circleRadius; % 表面间距
            if surfaceSpacing < minDist && surfaceSpacing > 0
                minDist = surfaceSpacing;
            end
        end
    end
    minSpacings(i) = minDist;
end

% 过滤掉无效值（如果有孤立颗粒）
validSpacings = minSpacings(minSpacings < inf & minSpacings > 0);

if isempty(validSpacings)
    avgSpacing = NaN;
else
    avgSpacing = mean(validSpacings);
end
end

function nmrConfig = loadNMRWorkflowConfig(rtmDir, userConfig)
addpath(rtmDir);
if exist('NMRSimulationConfig', 'file') == 2
    nmrConfig = NMRSimulationConfig();
else
    nmrConfig = struct('nmr_method', 'none', 'enablePNGSimulation', false, 'method', 'png_mesh');
end

if isstruct(userConfig) && isfield(userConfig, 'pngNMRConfig') && isstruct(userConfig.pngNMRConfig)
    nmrConfig = mergeConfigStructs(nmrConfig, userConfig.pngNMRConfig);
end

% Allow direct PNM_beauty3(config) callers to override the central file
% without constructing cfg.pngNMRConfig by hand.
nmrConfig = mergeTopLevelNMRFields(nmrConfig, userConfig);
if isstruct(userConfig) && ~isfield(userConfig, 'nmr_method') && ~isfield(userConfig, 'nmrMethod')
    if logical(cfgget(userConfig, 'enableNMRSimulation', false))
        nmrConfig.nmr_method = 'comsol';
    elseif logical(cfgget(userConfig, 'enableNMRSurrogate', false))
        nmrConfig.nmr_method = 'surrogate';
    elseif logical(cfgget(userConfig, 'enablePNGSimulation', false))
        nmrConfig.nmr_method = char(cfgget(userConfig, 'pngNMRMethod', cfgget(userConfig, 'pngMethod', cfgget(userConfig, 'method', 'png_mesh'))));
    elseif any(isfield(userConfig, {'enableNMRSimulation', 'enableNMRSurrogate', 'enablePNGSimulation'})) && ...
            ~logical(cfgget(userConfig, 'enableNMRSimulation', false)) && ...
            ~logical(cfgget(userConfig, 'enableNMRSurrogate', false)) && ...
            ~logical(cfgget(userConfig, 'enablePNGSimulation', false))
        nmrConfig.nmr_method = 'none';
    end
end
end

function out = mergeTopLevelNMRFields(base, userConfig)
out = base;
if ~isstruct(userConfig)
    return;
end
fields = { ...
    'nmr_method', 'nmrMethod', 'enableNMRSimulation', 'enableNMRSurrogate', ...
    'enablePNGSimulation', 'pngNMRMethod', 'pngMethod', ...
    'comsol_path', 'mph_file', 'comsol_python_exe', 'enable_comsol', 'enable_inversion', ...
    'overwrite_existing', 'scale_factor', 'export_mph', 'inversion_script', ...
    'nmrSurrogateModelPath', 'nmrSurrogateRoot', 'nmrSurrogatePythonExe', ...
    'nmrSurrogateDatasetPath', 'nmrSurrogateResolution', 'nmrSurrogateDevice'};
for iField = 1:numel(fields)
    fieldName = fields{iField};
    if isfield(userConfig, fieldName) && ~isempty(userConfig.(fieldName))
        out.(fieldName) = userConfig.(fieldName);
    end
end
end

function out = mergeConfigStructs(base, extra)
out = base;
if ~isstruct(extra)
    return;
end
names = fieldnames(extra);
for iName = 1:numel(names)
    name = names{iName};
    if ~isempty(extra.(name))
        out.(name) = extra.(name);
    end
end
end

function [pngConfig, outputDir, syncLogFile, runtimeConfigFile] = initializePNGNMRSync(rtmDir, reactiveRoot, resultsDir, pngConfig)
driverPath = fullfile(rtmDir, 'png_nmr_driver.py');
if ~exist(driverPath, 'file')
    error('PNG NMR driver 不存在: %s', driverPath);
end

pythonExe = char(cfgget(pngConfig, 'python_exe', 'python'));
advancedToolsDir = char(cfgget(pngConfig, 'advanced_tools_dir', ''));
if isempty(advancedToolsDir) || ~exist(advancedToolsDir, 'dir')
    error('PNG NMR advanced_tools_dir 不存在: %s', advancedToolsDir);
end

defaultT2Path = fullfile(reactiveRoot, 'T2_process');
t2ProcessPath = char(cfgget(pngConfig, 't2_process_path', defaultT2Path));
if isempty(t2ProcessPath) || ~exist(t2ProcessPath, 'dir')
    error('PNG NMR t2_process_path 不存在: %s', t2ProcessPath);
end

method = char(cfgget(pngConfig, 'method', 'png_mesh'));
outputRoot = char(cfgget(pngConfig, 'output_root', 'png_nmr_results'));
if isAbsolutePath(outputRoot)
    outputDir = outputRoot;
else
    outputDir = fullfile(resultsDir, outputRoot);
end
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

pngConfig.python_exe = pythonExe;
pngConfig.driver_path = driverPath;
pngConfig.advanced_tools_dir = advancedToolsDir;
pngConfig.t2_process_path = t2ProcessPath;
pngConfig.method = method;
pngConfig.output_root = outputDir;
runtimeConfigFile = fullfile(resultsDir, 'png_nmr_runtime_config.json');
writeJsonFile(runtimeConfigFile, pngConfig);

syncLogFile = fullfile(resultsDir, 'png_nmr_sync_log.csv');
if ~exist(syncLogFile, 'file')
    fid = fopen(syncLogFile, 'w');
    if fid ~= -1
        fprintf(fid, 'timestep,time_s,porosity,success,method,decay_csv,spectrum_mat,manifest,message\n');
        fclose(fid);
    else
        warning('MATLAB:PNMPNGNMRSync', '无法创建PNG NMR同步日志: %s', syncLogFile);
    end
end

fprintf('[PNG NMR] 已启用。方法: %s\n', method);
fprintf('[PNG NMR] Python: %s\n', pythonExe);
fprintf('[PNG NMR] advanced_tools: %s\n', advancedToolsDir);
fprintf('[PNG NMR] 输出: %s\n', outputDir);
end

function runSynchronousPNGNMRStep(pngConfig, runtimeConfigFile, resultsDir, interfacePngPath, ...
    timestep, timeSeconds, porosityValue, syncLogFile)

method = char(cfgget(pngConfig, 'method', 'png_mesh'));
success = false;
message = "OK";
decayCsv = "";
spectrumMat = "";
manifestPath = "";

try
    if ~exist(interfacePngPath, 'file')
        message = "interface PNG未生成";
        appendPNGNMRSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, success, ...
            method, decayCsv, spectrumMat, manifestPath, message);
        warning('MATLAB:PNMPNGNMRSync', 'PNG NMR跳过时间步 %04d: PNG不存在: %s', timestep, interfacePngPath);
        return;
    end

    fprintf('[PNG NMR] 时间步 %04d: %s + T2反演\n', timestep, method);
    cmd = sprintf('%s %s --exp-dir %s --png %s --timestep %d --time-s %.17g --porosity %.17g --config %s', ...
        quoteSystemArg(pngConfig.python_exe), ...
        quoteSystemArg(pngConfig.driver_path), ...
        quoteSystemArg(resultsDir), ...
        quoteSystemArg(interfacePngPath), ...
        timestep, ...
        timeSeconds, ...
        porosityValue, ...
        quoteSystemArg(runtimeConfigFile));

    [status, cmdout] = system(cmd);
    if ~isempty(strtrim(cmdout))
        fprintf('%s\n', strtrim(cmdout));
    end
    result = parseResultJson(cmdout);
    success = status == 0 && isfield(result, 'success') && logical(result.success);
    if isfield(result, 'decay') && isfield(result.decay, 'curve_csv')
        decayCsv = string(result.decay.curve_csv);
    end
    if isfield(result, 'inversion') && isfield(result.inversion, 'spectrum_mat')
        spectrumMat = string(result.inversion.spectrum_mat);
    end
    if isfield(result, 'manifest_path')
        manifestPath = string(result.manifest_path);
    end
    if ~success
        if isfield(result, 'error')
            message = string(result.error);
        else
            message = sprintf('PNG NMR退出码: %d', status);
        end
        warning('MATLAB:PNMPNGNMRSync', 'PNG NMR时间步 %04d 失败: %s', timestep, char(message));
        if logical(cfgget(pngConfig, 'fail_on_error', false))
            error('MATLAB:PNMPNGNMRSync', '%s', char(message));
        end
    end
catch ME
    message = string(ME.message);
    warning('MATLAB:PNMPNGNMRSync', 'PNG NMR时间步 %04d 异常: %s', timestep, ME.message);
    if logical(cfgget(pngConfig, 'fail_on_error', false))
        rethrow(ME);
    end
end

appendPNGNMRSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, success, ...
    method, decayCsv, spectrumMat, manifestPath, message);
end

function appendPNGNMRSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, success, ...
    method, decayCsv, spectrumMat, manifestPath, message)
fid = fopen(syncLogFile, 'a');
if fid == -1
    warning('MATLAB:PNMPNGNMRSync', '无法写入PNG NMR同步日志: %s', syncLogFile);
    return;
end
cleaner = onCleanup(@() fclose(fid));

method = strrep(char(method), '"', '""');
decayCsv = strrep(char(decayCsv), '"', '""');
spectrumMat = strrep(char(spectrumMat), '"', '""');
manifestPath = strrep(char(manifestPath), '"', '""');
message = strrep(char(message), '"', '""');
fprintf(fid, '%d,%.6f,%.8f,%d,"%s","%s","%s","%s","%s"\n', ...
    timestep, timeSeconds, porosityValue, logical(success), method, decayCsv, ...
    spectrumMat, manifestPath, message);
clear cleaner;
end

function tf = isAbsolutePath(pathValue)
pathValue = char(pathValue);
tf = numel(pathValue) >= 2 && pathValue(2) == ':';
if ~tf && startsWith(pathValue, filesep)
    tf = true;
end
end

function [nmrConfig, comsolOutputDir, inversionOutputDir, syncLogFile] = initializeNMRSync(reactiveRoot, resultsDir, userConfig)
automationDir = fullfile(reactiveRoot, 'automation');
if ~exist(automationDir, 'dir')
    error('NMR自动化目录不存在: %s', automationDir);
end
addpath(automationDir);

nmrConfig = AutomationConfig();
nmrConfig = applyAutomationConfigOverrides(nmrConfig, userConfig);
comsolOutputDir = fullfile(resultsDir, 'comsol_results');
inversionOutputDir = fullfile(resultsDir, 'inversion_results');
if ~exist(comsolOutputDir, 'dir'); mkdir(comsolOutputDir); end
if ~exist(inversionOutputDir, 'dir'); mkdir(inversionOutputDir); end

syncLogFile = fullfile(resultsDir, 'nmr_sync_log.csv');
if ~exist(syncLogFile, 'file')
    fid = fopen(syncLogFile, 'w');
    if fid ~= -1
        fprintf(fid, 'timestep,time_s,porosity,comsol_success,inversion_success,total_water,raw_spectrum_sum,calibration_factor,excel_output,message\n');
        fclose(fid);
    else
        warning('MATLAB:PNMNMRSync', '无法创建NMR同步日志: %s', syncLogFile);
    end
end

fprintf('[NMR同步] 已启用。参数来自: %s\n', fullfile(automationDir, 'AutomationConfig.m'));
fprintf('[NMR同步] 覆盖配置来自: %s\n', fullfile(fileparts(reactiveRoot), 'ReactiveTransport', 'RTM', 'NMRSimulationConfig.m'));
fprintf('[NMR同步] COMSOL输出: %s\n', comsolOutputDir);
fprintf('[NMR同步] 反演输出: %s\n', inversionOutputDir);
end

function configObj = applyAutomationConfigOverrides(configObj, userConfig)
if nargin < 2 || ~isstruct(userConfig)
    return;
end

configObj = setAutomationProperty(configObj, 'comsol_path', cfgget(userConfig, 'comsol_path', []));
configObj = setAutomationProperty(configObj, 'mph_file', cfgget(userConfig, 'mph_file', []));
configObj = setAutomationProperty(configObj, 'python_exe', cfgget(userConfig, 'comsol_python_exe', []));
configObj = setAutomationProperty(configObj, 'inversion_script', cfgget(userConfig, 'inversion_script', []));
configObj = setAutomationProperty(configObj, 'enable_comsol', cfgget(userConfig, 'enable_comsol', []));
configObj = setAutomationProperty(configObj, 'enable_inversion', cfgget(userConfig, 'enable_inversion', []));
configObj = setAutomationProperty(configObj, 'overwrite_existing', cfgget(userConfig, 'overwrite_existing', []));
configObj = setAutomationProperty(configObj, 'scale_factor', cfgget(userConfig, 'scale_factor', []));
configObj = setAutomationProperty(configObj, 'export_mph', cfgget(userConfig, 'export_mph', []));
configObj = setAutomationProperty(configObj, 'max_samples_per_folder', cfgget(userConfig, 'max_samples_per_folder', []));
end

function configObj = setAutomationProperty(configObj, propertyName, value)
if isempty(value)
    return;
end
try
    if isprop(configObj, propertyName)
        configObj.(propertyName) = value;
    end
catch ME
    warning('MATLAB:NMRConfigOverride', '无法设置 AutomationConfig.%s: %s', propertyName, ME.message);
end
end

function calibrationFactor = runSynchronousNMRStep(nmrConfig, poreDxfPath, solidDxfPath, ...
    lengthXAxis, lengthYAxis, timestep, timeSeconds, porosityValue, comsolOutputDir, ...
    inversionOutputDir, syncLogFile, calibrationFactor)

excelOutput = fullfile(comsolOutputDir, sprintf('T2_t%04d.xlsx', timestep));
comsolSuccess = false;
inversionSuccess = false;
totalWater = NaN;
rawSpectrumSum = NaN;
message = "OK";

try
    fprintf('[NMR同步] 时间步 %04d: COMSOL + T2反演\n', timestep);
    comsolSuccess = run_comsol_processing( ...
        nmrConfig.mph_file, ...
        poreDxfPath, ...
        solidDxfPath, ...
        lengthXAxis, ...
        lengthYAxis, ...
        excelOutput, ...
        nmrConfig);

    if ~comsolSuccess
        message = "COMSOL处理失败";
        appendNMRSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, comsolSuccess, ...
            inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message);
        warning('MATLAB:PNMNMRSync', 'NMR同步跳过时间步 %04d: COMSOL处理失败。', timestep);
        return;
    end

    if ~exist(excelOutput, 'file')
        message = "COMSOL Excel未生成";
        appendNMRSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, comsolSuccess, ...
            inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message);
        warning('MATLAB:PNMNMRSync', 'NMR同步跳过时间步 %04d: Excel未生成: %s', timestep, excelOutput);
        return;
    end

    shouldCalibrate = isempty(calibrationFactor);
    [inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor] = run_python_inversion( ...
        excelOutput, inversionOutputDir, nmrConfig, calibrationFactor);

    if inversionSuccess && shouldCalibrate && isfinite(porosityValue) && ...
            isfinite(rawSpectrumSum) && rawSpectrumSum > 0
        calibrationFactor = porosityValue / rawSpectrumSum;
        fprintf('[NMR同步] 基于时间步 %04d 孔隙率 %.6f 计算校准因子: %.6e\n', ...
            timestep, porosityValue, calibrationFactor);
        [inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor] = run_python_inversion( ...
            excelOutput, inversionOutputDir, nmrConfig, calibrationFactor);
    end

    if ~inversionSuccess
        message = "T2反演失败";
        warning('MATLAB:PNMNMRSync', 'NMR同步时间步 %04d 的T2反演失败。', timestep);
    end
catch ME
    message = string(ME.message);
    warning('MATLAB:PNMNMRSync', 'NMR同步时间步 %04d 失败: %s', timestep, ME.message);
end

appendNMRSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, comsolSuccess, ...
    inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message);
end

function [surrogateConfig, surrogateOutputDir, inversionOutputDir, maskDir, syncLogFile] = initializeNMRSurrogateSync(reactiveRoot, resultsDir, userConfig)
automationDir = fullfile(reactiveRoot, 'automation');
if ~exist(automationDir, 'dir')
    error('NMR自动化目录不存在: %s', automationDir);
end
addpath(automationDir);

inversionConfig = AutomationConfig();
inversionConfig = applyAutomationConfigOverrides(inversionConfig, userConfig);
defaultModelPath = 'C:\Users\imgw\Documents\Codex\NMR-agent\runs\IMGW_256_300_20260507-130311_3a583275\latest_model.pt';
modelPath = char(cfgget(userConfig, 'nmrSurrogateModelPath', defaultModelPath));
if isempty(modelPath) || ~exist(modelPath, 'file')
    error('NMR替代模型文件不存在: %s', modelPath);
end

scriptPath = fullfile(automationDir, 'run_nmr_surrogate_prediction.py');
if ~exist(scriptPath, 'file')
    error('NMR替代模型推理脚本不存在: %s', scriptPath);
end

nmrAgentRoot = char(cfgget(userConfig, 'nmrSurrogateRoot', inferNmrAgentRootFromModel(modelPath)));
datasetPath = char(cfgget(userConfig, 'nmrSurrogateDatasetPath', ''));
pythonExe = char(cfgget(userConfig, 'nmrSurrogatePythonExe', ''));
if isempty(pythonExe)
    candidatePython = fullfile(nmrAgentRoot, '.venv', 'Scripts', 'python.exe');
    if exist(candidatePython, 'file')
        pythonExe = candidatePython;
    elseif isprop(inversionConfig, 'python_exe') && exist(inversionConfig.python_exe, 'file')
        pythonExe = inversionConfig.python_exe;
    else
        pythonExe = 'python';
    end
end

surrogateConfig = struct();
surrogateConfig.inversionConfig = inversionConfig;
surrogateConfig.pythonExe = pythonExe;
surrogateConfig.scriptPath = scriptPath;
surrogateConfig.modelPath = modelPath;
surrogateConfig.datasetPath = datasetPath;
surrogateConfig.nmrAgentRoot = nmrAgentRoot;
surrogateConfig.device = char(cfgget(userConfig, 'nmrSurrogateDevice', 'auto'));
surrogateConfig.resolution = cfgget(userConfig, 'nmrSurrogateResolution', 256);

surrogateOutputDir = fullfile(resultsDir, 'surrogate_results');
inversionOutputDir = fullfile(resultsDir, 'surrogate_inversion_results');
maskDir = fullfile(surrogateOutputDir, 'interface_images');
if ~exist(surrogateOutputDir, 'dir'); mkdir(surrogateOutputDir); end
if ~exist(inversionOutputDir, 'dir'); mkdir(inversionOutputDir); end
if ~exist(maskDir, 'dir'); mkdir(maskDir); end

syncLogFile = fullfile(resultsDir, 'nmr_surrogate_sync_log.csv');
if ~exist(syncLogFile, 'file')
    fid = fopen(syncLogFile, 'w');
    if fid ~= -1
        fprintf(fid, 'timestep,time_s,porosity,surrogate_success,inversion_success,total_water,raw_spectrum_sum,calibration_factor,excel_output,message\n');
        fclose(fid);
    else
        warning('MATLAB:PNMNMRSurrogate', '无法创建NMR替代模型同步日志: %s', syncLogFile);
    end
end

fprintf('[NMR替代模型] 已启用。模型: %s\n', modelPath);
fprintf('[NMR替代模型] Python: %s\n', pythonExe);
fprintf('[NMR替代模型] 衰减曲线输出: %s\n', surrogateOutputDir);
fprintf('[NMR替代模型] 反演输出: %s\n', inversionOutputDir);
end

function calibrationFactor = runSynchronousNMRSurrogateStep(surrogateConfig, solidMask, poreMask, ...
    lengthXAxis, lengthYAxis, timestep, timeSeconds, porosityValue, surrogateOutputDir, ...
    inversionOutputDir, maskDir, syncLogFile, calibrationFactor)

excelOutput = fullfile(surrogateOutputDir, sprintf('T2_t%04d.xlsx', timestep));
csvOutput = fullfile(surrogateOutputDir, sprintf('T2_t%04d.csv', timestep));
interfaceImagePath = fullfile(maskDir, sprintf('interface_t%04d.png', timestep));
surrogateSuccess = false;
inversionSuccess = false;
totalWater = NaN;
rawSpectrumSum = NaN;
message = "OK";

try
    fprintf('[NMR替代模型] 时间步 %04d: U-Net预测 + T2反演\n', timestep);
    writeSurrogateInterfaceImage(interfaceImagePath, solidMask, poreMask);

    surrogateSuccess = runNMRSurrogatePrediction( ...
        surrogateConfig, ...
        interfaceImagePath, ...
        lengthXAxis, ...
        lengthYAxis, ...
        excelOutput, ...
        csvOutput);

    if ~surrogateSuccess
        message = "替代模型预测失败";
        appendNMRSurrogateSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, surrogateSuccess, ...
            inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message);
        warning('MATLAB:PNMNMRSurrogate', 'NMR替代模型跳过时间步 %04d: 预测失败。', timestep);
        return;
    end

    if ~exist(excelOutput, 'file')
        message = "替代模型Excel未生成";
        appendNMRSurrogateSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, surrogateSuccess, ...
            inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message);
        warning('MATLAB:PNMNMRSurrogate', 'NMR替代模型跳过时间步 %04d: Excel未生成: %s', timestep, excelOutput);
        return;
    end

    shouldCalibrate = isempty(calibrationFactor);
    [inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor] = run_python_inversion( ...
        excelOutput, inversionOutputDir, surrogateConfig.inversionConfig, calibrationFactor);

    if inversionSuccess && shouldCalibrate && isfinite(porosityValue) && ...
            isfinite(rawSpectrumSum) && rawSpectrumSum > 0
        calibrationFactor = porosityValue / rawSpectrumSum;
        fprintf('[NMR替代模型] 基于时间步 %04d 孔隙率 %.6f 计算校准因子: %.6e\n', ...
            timestep, porosityValue, calibrationFactor);
        [inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor] = run_python_inversion( ...
            excelOutput, inversionOutputDir, surrogateConfig.inversionConfig, calibrationFactor);
    end

    if ~inversionSuccess
        message = "T2反演失败";
        warning('MATLAB:PNMNMRSurrogate', 'NMR替代模型时间步 %04d 的T2反演失败。', timestep);
    end
catch ME
    message = string(ME.message);
    warning('MATLAB:PNMNMRSurrogate', 'NMR替代模型时间步 %04d 失败: %s', timestep, ME.message);
end

appendNMRSurrogateSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, surrogateSuccess, ...
    inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message);
end

function success = runNMRSurrogatePrediction(surrogateConfig, interfaceImagePath, ...
    lengthXAxis, lengthYAxis, excelOutput, csvOutput)
cmd = sprintf('%s %s --interface-image %s --length-x-axis %.17g --length-y-axis %.17g --model-path %s --output-excel %s --output-csv %s --resolution %d --device %s', ...
    quoteSystemArg(surrogateConfig.pythonExe), ...
    quoteSystemArg(surrogateConfig.scriptPath), ...
    quoteSystemArg(interfaceImagePath), ...
    lengthXAxis, ...
    lengthYAxis, ...
    quoteSystemArg(surrogateConfig.modelPath), ...
    quoteSystemArg(excelOutput), ...
    quoteSystemArg(csvOutput), ...
    surrogateConfig.resolution, ...
    quoteSystemArg(surrogateConfig.device));

if isfield(surrogateConfig, 'datasetPath') && ~isempty(surrogateConfig.datasetPath)
    cmd = sprintf('%s --dataset-path %s', cmd, quoteSystemArg(surrogateConfig.datasetPath));
end
if isfield(surrogateConfig, 'nmrAgentRoot') && ~isempty(surrogateConfig.nmrAgentRoot)
    cmd = sprintf('%s --nmr-agent-root %s', cmd, quoteSystemArg(surrogateConfig.nmrAgentRoot));
end

[status, cmdout] = system(cmd);
if ~isempty(strtrim(cmdout))
    fprintf('%s\n', strtrim(cmdout));
end

result = parseResultJson(cmdout);
success = status == 0 && isfield(result, 'success') && logical(result.success);
if ~success
    if isfield(result, 'error')
        fprintf('        x NMR替代模型预测失败: %s\n', result.error);
    else
        fprintf('        x NMR替代模型预测失败，退出码: %d\n', status);
    end
end
end

function appendNMRSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, comsolSuccess, ...
    inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message)
fid = fopen(syncLogFile, 'a');
if fid == -1
    warning('MATLAB:PNMNMRSync', '无法写入NMR同步日志: %s', syncLogFile);
    return;
end
cleaner = onCleanup(@() fclose(fid));

if isempty(calibrationFactor)
    calibrationValue = NaN;
else
    calibrationValue = calibrationFactor;
end

message = strrep(char(message), '"', '""');
excelOutput = strrep(char(excelOutput), '"', '""');
fprintf(fid, '%d,%.6f,%.8f,%d,%d,%.8g,%.8g,%.8g,"%s","%s"\n', ...
    timestep, timeSeconds, porosityValue, logical(comsolSuccess), logical(inversionSuccess), ...
    totalWater, rawSpectrumSum, calibrationValue, excelOutput, message);
clear cleaner;
end

function appendNMRSurrogateSyncLog(syncLogFile, timestep, timeSeconds, porosityValue, surrogateSuccess, ...
    inversionSuccess, totalWater, rawSpectrumSum, calibrationFactor, excelOutput, message)
fid = fopen(syncLogFile, 'a');
if fid == -1
    warning('MATLAB:PNMNMRSurrogate', '无法写入NMR替代模型同步日志: %s', syncLogFile);
    return;
end
cleaner = onCleanup(@() fclose(fid));

if isempty(calibrationFactor)
    calibrationValue = NaN;
else
    calibrationValue = calibrationFactor;
end

message = strrep(char(message), '"', '""');
excelOutput = strrep(char(excelOutput), '"', '""');
fprintf(fid, '%d,%.6f,%.8f,%d,%d,%.8g,%.8g,%.8g,"%s","%s"\n', ...
    timestep, timeSeconds, porosityValue, logical(surrogateSuccess), logical(inversionSuccess), ...
    totalWater, rawSpectrumSum, calibrationValue, excelOutput, message);
clear cleaner;
end

function writeSurrogateInterfaceImage(outputPath, solidMask, poreMask)
maskImage = zeros([size(solidMask), 3], 'uint8');
maskImage(:,:,1) = uint8((solidMask | poreMask) * 255);
maskImage(:,:,2) = uint8(solidMask * 255);
imwrite(maskImage, outputPath);
end

function nmrAgentRoot = inferNmrAgentRootFromModel(modelPath)
runDir = fileparts(modelPath);
runsDir = fileparts(runDir);
nmrAgentRoot = fileparts(runsDir);
end

function quoted = quoteSystemArg(value)
value = char(value);
value = strrep(value, '"', '\"');
quoted = ['"' value '"'];
end

function result = parseResultJson(cmdout)
result = struct();
marker = 'RESULT_JSON=';
lines = regexp(cmdout, '\r?\n', 'split');
jsonText = '';

for iLine = 1:length(lines)
    line = strtrim(lines{iLine});
    if startsWith(line, marker)
        jsonText = extractAfter(line, strlength(marker));
    end
end

if strlength(jsonText) == 0
    return;
end

try
    result = jsondecode(char(jsonText));
catch ME
    fprintf('        x JSON解析失败: %s\n', ME.message);
    result = struct();
end
end

function nextStepSize = adaptivePorosityStepSize(currentStepSize, porosityDelta, targetDelta, ...
    tolerance, upperFactor, growthFactor, shrinkSafety, minStepSize, maxStepSize)
targetDelta = max(eps, targetDelta);
tolerance = max(0, tolerance);
upperDelta = targetDelta * max(1, upperFactor);
growthFactor = max(1.0, growthFactor);
shrinkSafety = min(1.0, max(0.1, shrinkSafety));
lowerDelta = max(0, targetDelta - tolerance);

if isnan(porosityDelta)
    proposedStepSize = currentStepSize * growthFactor;
elseif porosityDelta <= eps
    proposedStepSize = currentStepSize * growthFactor;
elseif porosityDelta < lowerDelta
    targetGrowthFactor = max(1.0, targetDelta / porosityDelta);
    proposedStepSize = currentStepSize * min(growthFactor, targetGrowthFactor);
elseif porosityDelta <= upperDelta
    proposedStepSize = currentStepSize;
else
    shrinkFactor = max(0.05, shrinkSafety * upperDelta / porosityDelta);
    proposedStepSize = currentStepSize * min(1.0, shrinkFactor);
end

nextStepSize = min(maxStepSize, max(minStepSize, proposedStepSize));
end

function text = formatAdaptiveDelta(porosityDelta)
if isnan(porosityDelta)
    text = 'n/a';
else
    text = sprintf('%.4g', porosityDelta);
end
end

function extendTransportVariablesToStepper(transportProblem)
variableFields = {'Q', 'U', 'A', 'B', 'C', 'D', 'E', 'F', ...
    'uD', 'gF', 'uR', 'rob', 'L'};
for iField = 1:numel(variableFields)
    fieldName = variableFields{iField};
    if isprop(transportProblem, fieldName) && ~isempty(transportProblem.(fieldName))
        ExtendHyPHMVariableToStepper(transportProblem.(fieldName), 'copy_previous');
    end
end
end

function transportProblem = createPhreeqcTransport(grid, stepper, flow, name, initialConcentration, ...
    inletConcentration, diffusionCoefficient, inletVelocity, EPS, upwindMode)
concentration = Variable(grid, stepper, name, 'P0');
concentration.setdata(0, @(t, x) initialConcentration);
transportProblem = TransportLEVEL(grid, stepper, [name, ' Transport']);
transportProblem.id2N = {1, 2, 3};
transportProblem.id2F = {4};
transportProblem.U = concentration;
transportProblem.D.setdata(diffusionCoefficient * eye(2));
transportProblem.gF.setdata(@(t, x) -inletConcentration * inletVelocity * (x(1) < EPS));
transportProblem.A.setdata(0, @(t, x) 1);
transportProblem.C.setdata(0, flow.getdata(1));
transportProblem.isUpwind = upwindMode;
end

function setPhreeqcTransportStepData(transportProblem, timeStep, flowData)
transportProblem.A.setdata(timeStep, @(t, x) 1);
transportProblem.gF.setdata(timeStep, transportProblem.gF.getdata(0));
transportProblem.C.setdata(timeStep, flowData);
end

function geom = computePhreeqcGeometryState(grid, levels, molarVolume, thickness)
numTriangles = grid.numT;
interfaceArea = zeros(numTriangles, 1);
waterVolume = zeros(numTriangles, 1);
solidVolume = zeros(numTriangles, 1);

for kT = 1:numTriangles
    vertexIds = grid.V0T(kT, :);
    points = grid.coordV(vertexIds, :);
    values = levels(vertexIds);
    totalArea = grid.areaT(kT);
    waterArea = triangleAreaBelowLevelZero(points, values);
    solidArea = max(totalArea - waterArea, 0);
    interfaceLength = triangleInterfaceLength(points, values);
    waterVolume(kT) = max(waterArea * thickness, 0);
    solidVolume(kT) = max(solidArea * thickness, 0);
    interfaceArea(kT) = max(interfaceLength * thickness, 0);
end

geom.interface_area_cm2 = interfaceArea;
geom.water_volume_cm3 = waterVolume;
geom.calcite_moles = solidVolume ./ max(molarVolume, eps);
end

function area = triangleAreaBelowLevelZero(points, values)
polyPoints = points;
polyValues = values(:);
if isempty(polyPoints)
    area = 0;
    return;
end
inputPoints = polyPoints;
inputValues = polyValues;
polyPoints = zeros(0, 2);
polyValues = zeros(0, 1);
nInput = size(inputPoints, 1);
for iPoint = 1:nInput
    jPoint = mod(iPoint, nInput) + 1;
    p1 = inputPoints(iPoint, :);
    p2 = inputPoints(jPoint, :);
    v1 = inputValues(iPoint);
    v2 = inputValues(jPoint);
    inside1 = v1 <= 0;
    inside2 = v2 <= 0;
    if inside1 && inside2
        polyPoints(end + 1, :) = p2; %#ok<AGROW>
        polyValues(end + 1, 1) = v2; %#ok<AGROW>
    elseif inside1 && ~inside2
        [pCross, vCross] = interpolateLevelCrossing(p1, p2, v1, v2);
        polyPoints(end + 1, :) = pCross; %#ok<AGROW>
        polyValues(end + 1, 1) = vCross; %#ok<AGROW>
    elseif ~inside1 && inside2
        [pCross, vCross] = interpolateLevelCrossing(p1, p2, v1, v2);
        polyPoints(end + 1, :) = pCross; %#ok<AGROW>
        polyValues(end + 1, 1) = vCross; %#ok<AGROW>
        polyPoints(end + 1, :) = p2; %#ok<AGROW>
        polyValues(end + 1, 1) = v2; %#ok<AGROW>
    end
end

if size(polyPoints, 1) < 3
    area = 0;
else
    x = polyPoints(:, 1);
    y = polyPoints(:, 2);
    area = 0.5 * abs(sum(x .* y([2:end, 1]) - y .* x([2:end, 1])));
end
end

function lengthValue = triangleInterfaceLength(points, values)
crossings = zeros(0, 2);
edges = [1 2; 2 3; 3 1];
for iEdge = 1:3
    i1 = edges(iEdge, 1);
    i2 = edges(iEdge, 2);
    v1 = values(i1);
    v2 = values(i2);
    if v1 == 0
        crossings(end + 1, :) = points(i1, :); %#ok<AGROW>
    end
    if v1 * v2 < 0
        [pCross, ~] = interpolateLevelCrossing(points(i1, :), points(i2, :), v1, v2);
        crossings(end + 1, :) = pCross; %#ok<AGROW>
    end
end
if values(3) == 0
    crossings(end + 1, :) = points(3, :); %#ok<AGROW>
end

if size(crossings, 1) < 2
    lengthValue = 0;
    return;
end

crossings = unique(round(crossings, 14), 'rows');
if size(crossings, 1) < 2
    lengthValue = 0;
else
    distances = squareformDistance(crossings);
    lengthValue = max(distances(:));
end
end

function [point, value] = interpolateLevelCrossing(p1, p2, v1, v2)
denominator = v1 - v2;
if abs(denominator) < eps
    lambda = 0.5;
else
    lambda = v1 / denominator;
end
lambda = min(max(lambda, 0), 1);
point = p1 + lambda * (p2 - p1);
value = 0;
end

function distances = squareformDistance(points)
numPoints = size(points, 1);
distances = zeros(numPoints, numPoints);
for iPoint = 1:numPoints
    delta = points - points(iPoint, :);
    distances(iPoint, :) = sqrt(sum(delta.^2, 2));
end
end

function state = collectPhreeqcState(hTransport, caTransport, cTransport, naTransport, clTransport, timeStep, geometryState)
state = struct();
state.h_mol_cm3 = max(hTransport.U.getdata(timeStep), 0);
state.ca_mol_cm3 = max(caTransport.U.getdata(timeStep), 0);
state.c_mol_cm3 = max(cTransport.U.getdata(timeStep), 0);
state.na_mol_cm3 = max(naTransport.U.getdata(timeStep), 0);
state.cl_mol_cm3 = max(clTransport.U.getdata(timeStep), 0);
state.ca_total_mol_cm3 = state.ca_mol_cm3;
state.c_total_mol_cm3 = state.c_mol_cm3;
state.na_total_mol_cm3 = state.na_mol_cm3;
state.cl_total_mol_cm3 = state.cl_mol_cm3;
if nargin >= 7 && ~isempty(geometryState)
    state.interface_area_cm2 = geometryState.interface_area_cm2;
    state.water_volume_cm3 = geometryState.water_volume_cm3;
    state.calcite_moles = geometryState.calcite_moles;
end
end

function data = initializePhreeqcSpeciesDataFromTransports(hTransport, caTransport, cTransport, ...
    naTransport, clTransport, timeStep)
data = struct();
data.h_mol_cm3 = max(hTransport.U.getdata(timeStep), 0);
data.ca_total_mol_cm3 = max(caTransport.U.getdata(timeStep), 0);
data.c_total_mol_cm3 = max(cTransport.U.getdata(timeStep), 0);
data.na_total_mol_cm3 = max(naTransport.U.getdata(timeStep), 0);
data.cl_total_mol_cm3 = max(clTransport.U.getdata(timeStep), 0);
numCells = numel(data.h_mol_cm3);
data.ca_mol_cm3 = data.ca_total_mol_cm3;
data.hco3_mol_cm3 = zeros(numCells, 1);
data.co3_mol_cm3 = zeros(numCells, 1);
data.na_mol_cm3 = data.na_total_mol_cm3;
data.cl_mol_cm3 = data.cl_total_mol_cm3;
data.pH = -log10(max(data.h_mol_cm3 * 1000, 1e-7));
data.chargeBalance = zeros(numCells, 1);
data.calciteSI = NaN(numCells, 1);
data.calciteRate_mol_dm2_s = zeros(numCells, 1);
data.calciteRate_mol_s = zeros(numCells, 1);
data.calciteDeltaMoles = zeros(numCells, 1);
data.calciteDissolvedMoles = zeros(numCells, 1);
data.solutionNumber = (1:numCells)';
end

function state = updatePhreeqcStateFromResult(state, result)
fields = { ...
    'h_mol_cm3', 'h_activity_mol_cm3', 'ca_mol_cm3', 'hco3_mol_cm3', ...
    'co3_mol_cm3', 'cl_mol_cm3', 'na_mol_cm3', ...
    'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'cl_total_mol_cm3', 'na_total_mol_cm3', 'pH', 'chargeBalance'};
for iField = 1:numel(fields)
    fieldName = fields{iField};
    if isfield(result, fieldName) && ~isempty(result.(fieldName))
        state.(fieldName) = result.(fieldName)(:);
    end
end
if ~isfield(state, 'h_activity_mol_cm3') && isfield(state, 'pH')
    state.h_activity_mol_cm3 = 10 .^ (-state.pH(:)) ./ 1000;
end
end

function ratePerArea = computePhreeqcRatePerArea(result, interfaceAreaCm2, timeStepSize)
if isfield(result, 'calciteDissolvedMoles')
    dissolvedMoles = result.calciteDissolvedMoles;
else
    dissolvedMoles = max(result.calciteRate_mol_s * timeStepSize, 0);
end
ratePerArea = zeros(size(dissolvedMoles));
active = interfaceAreaCm2 > 0 & timeStepSize > 0;
ratePerArea(active) = dissolvedMoles(active) ./ max(timeStepSize, eps) ./ interfaceAreaCm2(active);
ratePerArea(~isfinite(ratePerArea)) = 0;
end

function ratePerArea = calciteMolesToInterfaceRate(dissolvedMoles, interfaceAreaCm2, timeStepSize)
ratePerArea = zeros(size(dissolvedMoles(:)));
active = interfaceAreaCm2(:) > 0 & timeStepSize > 0;
ratePerArea(active) = dissolvedMoles(active) ./ max(timeStepSize, eps) ...
    ./ max(interfaceAreaCm2(active), eps);
ratePerArea(~isfinite(ratePerArea)) = 0;
end

function applyPhreeqcResultToTransports(result, hTransport, caTransport, cTransport, ...
    naTransport, clTransport, timeStep)
hTransport.U.setdata(timeStep, max(result.h_mol_cm3, 0));
caTransport.U.setdata(timeStep, max(result.ca_total_mol_cm3, 0));
cTransport.U.setdata(timeStep, max(result.c_total_mol_cm3, 0));
naTransport.U.setdata(timeStep, max(result.na_total_mol_cm3, 0));
clTransport.U.setdata(timeStep, max(result.cl_total_mol_cm3, 0));
end

function result = applyAgglomeratedPhreeqcUpdate(result, baseState, reactionInputState, ...
    waterVolumeCm3, agglomerateWeightMatrix)
adjusted = ApplyPhreeqcAgglomeratedDelta(baseState, reactionInputState, result, ...
    waterVolumeCm3, agglomerateWeightMatrix);
memberCells = any(agglomerateWeightMatrix > 0, 2);
concentrationFields = { ...
    'h_mol_cm3', 'ca_total_mol_cm3', 'c_total_mol_cm3', ...
    'na_total_mol_cm3', 'cl_total_mol_cm3', ...
    'ca_mol_cm3', 'hco3_mol_cm3', 'co3_mol_cm3', ...
    'na_mol_cm3', 'cl_mol_cm3', 'pH'};
for iField = 1:numel(concentrationFields)
    fieldName = concentrationFields{iField};
    if isfield(result, fieldName) && isfield(adjusted, fieldName) && ...
            numel(result.(fieldName)) == numel(memberCells) && ...
            numel(adjusted.(fieldName)) == numel(memberCells)
        values = result.(fieldName)(:);
        adjustedValues = adjusted.(fieldName)(:);
        values(memberCells) = adjustedValues(memberCells);
        result.(fieldName) = values;
    end
end
if isfield(result, 'h_activity_mol_cm3') && isfield(result, 'pH')
    result.h_activity_mol_cm3 = 10 .^ (-result.pH(:)) ./ 1000;
end
end

function data = addPhreeqcWaterPhaseMassDiagnostics(data, baseState, waterVolumeCm3)
waterVolumeCm3 = max(waterVolumeCm3(:), 0);
data.water_phase_h_delta_moles_total = computeWaterPhaseDelta( ...
    data, baseState, waterVolumeCm3, 'h_mol_cm3');
data.water_phase_ca_delta_moles_total = computeWaterPhaseDelta( ...
    data, baseState, waterVolumeCm3, 'ca_total_mol_cm3');
data.water_phase_c_delta_moles_total = computeWaterPhaseDelta( ...
    data, baseState, waterVolumeCm3, 'c_total_mol_cm3');
data.water_phase_na_delta_moles_total = computeWaterPhaseDelta( ...
    data, baseState, waterVolumeCm3, 'na_total_mol_cm3');
data.water_phase_cl_delta_moles_total = computeWaterPhaseDelta( ...
    data, baseState, waterVolumeCm3, 'cl_total_mol_cm3');
end

function deltaMoles = computeWaterPhaseDelta(data, baseState, waterVolumeCm3, fieldName)
if ~isfield(data, fieldName) || ~isfield(baseState, fieldName)
    deltaMoles = NaN;
    return;
end
newValues = data.(fieldName)(:);
oldValues = baseState.(fieldName)(:);
if numel(newValues) ~= numel(waterVolumeCm3) || numel(oldValues) ~= numel(waterVolumeCm3)
    deltaMoles = NaN;
    return;
end
deltaMoles = sum((newValues - oldValues) .* waterVolumeCm3, 'omitnan');
end

function value = computeAgglomerateMaxRowOverlap(agglomerateWeightMatrix)
if isempty(agglomerateWeightMatrix) || nnz(agglomerateWeightMatrix) == 0
    value = 0;
    return;
end
rowCounts = full(sum(agglomerateWeightMatrix > 0, 2));
value = max(rowCounts, [], 'omitnan');
end

function value = computeAgglomerateMeanRowOverlap(agglomerateWeightMatrix)
if isempty(agglomerateWeightMatrix) || nnz(agglomerateWeightMatrix) == 0
    value = 0;
    return;
end
rowCounts = full(sum(agglomerateWeightMatrix > 0, 2));
rowCounts = rowCounts(rowCounts > 0);
if isempty(rowCounts)
    value = 0;
else
    value = mean(rowCounts, 'omitnan');
end
end

function value = computeAgglomerateMaxVolumeRatio(reactionWaterCm3, waterVolumeCm3, activeMoles)
if isempty(reactionWaterCm3) || isempty(waterVolumeCm3) || isempty(activeMoles)
    value = NaN;
    return;
end
active = activeMoles(:) > 0 & waterVolumeCm3(:) > 0;
if ~any(active)
    value = 1;
else
    value = max(reactionWaterCm3(active) ./ max(waterVolumeCm3(active), eps), [], 'omitnan');
end
end

function [nonlinearFunc, nonlinearJacFunc] = zeroPhreeqcNonlinearFunctions(numSpecies, ~)
nonlinearFunc = cell(numSpecies, 1);
nonlinearJacFunc = cell(numSpecies, numSpecies);
for iSpecies = 1:numSpecies
    nonlinearFunc{iSpecies} = @(varargin) zeros(size(varargin{1}));
    for jSpecies = 1:numSpecies
        nonlinearJacFunc{iSpecies, jSpecies} = @(varargin) zeros(size(varargin{1}));
    end
end
end

function applyTransportCorrector(speciesCells, Corrector, timeStep)
for iSpecies = 1:numel(speciesCells)
    if size(Corrector, 2) >= iSpecies
        speciesCells{iSpecies}.U.setdata(timeStep - 1, ...
            speciesCells{iSpecies}.U.getdata(timeStep - 1) - Corrector(:, iSpecies));
    end
end
end

function speciesCells = limitPhreeqcTransportBounds(speciesCells, upperBounds, speciesNames, timeStep)
for iSpecies = 1:min(numel(speciesCells), numel(upperBounds))
    data = speciesCells{iSpecies}.U.getdata(timeStep);
    data = max(data, 0);
    upperBound = upperBounds(iSpecies);
    if isfinite(upperBound)
        tolerance = max(1e-15, 1e-9 * max(abs(upperBound), 1));
        oldMax = max(data, [], 'omitnan');
        clipped = data > upperBound;
        if any(clipped) && oldMax > upperBound + tolerance
            fprintf('  [PHREEQC limiter] %s step %d: clipped %d cells from max %.3e to %.3e mol/cm^3\n', ...
                speciesNames{iSpecies}, timeStep, nnz(clipped), oldMax, upperBound);
        end
        data(clipped) = upperBound;
    end
    speciesCells{iSpecies}.U.setdata(timeStep, data);
end
end

function vertexData = averageTriangleDataToVertices(triangleData, vertexTriMatrix, numVertices)
vertexData = zeros(numVertices, 1);
for iVertex = 1:numVertices
    triIndices = vertexTriMatrix{1, iVertex};
    values = triangleData(triIndices);
    values = values(isfinite(values));
    if isempty(values)
        vertexData(iVertex) = 0;
    else
        vertexData(iVertex) = mean(values);
    end
end
end

function vertexData = averagePositiveTriangleDataToVertices(triangleData, vertexTriMatrix, numVertices)
vertexData = zeros(numVertices, 1);
for iVertex = 1:numVertices
    triIndices = vertexTriMatrix{1, iVertex};
    values = triangleData(triIndices);
    values = values(isfinite(values) & values > 0);
    if isempty(values)
        vertexData(iVertex) = 0;
    else
        vertexData(iVertex) = mean(values);
    end
end
end

function [projectedMoles, reactionWaterCm3] = projectInterfaceReactionToPoreCells( ...
    interfaceMoles, grid, vertexTriMatrix, levelSetData, waterVolumeCm3, ...
    minCells, maxCells, maxRings, triangleVelocityData, flowBias)
numTriangles = numel(interfaceMoles);
projectedMoles = zeros(numTriangles, 1);
reactionWaterCm3 = max(waterVolumeCm3(:), 0);
if numTriangles == 0
    return;
end

if nargin < 9 || isempty(triangleVelocityData)
    triangleVelocityData = zeros(numTriangles, 2);
end
if nargin < 10 || isempty(flowBias)
    flowBias = 0;
end
minCells = max(1, round(minCells));
maxCells = max(minCells, round(maxCells));
maxRings = max(0, round(maxRings));
flowBias = max(flowBias, 0);
levelSetOnTriangles = levelSetData(grid.V0T);
fullyPoreTriangles = all(levelSetOnTriangles < 0, 2) & reactionWaterCm3 > 0;
positiveSources = find(interfaceMoles(:) > 0 & isfinite(interfaceMoles(:)));
for iSource = 1:numel(positiveSources)
    sourceTri = positiveSources(iSource);
    if fullyPoreTriangles(sourceTri)
        candidates = sourceTri;
    else
        candidates = findNearbyPoreTriangles(sourceTri, grid, vertexTriMatrix, ...
            fullyPoreTriangles, minCells, maxCells, maxRings);
    end
    if isempty(candidates)
        candidates = sourceTri;
    end

    candidateWater = max(reactionWaterCm3(candidates), 0);
    delta = grid.baryT(candidates, :) - grid.baryT(sourceTri, :);
    distances = sqrt(sum(delta.^2, 2));
    distanceScale = max(sqrt(max(grid.areaT(sourceTri), eps)), eps);
    weights = candidateWater ./ max(distances + distanceScale, eps);
    localVelocity = triangleVelocityData(sourceTri, :);
    localSpeed = sqrt(sum(localVelocity.^2));
    if localSpeed > eps && flowBias > 0
        directionScore = (delta * localVelocity(:)) ./ max(distances * localSpeed, eps);
        directionScore = max(min(directionScore, 1), -1);
        weights = weights .* exp(flowBias .* directionScore);
    end
    if ~any(weights > 0) || any(~isfinite(weights))
        weights = ones(numel(candidates), 1);
    end
    weights = weights(:) ./ sum(weights(:));
    projectedMoles(candidates) = projectedMoles(candidates) + interfaceMoles(sourceTri) .* weights;
end
end

function [reactionWaterCm3, agglomerateWeightMatrix] = computeAgglomeratedReactionWaterVolumes( ...
    interfaceMoles, grid, vertexTriMatrix, levelSetData, waterVolumeCm3, ...
    minCells, maxCells, maxRings)
% Use a local cut-cell agglomerate as the PHREEQC mixing volume while the
% prescribed CaCO3 amount remains attached to the interface cell. This avoids
% singular mol/kgw reactions caused by arbitrarily small water slivers.
numTriangles = numel(interfaceMoles);
reactionWaterCm3 = max(waterVolumeCm3(:), 0);
agglomerateWeightMatrix = sparse(numTriangles, numTriangles);
if numTriangles == 0
    return;
end

minCells = max(1, round(minCells));
maxCells = max(minCells, round(maxCells));
maxRings = max(0, round(maxRings));
levelSetOnTriangles = levelSetData(grid.V0T);
fullyPoreTriangles = all(levelSetOnTriangles < 0, 2) & reactionWaterCm3 > 0;
positiveSources = find(interfaceMoles(:) > 0 & isfinite(interfaceMoles(:)));

for iSource = 1:numel(positiveSources)
    sourceTri = positiveSources(iSource);
    candidates = sourceTri;
    if ~fullyPoreTriangles(sourceTri)
        nearbyPores = findNearbyPoreTriangles(sourceTri, grid, vertexTriMatrix, ...
            fullyPoreTriangles, minCells, maxCells, maxRings);
        candidates = unique([sourceTri; nearbyPores(:)]);
    end
    candidateWater = max(reactionWaterCm3(candidates), 0);
    positiveWater = candidateWater > 0;
    candidates = candidates(positiveWater);
    candidateWater = candidateWater(positiveWater);
    localWater = sum(candidateWater);
    if isfinite(localWater) && localWater > reactionWaterCm3(sourceTri)
        reactionWaterCm3(sourceTri) = localWater;
    end
    if isfinite(localWater) && localWater > 0
        agglomerateWeightMatrix(candidates, sourceTri) = candidateWater ./ localWater;
    end
end
end

function candidates = findNearbyPoreTriangles(sourceTri, grid, vertexTriMatrix, ...
    fullyPoreTriangles, minCells, maxCells, maxRings)
numTriangles = size(grid.V0T, 1);
visited = false(numTriangles, 1);
frontier = sourceTri;
candidates = [];
for iRing = 0:maxRings
    frontier = unique(frontier(:));
    frontier = frontier(frontier >= 1 & frontier <= numTriangles & ~visited(frontier));
    if isempty(frontier)
        break;
    end
    visited(frontier) = true;
    poreNow = frontier(fullyPoreTriangles(frontier));
    if ~isempty(poreNow)
        candidates = unique([candidates; poreNow(:)]); %#ok<AGROW>
        if numel(candidates) >= minCells || iRing == maxRings
            break;
        end
    end

    nextFrontier = [];
    for iTri = 1:numel(frontier)
        vertices = grid.V0T(frontier(iTri), :);
        for iVertex = 1:numel(vertices)
            nextFrontier = [nextFrontier, vertexTriMatrix{1, vertices(iVertex)}]; %#ok<AGROW>
        end
    end
    frontier = setdiff(unique(nextFrontier(:)), find(visited));
end

if numel(candidates) > maxCells
    distances = sqrt(sum((grid.baryT(candidates, :) - grid.baryT(sourceTri, :)).^2, 2));
    [~, order] = sort(distances, 'ascend');
    candidates = candidates(order(1:maxCells));
end
end

function triangleVelocityData = stokesVelocityToTriangleVelocity(grid, stokesVelocityData)
triangleVelocityData = zeros(grid.numT, 2);
nP2 = grid.numV + grid.numE;
if isempty(stokesVelocityData)
    return;
end
if size(stokesVelocityData, 2) == 2 && size(stokesVelocityData, 1) >= grid.numV
    vertexVelocity = stokesVelocityData(1:grid.numV, :);
elseif size(stokesVelocityData, 2) == 1 && numel(stokesVelocityData) >= 2 * nP2
    vertexVelocity = [stokesVelocityData(1:grid.numV), ...
        stokesVelocityData(nP2 + (1:grid.numV))];
else
    return;
end
vertexVelocityX = vertexVelocity(:, 1);
vertexVelocityY = vertexVelocity(:, 2);
triangleVelocityData(:, 1) = mean(vertexVelocityX(grid.V0T), 2);
triangleVelocityData(:, 2) = mean(vertexVelocityY(grid.V0T), 2);
triangleVelocityData(~isfinite(triangleVelocityData)) = 0;
end

function normalSpeed = computeLegacyTstNormalSpeedFromTriangleHydrogen( ...
    triangleHydrogen, vertexTriMatrix, numVertices, molarVolume, rateCoefficientTST)
hydrogenDataV = zeros(numVertices, 1);
for iVertex = 1:numVertices
    triIndices = vertexTriMatrix{1, iVertex};
    values = triangleHydrogen(triIndices);
    hydrogenDataV(iVertex) = min(values);
end
normalSpeed = molarVolume .* max(hydrogenDataV, 0) .* 1000 .* rateCoefficientTST;
normalSpeed(~isfinite(normalSpeed)) = 0;
end

function rateDriver = getPhreeqcTstRateDriver(phreeqcSpeciesData)
if isfield(phreeqcSpeciesData, 'h_activity_mol_cm3') && ...
        ~isempty(phreeqcSpeciesData.h_activity_mol_cm3)
    rateDriver = phreeqcSpeciesData.h_activity_mol_cm3;
else
    rateDriver = phreeqcSpeciesData.h_mol_cm3;
end
rateDriver = rateDriver(:);
end

function rateValue = computeHydrogenFluxApparentRate(hTransport, timeStep, ...
    inletConcentration, inletVelocity, flowConfig, grid, surfaceAreaCm2)
% Historical RTSPHEM rate diagnostic from outlet H+ flux deficit.
% For PHREEQC multi-component chemistry this remains an apparent acid
% consumption diagnostic, not a strict CaCO3 molar dissolution balance.
outletEdges = flowConfig.outletEdgePredicate(grid.baryE);
hydrogenFlux = hTransport.Q.getdata(timeStep);
rateValue = -((hydrogenFlux - inletConcentration * inletVelocity) ...
    .* outletEdges)' * grid.areaE ./ max(surfaceAreaCm2, eps);
if ~isfinite(rateValue)
    rateValue = NaN;
end
end

function initializePhreeqcSummaryLog(csvFile)
fid = fopen(csvFile, 'w');
if fid ~= -1
    fprintf(fid, ['timestep,time_s,pH_mean,pH_min,pH_max,H_mean_mol_cm3,', ...
        'Ca_mean_mol_cm3,HCO3_mean_mol_cm3,CO3_mean_mol_cm3,Cl_mean_mol_cm3,', ...
        'Na_mean_mol_cm3,calcite_rate_mean_mol_cm2_s,calcite_dissolved_moles_mean,', ...
        'calcite_si_mean,charge_balance_mean,interface_prescribed_moles_total,', ...
        'projected_prescribed_moles_total,result_dissolved_moles_total,', ...
        'interface_prescribed_active_cells,projected_prescribed_active_cells,', ...
        'interface_projected_overlap_cells,water_phase_h_delta_moles_total,', ...
        'water_phase_ca_delta_moles_total,water_phase_c_delta_moles_total,', ...
        'water_phase_na_delta_moles_total,water_phase_cl_delta_moles_total,', ...
        'agglomerate_max_row_overlap,agglomerate_mean_row_overlap,', ...
        'agglomerate_max_volume_ratio\n']);
    fclose(fid);
end
end

function appendPhreeqcSummaryLog(csvFile, stepIndex, timeSeconds, data)
fid = fopen(csvFile, 'a');
if fid == -1
    warning('RTSPHEM:Phreeqc:SummaryLogOpenFailed', 'Cannot open PHREEQC summary log: %s', csvFile);
    return;
end
if isfield(data, 'calciteRatePerArea_mol_cm2_s')
    rateForLog = data.calciteRatePerArea_mol_cm2_s;
else
    rateForLog = data.calciteRate_mol_dm2_s;
end
if isfield(data, 'calciteDissolvedMoles')
    dissolvedMolesForLog = data.calciteDissolvedMoles;
else
    dissolvedMolesForLog = NaN(size(rateForLog));
end
[interfaceMolesTotal, projectedMolesTotal, resultMolesTotal, ...
    interfaceActiveCells, projectedActiveCells, overlapCells] = ...
    computePhreeqcMoleSupportDiagnostics(data, dissolvedMolesForLog);
fprintf(fid, ['%d,%.12g,%.12g,%.12g,%.12g,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,', ...
    '%.12e,%.12e,%.12g,%.12g,%.12e,%.12e,%.12e,%d,%d,%d,', ...
    '%.12e,%.12e,%.12e,%.12e,%.12e,%.12g,%.12g,%.12g\n'], ...
    stepIndex, timeSeconds, mean(data.pH, 'omitnan'), min(data.pH, [], 'omitnan'), ...
    max(data.pH, [], 'omitnan'), mean(data.h_mol_cm3, 'omitnan'), ...
    mean(data.ca_mol_cm3, 'omitnan'), mean(data.hco3_mol_cm3, 'omitnan'), ...
    mean(data.co3_mol_cm3, 'omitnan'), mean(data.cl_mol_cm3, 'omitnan'), ...
    mean(data.na_mol_cm3, 'omitnan'), mean(rateForLog, 'omitnan'), ...
    mean(dissolvedMolesForLog, 'omitnan'), ...
    mean(data.calciteSI, 'omitnan'), mean(data.chargeBalance, 'omitnan'), ...
    interfaceMolesTotal, projectedMolesTotal, resultMolesTotal, ...
    interfaceActiveCells, projectedActiveCells, overlapCells, ...
    optionalScalarField(data, 'water_phase_h_delta_moles_total'), ...
    optionalScalarField(data, 'water_phase_ca_delta_moles_total'), ...
    optionalScalarField(data, 'water_phase_c_delta_moles_total'), ...
    optionalScalarField(data, 'water_phase_na_delta_moles_total'), ...
    optionalScalarField(data, 'water_phase_cl_delta_moles_total'), ...
    optionalScalarField(data, 'agglomerateMaxRowOverlap'), ...
    optionalScalarField(data, 'agglomerateMeanRowOverlap'), ...
    optionalScalarField(data, 'agglomerateMaxVolumeRatio'));
fclose(fid);
end

function value = optionalScalarField(data, fieldName)
if isfield(data, fieldName) && ~isempty(data.(fieldName))
    value = data.(fieldName);
else
    value = NaN;
end
end

function [interfaceMolesTotal, projectedMolesTotal, resultMolesTotal, ...
    interfaceActiveCells, projectedActiveCells, overlapCells] = ...
    computePhreeqcMoleSupportDiagnostics(data, dissolvedMolesForLog)
if isfield(data, 'interfacePrescribedCalciteMoles') && ...
        ~isempty(data.interfacePrescribedCalciteMoles)
    interfaceMoles = data.interfacePrescribedCalciteMoles(:);
else
    interfaceMoles = NaN(size(dissolvedMolesForLog(:)));
end
if isfield(data, 'projectedPrescribedCalciteMoles') && ...
        ~isempty(data.projectedPrescribedCalciteMoles)
    projectedMoles = data.projectedPrescribedCalciteMoles(:);
else
    projectedMoles = NaN(size(dissolvedMolesForLog(:)));
end
resultMoles = dissolvedMolesForLog(:);

interfaceMolesTotal = sum(interfaceMoles, 'omitnan');
projectedMolesTotal = sum(projectedMoles, 'omitnan');
resultMolesTotal = sum(resultMoles, 'omitnan');

interfaceActive = isfinite(interfaceMoles) & interfaceMoles > 0;
projectedActive = isfinite(projectedMoles) & projectedMoles > 0;
interfaceActiveCells = nnz(interfaceActive);
projectedActiveCells = nnz(projectedActive);
if numel(interfaceActive) == numel(projectedActive)
    overlapCells = nnz(interfaceActive & projectedActive);
else
    overlapCells = 0;
end
end

function databasePath = defaultPhreeqcDatabasePath()
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
repoCandidate = fullfile(projectRoot, 'ReactiveTransport', 'RTM', 'couplePhreeqc', 'phreeqc-m.dat');
referenceCandidate = 'C:\Users\imgw\Downloads\RTSPHEM-P-main (1)\RTSPHEM-P-main\SourceCode\phreeqc-m.dat';
usgsCandidate = 'C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat';
if exist(repoCandidate, 'file')
    databasePath = repoCandidate;
elseif exist(referenceCandidate, 'file')
    databasePath = referenceCandidate;
elseif exist(usgsCandidate, 'file')
    databasePath = usgsCandidate;
else
    databasePath = 'phreeqc.dat';
end
end

function rateLaw = defaultPhreeqcRateLaw(runGroup)
normalized = lower(strrep(strtrim(char(runGroup)), '-', '_'));
switch normalized
    case {'phreeqc_database_calcite', 'database_calcite', 'database', 'calcite'}
        rateLaw = 'database_calcite';
    case {'phreeqc_tst_match', 'tst_match', 'calcite_tst_match'}
        rateLaw = 'tst_match';
    otherwise
        error('RTSPHEM:Phreeqc:UnknownRunGroup', ...
            'Unknown PHREEQC run group: %s.', char(runGroup));
end
end

function rateLaw = normalizePhreeqcRateLaw(rateLaw)
normalized = lower(strrep(strtrim(char(rateLaw)), '-', '_'));
switch normalized
    case {'phreeqc_database_calcite', 'database_calcite', 'database', 'calcite'}
        rateLaw = 'database_calcite';
    case {'phreeqc_tst_match', 'tst_match', 'calcite_tst_match'}
        rateLaw = 'tst_match';
    otherwise
        error('RTSPHEM:Phreeqc:UnknownRateLaw', ...
            'Unknown PHREEQC rate law: %s.', char(rateLaw));
end
end

function value = cfgget(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end

function tf = cfgHasNonEmptyField(config, fieldName)
tf = isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName));
end

function setupPNMPaths(reactiveRoot, rtmDir)
addpath(rtmDir);
addpath(fullfile(rtmDir, 'couplePhreeqc'));
addpath(genpath(fullfile(reactiveRoot, 'src')));
addpath(fullfile(reactiveRoot, 'HyPHM'));
addpath(fullfile(reactiveRoot, 'HyPHM', 'tools'));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'classes')));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'domains')));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'opt')));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'symbolic')));
end

function writeJsonFile(outputPath, data)
try
    fid = fopen(outputPath, 'w');
    if fid == -1
        warning('MATLAB:PNMJsonWrite', 'Cannot open metadata file: %s', outputPath);
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', jsonencode(data));
    clear cleaner;
catch ME
    warning('MATLAB:PNMJsonWrite', 'Failed to write metadata JSON: %s', ME.message);
end
end

function tf = shouldExportStep(stepIndex, lastStepIndex, exportEvery, forceExport)
tf = forceExport || stepIndex == 1 || stepIndex == lastStepIndex || mod(stepIndex, exportEvery) == 0;
end

function limitedTimeSteps = limitTimeGridByMaxStep(timeSteps, maxStep)
if ~isfinite(maxStep) || maxStep <= 0
    limitedTimeSteps = timeSteps;
    return;
end

limitedTimeSteps = timeSteps(1);
for iStep = 1:(numel(timeSteps) - 1)
    t0 = timeSteps(iStep);
    t1 = timeSteps(iStep + 1);
    dt = t1 - t0;
    if dt <= maxStep * (1 + 10 * eps)
        limitedTimeSteps(end + 1) = t1; %#ok<AGROW>
    else
        numSubSteps = max(1, ceil(dt / maxStep));
        subSteps = linspace(t0, t1, numSubSteps + 1);
        limitedTimeSteps = [limitedTimeSteps, subSteps(2:end)]; %#ok<AGROW>
    end
end
limitedTimeSteps = unique(limitedTimeSteps, 'stable');
end

function numSteps = countLimitedTimeGridSteps(timeSteps, maxStep)
if ~isfinite(maxStep) || maxStep <= 0
    numSteps = numel(timeSteps) - 1;
    return;
end

numSteps = 0;
for iStep = 1:(numel(timeSteps) - 1)
    dt = timeSteps(iStep + 1) - timeSteps(iStep);
    numSteps = numSteps + max(1, ceil(dt / maxStep));
end
end

function diagnostics = computeStabilityDiagnostics(grid, levels, previousConcentration, currentConcentration, ...
    flowData, stokesVelocityData, diffusionCoefficient, hMin, dt, sourceScale, ...
    rateCoefficientTST, thickness, inletConcentration, inletVelocity, flowConfig, ...
    normalSpeedMax, porosityValue, permeabilityValue, permeabilityRatio, upwindMode)

fluidTriangles = any(levels(grid.V0T) < -eps, 2);
cCurrent = currentConcentration(:);
cPrevious = previousConcentration(:);

diagnostics = struct();
diagnostics.cMin = min(cCurrent, [], 'omitnan');
diagnostics.cMax = max(cCurrent, [], 'omitnan');
diagnostics.totalMass = computeTotalSoluteMass(grid, fluidTriangles, cCurrent, thickness);
previousMass = computeTotalSoluteMass(grid, fluidTriangles, cPrevious, thickness);

outletTriangles = flowConfig.outletTrianglePredicate(grid.baryT);
outletConcentration = ComputeOutletConcentration(cCurrent, outletTriangles);
outletEdges = flowConfig.outletEdgePredicate(grid.baryE);
outletWaterFlux = sum(abs(flowData(outletEdges)) .* grid.areaE(outletEdges)) * thickness;
inletMassFlux = inletConcentration * abs(inletVelocity) * flowConfig.crossSectionLength * thickness;
outletMassFlux = max(outletConcentration, 0) * outletWaterFlux;

reactionSink = 0;
if numel(sourceScale) == grid.numT
    reactionRate = max(sourceScale(:), 0) .* (rateCoefficientTST * 1000) .* max(cCurrent, 0);
    reactionSink = sum(reactionRate(fluidTriangles) .* grid.areaT(fluidTriangles)) * thickness;
end

diagnostics.massBalanceResidual = diagnostics.totalMass - previousMass - ...
    dt * (inletMassFlux - outletMassFlux - reactionSink);
massScale = max([abs(previousMass), abs(diagnostics.totalMass), abs(dt * inletMassFlux), eps]);
diagnostics.massBalanceRelative = abs(diagnostics.massBalanceResidual) / massScale;

[diagnostics.localVelocityMax, diagnostics.localPeMax] = computeLocalVelocityStats( ...
    grid, levels, flowData, stokesVelocityData, diffusionCoefficient, hMin);
diagnostics.advectiveCfl = diagnostics.localVelocityMax * dt / max(hMin, eps);
diagnostics.reactionCfl = normalSpeedMax * dt / max(hMin, eps);
diagnostics.normalSpeedMax = normalSpeedMax;
[diagnostics.poreComponents, diagnostics.largestPoreComponentFraction] = ...
    computePoreConnectivity(grid, fluidTriangles);
diagnostics.porosity = porosityValue;
diagnostics.permeability = permeabilityValue;
diagnostics.permeabilityRatio = permeabilityRatio;
diagnostics.upwindMode = char(upwindMode);
diagnostics.flag = joinDiagnosticFlags(diagnostics, inletConcentration);
end

function totalMass = computeTotalSoluteMass(grid, fluidTriangles, concentration, thickness)
if isempty(concentration)
    totalMass = NaN;
    return;
end
valid = fluidTriangles & isfinite(concentration);
totalMass = sum(concentration(valid) .* grid.areaT(valid)) * thickness;
end

function [localVelocityMax, localPeMax] = computeLocalVelocityStats(grid, levels, flowData, ...
    stokesVelocityData, diffusionCoefficient, hMin)
fluidTriangles = any(levels(grid.V0T) < -eps, 2);
localSpeed = NaN(grid.numT, 1);
try
    flowP0P0 = RT0.RT0toP0P0slice(grid, flowData);
    localSpeed = sqrt(sum(flowP0P0.^2, 2));
catch
end

stokesVelocityMax = NaN;
nP2 = grid.numV + grid.numE;
if size(stokesVelocityData, 2) == 2 && size(stokesVelocityData, 1) == nP2
    vx = stokesVelocityData(:, 1);
    vy = stokesVelocityData(:, 2);
elseif size(stokesVelocityData, 2) == 1 && length(stokesVelocityData) == 2 * nP2
    vx = stokesVelocityData(1:nP2);
    vy = stokesVelocityData(nP2+1:end);
else
    vx = [];
    vy = [];
end
if ~isempty(vx)
    lsEdgeMid = 0.5 * (levels(grid.V0E(:,1)) + levels(grid.V0E(:,2)));
    poreP2 = [levels; lsEdgeMid] < -eps;
    p2Speed = sqrt(vx.^2 + vy.^2);
    if any(poreP2)
        stokesVelocityMax = max(p2Speed(poreP2), [], 'omitnan');
    end
end

if any(fluidTriangles) && any(isfinite(localSpeed(fluidTriangles)))
    rt0VelocityMax = max(localSpeed(fluidTriangles), [], 'omitnan');
else
    rt0VelocityMax = NaN;
end
localVelocityMax = max([rt0VelocityMax, stokesVelocityMax], [], 'omitnan');
if isempty(localVelocityMax) || ~isfinite(localVelocityMax)
    localVelocityMax = 0;
end

if any(fluidTriangles) && any(isfinite(localSpeed(fluidTriangles)))
    localPeMax = max(localSpeed(fluidTriangles) * hMin ./ max(2 * diffusionCoefficient, eps), [], 'omitnan');
else
    localPeMax = localVelocityMax * hMin / max(2 * diffusionCoefficient, eps);
end
end

function [numComponents, largestFraction] = computePoreConnectivity(grid, fluidTriangles)
fluidIds = find(fluidTriangles);
numFluid = numel(fluidIds);
if numFluid == 0
    numComponents = 0;
    largestFraction = 0;
    return;
end

fluidMap = zeros(grid.numT, 1);
fluidMap(fluidIds) = 1:numFluid;
neighbors = cell(numFluid, 1);
for iEdge = 1:grid.numE
    adjacent = grid.T0E(iEdge, :);
    adjacent = adjacent(adjacent > 0);
    if numel(adjacent) == 2 && fluidTriangles(adjacent(1)) && fluidTriangles(adjacent(2))
        a = fluidMap(adjacent(1));
        b = fluidMap(adjacent(2));
        neighbors{a}(end + 1) = b; %#ok<AGROW>
        neighbors{b}(end + 1) = a; %#ok<AGROW>
    end
end

visited = false(numFluid, 1);
componentSizes = [];
for iFluid = 1:numFluid
    if visited(iFluid)
        continue;
    end
    stack = iFluid;
    visited(iFluid) = true;
    count = 0;
    while ~isempty(stack)
        node = stack(end);
        stack(end) = [];
        count = count + 1;
        nextNodes = neighbors{node};
        for iNext = 1:numel(nextNodes)
            nextNode = nextNodes(iNext);
            if ~visited(nextNode)
                visited(nextNode) = true;
                stack(end + 1) = nextNode; %#ok<AGROW>
            end
        end
    end
    componentSizes(end + 1) = count; %#ok<AGROW>
end

numComponents = numel(componentSizes);
largestFraction = max(componentSizes) / numFluid;
end

function nextStepSize = computeHardCflStepSize(diagnostics, hMin, advectiveSafety, reactiveSafety, minStep)
advectiveLimit = Inf;
reactiveLimit = Inf;
if diagnostics.localVelocityMax > 0
    advectiveLimit = advectiveSafety * hMin / diagnostics.localVelocityMax;
end
if diagnostics.normalSpeedMax > 0
    reactiveLimit = reactiveSafety * hMin / diagnostics.normalSpeedMax;
end
nextStepSize = max(minStep, min(advectiveLimit, reactiveLimit));
end

function flagText = joinDiagnosticFlags(diagnostics, inletConcentration)
flags = {};
if ~isfinite(diagnostics.cMin) || ~isfinite(diagnostics.cMax)
    flags{end + 1} = 'nonfinite_c'; %#ok<AGROW>
end
if diagnostics.cMin < -max(1e-12, 1e-6 * max(abs(inletConcentration), eps))
    flags{end + 1} = 'negative_c'; %#ok<AGROW>
end
if diagnostics.cMax > max(10 * abs(inletConcentration), abs(inletConcentration) + 1e-12)
    flags{end + 1} = 'overshoot_c'; %#ok<AGROW>
end
if diagnostics.advectiveCfl > 1
    flags{end + 1} = 'advective_cfl_gt_1'; %#ok<AGROW>
end
if diagnostics.reactionCfl > 1
    flags{end + 1} = 'reaction_cfl_gt_1'; %#ok<AGROW>
end
if diagnostics.poreComponents > 1
    flags{end + 1} = 'disconnected_pores'; %#ok<AGROW>
end
if diagnostics.massBalanceRelative > 0.25
    flags{end + 1} = 'mass_balance_drift'; %#ok<AGROW>
end
if isempty(flags)
    flagText = 'ok';
else
    flagText = strjoin(flags, ';');
end
end

function plotEvolutionSummary(timeValues, rateValues, porosityValues, permeabilityValues, ...
    surfaceAreaValues, volumeValues, pvValues, tortuosityValues, tortuositySegments, outputPath)

finalFig = figure('Visible', 'off', 'Position', [100, 100, 1200, 2100]);
fontSize = 14;
titleFontSize = 16;
fontName = 'Helvetica';

subplot(7,1,1);
plot(timeValues, rateValues, 'b-', 'LineWidth', 1.5);
xlabel('{\itt} (s)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Avg. Dissolution Rate', 'FontSize', fontSize, 'FontName', fontName);
title('Dissolution Rate Evolution', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
set(gca, 'FontSize', fontSize, 'FontName', fontName); grid on;

subplot(7,1,2);
plot(timeValues, porosityValues, 'k-', 'LineWidth', 1.5);
xlabel('{\itt} (s)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Porosity [-]', 'FontSize', fontSize, 'FontName', fontName);
title('Porosity Evolution', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
ylim([0 1]); set(gca, 'FontSize', fontSize, 'FontName', fontName); grid on;

subplot(7,1,3);
plot(timeValues, permeabilityValues, 'm-', 'LineWidth', 1.5);
xlabel('{\itt} (s)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Permeability (mD)', 'FontSize', fontSize, 'FontName', fontName);
title('Permeability Evolution', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
if any(permeabilityValues > 0 & ~isnan(permeabilityValues))
    set(gca, 'YScale', 'log');
end
set(gca, 'FontSize', fontSize, 'FontName', fontName); grid on;

subplot(7,1,4);
plot(timeValues, surfaceAreaValues, 'r-', 'LineWidth', 1.5);
xlabel('{\itt} (s)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Surface Area (cm^2)', 'FontSize', fontSize, 'FontName', fontName);
title('Surface Area Evolution', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
set(gca, 'FontSize', fontSize, 'FontName', fontName); grid on;

subplot(7,1,5);
plot(timeValues, volumeValues, 'g-', 'LineWidth', 1.5);
xlabel('{\itt} (s)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Grain Volume (cm^3)', 'FontSize', fontSize, 'FontName', fontName);
title('Grain Volume Evolution', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
set(gca, 'FontSize', fontSize, 'FontName', fontName); grid on;

subplot(7,1,6);
plot(timeValues, pvValues, 'c-', 'LineWidth', 1.5);
xlabel('{\itt} (s)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Injected PV', 'FontSize', fontSize, 'FontName', fontName);
title('Injected Pore Volumes', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
set(gca, 'FontSize', fontSize, 'FontName', fontName); grid on;

subplot(7,1,7);
hold on;
plot(timeValues, tortuosityValues, 'k-', 'LineWidth', 1.8, 'DisplayName', '\tau all');
segmentColors = {'r', 'g', 'b', 'm', 'c'};
for iSeg = 1:min(5, size(tortuositySegments, 1))
    plot(timeValues, tortuositySegments(iSeg,:), [segmentColors{iSeg} '-'], ...
        'LineWidth', 1.2, 'DisplayName', sprintf('\\tau seg%d', iSeg));
end
xlabel('{\itt} (s)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Tortuosity {\tau} [-]', 'FontSize', fontSize, 'FontName', fontName);
title('Tortuosity Evolution (all + 5 X-segments)', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
legend('Location', 'best'); hold off;
set(gca, 'FontSize', fontSize, 'FontName', fontName); grid on;

safePrintPng(finalFig, outputPath, 'final_global_evolution_png');
safeCloseFigure(finalFig);
end
