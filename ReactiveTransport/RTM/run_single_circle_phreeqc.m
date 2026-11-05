function result = run_single_circle_phreeqc(overrides)
%% run_single_circle_phreeqc - 单个方解石圆粒 + PHREEQC 的简单溶解算例
%
% 适用场景：
%   快速检查 PHREEQC 耦合的 RTM 溶解流程（界面推进、水化学导出、DXF 导出），
%   几何为通道中心的一个圆形方解石颗粒（Molins Part I/II 尺寸：0.1 x 0.05 cm，R = 0.01 cm）。
%
% 输入：
%   overrides  (可选) struct，覆盖下面任意字段，再传给 run_phreeqc_se。
% 输出：
%   outputs/rtm_runs/<runName>/  （结构同 run_phreeqc_se：global_evolution.xlsx、
%   dxf_pore/、dxf_solid/、phreeqc_results/、run_metadata.json 等）
% 依赖：
%   run_phreeqc_se.m -> PNM_beauty3.m（layoutType = 'single_circle'），PHREEQC COM/IPhreeqc。
%
% 用法：
%   result = run_single_circle_phreeqc();
%   result = run_single_circle_phreeqc(struct('acidConcentration_mol_L', 0.01));

if nargin < 1 || isempty(overrides)
    overrides = struct();
end

cfg = struct();
cfg.runName = 'single_circle_R100um_0p1x0p05cm_phreeqc_100x50';

% ---- 几何：一个居中的圆形颗粒 ----
cfg.layoutType = 'single_circle';
cfg.loadExistingGeometry = false;
cfg.useExternalGeometry = false;
cfg.useExternalDxfGeometry = false;
cfg.circleRadius = 0.01;               % [cm] = 100 um
cfg.targetLengthYAxis = 0.05;          % [cm]
cfg.targetAspectRatio = 2;             % X = 0.1 cm
cfg.singleCircleCenter = [0.05, 0.025];
cfg.circleSpacing = 0;
cfg.circleSpacingXLeft = 0;
cfg.characteristicLength = 0.02;       % 颗粒直径 [cm]，用于命名/Pe

% ---- 网格：简单算例用较粗网格（h = 10 um）----
cfg.meshTargetElementSizeCm = [];
cfg.meshNumPartitionsX = 100;
cfg.meshNumPartitionsY = 50;
cfg.numPartitionsMicroscale = 100;
cfg.dxfResolutionX = 200;
cfg.dxfResolutionY = 100;

% ---- 导出：从初始到完全溶解取 20 个切片 ----
cfg.targetDissolutionSlices = 20;
cfg.saveMeshDiagnostics = false;

names = fieldnames(overrides);
for i = 1:numel(names)
    cfg.(names{i}) = overrides.(names{i});
end

result = run_phreeqc_se(cfg);
end
