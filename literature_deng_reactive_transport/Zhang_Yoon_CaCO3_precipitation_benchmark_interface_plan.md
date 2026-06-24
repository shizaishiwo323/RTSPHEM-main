# Zhang-Yoon 微流控 CaCO3 沉淀 Benchmark 与 RTSPHEM 接口构建方案

日期：2026-06-23

目标：在现有 RTSPHEM 碳酸钙溶解代码体系上，构建一个可验证的碳酸钙沉淀模拟接口。第一阶段以 Zhang et al. (2010) 的微流控 mixing-induced CaCO3 沉淀实验为 benchmark，并以 Yoon et al. (2012) 的孔尺度数值复现作为模拟对照。方案重点是实验几何、边界条件、流速、浓度、沉淀动力学接口、输出指标和与实验/已有模拟的定量对比。

## 0. 本次修订后的实现边界

本计划的代码边界固定为：

```text
C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate
```

实施原则：

1. 所有新增 MATLAB 代码、benchmark runner、沉淀适配器、测试脚本和结果比较脚本均放在 `ReactiveTransport/RTM/precipitate`。
2. 首轮不修改原始 `ReactiveTransport/RTM/PNM_beauty3.m`、`run_single_pnm_mine.m`、`couplePhreeqc/*.m` 或已有 batch/automation 入口。
3. 允许将 `ReactiveTransport/RTM/PNM_beauty3.m` 复制为 `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`，作为沉淀模拟专用主求解器；后续所有 split inlet、signed precipitation、benchmark 输出改动都在这份副本中完成。
4. 可以复用已有公共接口，例如 `BuildCalcitePhreeqcInput.m`、`ExportBinaryMaskToDXF.m`、`ReadLayeredDxfGeometry.m`、HyPHM/Transport/Stokes/level-set 类。`PNM_beauty3.m` 的 local functions 通过复制后的 `precip_PNM_beauty3.m` 保留下来，并在副本内部按沉淀需求改造。
5. 需要从旧代码拆出的通用沉淀逻辑，在 `precipitate` 目录中用 `precip_` 前缀实现独立 helper，避免函数名覆盖或改变原有溶解模拟行为。

可行性审核结论：

```text
物理 benchmark 路线：可行
复制 PNM_beauty3 为沉淀专用脚本的路线：最可行，能保留 local helper 且不影响原溶解主线
完全从零重构独立 runner 的路线：可行但成本更高，不作为首选
直接调用原始 PNM_beauty3 并实现 split inlet + signed precipitation：不可行，除非修改原始 PNM_beauty3
严格复现 Zhang/Yoon 全部后期沉淀下降：首轮不可保证，应作为第二阶段诊断目标
```

## 1. 核心结论

最适合作为当前代码沉淀模块验证对象的是：

```text
Zhang et al. 2010 微流控 CaCO3 沉淀实验
+ Yoon et al. 2012 孔尺度 LBM/FVM 复现
```

推荐首版 benchmark：

```text
25 mM CaCl2 + 25 mM Na2CO3
Darcy velocity = 1.25 cm/min = 0.020833 cm/s
geometry = Zhang 2010 micromodel pore network
comparison times = 13 min, 18 min, 118 min
comparison metric = precipitated area in first pore / first three pore bodies
```

首版沉淀模块不应先追求复杂均质成核，而应先完成：

1. 双入口横向混合边界。
2. PHREEQC `KIN_DELTA("Calcite")` 的 signed mineral reaction 接口。
3. 表面生长/已有表面沉淀的 level-set 几何推进。
4. 沉淀面积-时间曲线与 Zhang 实验、Yoon Case 1/Case 5 的对比。

## 2. 文献依据与用途

### 2.1 Zhang et al. 2010：实验 benchmark

Zhang, C., Dehoff, K., Hess, N., Oostrom, M., Wietsma, T. W., Valocchi, A. J., Fouke, B. W., and Werth, C. J. (2010). Pore-Scale Study of Transverse Mixing Induced CaCO3 Precipitation and Permeability Reduction in a Model Subsurface Sedimentary System. Environmental Science & Technology, 44, 7833-7838. DOI: https://doi.org/10.1021/es1019788

用途：

- 提供微流控孔隙网络几何尺寸。
- 提供 CaCl2/Na2CO3 浓度、pH、饱和度和流速。
- 提供沉淀图像、沉淀面积随时间变化、堵塞和后期溶解/重构证据。

### 2.2 Yoon et al. 2012：已有模拟复现

Yoon, H., Valocchi, A. J., Werth, C. J., and Dewers, T. (2012). Pore-scale simulation of mixing-induced calcium carbonate precipitation and dissolution in a microfluidic pore network. Water Resources Research, 48, W02524. DOI: https://doi.org/10.1029/2011WR011192

用途：

- 提供孔尺度模拟复现框架。
- 提供网格尺寸、流场求解、反应输运和沉淀面积统计方法。
- 提供 Case 1 和 Case 5 的结果，用于判断当前代码能否复现早期沉淀和后期面积下降。

### 2.3 PHREEQC：水化学与矿物动力学接口

Parkhurst, D. L., and Appelo, C. A. J. (2013). Description of input and examples for PHREEQC version 3. U.S. Geological Survey Techniques and Methods 6-A43. DOI: https://doi.org/10.3133/tm6A43

用途：

- 计算 pH、Ca/C/Na/Cl speciation、Calcite saturation index。
- 通过 KINETICS 和 `KIN_DELTA("Calcite")` 提供矿物沉淀/溶解净反应量。

### 2.4 Deng 相关沉淀理论文献

Deng et al. (2021), Water Resources Research, DOI: https://doi.org/10.1029/2020WR028483

- 支持区分 surface growth 和 homogeneous nucleation 两类沉淀模式。
- 支持沉淀位置/方式会显著改变扩散率和堵塞路径。

Masoudi, Nooraiepour, Deng, and Hellevang (2024), Energy & Fuels, DOI: https://doi.org/10.1021/acs.energyfuels.4c01432

- 支持沉淀导致的 k-phi 关系具有强路径依赖。
- 支持不能只用总孔隙率下降评价沉淀模拟。

Jiang et al. (2025), Langmuir, DOI: https://doi.org/10.1021/acs.langmuir.4c03532

- 支持用成核/沉淀速率与对流速率竞争来解释沉淀位置。
- 可作为后续 `Da_nuc`、`Da_precip` 诊断指标依据。

Yang et al. (2024), PNAS, DOI: https://doi.org/10.1073/pnas.2401318121

- 可作为后续 COMSOL 3D 惯性沉淀分支参考。
- 不建议直接并入首版 2D Stokes/PNM 主线。

## 3. Benchmark 实验设置

### 3.1 几何尺寸

Zhang 实验的微流控模型为 2D model subsurface sedimentary system：

| 参数 | 实验值 | RTSPHEM 单位 |
|---|---:|---:|
| 孔隙网络长度 | 2 cm | 2.0 cm |
| 孔隙网络宽度 | 1 cm | 1.0 cm |
| 深度 | 17-21 um | 建议 0.002 cm |
| 圆柱颗粒直径 | 300 um | 0.03 cm |
| 孔体宽度 | 180 um | 0.018 cm |
| 孔喉宽度 | 40 um | 0.004 cm |
| 孔隙率 | 0.39 | 0.39 |

Yoon 复现中主要关注靠近入口的 first pore body 和 first three pore bodies。因此推荐两级几何：

### 3.2 几何级别 A：入口局部 benchmark

用途：快速验证沉淀面积-时间曲线和图像对比。

```matlab
cfg.benchmarkName = 'zhang_yoon_2010_2012_local';
cfg.layoutType = 'zhang2010_micromodel_local';
cfg.lengthXAxis = 0.30;      % cm, 约覆盖入口附近前三个孔体
cfg.lengthYAxis = 0.30;      % cm, 覆盖上下双入口混合区
cfg.thickness = 0.002;       % cm
cfg.postDiameter = 0.03;     % cm
cfg.poreBody = 0.018;        % cm
cfg.poreThroat = 0.004;      % cm
cfg.targetPorosity = 0.39;
```

优点：

- 计算量小。
- 最贴近 Yoon 图像和面积曲线比较。
- 适合作为单元测试和接口验收 benchmark。

### 3.3 几何级别 B：完整实验域

用途：复现 Zhang 实验中更长距离的混合带、堵塞和通道重分配。

```matlab
cfg.benchmarkName = 'zhang_yoon_2010_2012_full';
cfg.layoutType = 'zhang2010_micromodel_full';
cfg.lengthXAxis = 2.0;       % cm
cfg.lengthYAxis = 1.0;       % cm
cfg.thickness = 0.002;       % cm
cfg.postDiameter = 0.03;     % cm
cfg.poreBody = 0.018;        % cm
cfg.poreThroat = 0.004;      % cm
cfg.targetPorosity = 0.39;
```

推荐优先级：

1. 先做几何级别 A。
2. 验证 signed reaction、沉淀面积、图像输出后再做完整域。
3. 若周期阵列无法准确表达实验孔喉，应从 Zhang/Yoon 图像数字化生成外部 DXF 几何。

## 4. 流场边界条件

### 4.1 实验边界

实验中 CaCl2 和 Na2CO3 分别从两个入口进入，等流量并排流入孔隙网络，在中心混合带诱导 CaCO3 沉淀。Darcy velocity 为：

```text
uD = 1.25 cm/min = 0.020833333 cm/s
```

第三股 HCl 用于防止出口管路堵塞，不作为首版模拟域内主反应条件。若后续要复现实验出口管路行为，可在出口下游增加 acid guard channel；首版不建议加入。

### 4.2 RTSPHEM 双入口设置

新增边界模式：

```matlab
cfg.flowBoundaryMode = 'split_left_inlet';
cfg.flowDirection = 'left_to_right';
cfg.splitInletY = 0.5 * cfg.lengthYAxis;
cfg.inletVelocity = 0.020833333;   % cm/s
```

流场：

```text
left lower inlet  : u = (uD, 0)
left upper inlet  : u = (uD, 0)
right outlet      : pressure outlet / reference pressure
top wall          : no-slip
bottom wall       : no-slip
solid posts       : no-slip
```

如果两个入口各占左边界一半高度，且两个入口都赋同一 `uD`，则整个左边界平均 Darcy velocity 仍为 `uD`。

## 5. 传质边界条件与浓度设置

### 5.1 25 mM 主 benchmark

Zhang/Yoon 主复现组：

```text
CaCl2 inlet:  Ca = 25 mM, Cl = 50 mM, pH = 6.1
Na2CO3 inlet: C_total = 25 mM, Na = 50 mM, pH = 10.9
```

RTSPHEM 使用 `mol/cm3`：

```text
25 mM = 25e-3 mol/L = 25e-6 mol/cm3
50 mM = 50e-6 mol/cm3
```

推荐配置：

```matlab
cfg.inletA.name = 'CaCl2';
cfg.inletA.yRange = [0, cfg.splitInletY];
cfg.inletA.H_total  = 10^(-6.1) * 1e-3;  % mol/cm3
cfg.inletA.Ca_total = 25e-6;
cfg.inletA.C_total  = 0;
cfg.inletA.Na_total = 0;
cfg.inletA.Cl_total = 50e-6;

cfg.inletB.name = 'Na2CO3';
cfg.inletB.yRange = [cfg.splitInletY, cfg.lengthYAxis];
cfg.inletB.H_total  = 10^(-10.9) * 1e-3; % mol/cm3
cfg.inletB.Ca_total = 0;
cfg.inletB.C_total  = 25e-6;
cfg.inletB.Na_total = 50e-6;
cfg.inletB.Cl_total = 0;
```

扩散系数：

```matlab
cfg.diffusionCoefficient = 0.9e-5;  % cm2/s, Yoon 2012
```

如果为了和现有代码统一，也可以先用 `1e-5 cm2/s`，但 benchmark 文档中应记录差异。

### 5.2 饱和度系列扩展

Zhang 实验还包含 6.5、10、50 mM。用于扩展验证：

| 组别 | CaCl2 / Na2CO3 浓度 | 用途 |
|---|---:|---|
| low | 6.5 mM | 宽沉淀带、低速增长 |
| mid-low | 10 mM | 饱和度敏感性 |
| main | 25 mM | Yoon 主复现组 |
| high | 50 mM | 快速沉淀、后期重构明显 |

首版只实现 25 mM；第二阶段加入浓度扫描。

## 6. 当前代码基础与缺口

### 6.1 已具备能力

现有主入口：

- `ReactiveTransport/RTM/run_single_pnm_mine.m`
- `ReactiveTransport/RTM/PNM_beauty3.m`
- `ReactiveTransport/RTM/couplePhreeqc/`

已具备：

- PHREEQC calcite chemistry。
- `H`、`Ca_total`、`C_total`、`Na_total`、`Cl_total` 多组分输运。
- 每单元 `interface_area_cm2`、`water_volume_cm3`、`calcite_moles`。
- PHREEQC 输出 pH、Ca/C/Na/Cl、`calciteSI`、`KIN_DELTA_Calcite`、`RATE_Calcite`。
- level-set 几何演化。
- 孔隙率、渗透率、DXF、NMR 后处理链路。

### 6.2 主要缺口

当前代码仍是溶解语义：

- `ParsePhreeqcSelectedOutput.m` 读入 `calciteDeltaMoles` 后，只保留 `calciteDissolvedMoles = max(-delta, 0)`。
- `RunPhreeqcCalciteBatch.m` 中 `scaleKineticDissolutionToCellInventory` 把 PHREEQC delta 强制变成负的溶解量。
- `ComputePhreeqcInterfaceRatePerArea.m` 返回非负溶解速率。
- `PNM_beauty3.m` 中 `normalSpeed = molarVolume * rate` 默认正值表示溶解。
- `createPhreeqcTransport` 只支持单入口统一浓度，不支持上下 split inlet。

因此，沉淀接口需要同时解决两个问题：

1. 多入口混合边界。
2. signed mineral reaction 与几何推进方向。

## 7. 沉淀反应理论与接口公式

### 7.1 饱和度

PHREEQC 输出：

```text
SI_Calcite = log10(Omega)
Omega = IAP / Ksp
```

判断：

```text
SI > 0  : calcite thermodynamically supersaturated
SI = 0  : equilibrium
SI < 0  : undersaturated
```

### 7.2 Signed mineral reaction

统一约定：

```text
calciteDeltaMoles > 0 : CaCO3 沉淀，矿物库存增加
calciteDeltaMoles < 0 : CaCO3 溶解，矿物库存减少
```

输出字段：

```text
calciteDeltaMoles
calcitePrecipitatedMoles = max(calciteDeltaMoles, 0)
calciteDissolvedMoles    = max(-calciteDeltaMoles, 0)
calciteNetRate_mol_s     = calciteDeltaMoles / dt
```

### 7.3 几何速度

当前 level-set 约定中，正 `normalSpeed` 对应溶解时固体后退。因此沉淀必须给负速度。

推荐内部使用 dissolution-positive speed：

```text
r_diss_pos = -Delta_n_calcite / (A_reactive Delta_t)
normalSpeed = Vm r_diss_pos
```

其中：

```text
Vm = 36.9 cm3/mol for calcite
A_reactive = interface_area_cm2
```

等价写法：

```text
r_signed = Delta_n_calcite / (A_reactive Delta_t)
normalSpeed = -Vm r_signed
```

### 7.4 表面生长模型

首版使用 PHREEQC signed delta，不另写速率律：

```text
Delta V_solid = Vm Delta n_calcite
```

适用范围：

- 已有 calcite 表面。
- 表面沉淀和再溶解。
- Zhang/Yoon benchmark 的首版复现。

限制：

- 不在无矿物表面的位置凭空创建晶核。
- 不能表达孔隙水体内部均质成核形成的孤立晶体。

### 7.5 异质成核接口

第二阶段引入：

```text
P_het = 1 - exp(-J_het A_available Delta_t)
J_het = J0_het exp(-f(theta) DeltaG_hom / kB T)
f(theta) = (2 + cos theta)(1 - cos theta)^2 / 4
```

推荐状态变量：

```text
nucleiMask
secondaryCalciteMoles
secondarySurfaceAreaCm2
nucleationEventCount
```

首版参数只作为有效参数，不作绝对 CNT 预测：

```matlab
cfg.precipitation.enableHeterogeneousNucleation = false;
cfg.precipitation.siNucleationThresholdHet = 0.5;
cfg.precipitation.contactAngleDeg = 60;
cfg.precipitation.nucleationSiteDensity_m2 = 1e12;
```

### 7.6 均质成核接口

第三阶段引入：

```text
P_hom = 1 - exp(-J_hom V_water Delta_t)
```

推荐采用 micro-continuum pore filling，而不是立即生成大量显式 level-set 小岛：

```text
precipitateVolumeFraction
immobileParticleMoles
cloggedCellFraction
```

局部孔隙率：

```text
phi_local = phi_pore - precipitateVolumeFraction
```

当 `precipitateVolumeFraction` 超过阈值后，再转成显式固体 mask 或堵塞单元。

## 8. `precipitate` 独立模块配置接口

### 8.1 Benchmark 配置

新增入口：

```text
ReactiveTransport/RTM/precipitate/run_zhang_yoon_caco3_precipitation_benchmark.m
```

核心配置：

```matlab
cfg.benchmarkName = 'zhang_yoon_caco3_25mM';
cfg.reactionModel = 'phreeqc';
cfg.phreeqcRunGroup = 'phreeqc_database_calcite';
cfg.mineralEvolutionMode = 'signed_calcite_surface';

cfg.flowBoundaryMode = 'split_left_inlet';
cfg.flowDirection = 'left_to_right';
cfg.inletVelocity = 0.020833333;
cfg.diffusionCoefficient = 0.9e-5;
cfg.molarVolume = 36.9;

cfg.precipitation.enableSurfaceGrowth = true;
cfg.precipitation.enableHeterogeneousNucleation = false;
cfg.precipitation.enableHomogeneousNucleation = false;
cfg.precipitation.usePhreeqcSignedDelta = true;
cfg.precipitation.allowRedissolution = true;
```

### 8.2 对比时间

```matlab
cfg.benchmarkComparisonTimes_s = [13, 18, 118] * 60;
cfg.endTime = 4 * 3600;
cfg.phreeqcExportEvery = 1;
cfg.exportEvery = 1;
```

### 8.3 输出目录

```text
outputs/rtm_runs/zhang_yoon_caco3_25mM_<timestamp>/
  run_metadata.json
  global_evolution.xlsx
  precipitation_metrics.csv
  benchmark_comparison/
    area_time_curve.png
    snapshots_13_18_118min.png
    yoon_case1_case5_comparison.png
    precipitate_masks/
  phreeqc_results/
  dxf_pore/
  dxf_solid/
```

## 9. `precipitate` 内拟新增函数

### 9.1 沉淀专用模块目录结构

```text
ReactiveTransport/RTM/precipitate/
  run_zhang_yoon_caco3_precipitation_benchmark.m
  precip_ConfigureZhangYoonBenchmark.m
  precip_PNM_beauty3.m
  precip_GenerateZhang2010MicromodelGeometry.m
  precip_GetSplitInletBoundaryConfig.m
  precip_CreateTransportMultiInlet.m
  precip_RunPhreeqcCalciteBatchSigned.m
  precip_ParsePhreeqcSelectedOutputSigned.m
  precip_ScaleSignedCalciteDeltaToCellInventory.m
  precip_ComputeSignedCalciteInterfaceRatePerArea.m
  precip_ComputePrecipitationMetrics.m
  precip_CompareZhangYoonBenchmark.m
  precip_WriteSpeciesTable.m
  precip_WriteRunMetadata.m
  tests/
    test_precip_signed_phreeqc.m
    test_precip_levelset_sign.m
    test_precip_split_inlet_transport.m
    test_precip_metrics.m
```

命名规则：

- 新文件统一使用 `precip_` 前缀，避免与主项目函数重名。
- `precip_PNM_beauty3.m` 由原始 `PNM_beauty3.m` 复制而来，只在 `precipitate` 目录中改造。
- runner 可以加入 `precipitate` 目录和现有 `RTM`、`couplePhreeqc`、`HyPHM`、`src` 路径；路径顺序应让 `precipitate` 优先于 `RTM`，避免误调用旧的同名 helper。
- `precipitate` 目录不应向上写入代码文件，只在自己的输出目录写结果。

### 9.2 新增函数职责

```text
precip_ConfigureZhangYoonBenchmark.m
```

生成 benchmark cfg，包括几何、双入口边界、PHREEQC、输出时间和对比指标。配置中应显式记录 `moduleRoot`、`rtmRoot`、`couplePhreeqcRoot` 和输出目录。

```text
precip_PNM_beauty3.m
```

沉淀模拟专用主求解器。该文件从 `ReactiveTransport/RTM/PNM_beauty3.m` 复制而来，在副本中完成以下改造：

```text
1. 将函数名改为 precip_PNM_beauty3
2. 增加 split inlet 边界和多入口 transport
3. 将 PHREEQC delta 改为 signed calcite reaction
4. 沉淀时允许 normalSpeed 为负
5. 增加 Zhang/Yoon benchmark 输出时间和沉淀面积指标
6. 输出 precipitate 专用 species table、metrics 和 comparison figures
```

保留 local helper 的好处是：`computePhreeqcGeometryState`、`collectPhreeqcState`、`applyPhreeqcResultToTransports`、`createPhreeqcTransport` 等原本只能在 `PNM_beauty3.m` 内部访问的函数，可以在副本中直接改造成沉淀版本，避免从零重写主循环。

```text
precip_GenerateZhang2010MicromodelGeometry.m
```

生成近似 Zhang 微流控 post array。支持 local/full 两种域。

```text
precip_GetSplitInletBoundaryConfig.m
```

返回 split inlet 的 Stokes 和 transport 边界谓词。

```text
precip_CreateTransportMultiInlet.m
```

替代单入口 `createPhreeqcTransport`，支持每个 species 的 piecewise inlet concentration。

```text
precip_RunPhreeqcCalciteBatchSigned.m
```

沉淀专用 PHREEQC runner。`precip_PNM_beauty3.m` 调用该函数替代旧的 `RunPhreeqcCalciteBatch.m`。该函数可复用 `BuildCalcitePhreeqcInput.m` 生成 PHREEQC 输入，但必须用 `precip_ParsePhreeqcSelectedOutputSigned.m` 和 `precip_ScaleSignedCalciteDeltaToCellInventory.m` 保留沉淀正值。

```text
precip_ParsePhreeqcSelectedOutputSigned.m
```

解析 `KIN_DELTA_Calcite`，保留：

```text
calciteDeltaMoles
calcitePrecipitatedMoles
calciteDissolvedMoles
calciteSignedRate_mol_s
```

```text
precip_ScaleSignedCalciteDeltaToCellInventory.m
```

替代当前 `RunPhreeqcCalciteBatch.m` 内部的 dissolution-only 缩放逻辑，保留沉淀和溶解方向。

```text
precip_ComputeSignedCalciteInterfaceRatePerArea.m
```

从 signed delta 得到 level-set 几何推进需要的 dissolution-positive rate。

`computePhreeqcGeometryState`、`collectPhreeqcState`、`applyPhreeqcResultToTransports` 等 local helper 不再强制拆成外部文件。首轮优先在 `precip_PNM_beauty3.m` 副本内直接改造，只有当某段逻辑需要被单元测试或多个 runner 复用时，再抽出为独立 `precip_*.m` helper。

```text
precip_ComputePrecipitationMetrics.m
precip_CompareZhangYoonBenchmark.m
```

从 level-set、PHREEQC 输出和 benchmark mask 中计算沉淀面积、峰值、后期下降、堵塞指标，并与 Zhang/Yoon 数据对比。

### 9.3 不修改旧函数的处理策略

首轮不直接修改这些文件：

```text
ReactiveTransport/RTM/PNM_beauty3.m
ReactiveTransport/RTM/run_single_pnm_mine.m
ReactiveTransport/RTM/couplePhreeqc/ParsePhreeqcSelectedOutput.m
ReactiveTransport/RTM/couplePhreeqc/RunPhreeqcCalciteBatch.m
ReactiveTransport/RTM/couplePhreeqc/ComputePhreeqcInterfaceRatePerArea.m
ReactiveTransport/RTM/couplePhreeqc/WritePhreeqcSpeciesTable.m
```

原因：

- 这些文件服务于既有溶解模拟和 PHREEQC-TST benchmark。
- 直接改动会改变已有数据生成链路的行为。
- 当前任务要求所有沉淀代码放在 `precipitate`。

对应替代方式：

```text
ParsePhreeqcSelectedOutput.m              -> precip_ParsePhreeqcSelectedOutputSigned.m
RunPhreeqcCalciteBatch.m                  -> precip_RunPhreeqcCalciteBatchSigned.m
ComputePhreeqcInterfaceRatePerArea.m      -> precip_ComputeSignedCalciteInterfaceRatePerArea.m
WritePhreeqcSpeciesTable.m                -> precip_WriteSpeciesTable.m
PNM_beauty3.m                             -> precip_PNM_beauty3.m
createPhreeqcTransport local function     -> precip_PNM_beauty3.m 内改造，必要时抽出 precip_CreateTransportMultiInlet.m
computePhreeqcGeometryState local function -> precip_PNM_beauty3.m 内保留
collectPhreeqcState local function        -> precip_PNM_beauty3.m 内保留
applyPhreeqcResultToTransports local function -> precip_PNM_beauty3.m 内保留
```

### 9.4 signed delta 缩放规则

```text
rawDeltaReference = KIN_DELTA from PHREEQC reference water
scale = water_volume_cm3 * 1e-3 / solutionWaterKg
deltaCell = rawDeltaReference * scale
```

限幅：

```text
if deltaCell < 0:
    abs(deltaCell) <= calcite_moles
if deltaCell > 0:
    deltaCell <= min(available_Ca_moles, available_C_moles)
```

几何速度：

```text
rateDissolutionPositive = -deltaCell / (interfaceAreaCm2 * dt)
normalSpeed = molarVolume * rateDissolutionPositive
```

其中沉淀时 `deltaCell > 0`，因此 `normalSpeed < 0`。

## 10. 双入口输运实现细节

现有 `createPhreeqcTransport` 中入口通量类似：

```matlab
gF = -inletConcentration * inletVelocity * (x(1) < EPS)
```

新函数应支持：

```matlab
function c = splitInletConcentration(x, speciesName, cfg)
    isLeft = x(1) < EPS;
    isLower = x(2) >= cfg.inletA.yRange(1) && x(2) < cfg.inletA.yRange(2);
    isUpper = x(2) >= cfg.inletB.yRange(1) && x(2) <= cfg.inletB.yRange(2);

    c = 0;
    if isLeft && isLower
        c = cfg.inletA.(speciesName);
    elseif isLeft && isUpper
        c = cfg.inletB.(speciesName);
    end
end
```

说明：CaCl2 位于 lower-y inlet、Na2CO3 位于 upper-y inlet 是本 RTSPHEM 局部坐标约定，用于实现 Zhang/Yoon 的双入口并排注入；不要把上下方向写成论文图像的绝对方向。共享分界点归 upper inlet，避免中线双计数。实际 benchmark 运行前需确认 split/yRange 端点落在左边界网格节点上，否则跨 split 的边界边会产生点值归类误差。

通量：

```matlab
transportProblem.gF.setdata( ...
    @(t, x) -splitInletConcentration(x, speciesName, cfg) .* inletVelocity .* (x(1) < EPS));
```

推荐 species map：

```text
H_total  -> hydrogenTransport
Ca_total -> calciumTransport
C_total  -> carbonTransport
Na_total -> sodiumTransport
Cl_total -> chlorideTransport
```

## 11. 输出指标

### 11.1 全局指标

`global_evolution.xlsx` 增加：

```text
calcite_dissolved_mol_step
calcite_precipitated_mol_step
calcite_net_delta_mol_step
calcite_dissolved_mol_cumulative
calcite_precipitated_mol_cumulative
calcite_net_delta_mol_cumulative
precipitation_volume_cm3
dissolution_volume_cm3
net_solid_volume_change_cm3
mean_calcite_SI
max_calcite_SI
min_calcite_SI
mean_precipitation_rate_mol_cm2_s
max_precipitation_rate_mol_cm2_s
clogged_cell_fraction
```

### 11.2 Benchmark 指标

`precipitation_metrics.csv`：

```text
timestep
time_s
time_min
precipitated_area_first_pore_cm2
precipitated_area_first_three_pores_cm2
precipitated_area_total_cm2
normalized_area_first_three_pores
peak_area_so_far
time_to_peak_min
late_stage_decline_fraction
```

### 11.3 图像指标

导出：

```text
precipitate_mask_0013min.png
precipitate_mask_0018min.png
precipitate_mask_0118min.png
si_calcite_0013min.png
pH_0013min.png
velocity_magnitude_0013min.png
```

图像对比：

```text
mask IoU
centroid offset
precipitation band width
precipitation band length
area ratio against experiment
area ratio against Yoon Case 1 / Case 5
```

## 12. 与 Zhang/Yoon 的对比方案

### 12.1 对比对象

三层对比：

1. Zhang 实验图像和面积曲线。
2. Yoon Case 1 普通模型。
3. Yoon Case 5 增强溶解模型。

### 12.2 首版验收标准

首版 surface growth 模型通过标准：

```text
1. 沉淀主要出现在两股流体混合界面附近。
2. 13 min 时已有明显沉淀带。
3. 18 min 沉淀面积不小于 13 min。
4. 118 min 能输出沉淀图像、pH、SI 和速度场。
5. calciteDeltaMoles 出现正值，且 Ca/C 水相减少与沉淀量方向一致。
6. 孔隙率下降，局部流速发生重分配。
```

不要求首版完美复现实验后期下降。Yoon 2012 已显示普通模型难以自然复现实验中快速下降，需要增强溶解或额外重构机制。

### 12.3 第二阶段验收标准

增强溶解/再溶解模型通过标准：

```text
1. 沉淀面积先升高后下降。
2. peak time 接近实验量级。
3. 30 min 前面积下降趋势接近 Yoon Case 5。
4. 118 min 沉淀面积不再持续无限增长。
```

参数扫描：

```text
kdiss_multiplier = [1, 10, 100, 300]
```

其中 `300` 对应 Yoon Case 5 的经验增强溶解思想，不应解释为普适物理常数。

## 13. 实施阶段

### Phase 0：单 cell signed PHREEQC 验证

目标：

- 确认 PHREEQC `KIN_DELTA("Calcite")` 符号。
- 确认过饱和输入产生 `calciteDeltaMoles > 0`。
- 确认欠饱和输入产生 `calciteDeltaMoles < 0`。

验证：

```text
single cell, calcite surface area > 0
case A: undersaturated acidic water
case B: supersaturated Ca-HCO3/CO3 water
```

### Phase 1：level-set 符号验证

目标：

- 确认 `normalSpeed > 0` 固体后退。
- 确认 `normalSpeed < 0` 固体增长。

验证：

```text
single circular grain
constant normalSpeed = +v
constant normalSpeed = -v
compare solid area change
```

### Phase 2：单入口均一过饱和沉淀

目标：

- 不处理混合线，只验证 signed surface growth。

设置：

```text
uniform supersaturated inlet
existing calcite surface
PHREEQC signed delta
```

通过标准：

```text
calciteDeltaMoles > 0
porosity decreases
permeability decreases
solid volume increases by Vm * precipitated moles
```

### Phase 3：Zhang/Yoon 双入口 local benchmark

目标：

- 复现入口附近 first three pore bodies 的沉淀。

设置：

```text
local Zhang geometry
split inlet
25 mM CaCl2 / 25 mM Na2CO3
uD = 0.020833 cm/s
D = 0.9e-5 cm2/s
```

输出：

```text
precipitated area vs time
13/18/118 min snapshots
pH/SI/velocity maps
```

### Phase 4：后期下降与 Yoon Case 5 对比

目标：

- 测试是否能复现实验中沉淀面积下降。

设置：

```text
allowRedissolution = true
kdiss_multiplier = [1, 10, 100, 300]
```

输出：

```text
area_time_curve_with_yoon_case1_case5.png
late_stage_decline_fraction
```

### Phase 5：异质成核

目标：

- 在无原始 calcite 或沉淀覆盖不足的壁面生成新沉淀位点。

新增状态：

```text
nucleiMask
secondaryCalciteMoles
secondarySurfaceAreaCm2
```

### Phase 6：均质成核与 pore filling

目标：

- 表达孔隙体内晶核和微连续体堵塞。

新增状态：

```text
precipitateVolumeFraction
immobileParticleMoles
cloggedCellFraction
```

## 14. 风险与处理

### 14.1 几何不确定性

实验几何不是简单规则阵列，`300 um` posts、`180 um` pore body 和 `40 um` pore throat 可能来自特定排列。若规则阵列无法复现沉淀带，应优先数字化 Zhang/Yoon 图像，生成外部 DXF。

### 14.2 反应符号风险

这是最高优先级风险。所有输出必须同时记录：

```text
calciteDeltaMoles
calcitePrecipitatedMoles
calciteDissolvedMoles
normalSpeed
```

### 14.3 沉淀相风险

Yoon 2012 主要考虑 vaterite formation。当前代码使用 Calcite。首版可继续用 Calcite，但文档中必须说明：

```text
This benchmark uses Calcite as an effective CaCO3 phase.
Experimental/Yoon dominant phase may be vaterite under the reported conditions.
```

若要严格复现 Yoon，需要在 PHREEQC 中增加 Vaterite phase 或用有效动力学参数近似。

### 14.4 后期溶解/重构风险

实验出现先沉淀后下降。普通表面沉淀模型可能不能自然复现。应把这作为模型诊断结果，而不是接口失败。

### 14.5 3D 效应风险

Zhang/Yoon benchmark 是 2D 微流控近似，PNAS 2024 的惯性沉淀堵塞机制涉及 3D Dean flow 和壁面可达性。当前首版不纳入 3D 惯性。

## 15. 可行性审核

### 15.1 总体判断

本 benchmark 作为独立沉淀模块是可行的，但它不是对原始 `PNM_beauty3.m` 加几个参数就能完成的轻量改造。原因是当前主函数把流场、输运初始化、PHREEQC 几何状态、PHREEQC 结果应用和 level-set 更新写在同一个长函数内，其中若干关键函数是 local function，外部无法直接复用。

因此，本计划的可行实现应是：

```text
复制 PNM_beauty3.m -> precipitate/precip_PNM_beauty3.m
+ 复用 HyPHM/level-set/公共 DXF/PHREEQC input builder
+ 在 precip_PNM_beauty3.m 副本中改造 split inlet、signed PHREEQC、negative normalSpeed 和 benchmark 输出
+ 将需要单元测试或复用的沉淀逻辑抽成 precip_*.m helper
```

不推荐的实现是：

```text
从 precipitate 目录直接调用 PNM_beauty3.m，期望通过 cfg 自动获得 split inlet 和 signed precipitation
```

因为当前原始 `PNM_beauty3.m` 尚未暴露这些 hook。复制副本后，hook 可以直接在 `precip_PNM_beauty3.m` 内部实现，而不会影响原溶解主线。

### 15.2 分项可行性

| 模块 | 可行性 | 审核结论 |
|---|---|---|
| 复制 `PNM_beauty3.m` 为 `precip_PNM_beauty3.m` | 高 | 最推荐路线，保留 local helper，同时隔离原项目功能 |
| Zhang/Yoon 25 mM local benchmark | 高 | 几何、流速、浓度、对比时间明确，适合首版 |
| split inlet 多组分输运 | 中-高 | 在 `precip_PNM_beauty3.m` 副本内改造 `createPhreeqcTransport` 或抽出 `precip_CreateTransportMultiInlet` |
| PHREEQC signed calcite delta | 高 | 可复用 `BuildCalcitePhreeqcInput.m`，但 runner/parser/scaling 使用沉淀专用版本 |
| level-set 沉淀几何推进 | 中 | 需要先做符号测试；沉淀时 `normalSpeed` 应为负 |
| 沉淀面积 benchmark 指标 | 高 | 可从 level-set/mask 直接计算 |
| 后期沉淀面积下降 | 中-低 | 普通 surface growth 未必自然复现，需要再溶解增强或相变/重构机制 |
| 严格相组成复现 | 中-低 | Zhang/Yoon 涉及 vaterite；当前首版用 calcite 作为有效 CaCO3 phase |
| 不影响已有溶解功能 | 高 | 前提是只复制旧文件，不修改旧文件；所有新逻辑留在 `precipitate` |

### 15.3 通过/停止门槛

继续进入 Zhang/Yoon local benchmark 前，必须先通过：

```text
1. single-cell PHREEQC signed delta 测试
2. level-set 正负 normalSpeed 符号测试
3. split inlet transport 保守性测试
4. 过饱和单入口 calcite surface growth 测试
```

若任一失败，不应直接跑 13/18/118 min benchmark。

### 15.4 对旧项目功能的保护

实现时应遵守：

```text
不改原始 PNM_beauty3.m
不改 run_single_pnm_mine.m
不改 couplePhreeqc/*.m
不改 automation/*.m
不覆盖已有 results、DXF、Excel、MAT、mph 文件
```

允许的操作：

```text
复制 ReactiveTransport/RTM/PNM_beauty3.m
生成 ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m
只修改 precipitate 目录下的副本和 helper
```

如果后续必须接入主流程，先形成一个独立 PR/补丁方案，只加入默认关闭的 hook，例如：

```matlab
cfg.experimentalPrecipitationAdapter = [];
```

默认值为空时，旧溶解流程输出应逐字节或数值等价于原流程。

## 16. 推荐交付物

首轮实现完成后，应生成：

```text
1. ReactiveTransport/RTM/precipitate/run_zhang_yoon_caco3_precipitation_benchmark.m
2. ReactiveTransport/RTM/precipitate/precip_ConfigureZhangYoonBenchmark.m
3. ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m
4. ReactiveTransport/RTM/precipitate/precip_GenerateZhang2010MicromodelGeometry.m
5. ReactiveTransport/RTM/precipitate/precip_GetSplitInletBoundaryConfig.m
6. ReactiveTransport/RTM/precipitate/precip_CreateTransportMultiInlet.m
7. ReactiveTransport/RTM/precipitate/precip_RunPhreeqcCalciteBatchSigned.m
8. ReactiveTransport/RTM/precipitate/precip_ComputeSignedCalciteInterfaceRatePerArea.m
9. ReactiveTransport/RTM/precipitate/precip_ScaleSignedCalciteDeltaToCellInventory.m
10. ReactiveTransport/RTM/precipitate/precip_ComputePrecipitationMetrics.m
11. ReactiveTransport/RTM/precipitate/precip_CompareZhangYoonBenchmark.m
12. benchmark result package with figures and metrics
```

## 17. 最小验收清单

代码层面：

```text
PHREEQC signed delta preserved
precipitated/dissolved/net moles exported
normalSpeed sign verified
split inlet transport verified
local Zhang geometry generated or imported
```

物理层面：

```text
CaCl2 and Na2CO3 mix near central transverse zone
Calcite SI becomes positive near mixing zone
Ca/C decrease where calciteDeltaMoles > 0
solid volume increases during precipitation
flow field updates after geometry change
```

benchmark 层面：

```text
area-time curve exported
13/18/118 min snapshots exported
comparison with Zhang/Yoon figures exported
case1-like and case5-like parameter sets documented
```

## 18. 建议下一步

优先执行顺序：

1. 复制 `ReactiveTransport/RTM/PNM_beauty3.m` 到 `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`。
2. 在副本中将主函数名改为 `precip_PNM_beauty3`，并确认 runner 调用副本而不是原始 `PNM_beauty3.m`。
3. 新建 `run_zhang_yoon_caco3_precipitation_benchmark.m` 和 `precip_ConfigureZhangYoonBenchmark.m`。
4. 先实现 `precip_RunPhreeqcCalciteBatchSigned.m`，只跑 single-cell signed PHREEQC。
5. 新建 `precip_ComputeSignedCalciteInterfaceRatePerArea.m` 和 `precip_ScaleSignedCalciteDeltaToCellInventory.m`。
6. 在 `precip_PNM_beauty3.m` 副本中接入 signed PHREEQC，并做 level-set 正负速度单元测试。
7. 在 `precip_PNM_beauty3.m` 副本中实现 split inlet transport，必要时抽出 `precip_CreateTransportMultiInlet.m`。
8. 跑 local Zhang/Yoon geometry。
9. 导出沉淀面积-时间曲线和 13/18/118 min 图像。
10. 与 Zhang 实验和 Yoon Case 1/Case 5 做图像与曲线对比。

首版成功的定义不是“完全拟合实验”，而是形成一个可重复、可诊断、可扩展的沉淀 benchmark：能清楚说明哪些现象由当前 surface growth + PHREEQC signed delta 已经解释，哪些现象需要再溶解增强、异质成核、均质成核或 3D 惯性机制。
