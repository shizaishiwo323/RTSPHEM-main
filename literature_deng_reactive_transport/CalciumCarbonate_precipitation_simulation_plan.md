# 碳酸钙沉淀模拟调研与模块接入方案

日期：2026-06-23

目标：在现有碳酸钙溶解 RTM/PNM 代码体系上，逐步加入可验证、可扩展的碳酸钙沉淀模拟能力。优先复用当前 PHREEQC 碳酸盐化学、level-set 几何演化、孔隙率/渗透率/扩散率后处理与 DXF/NMR 数据链路。

## 1. 快速结论

现有代码已经具备沉淀模块的基础：`run_single_pnm_mine.m` 已经可以通过 PHREEQC 处理 Calcite 碳酸盐化学，`PNM_beauty3.m` 已经保存每个三角单元的水体积、界面面积、方解石摩尔量、SI、pH、Ca/C/Na/Cl 等信息。因此，最稳妥的路线不是重写 RTM，而是把当前“只认溶解”的 PHREEQC 接口改成“有符号矿物反应”。

推荐分三层实现：

1. 第一层：表面生长/异质成核近似。把已有 calcite interface 上的 PHREEQC `KIN_DELTA("Calcite")` 作为有符号反应量。溶解时界面后退，沉淀时界面向孔隙推进。这一层最容易接入现有 level-set。
2. 第二层：壁面异质成核。对没有初始 calcite 或沉淀覆盖不足的孔壁，基于 SI、接触角和成核位点密度生成新沉淀位点，再让其表面生长。这一层需要新增每单元沉淀库存和成核 mask。
3. 第三层：孔隙内均质成核。对高 SI 孔隙体内生成晶核/颗粒，采用局部沉淀体积分数或微连续体孔隙率变量。它不能只靠现有 level-set 边界推进表达，因为会出现孤立固体岛、堵塞孔喉和概率性空间分布。

PNAS 的 `Fluid inertia controls mineral precipitation and clogging...` 不应直接并入现有 2D Stokes PNM 主线。它的关键机制是 3D 惯性二次流增强壁面可达性和 `Ω_surf`，现有 2D Stokes/level-set 框架不能自然复现 Dean flow 和 z 向表面输运。建议把 `even flow.mph` 作为 COMSOL 3D 标定分支：先复现表面反应热点和 reaction area vs Re，再决定是否把 Navier-Stokes 或惯性修正引入主代码。

## 2. 现有代码基础和主要缺口

### 2.1 已有能力

主要入口：

- `ReactiveTransport/RTM/run_single_pnm_mine.m`
- `ReactiveTransport/RTM/PNM_beauty3.m`
- `ReactiveTransport/RTM/couplePhreeqc/`

当前配置已经支持：

- `cfg.reactionModel = 'phreeqc'`
- `cfg.phreeqcRunGroup = 'phreeqc_database_calcite'`
- PHREEQC 数据库 `phreeqc-m.dat`
- 传输组分：`H`、`Ca_total`、`C_total`、`Na_total`、`Cl_total`
- 每单元几何状态：`interface_area_cm2`、`water_volume_cm3`、`calcite_moles`
- 每步 PHREEQC 输出：pH、Ca/C/Na/Cl、`calciteSI`、`KIN_DELTA_Calcite`、`RATE_Calcite`

现有 level-set 约定为：

- `phi > 0`：固体
- `phi < 0`：孔隙
- 更新式为 `phi_new = phi_old - dt * Vn * |grad phi|`

因此在当前符号约定下，正的 `Vn` 会使固体区缩小，适合表示溶解；沉淀使固体向孔隙推进时，几何速度应取负号。

### 2.2 关键缺口

当前 PHREEQC 接口在多个位置把反应硬编码为“溶解”：

- `ParsePhreeqcSelectedOutput.m` 中 `calciteDissolvedMoles = max(-calciteDeltaMoles, 0)`，正向沉淀量没有被保留。
- `RunPhreeqcCalciteBatch.m` 中 `scaleKineticDissolutionToCellInventory` 只对溶解做库存限幅。
- `ComputePhreeqcInterfaceRatePerArea.m` 用 `calciteDissolvedMoles` 推界面，返回值非负。
- `WritePhreeqcSpeciesTable.m` 只写出 `calciteDissolvedMoles`，没有写 `calcitePrecipitatedMoles` 和 signed reaction。
- 日志、图表、变量名多数是 dissolution 语义。

这意味着：即使 PHREEQC 算出了过饱和和沉淀，当前几何和输出也基本不会响应沉淀。

## 3. 文献机制总结

### 3.1 邓航老师 2026 综述：沉淀的基本理论

来源：

- Hang Deng and Jenna Poonoosamy, `Mineral precipitation in porous media systems: Controlling factors and impacts on porous media evolution`, Advances in Colloid and Interface Science, 2026. DOI: `10.1016/j.cis.2025.103745`
- 本地 PDF：`literature_deng_reactive_transport/pdfs/agent_recent/2026_Mineral_precipitation_porous_media_systems.pdf`

核心思想：

- 沉淀由热力学是否过饱和与动力学是否足够快共同控制。
- 饱和度通常写为 `Ω = IAP / Ksp`。PHREEQC 中更常用 `SI = log10(Ω)`，而部分理论推导用 `ln(Ω)`。
- 沉淀包含两个过程：成核和晶体生长。宏观模型中二者可以同时发生。
- 表面性质、孔径、粗糙度、已有同种矿物、局部流场和混合线都会改变沉淀位置和速率。

对碳酸钙模拟的启示：

- Calcite 沉淀不应只写成 dissolution 的反号。需要明确区分表面生长、异质成核和均质成核。
- CaCO3 可能有 ACC、aragonite、vaterite、calcite 等路径，若第一版只模拟 calcite，应在模型说明中写清楚“只追踪 calcite phase”。
- 小孔中的沉淀不仅受 SI 控制，还受成核概率、孔体积、孔径限制和局部表面功能团控制。

### 3.2 成核公式

对均质成核，可采用经典成核理论 CNT 作为第一版可调参数化：

```text
ΔG(r) = -(4πr^3 / 3v0) kB T ln(Ω) + 4πr^2 σ
rc = 2σv0 / (kB T ln(Ω))
ΔGc = B v0^2 σ^3 / (kB T ln(Ω))^2
Jhom = Γ exp(-ΔGc / kB T)
```

其中：

- `r`：晶核半径
- `v0`：单个分子或公式单元体积
- `σ`：晶体与溶液界面能
- `Ω`：饱和比
- `Jhom`：均质成核率，单位通常为 `#/m3/s`

对异质成核，可以用接触角降低能垒：

```text
f(θ) = (2 + cosθ)(1 - cosθ)^2 / 4
ΔGhet = f(θ) ΔGhom
Jhet = Γhet exp(-ΔGhet / kB T)
```

壁面成核的概率表达：

```text
Phet = 1 - exp(-Jhet Areactive Δt)
```

孔隙内均质成核的概率表达：

```text
Phom = 1 - exp(-Jhom Vwater Δt)
```

建议第一版不要直接追求绝对成核率精确预测，而是把 `N0`、`σ`、`θ`、`SIcrit` 作为校准参数。碳酸钙的非经典成核会让 CNT 只能作为有效模型，不应过度解释参数物理意义。

### 3.3 晶体生长和 TST 速率

常用的表面生长速率可写成：

```text
rA = kp (Ω - 1)^n, Ω > 1
rA = -kd (1 - Ω)^n, Ω < 1
```

其中 `rA` 是单位表面积矿物摩尔变化率。若定义 `rA > 0` 为沉淀，则几何速度为：

```text
Vn = -Vm rA
```

这里 `Vm` 是 calcite 摩尔体积。负号来自当前代码 level-set 的符号约定：沉淀要让固体区扩大。

如果沿用 PHREEQC database `RATES Calcite`，更推荐将 PHREEQC 的 signed `KIN_DELTA("Calcite")` 作为源项：

```text
Δn_calcite > 0  表示沉淀，矿物库存增加
Δn_calcite < 0  表示溶解，矿物库存减少
rA_signed = Δn_calcite / (Ainterface Δt)
Vn = -Vm rA_signed
```

为了兼容当前代码中“正速度为溶解”的写法，也可以内部定义：

```text
rA_dissolution_positive = -Δn_calcite / (Ainterface Δt)
Vn = Vm rA_dissolution_positive
```

这样溶解时 `Δn_calcite < 0`，`Vn > 0`；沉淀时 `Δn_calcite > 0`，`Vn < 0`。

### 3.4 邓航老师 2021 WRR：均质成核 vs 表面生长对扩散率的影响

来源：

- Deng et al., `A Pore-Scale Investigation of Mineral Precipitation Driven Diffusivity Change at the Column-Scale`, Water Resources Research, 2021. DOI: `10.1029/2020WR028483`
- 本地 PDF：`literature_deng_reactive_transport/pdfs/agent_precip_pore/2021_Pore_scale_precipitation_diffusivity.pdf`

该文把沉淀概念化为两类：

- HN：孔隙内均质成核提供初始表面积，随后晶体生长。
- SG：表面生长，等价于异质成核很快，晶体可立即在已有表面上生长。

重要启示：

- 均质成核和表面生长即使沉淀量相同，也会造成不同的堵塞路径和有效扩散率下降。
- 传统 Archie 关系难以捕捉临界堵塞，建议引入临界孔隙率 `φc`：

```text
Deff = Deff,c + a Dm (φ - φc)^n
```

这对本项目很有用：当前代码每步都显式计算几何和扩散/渗透率，短期不必强行使用经验关系；但当后续要上尺度或生成快速代理数据时，应输出并拟合 `φc`、`Deff,c`、`n`。

### 3.5 邓航老师 2024 Energy & Fuels：沉淀导致 k-φ 关系高度路径依赖

来源：

- Masoudi, Nooraiepour, Deng and Hellevang, `Mineral Precipitation and Geometry Alteration in Porous Structures: How to Upscale Variations in Permeability-Porosity Relationship?`, Energy & Fuels, 2024. DOI: `10.1021/acs.energyfuels.4c01432`
- 本地 PDF：`literature_deng_reactive_transport/pdfs/agent_recent/2024_mineral_precipitation_geometry_alteration.pdf`

关键点：

- 沉淀不是均匀缩孔。晶核数量、晶核位置和初始孔隙结构共同控制堵塞。
- 同样的孔隙率损失会导致非常不同的渗透率下降。
- `k/k0 = (φ/φ0)^n` 中的 `n` 可变化很大，不能简单沿用溶解问题中的 k-φ 经验式。
- 这篇文章适合作为沉淀训练数据设计的提醒：必须记录晶核数量、空间分布、局部 SI、局部速度和堵塞路径，而不只是记录总孔隙率。

### 3.6 邓航老师 2025 Langmuir：成核速率和对流速率的竞争

来源：

- Jiang et al., `Controls of the Nucleation Rate and Advection Rate on Barite Precipitation in Fractured Porous Media`, Langmuir, 2025. DOI: `10.1021/acs.langmuir.4c03532`
- 本地 PDF：`literature_deng_reactive_transport/pdfs/agent_recent/2025_Nucleation_advection_barite_precipitation.pdf`

该文提出一个简单但很适合做模型诊断的无量纲数：

```text
Da = J' / v
v = Pe D / b
Da = J' b / (Pe D)
```

其中：

- `J'` 是换算成速度量纲的成核/沉淀速率。
- `v` 是对流速率。
- `b` 是裂隙宽度。
- `D` 是扩散系数。

观测到两种 regime：

- Regime I：大 `Da`，沉淀主要在裂隙/岩石表面。
- Regime II：小 `Da`，沉淀偏向蚀变带或基质内部。

对本项目的启示：

- 可以把 `Da_nuc = J' / v` 或 `Da_precip = rA Vm / u` 作为每步输出指标。
- 对碳酸钙沉淀，pH、碱度和 Ca/C 比值改变 `Ω` 与成核率；流速改变成核前被带走的程度。
- 对已有碳酸钙溶解模型，可自然扩展到“先溶解释放 Ca，再在局部 pH/碳酸盐条件改变后二次沉淀”的耦合问题。

### 3.7 PNAS 2024：惯性控制混合诱导沉淀和堵塞

来源：

- Yang et al., `Fluid inertia controls mineral precipitation and clogging in pore to network-scale flows`, PNAS, 2024. DOI: `10.1073/pnas.2401318121`
- 本地 PDF：`C:\Users\imgw\Downloads\Fluid inertia controls mineral precipitation and clogging in pore to network-scale flows.pdf`
- 补充材料：`C:\Users\imgw\Downloads\pnas.2401318121.sapp (1).pdf`
- COMSOL 模型：`C:\Users\imgw\Documents\COMSOL\precipitation\even flow.mph`

该文的模拟设置：

- 3D COMSOL, steady reactive transport。
- 流场：不可压 Navier-Stokes。
- 传质：对流扩散。
- 反应：壁面异质沉淀，速率近似为 `R = k(Ω - 1)`。
- 关键输出：体相 `Ω_soln` 和壁面 `Ω_surf`。
- 几何不随沉淀动态更新，主要用反应区域预测实验沉淀图案。

主要机制：

- 低 Re 下，沉淀形成薄的混合屏障，阻断进一步混合，水力传导率下降有限。
- 高 Re 下，3D Dean flow 和回流把反应物带到壁面，显著提高 `Ω_surf`，交叉口快速堵塞。
- 惯性对表面反应的增强大于对体相反应的增强。

对本项目的边界条件：

- 现有 2D Stokes/PNM 很难复现这篇的核心 3D 惯性机制。
- 可以在 COMSOL 分支复现 `Ω_surf` 和 reaction area，再把经验的 `surface_accessibility_factor(Re, AR)` 作为后续近似修正。
- 若要真正并入主代码，需要从 Stokes 过渡到 Navier-Stokes 或至少引入惯性修正。HyPHM 目录中已有 NavierStokes 类，但当前 PNM 主线没有使用它。

## 4. 推荐的代码接入架构

### 4.1 新增配置

建议在 `run_single_pnm_mine.m` 或一个新入口 `run_single_pnm_precipitation.m` 中增加：

```matlab
cfg.reactionModel = 'phreeqc';
cfg.phreeqcRunGroup = 'phreeqc_database_calcite';
cfg.mineralEvolutionMode = 'signed_calcite_surface';

cfg.precipitation.enableSurfaceGrowth = true;
cfg.precipitation.enableHeterogeneousNucleation = false;
cfg.precipitation.enableHomogeneousNucleation = false;
cfg.precipitation.usePhreeqcSignedDelta = true;

cfg.precipitation.siThreshold = 0;
cfg.precipitation.siNucleationThresholdHet = 0.5;
cfg.precipitation.siNucleationThresholdHom = 1.5;
cfg.precipitation.contactAngleDeg = 60;
cfg.precipitation.nucleationSiteDensity_m2 = 1e12;
cfg.precipitation.nucleationSiteDensity_m3 = 1e18;
cfg.precipitation.interfacialEnergy_J_m2 = [];
cfg.precipitation.randomSeed = 323;
```

第一版建议关闭显式成核，只使用 PHREEQC signed surface reaction。这样能最快检验沉淀方向、质量守恒和几何演化。

### 4.2 PHREEQC 有符号矿物反应

建议新增或替换函数：

- `ParsePhreeqcSelectedOutput.m`
  - 保留 `calciteDeltaMoles`
  - 新增 `calcitePrecipitatedMoles = max(calciteDeltaMoles, 0)`
  - 保留 `calciteDissolvedMoles = max(-calciteDeltaMoles, 0)`
  - 新增 `calciteSignedMoles = calciteDeltaMoles`

- `RunPhreeqcCalciteBatch.m`
  - 将 `scaleKineticDissolutionToCellInventory` 改为 signed scaling。
  - 溶解限幅：`-Δn <= calcite_moles`
  - 沉淀限幅：`Δn` 不应超过水相可提供的 Ca/C 摩尔数，PHREEQC 已更新溶液浓度，但几何库存仍需一致。
  - 对没有初始 calcite 的孔隙单元，第一版不允许沉淀；第二版异质成核打开后才允许创建新 calcite 库存。

- `ComputePhreeqcInterfaceRatePerArea.m`
  - 建议改名为 `ComputeSignedCalciteInterfaceRatePerArea.m`
  - 返回 signed rate 或 dissolution-positive rate，必须在函数名和注释中固定约定。

推荐几何速度公式：

```matlab
signedDelta = phreeqcSpeciesData.calciteDeltaMoles; % + precip, - dissolve
rateDissPositive = -signedDelta ./ max(dt, eps) ./ max(interfaceArea, eps);
normalSpeed = molarVolume * rateDissPositive;
```

这样沉淀时 `normalSpeed < 0`，固体向孔隙增长。

### 4.3 表面生长模型

第一版表面生长只在已有 calcite interface 上发生：

```text
Areactive = interface_area_cm2
ΔVsolid = Vm Δn_calcite
Vn = -Vm Δn_calcite / (Areactive Δt)
```

优点：

- 直接复用 `computePhreeqcGeometryState`。
- 直接复用 level-set 和后续孔隙率/渗透率/扩散率计算。
- 能先做质量守恒、符号方向和堵塞行为验证。

限制：

- 无法在孔隙体内凭空生成晶核。
- 对非 calcite 基底的异质成核还没有表达。
- 如果沉淀导致孤立晶体或桥联孔喉，level-set 单界面推进可能不足。

### 4.4 异质成核模型

第二版引入壁面成核：

状态变量：

```text
nucleiMask(T)                 是否已有沉淀晶核
secondaryCalciteMoles(T)      二次沉淀 calcite 摩尔量
secondarySurfaceAreaCm2(T)    二次沉淀表面积
nucleationEventCount(T)       成核事件数
```

触发条件：

```text
SI > SIcrit_het
water_volume_cm3 > water_min
available_wall_area_cm2 > area_min
rand < 1 - exp(-Jhet Aavailable Δt)
```

壁面可用面积可以先近似为当前 solid-pore interface 面积。若需要区分 calcite 原生表面与非反应基底，可从 DXF layer 或 level-set 初始矿物类型扩展。

异质成核后，沉淀库存应进入 `calcite_moles_total = primaryCalciteMoles + secondaryCalciteMoles`，但输出中要分开记录 primary 与 secondary，便于分析溶解-再沉淀。

### 4.5 均质成核模型

第三版引入孔隙内均质成核。建议先用微连续体体积分数，而不是立即把每个晶核变成显式 level-set 小岛：

状态变量：

```text
precipitateVolumeFraction(T)  孔隙单元内沉淀占据体积分数
mobileParticleMoles(T)        可选，若考虑悬浮颗粒迁移
immobileParticleMoles(T)      默认，沉积/堵塞颗粒
```

触发条件：

```text
SI > SIcrit_hom
rand < 1 - exp(-Jhom Vwater Δt)
```

几何/输运影响：

- 局部孔隙率：`φ_local = φ_pore - precipitateVolumeFraction`
- 局部扩散：可用 Deng 2021 的临界孔隙率修正式。
- 局部渗透率：可记录用于后处理，主代码短期仍可由显式几何计算；若沉淀体积分数尚未显式转成几何，则需引入微连续体修正。

当 `precipitateVolumeFraction` 超过阈值时，再把该区域转换成 level-set 固体或导出到二值掩膜。这样可避免一开始就处理大量孤立小晶体拓扑变化。

## 5. 碳酸钙沉淀场景设计

### 5.1 最小可验证场景：均一过饱和流入

目的：验证 PHREEQC signed delta、界面运动方向、质量守恒。

设置：

- 初始几何：简单圆形 calcite grain 或当前 square/random 几何。
- 流体：均一 supersaturated Ca-HCO3/CO3 溶液。
- 只启用 `signed_calcite_surface`。
- 不启用显式成核。

预期：

- `calciteSI > 0`
- `calciteDeltaMoles > 0`
- `normalSpeed < 0`
- 固体体积增加，孔隙率和渗透率下降。

### 5.2 溶解后再沉淀场景

目的：贴近现有碳酸钙溶解体系。

设置：

- 第一阶段：酸性流体溶解 calcite，释放 Ca 和 C。
- 第二阶段：提高 pH 或注入碳酸盐/碱性流体，使局部 SI 转正。
- 观察沉淀是否在下游、低速区、孔喉或原 calcite 表面发生。

预期输出：

- 同一几何中出现先孔隙扩大、后局部缩孔。
- k-φ 曲线出现路径依赖或 hysteresis。
- 沉淀区与局部 `SI`、`Da_precip`、速度场相关。

### 5.3 混合诱导 CaCO3 沉淀

目的：复现 PNAS/经典 mixing-induced precipitation 的思想，但改成碳酸钙。

候选反应流体：

- Inlet A：CaCl2 溶液。
- Inlet B：Na2CO3 或 NaHCO3 溶液。
- 用 PHREEQC 计算 pH、碳酸盐配分和 calcite SI。

现有主代码可能只有单入口边界，若要做混合线，需要新增多入口边界或先在 COMSOL 中做标定。

## 6. COMSOL `even flow.mph` 的复用路线

当前文件：

```text
C:\Users\imgw\Documents\COMSOL\precipitation\even flow.mph
```

当前 MATLAB 探测状态：

- MATLAB R2025b 可用。
- `C:\Program Files\COMSOL\COMSOL63\Multiphysics\mli` 存在。
- 直接 `mphload` 需要 COMSOL server，当前未连接到 server，因此本次没有读取模型内部 physics。

建议后续只读抽取步骤：

```matlab
addpath('C:\Program Files\COMSOL\COMSOL63\Multiphysics\mli')
mphstart(2036)
model = mphload('C:\Users\imgw\Documents\COMSOL\precipitation\even flow.mph')
model.component.tags
model.study.tags
model.param.varnames
```

需要抽取的信息：

- 几何尺寸、深度、交叉口宽度、aspect ratio。
- physics tags：Laminar Flow / Transport of Diluted Species / boundary reaction。
- 入口速度或流量参数与 Re 的对应关系。
- 反应速率常数 `k`、扩散系数 `D`、`Ksp` 和 `Ω` 定义。
- 输出变量：`Ω_surf`、surface reaction rate、reaction area。

推荐复现实验：

1. 保持 `even flow.mph` 原样，扫 Re = 10, 30, 50, 100, 150。
2. 导出壁面 `Ω_surf` 和 reaction rate map。
3. 用二值化 reaction area 对比 PNAS 图 1 和 Fig. S4。
4. 若改成 CaCO3，把反应物从 Ba/SO4 改成 Ca/CO3，先不更新几何，只比较沉淀热点。

不建议立即做的事：

- 不建议直接用 COMSOL 替代 RTSPHEM 主流程。
- 不建议把 3D 惯性机制硬塞进当前 2D Stokes PNM。
- 不建议在未验证符号和质量守恒前让几何随沉淀动态更新。

## 7. 输出和数据结构建议

`global_evolution.xlsx` 建议新增字段：

```text
mean_calcite_SI
max_calcite_SI
min_calcite_SI
calcite_dissolved_mol_step
calcite_precipitated_mol_step
calcite_net_delta_mol_step
calcite_dissolved_mol_cumulative
calcite_precipitated_mol_cumulative
calcite_net_delta_mol_cumulative
precipitation_volume_cm3
dissolution_volume_cm3
net_solid_volume_change_cm3
mean_precipitation_rate_mol_cm2_s
max_precipitation_rate_mol_cm2_s
nucleation_event_count
active_nucleation_cells
clogged_cell_fraction
critical_throat_width_cm
Da_precip_mean
Da_precip_max
```

`phreeqc_species_*.csv` 建议新增字段：

```text
calciteDeltaMoles
calciteDissolvedMoles
calcitePrecipitatedMoles
calciteSignedRate_mol_cm2_s
calciteSurfaceGrowthRate_mol_cm2_s
nucleationProbabilityHet
nucleationProbabilityHom
secondaryCalciteMoles
precipitateVolumeFraction
```

这些字段对后续 NMR-agent 训练很重要，因为沉淀导致的孔隙率/渗透率变化比溶解更强烈依赖空间位置和局部堵塞。

## 8. 分阶段实施计划

### Phase 0：符号和 0D 化学验证

目标：

- 确认 PHREEQC 对 calcite 沉淀时 `KIN_DELTA("Calcite")` 的符号。
- 确认当前 level-set 中 `normalSpeed < 0` 会让固体增长。

验证：

- 单 cell PHREEQC：欠饱和水溶解 calcite，过饱和水沉淀 calcite。
- 单圆孔/单 calcite grain：手动给常数 `normalSpeed > 0` 和 `< 0`，检查固体面积变化方向。

### Phase 1：signed surface growth

目标：

- 完成最小可用沉淀模块。

改动范围：

- `ParsePhreeqcSelectedOutput.m`
- `RunPhreeqcCalciteBatch.m`
- `ComputePhreeqcInterfaceRatePerArea.m`
- `WritePhreeqcSpeciesTable.m`
- `PNM_beauty3.m` 中 PHREEQC 几何速度与日志字段
- `couplePhreeqc/tests/test_CalcitePhreeqcHelpers.m`

验收标准：

- 过饱和入口导致孔隙率下降。
- `calcitePrecipitatedMoles` 与水相 Ca/C 减少量质量守恒。
- 欠饱和入口仍能复现原溶解行为。
- 沉淀和溶解不会在同一步被错误截断为非负溶解量。

### Phase 2：heterogeneous nucleation

目标：

- 允许在孔壁或指定基底上产生新沉淀位点。

新增：

- `ComputeCalciteNucleationProbability.m`
- `UpdateCalciteNucleationState.m`
- per-cell `secondaryCalciteMoles` 和 `nucleiMask`

验收标准：

- 在 `SI > SIcrit_het` 且有可用壁面面积时出现概率性成核。
- 固定 random seed 后结果可复现。
- 晶核数量增加会导致更快堵塞和更陡的 k-φ 下降。

### Phase 3：homogeneous nucleation / micro-continuum pore filling

目标：

- 表达孔隙体内晶核与堵塞，而不要求每个晶体都显式几何化。

新增：

- `precipitateVolumeFraction`
- 局部孔隙率/扩散率修正
- `Phom = 1 - exp(-Jhom Vwater Δt)`

验收标准：

- 高 SI 下孔隙内出现沉淀体积分数。
- 能复现 Deng 2021 中 HN 比 SG 更容易造成扩散率临界下降的趋势。

### Phase 4：COMSOL/PNAS 标定分支

目标：

- 从 `even flow.mph` 复现惯性控制的壁面反应热点。

步骤：

- 启动 COMSOL server。
- 只读加载 mph，导出参数和变量。
- 扫 Re/AR，输出 `Ω_surf`、reaction area。
- 若做 CaCO3 版本，使用相同几何和流场，只替换反应物与 PHREEQC 或等效 `Ω`。

验收标准：

- 低 Re 反应区沿混合线集中。
- 高 Re 壁面反应面积显著扩大。
- 可得到一个经验函数 `F_access(Re, AR)`，用于后续 2D 模型的表面可达性修正。

## 9. 风险和注意事项

1. 符号风险最高。当前正 `normalSpeed` 对应溶解，沉淀必须反号；实施前必须做单步几何测试。
2. PHREEQC 的 `KIN_DELTA` 与溶液组分更新要保持一致，不能只改几何不改水相，也不能只改水相不改几何。
3. Calcite 沉淀可能经过 ACC 或其他多晶型，第一版只模拟 calcite 时要限定适用范围。
4. 成核参数强依赖表面、杂质、pH 和离子强度。CNT 应作为有效参数化，而不是绝对预测。
5. 当前 2D Stokes 主线不能解释 PNAS 的 3D 惯性二次流。PNAS 机制应通过 COMSOL 分支标定。
6. 沉淀导致 k-φ 关系强路径依赖，不能沿用溶解数据集中的简单单调关系作为唯一标签。

## 10. 建议优先级

最高优先级：

1. 做 PHREEQC signed calcite delta。
2. 做 level-set 符号测试。
3. 做均一过饱和 calcite surface growth 验证。

第二优先级：

4. 做 pH/碱度/Ca-C 输入场配置，形成沉淀专用 run script。
5. 输出沉淀字段和质量守恒日志。
6. 加入异质成核概率模型。

第三优先级：

7. 均质成核与微连续体孔隙内沉淀。
8. COMSOL `even flow.mph` Re/AR 标定。
9. 将 COMSOL 标定得到的表面可达性修正用于 2D PNM 近似模型。

## 11. 参考文献与本地资料

- Deng, H., and Poonoosamy, J. Mineral precipitation in porous media systems: Controlling factors and impacts on porous media evolution. Advances in Colloid and Interface Science, 2026. DOI: `10.1016/j.cis.2025.103745`.
- Deng, H., Tournassat, C., Molins, S., Claret, F., and Steefel, C. I. A Pore-Scale Investigation of Mineral Precipitation Driven Diffusivity Change at the Column-Scale. Water Resources Research, 2021. DOI: `10.1029/2020WR028483`.
- Masoudi, M., Nooraiepour, M., Deng, H., and Hellevang, H. Mineral Precipitation and Geometry Alteration in Porous Structures: How to Upscale Variations in Permeability-Porosity Relationship? Energy & Fuels, 2024. DOI: `10.1021/acs.energyfuels.4c01432`.
- Jiang, Q.-R., Hu, R., Deng, H., Ling, B., Yang, Z., and Chen, Y.-F. Controls of the Nucleation Rate and Advection Rate on Barite Precipitation in Fractured Porous Media. Langmuir, 2025. DOI: `10.1021/acs.langmuir.4c03532`.
- Yang, W., Chen, M. A., Lee, S. H., and Kang, P. K. Fluid inertia controls mineral precipitation and clogging in pore to network-scale flows. PNAS, 2024. DOI: `10.1073/pnas.2401318121`.
- COMSOL model: `C:\Users\imgw\Documents\COMSOL\precipitation\even flow.mph`.
