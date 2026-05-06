# AGENTS.md

## 项目定位

本仓库的主要任务是为 `shizaishiwo323/NMR-agent` 的 `IMGW-local-gpu-training` 分支准备训练、验证和测试数据集。

目标数据来自碳酸钙反应输运模拟，并结合 NMR 数值模拟与 T2 反演结果。这里不是优先开发通用软件功能，而是优先保证数据生成流程、数据质量、元数据记录和可复现实验。

## 领域背景

- 反应对象：以碳酸钙溶解/沉淀、孔隙结构演化、渗透率和孔隙率变化为核心。
- 数值流程：RTM/PNM 反应输运模拟生成随时间变化的孔隙/固体几何，COMSOL NMR 模型生成信号衰减数据，T2 反演得到谱分布和含水率/孔隙率相关标签。
- 下游用途：为 NMR-agent 提供可用于本地 GPU 训练的数据样本、监督标签、参数元数据和可追溯的模拟来源。

## 重点目录

- `ReactiveTransport/RTM/`：反应输运与孔隙网络模拟相关脚本。
- `ReactiveTransport/RTM/Batch_Simu/`：批量实验设计和批量模拟入口，通常生成 `BatchResults_*` 目录。
- `ReactiveTransport/NMR/`：NMR/COMSOL 模型文件，例如 `CT-simulation.mph`。
- `ReactiveTransport/automation/`：COMSOL 批处理、T2 反演、结果对比图和报告生成脚本。
- `Data/`：原始 RTSPHEM 示例数据和 CNN 训练数据，不要随意改写。

## 数据集产物约定

每个有效样本应尽量保留以下信息：

- 模拟参数：`Da`、`Pe`、`L`、`lengthXAxis`、`lengthYAxis`、几何布局类型、时间步、物理单位。
- 几何数据：`dxf_pore/`、`dxf_solid/` 中成对的 pore/solid DXF 文件，时间步必须可对应。
- 反应输运标签：`global_evolution.xlsx` 或等价表格中的 `time_s`、`porosity`、`permeability_mD`、`k_k0`、`surface_area_cm2`、`tortuosity` 等列。
- NMR 输入/输出：`comsol_results/` 中的 T2 信号表格，以及 `inversion_results/` 中的 `*_T2.mat`、`*_T2.png`、`porosity_comparison.xlsx/png`。
- 元数据：为每批数据生成或维护 manifest，至少记录样本路径、参数、时间步数量、处理状态、失败原因、生成脚本版本和生成日期。

推荐沿用现有文件夹命名：

```text
dissolution_results-Da_{Da}_Pe_{Pe}_L_{L}_lengthXAxis_{X}_lengthYAxis_{Y}_{layout}
```

推荐沿用现有子目录结构：

```text
dissolution_results-.../
  dxf_pore/
  dxf_solid/
  comsol_results/
  inversion_results/
  batch_logs/
```

## 工作原则

- 优先保护已有数据。不要覆盖 `.mat`、`.mph`、DXF、Excel、反演结果或批处理输出，除非用户明确要求。
- 修改 MATLAB/Python 脚本前，先确认输入目录、输出目录和路径配置，尤其是 `AutomationConfig.m` 中的 `project_root`、`data_root`、`mph_file`、`python_exe`。
- Windows 长路径问题已经在 `ReactiveTransport/automation/` 中有相关处理说明，新增脚本时要考虑路径长度和中文路径。
- 数据生成脚本应尽量支持断点续跑、跳过已完成样本、记录失败样本，不要因为单个样本失败中断整批任务。
- 输出给 NMR-agent 的数据应保持机器可读，优先使用 `.csv`、`.xlsx`、`.mat`、`.json` 或 `.npz`，并附带字段说明。
- 对每次批量生成或转换，记录运行环境：MATLAB 版本、COMSOL 版本、Python 环境、关键依赖、GPU/CPU 设置。

## 验证清单

在认为一批数据可交付前，至少检查：

- pore/solid DXF 文件数量一致，时间步编号一致。
- 文件夹名中的 `Da`、`Pe`、`L`、尺寸和 layout 能被 `parse_folder_name.m` 正确解析。
- `global_evolution.xlsx` 或等价文件存在，关键列非空且数值范围合理。
- COMSOL 输出的 T2 Excel 文件能被反演脚本读取。
- `*_T2.mat` 中包含 T2 bins、T2 log、combined spectrum、谱峰或累计含水率等下游训练需要的变量。
- `porosity_comparison` 图表/表格能反映 NMR 反演孔隙率与原始模拟孔隙率的对应关系。
- 随机抽样检查若干样本的图像、谱图和参数，排除明显错配。

## 代码修改偏好

- 保持现有 MATLAB 代码风格，优先做小范围、可回溯的修改。
- 新增脚本时，在文件头写清楚输入、输出、依赖和适用场景。
- 尽量不要引入新的大型依赖；如必须引入，说明用途和安装方式。
- 不要把绝对路径硬编码为唯一可用路径；需要默认值时，允许用户在配置文件中覆盖。
- 不要重构与当前数据准备无关的 RTSPHEM 原始算法代码。
- 处理数据转换时，优先保留原始文件，再生成派生文件。

## 禁止批量删除

禁止批量删除文件或目录。

不要使用：

- `del /s`
- `rd /s`
- `rmdir /s`
- `Remove-Item -Recurse`
- `rm -rf`

需要删除文件时，只能一次删除一个明确路径的文件，例如：

```powershell
Remove-Item "C:\path\to\file.txt"
```

如果需要批量删除文件，应停止操作，并请求用户手动删除。

## 常用入口

- 批量模拟：查看 `ReactiveTransport/RTM/Batch_Simu/RunBatchExperiments.m`。
- 单个/批量 NMR 后处理：查看 `ReactiveTransport/automation/run_single_folder.m` 和 `ReactiveTransport/automation/run_automation.m`。
- 自动化配置：查看 `ReactiveTransport/automation/AutomationConfig.m`。
- T2 反演：查看 `ReactiveTransport/automation/run_python_inversion.m` 和 `ReactiveTransport/automation/inversionBatch_matlab.py`。

## 交付给用户时

说明本次改动影响了哪些数据准备阶段，例如：

- RTM 参数设计
- DXF 几何导出
- COMSOL NMR 模拟
- T2 反演
- 数据清洗/manifest 生成
- 面向 NMR-agent 的训练数据打包

如果没有实际运行 MATLAB、COMSOL 或 Python 批处理，要明确说明未运行，以及原因。
