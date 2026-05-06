# NMR T2 反演自动化框架

## 概述

该自动化框架用于批量处理 `ALLDissolutionResults` 目录中的溶解模拟结果，自动执行以下流程：

1. **COMSOL 几何处理**：读取 DXF 文件，构建几何，求解模型，导出结果表格
2. **Python T2 反演**：读取结果表格，进行 NMR T2 反演，生成分析图像

## 目录结构

```
automation/
├── run_automation.m          # 主批处理脚本
├── run_single_folder.m       # 单文件夹处理脚本
├── AutomationConfig.m        # 配置类
├── scan_dissolution_folders.m # 扫描文件夹函数
├── parse_folder_name.m       # 解析文件夹名称函数
├── get_dxf_files.m          # 获取DXF文件函数
├── extract_timestep.m       # 提取时间步函数
├── run_comsol_processing.m  # COMSOL处理函数
├── run_python_inversion.m   # Python反演函数
├── generate_summary_report.m # 生成报告函数
├── run_t2_process_inversion.py # T2_process反演桥接脚本
└── QUICKSTART.md           # 本文档
```

## 快速开始

### 1. 环境准备

确保以下软件已安装并配置：

- **MATLAB** (已安装并可运行)
- **COMSOL Multiphysics 6.3** (默认路径: `C:\Program Files\COMSOL\COMSOL63\Multiphysics`)
- **Python** (Anaconda, 默认路径: `C:\ProgramData\anaconda3\python.exe`)
- Python 依赖库：`pandas`, `numpy`, `matplotlib`, `scipy`

### 2. 启动 COMSOL 服务器

在运行自动化脚本之前，需要先启动 COMSOL 服务器：

```
方法1: COMSOL GUI 中选择 "文件" -> "客户端/服务器" -> "COMSOL Multiphysics 服务器"
方法2: 命令行运行 comsolmphserver -port 2036
```

### 3. 修改配置 (如需要)

打开 `AutomationConfig.m`，根据您的环境修改路径：

```matlab
% 项目根目录
obj.project_root = 'C:\Users\imgw\Documents\MATLAB\RTSPHEM-main';

% COMSOL配置
obj.comsol_path = 'C:\Program Files\COMSOL\COMSOL63\Multiphysics';
obj.mph_file = fullfile(obj.project_root, 'NMR', 'CT-simulation.mph');

% Python配置
obj.python_exe = 'C:\ProgramData\anaconda3\python.exe';
obj.inversion_script = fullfile(obj.project_root, 'automation', 'run_t2_process_inversion.py');
```

### 4. 运行批处理

#### 方式A: 处理所有文件夹

```matlab
cd('C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\automation')
run_automation
```

#### 方式B: 处理单个文件夹

```matlab
cd('C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\automation')
run_single_folder
```

然后按提示输入目标文件夹名称。

## 输入数据格式

### 文件夹命名规则

```
dissolution_results-Da_{Da值}_Pe_{Pe值}_L_{L值}_lengthXAxis_{X值}_lengthYAxis_{Y值}_{布局类型}
```

例如：
```
dissolution_results-Da_0.04_Pe_1.00_L_0.001_lengthXAxis_0.055_lengthYAxis_0.045_square
```

### 子文件夹结构

```
dissolution_results-xxx/
├── dxf_pore/
│   ├── pore_t0001.dxf
│   ├── pore_t0002.dxf
│   └── ...
├── dxf_solid/
│   ├── solid_t0001.dxf
│   ├── solid_t0002.dxf
│   └── ...
├── comsol_results/     # 自动创建
│   ├── T2_xxx_t0001.xlsx
│   └── ...
└── inversion_results/  # 自动创建
    ├── T2_xxx_t0001.png
    └── ...
```

## 输出结果

### COMSOL 结果 (`comsol_results/`)

Excel 表格包含 NMR 信号衰减数据，命名格式：
```
T2_Da{Da}_Pe{Pe}_X{X}_Y{Y}_t{时间步}.xlsx
```

### 反演结果 (`inversion_results/`)

反演由 `ReactiveTransport/T2_process/nmr_t2` 工具包执行，数值模拟统一使用固定平滑/正则化因子 `0.01`，不使用 L-curve 搜索。输出包括：

- `*_T2.png`：采用 T2_process 绘图风格的衰减数据与 T2 谱配对图
- `*_T2.mat`：保留自动化后续对比所需的 `total_water`、`raw_spectrum_sum`、`calibration_factor` 等字段

### 日志文件 (`batch_logs/`)

- `batch_log_YYYYMMDD_HHMMSS.txt` - 处理日志
- `summary_report_YYYYMMDD_HHMMSS.txt` - 总结报告

## 常见问题

### Q: COMSOL 连接失败

确保 COMSOL 服务器正在运行，端口为 2036。

### Q: Python 反演失败

1. 检查 Python 路径是否正确
2. 确保已安装依赖库：
   ```bash
   pip install pandas numpy matplotlib scipy openpyxl
   ```

### Q: DXF 文件未识别

确保 DXF 文件命名符合 `pore_t{数字}.dxf` 和 `solid_t{数字}.dxf` 格式。

## 高级配置

### 禁用 COMSOL 处理 (仅运行反演)

```matlab
config = AutomationConfig();
config.enable_comsol = false;
```

### 禁用 Python 反演 (仅运行 COMSOL)

```matlab
config = AutomationConfig();
config.enable_inversion = false;
```

### 覆盖已存在的结果

```matlab
config = AutomationConfig();
config.overwrite_existing = true;
```
