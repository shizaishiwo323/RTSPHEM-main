# RTSPHEM-main

碳酸钙反应输运（RTM/PNM）+ PHREEQC 水化学 + NMR 数值模拟的数据生成仓库，为
`shizaishiwo323/NMR-agent` 准备训练数据。项目约定和数据规范见 [AGENTS.md](AGENTS.md)，
工作流程见 [PROJECT_WORKFLOW.md](PROJECT_WORKFLOW.md)。

## 环境配置（新电脑按此顺序操作）

> 给 AI agent 的说明：在新电脑上打开本仓库时，按下面 1 → 4 逐步检查和配置。
> 第 2 步的 MSI 安装需要管理员权限，若无法自动完成，请把命令交给用户手动执行。
> 最后用第 3 步的检查脚本确认全部 PASS。

### 参考环境（作者本机，结果在此环境下得到）

| 组件 | 版本 | 是否必需 | 获取方式 |
|------|------|----------|----------|
| 操作系统 | Windows 10/11 x64 | 必需 | IPhreeqcCOM 是 Windows COM 组件，Linux/macOS 不能运行 PHREEQC 耦合 |
| MATLAB | R2025b（25.2），win64 | 必需 | MathWorks 账号下载安装 |
| Parallel Computing Toolbox | 25.2 | 可选 | 没有时 `parfor` 自动串行，结果相同只是更慢 |
| USGS IPhreeqcCOM | **3.8.6-17100**，x64 | 必需 | 见下方第 2 步 |
| PHREEQC 数据库 | 仓库自带 | 必需 | `ReactiveTransport/RTM/phreeqc/database/backup_legacy_20260626/phreeqc.dat`，无需另装 |

作者本机还装有 Image Processing / PDE / Signal Processing Toolbox（R2025b），
RTM + PHREEQC 主流程不依赖它们。COMSOL NMR、PNG NMR、surrogate NMR 的 Python 环境见 AGENTS.md
的“NMR 运行环境选择”一节，单纯跑 RTM + PHREEQC 不需要。

### 1. 获取代码

```powershell
git clone https://github.com/shizaishiwo323/RTSPHEM-main.git
cd RTSPHEM-main
```

### 2. 安装 IPhreeqcCOM 3.8.6-17100

官方页面：<https://water.usgs.gov/water-resources/software/PHREEQC/index.html>
（PHREEQC 总页面：<https://www.usgs.gov/software/phreeqc-version-3>）

- 64 位安装包（MATLAB win64 必须用这个）：
  <https://water.usgs.gov/water-resources/software/PHREEQC/IPhreeqcCOM-3.8.6-17100-x64.msi>
- 32 位安装包（USGS 建议 64 位 Windows 上两个都装，可选）：
  <https://water.usgs.gov/water-resources/software/PHREEQC/IPhreeqcCOM-3.8.6-17100-win32.msi>

在**管理员** PowerShell 中下载并安装：

```powershell
$msi = "$env:TEMP\IPhreeqcCOM-3.8.6-17100-x64.msi"
Invoke-WebRequest "https://water.usgs.gov/water-resources/software/PHREEQC/IPhreeqcCOM-3.8.6-17100-x64.msi" -OutFile $msi
Start-Process msiexec.exe -ArgumentList "/i `"$msi`" /qn" -Wait -Verb RunAs
```

也可以直接双击 MSI，一路默认安装。安装完成后应满足：

- 目录存在：`C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\bin\IPhreeqcCOM.dll`
- COM 已注册：

  ```powershell
  reg query "HKCR\IPhreeqcCOM.Object\CurVer"   # 应显示 IPhreeqcCOM.Object.4
  ```

要求版本必须是 3.8.6-17100，才能和作者本机的结果完全一致。如果 USGS 已发布更新版本并下架旧链接，
可以先装新版，但要注意 PHREEQC 求解器的改动可能让结果产生微小差异。

### 3. 运行环境检查

在 MATLAB R2025b 中：

```matlab
cd ReactiveTransport/RTM
report = check_rtsphem_environment();
```

或在 PowerShell 中：

```powershell
matlab -batch "cd('ReactiveTransport/RTM'); check_rtsphem_environment();"
```

全部为 `[PASS]` 即配置完成。脚本会检查 MATLAB 版本、`actxserver('IPhreeqcCOM.Object')`，
用仓库自带数据库做一次方解石平衡试算（参考结果约 pH = 8.007，Ca ≈ 9.92e-4 mol/kgw），
并检查 RTM 代码路径。脚本不修改任何文件。

常见问题：

| 现象 | 原因 / 处理 |
|------|-------------|
| `actxserver` 报 `Invalid class string` | IPhreeqcCOM 未安装或未注册；重装 x64 MSI |
| COM 能创建但 `LoadDatabase` 失败 | 路径含特殊字符或数据库被改动；确认 `git status` 干净 |
| 只装了 win32 MSI | win64 MATLAB 无法加载 32 位 COM，必须装 x64 |

### 4. 运行示例：单圆粒 PHREEQC 溶解

```matlab
cd ReactiveTransport/RTM
r = run_single_circle_phreeqc();
```

- 几何：0.1 × 0.05 cm 通道，中心一个 R = 100 µm 方解石圆粒；网格 100 × 50。
- 化学：0.1 mol/L HCl，u = 0.12 cm/s，D = 1e-5 cm²/s，`phreeqc_database_calcite`。
- 输出：`outputs/rtm_runs/single_circle_R100um_0p1x0p05cm_phreeqc_100x50/`，包括
  `interface_images/`（几何切片 PNG）、`dxf_pore/`、`dxf_solid/`、`global_evolution.xlsx`、
  `phreeqc_results/`、`run_metadata.json`。
- 参考结果（作者本机，约 40 分钟）：第 53 步、t = 6.67 s 时颗粒完全溶解，最终孔隙率 1.0，k/k0 ≈ 17.08，
  导出 20 张几何切片。

单圆粒几何没有随机数；在相同 MATLAB 和 IPhreeqcCOM 版本下应得到相同的时间步和切片。
换 CPU/BLAS 时可能有 1e-12 量级的浮点差异。`random` 布局还需要把对应的几何 `.mat`
（`geometryLoadFile`）一起拷贝过去。

## 目录速览

| 目录 | 内容 |
|------|------|
| `ReactiveTransport/RTM/` | RTM 主程序 `PNM_beauty3.m`、PHREEQC 入口 `run_phreeqc_se.m`、示例 `run_single_circle_phreeqc.m` |
| `ReactiveTransport/RTM/couplePhreeqc/` | PHREEQC 耦合（IPhreeqcCOM 调用、输入生成、结果解析） |
| `ReactiveTransport/RTM/phreeqc/database/` | 随仓库分发的 PHREEQC 数据库 |
| `ReactiveTransport/HyPHM/`、`ReactiveTransport/src/` | 有限元流动/输运与渗透率胞问题求解器 |
| `ReactiveTransport/automation/` | COMSOL NMR 批处理与 T2 反演 |
| `outputs/` | 运行结果（已在 `.gitignore` 中忽略，不随仓库分发） |
