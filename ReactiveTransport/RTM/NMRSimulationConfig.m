function cfg = NMRSimulationConfig()
% NMRSimulationConfig - RTM/PNM 中 NMR 模拟路径的集中配置文件。
%
% 这个配置文件统一控制三类 NMR 路径：
%   1. COMSOL NMR:
%        RTM/PNM 导出 pore/solid DXF
%        -> ReactiveTransport/automation/run_comsol_processing.m
%        -> COMSOL 输出 T2 Excel
%        -> run_python_inversion.m 做 T2 反演。
%
%   2. NMR-agent 机器学习替代模型:
%        RTM/PNM 生成 interface PNG
%        -> ReactiveTransport/automation/run_nmr_surrogate_prediction.py
%        -> 预测归一化 T2 衰减曲线
%        -> run_python_inversion.m 做 T2 反演。
%
%   3. PNG 直接数值模拟:
%        RTM/PNM 生成 interface PNG
%        -> png_nmr_driver.py 调用 advanced_tools 中已有的 NMR 求解器
%        -> 生成 T2 衰减曲线
%        -> ReactiveTransport/T2_process 做 NNLS T2 反演。
%
% COMSOL 仍会复用 AutomationConfig.m 的类结构；这里暴露的是同步 RTM
% 运行时最常需要调的 COMSOL/反演参数覆盖项，避免入口脚本里散落配置。
%
% 输入:
%   无输入参数。用户直接修改本文件中的 cfg 字段。
%
% 输出:
%   cfg: MATLAB struct，由 PNM_beauty3.m 读取，并写成
%        PNM_beauty3.m 会读取它，并在需要时写成
%        png_nmr_runtime_config.json 传给 png_nmr_driver.py。
%
% 常用修改方式:
%   1. 选择 NMR 方法:
%        cfg.nmr_method = 'none';          % 不同步跑 NMR
%        cfg.nmr_method = 'comsol';        % COMSOL + T2 反演
%        cfg.nmr_method = 'surrogate';     % NMR-agent 机器学习替代模型
%        cfg.nmr_method = 'png_pixel_cpu'; % PNG pixel CPU
%        cfg.nmr_method = 'png_pixel_gpu'; % PNG pixel GPU
%        cfg.nmr_method = 'png_mesh';      % PNG triangular mesh
%   2. 调 PNG 网格/速度:
%        cfg.max_grid_size、cfg.mesh_bulk_size_um、
%        cfg.mesh_boundary_size_um、cfg.mesh_max_points。
%   3. 调 T2 反演:
%        cfg.nnls_alpha、cfg.t2_min_ms、cfg.t2_max_ms、
%        cfg.t2_num_bins。
%
% 重要约定:
%   - PNG 中红色代表水/孔隙液相，黄色代表固体，白色代表样品外部空白。
%   - 默认使用 PNG 的非白色 bbox 计算像素物理尺寸，不把外缘白边算入
%     lengthXAxis_cm / lengthYAxis_cm。
%   - 每个 RTM 时间步的 PNG NMR 输出会写到当前实验目录下的
%     png_nmr_results/<method>/。
%   - cfg.skip_existing=true 时，已有 decay/T2 结果会被跳过，便于断点续跑。
%   - cfg.fail_on_error=false 时，单个时间步 NMR 失败只记日志，不中断
%     RTM 主流程；批量数据生成时通常建议保持 false。

rtmDir = fileparts(mfilename('fullpath'));
reactiveRoot = fileparts(rtmDir);
projectRoot = fileparts(reactiveRoot);

cfg = struct();

%% ===================== 总开关与方法选择 =====================
% 统一 NMR 方法选择。建议以后只改这一项，不再分别改
% enableNMRSimulation / enableNMRSurrogate / enablePNGSimulation。
%
% 可选值：
%   'none':
%       不在 RTM 过程中同步运行任何 NMR。只输出 RTM、DXF、PNG、Excel 等。
%   'comsol':
%       使用 COMSOL 模型读取 pore/solid DXF，求解 T2 信号，再做 T2 反演。
%       这是物理模型最完整但最慢、环境依赖最重的路径。
%   'surrogate':
%       使用 NMR-agent 训练好的图像模型，从 interface PNG 预测 T2 衰减曲线，
%       再做 T2 反演。速度快，适合批量数据扩展，但依赖训练模型质量。
%   'png_pixel_cpu':
%       直接在 PNG 像素网格上运行 CPU 数值 NMR。
%   'png_pixel_gpu':
%       直接在 PNG 像素网格上运行 GPU 数值 NMR。
%   'png_mesh':
%       从 PNG 生成三角网格后运行数值 NMR。
cfg.nmr_method = 'png_mesh';

% 下面三个 enable* 字段由 cfg.nmr_method 派生，保留是为了兼容
% PNM_beauty3.m 旧接口和 run_metadata.json 记录。通常不需要手动修改。
cfg.enableNMRSimulation = strcmpi(cfg.nmr_method, 'comsol');
cfg.enableNMRSurrogate = any(strcmpi(cfg.nmr_method, {'surrogate', 'ml', 'machine_learning', 'nmr_agent'}));
cfg.enablePNGSimulation = any(strcmpi(cfg.nmr_method, {'png_pixel_cpu', 'png_pixel_gpu', 'png_mesh'}));

% PNG NMR 求解方法。只有 cfg.nmr_method 为 png_* 时这个字段才会被使用。
% 保留 cfg.method 是为了兼容 png_nmr_driver.py 既有参数名。
%   'png_pixel_cpu':
%       直接在 PNG 像素网格上做有限差分/隐式 Euler T2 衰减求解。
%       优点是依赖少、最稳、适合快速验证；缺点是大图会慢。
%   'png_pixel_gpu':
%       与 pixel CPU 相同的离散模型，但用 CuPy/CUDA 解稀疏线性系统。
%       需要当前 Python 环境安装 CuPy 并能访问 CUDA GPU。
%       若 GPU 环境不稳定，建议先用 'png_pixel_cpu' 做基准。
%   'png_mesh':
%       将 PNG 相边界转换为三角网格后做有限元 T2 衰减求解。
%       更接近连续几何边界，适合正式数据生成；依赖 pyGIMLi/Triangle，
%       对 mesh_bulk_size_um、mesh_boundary_size_um 和 mesh_max_points 更敏感。
if cfg.enablePNGSimulation
    cfg.method = cfg.nmr_method;  % png_pixel_cpu | png_pixel_gpu | png_mesh
else
    cfg.method = 'png_mesh';
end
cfg.pngNMRMethod = cfg.method;

% PNG NMR 输出根目录。
%   相对路径: 写到当前 RTM 实验目录下，例如 <exp>/png_nmr_results。
%   绝对路径: 写到指定目录。
% 推荐保持默认相对路径，方便每个样本自包含。
cfg.output_root = 'png_nmr_results';

% 是否跳过已有输出，支持断点续跑。
% true 时，如果某时间步的 decay CSV 或 T2 spectrum 已存在，就尽量复用。
% false 时，会重新运行并覆盖同名派生结果；原始 RTM PNG 不会被改写。
cfg.skip_existing = true;

% NMR 出错时是否中断 RTM 主流程。
%   false: 记录到 png_nmr_sync_log.csv，RTM 继续跑；适合批量数据生成。
%   true : 立即报错停止；适合调试配置或定位某个时间步的问题。
cfg.fail_on_error = false;

%% ===================== COMSOL 与自动化反演参数 =====================
% 这些字段会覆盖 AutomationConfig.m 里的同名/对应字段。
% 如果留空或不修改，就沿用 AutomationConfig.m 的默认值。
%
% COMSOL 安装目录。run_comsol_processing.m 当前主要依赖 MATLAB LiveLink
% 的 mphload/mphsave 等函数；这个路径更多用于配置记录和 validate 检查。
cfg.comsol_path = 'C:\Program Files\COMSOL\COMSOL63\Multiphysics';

% COMSOL mph 模型文件。默认使用 ReactiveTransport/NMR/CT-simulation.mph。
% 若换了 NMR 模型，优先在这里改，RTM 同步 COMSOL 会读取此路径。
cfg.mph_file = fullfile(projectRoot, 'ReactiveTransport', 'NMR', 'CT-simulation.mph');

% COMSOL 路径是否实际执行 COMSOL 求解。通常 cfg.nmr_method='comsol' 时保持 true。
% 调试 T2 反演或复用已有 Excel 时可设为 false。
cfg.enable_comsol = true;

% COMSOL/Surrogate 输出 Excel 后是否继续做 Python T2 反演。
cfg.enable_inversion = true;

% COMSOL/Surrogate 路径中执行 T2_process 反演的 Python。
% 这和 cfg.python_exe 分开：cfg.python_exe 只给 PNG NMR driver 使用；
% cfg.comsol_python_exe 给 run_python_inversion.m / T2_process 使用。
cfg.comsol_python_exe = 'C:\ProgramData\anaconda3\python.exe';

% COMSOL/Surrogate T2 反演是否覆盖已有 *_T2.mat 和 *_T2.png。
% false 更适合断点续跑，true 适合调参后重算派生结果。
cfg.overwrite_existing = false;

% DXF 导入 COMSOL 时的几何缩放系数，沿用 AutomationConfig.m 约定。
cfg.scale_factor = 10000;

% 是否为每个时间步导出求解后的 mph 文件。会明显增加磁盘占用。
cfg.export_mph = false;

% COMSOL/Surrogate 反演脚本路径。默认使用 automation/run_t2_process_inversion.py。
cfg.inversion_script = fullfile(projectRoot, 'ReactiveTransport', 'automation', 'run_t2_process_inversion.py');

%% ===================== NMR-agent 机器学习替代模型参数 =====================
% 训练好的 NMR-agent 模型权重。
cfg.nmrSurrogateModelPath = 'C:\Users\imgw\Documents\Codex\NMR-agent\runs\IMGW_256_300_20260507-130311_3a583275\latest_model.pt';

% NMR-agent 项目根目录。用于推理脚本定位模型代码和数据处理逻辑。
cfg.nmrSurrogateRoot = 'C:\Users\imgw\Documents\Codex\NMR-agent';

% 运行 NMR-agent 推理的 Python。通常应指向 NMR-agent 自己的虚拟环境。
cfg.nmrSurrogatePythonExe = 'C:\Users\imgw\Documents\Codex\NMR-agent\.venv\Scripts\python.exe';

% 可选数据集路径。留空时由推理脚本使用模型或项目默认设置。
cfg.nmrSurrogateDatasetPath = '';

% 输入给替代模型的图像分辨率。
cfg.nmrSurrogateResolution = 256;

% 推理设备：'auto'、'cpu'、'cuda' 等，具体取决于 NMR-agent 脚本支持。
cfg.nmrSurrogateDevice = 'auto';

%% ===================== PNG Python 与模块路径 =====================
% 用于运行 png_nmr_driver.py 的 Python。
% 可写 'python' 使用当前 PATH，也可写绝对路径，例如：
%   C:\Users\imgw\miniconda3\envs\ml\python.exe
%   C:\Python314\python.exe
% 如果使用 png_mesh，需要这个 Python 能导入 pyGIMLi；
% 如果使用 png_pixel_gpu，需要这个 Python 能导入 cupy。
cfg.python_exe = 'python';

% 已有 PNG NMR 工具目录。driver 会把它的上级目录加入 sys.path，
% 然后以模块方式导入 advanced_tools.png_phase_nmr_decay 等函数。
% 不建议复制这些算法代码到 RTSPHEM；保持这里指向成熟工具目录即可。
cfg.advanced_tools_dir = 'C:\Users\imgw\Documents\Codex\NMR模拟\advanced_tools';

% RTSPHEM 内部 T2 反演包路径。
% driver 会导入 nmr_t2.config 和 nmr_t2.nnls，执行固定 alpha 的 NNLS 反演。
cfg.t2_process_path = fullfile(projectRoot, 'ReactiveTransport', 'T2_process');

%% ===================== PNG 几何尺寸解释 =====================
% 是否用非白色 bbox 计算像素物理尺寸。
%   true:
%       lengthXAxis_cm / lengthYAxis_cm 对应红+黄样品区域的外接 bbox。
%       PNG 外缘白边不参与尺寸换算。适合当前 RTM interface_images 输出。
%   false:
%       lengthXAxis_cm / lengthYAxis_cm 对应整张 PNG 的宽高。
%       只有在 PNG 没有白边或白边本身代表物理域时才使用。
cfg.use_nonwhite_bbox_for_geometry_size = true;

% PNG 求解前的最大图像边长限制。
% 若原图最大边长超过该值，会按最近邻下采样，像素物理尺寸同步放大。
%   数值更小: 更快、内存更省，但几何边界更粗。
%   数值更大: 更接近原图，但 pixel/mesh 求解都更慢。
%   []      : 尽量使用原图分辨率；大批量或大图时要谨慎。
cfg.max_grid_size = [];

%% ===================== NMR 物理参数 =====================
% 扩散系数 D，单位 um^2/ms。
% 例如水中 2e-9 m^2/s = 2 um^2/ms。
% 这个参数控制磁化强度在孔隙水域中的扩散速度。
cfg.diffusion_um2_per_ms = 2.0;

% 体弛豫时间 T2B，单位 ms。
% 仅代表体相弛豫，不包括水-固/水-气界面表面弛豫。
% 数值越小，整体衰减越快。
cfg.bulk_t2_ms = 3000.0;

% 水-固界面 surface relaxivity，单位 um/ms。
% 对 T2 峰位置非常敏感；数值越大，靠近固体边界的水衰减越快，
% T2 谱整体更容易向短 T2 移动。
% COMSOL rho=5e-5 m/s 时，换算为 0.05 um/ms。
cfg.rho_solid_um_per_ms = 0.05;

% 水-气界面 surface relaxivity，单位 um/ms。
% 当前 RTSPHEM phase map 中白色外部区域的左右边界按 gas-liquid 处理，
% 上下边界按 solid-liquid 处理；若不希望水-气界面额外弛豫，保持 0。
cfg.rho_gas_um_per_ms = 0.0;

%% ===================== PNG 外边缘边界条件 =====================
% 这里控制 PNG 样品外部白色区域与水相接触时使用哪类表面弛豫。
% 注意：红-黄内部接触永远按 water-solid/solid-liquid 处理；
% 本节只控制红色水相碰到外部白色区域或图像外边缘时的边界类型。
%
% 可选模式：
%   'left_right_gas_top_bottom_solid':
%       默认 RTSPHEM 约定：左右外边缘为 gas-liquid，上下外边缘为 solid-liquid。
%       这对应 png_phase_nmr_decay.py 原始 notebook 的 boundary convention。
%   'all_solid':
%       四个外边缘全部作为 solid-liquid，全部使用 rho_solid_um_per_ms。
%   'all_gas':
%       四个外边缘全部作为 gas-liquid，全部使用 rho_gas_um_per_ms。
%   'custom':
%       使用下面 outside_boundary_left/right/top/bottom 四个字段逐边指定。
%
% 边界类型字段可写：
%   'solid' 或 'solid_liquid'：使用 rho_solid_um_per_ms。
%   'gas'   或 'gas_liquid'  ：使用 rho_gas_um_per_ms。
cfg.outside_boundary_mode = 'left_right_gas_top_bottom_solid';

% 仅当 cfg.outside_boundary_mode='custom' 时通常需要修改这些逐边设置。
% 图像坐标约定：top 是 PNG 上边，bottom 是 PNG 下边；
% left/right 是 PNG 左右边。角点按 left/right 优先，与原脚本逻辑一致。
cfg.outside_boundary_left = 'gas';
cfg.outside_boundary_right = 'gas';
cfg.outside_boundary_top = 'solid';
cfg.outside_boundary_bottom = 'solid';

%% ===================== T2 衰减求解时间参数 =====================
% 衰减模拟时间步长，单位 ms。
% 数值更小: 时间积分更细，结果更平滑，但每个时间步线性求解次数更多。
% 数值更大: 更快，但可能丢失短 T2 的早期衰减细节。
cfg.decay_dt_ms = 5.0;

% 衰减模拟最大时间，单位 ms。
% 应覆盖你关心的最长 T2 区间；若 t2_max_ms 很大但 decay_t_max_ms 太短，
% 长 T2 部分反演约束会偏弱。
cfg.decay_t_max_ms = 5500.0;

%% ===================== PNG pixel GPU 参数 =====================
% GPU 稀疏线性求解器，仅 cfg.method='png_pixel_gpu' 时使用。
%   'cg':
%       共轭梯度法，推荐默认值；隐式 Euler 矩阵通常是对称正定的。
%   'factorized':
%       GPU 稀疏分解，某些 CUDA/cuSPARSE 组合上可能不稳定。
%   'spsolve':
%       直接稀疏求解，适合调试或小问题。
cfg.pixel_gpu_solver = 'cg';  % cg | factorized | spsolve

% GPU CG 收敛容差。越小越接近 CPU 结果，但迭代次数可能增加。
cfg.pixel_gpu_cg_tol = 1e-8;

% GPU CG 最大迭代次数。若出现不收敛，可适当增大；
% 若经常达到上限，优先检查网格尺寸、时间步长和 GPU 环境。
cfg.pixel_gpu_cg_maxiter = 2000;

%% ===================== PNG triangular mesh 参数 =====================
% 三角网格内部目标点距，单位 um，仅 cfg.method='png_mesh' 时使用。
% 数值越小，水域内部网格越密，几何/扩散解更细，但速度更慢。
cfg.mesh_bulk_size_um = 16.0;

% 三角网格边界目标点距，单位 um，仅 cfg.method='png_mesh' 时使用。
% 控制水-固和水-气边界采样密度。通常边界需要比内部更细，
% 因为表面弛豫项直接作用在边界附近。
cfg.mesh_boundary_size_um = 5.0;

% 三角网格节点数安全上限。
% 防止复杂 PNG 或过细边界参数导致节点数爆炸、内存占用过大。
% 如果 mesh 求解被这个上限截断或报错，优先增大 mesh_bulk_size_um /
% mesh_boundary_size_um，必要时再提高本上限。
cfg.mesh_max_points = 50000;

%% ===================== T2 反演参数 =====================
% 是否对 PNG 求解得到的 decay 曲线立即做 T2 反演。
% true 会生成 spectrum CSV、MAT 和 PNG；false 只生成 decay 曲线。
cfg.run_inversion = true;

% 固定 alpha 的 NNLS 正则化强度。
% 数值越大，T2 谱越平滑；数值越小，谱峰可能更尖锐但更容易受噪声/数值误差影响。
% 当前 driver 使用固定 alpha，不做 L-curve 自动选择。
cfg.nnls_alpha = 0.1;

% T2 反演下限，单位 ms。
% 若关心强表面弛豫导致的短 T2 峰，可适当降低。
cfg.t2_min_ms = 1.0;

% T2 反演上限，单位 ms。
% 应与 decay_t_max_ms 协同设置；上限过大而衰减时长不足时，
% 长 T2 端会缺少充分约束。
cfg.t2_max_ms = 100000.0;

% T2 对数 bins 数量。
% 数值越大，谱分辨率越高，但 NNLS 规模增大、结果也可能更敏感。
cfg.t2_num_bins = 240;

% 从 decay 峰值之后用于反演的最少采样点数量。
% 若 decay_t_max_ms / decay_dt_ms 太小导致点数不足，反演会失败。
cfg.min_points_after_trim = 10;
end
