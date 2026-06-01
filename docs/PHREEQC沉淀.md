下面我按“**模型思想 → 控制方程 → 上尺度 → 数值耦合 → 两个算例设置 → 误差/假设**”的逻辑，把这篇 *An improved micro–macro model of multicomponent reactive transport for multimineral systems* 里模拟全过程涉及的方程、边界条件和关键假设梳理出来。核心上，这篇文章做的是：**宏观尺度求达西流和溶质运移，微观 REV 单元里追踪矿物溶解/沉淀导致的孔隙结构演化，再把孔隙率、渗透率、有效扩散张量、反应比表面积反馈给宏观方程**。图 2 给出的流程就是“宏观输运/反应 ↔ 微观几何演化 ↔ 有效参数更新”的循环。

---

# 1. 几何与尺度设定

## 1.1 宏观区域与周期性微观单元

论文假设宏观区域为二维矩形区域：

[
\Omega \subset \mathbb{R}^2
]

其边界为：

[
\partial \Omega
]

每一个宏观网格点都附带一个微观代表性体积单元，即 REV。论文假设孔隙介质具有周期性微结构，周期长度为：

[
\varepsilon \ll 1
]

其中 (\varepsilon) 表示微观尺度与宏观尺度的比值。

整体区域被分解为许多尺度化单元：

[
\Omega = \bigcup_{i,j} Y_{\varepsilon}^{i,j}
]

流体区域为：

[
\Omega_\varepsilon = \bigcup_{i,j} (\chi_\varepsilon^1)^{i,j}
]

固-液界面为：

[
\Gamma_\varepsilon = \bigcup_{i,j} \Gamma_\varepsilon^{i,j}
]

这里 (\chi_\varepsilon^1) 表示流体相，(\chi_\varepsilon^2,\chi_\varepsilon^3,\dots) 表示不同矿物固相。化学反应发生在固-液界面上，所以 (\Gamma_\varepsilon) 是控制孔隙结构演化的关键。

---

# 2. 孔隙尺度模型：流体、溶质、界面演化

这一部分是整个模型的“物理源头”。

## 2.1 Level-set 描述孔隙/固体/界面

论文用 level-set 函数 (L_\varepsilon(t,x)) 描述孔隙结构：

[
L_\varepsilon(t,x)=
\begin{cases}
<0, & \text{fluid phase / pore space} \
=0, & \text{interface} \

> 0, & \text{solid phase}
> \end{cases}
> ]

也就是：

| (L_\varepsilon) 的符号 | 含义     |
| ------------------- | ------ |
| (L_\varepsilon<0)   | 孔隙流体区域 |
| (L_\varepsilon=0)   | 固-液界面  |
| (L_\varepsilon>0)   | 固体矿物区域 |

界面演化方程为：

[
\frac{\partial L_\varepsilon}{\partial t}
+
v_{n,\varepsilon} |\nabla L_\varepsilon|
=0
]

其中 (v_{n,\varepsilon}) 是固-液界面的法向移动速度。

矿物溶解/沉淀导致界面移动，界面速度由表面反应速率控制：

[
v_{n,\varepsilon}=V_m r_\varepsilon
]

其中：

| 符号              | 含义     |
| --------------- | ------ |
| (V_m)           | 矿物摩尔体积 |
| (r_\varepsilon) | 界面反应速率 |

溶解时固体减少，孔隙扩大；沉淀时固体增加，孔隙缩小。

---

## 2.2 孔隙尺度溶质运移方程

对组分浓度 (c_\varepsilon)，孔隙尺度控制方程为：

[
\frac{\partial c_\varepsilon}{\partial t}
+
\nabla \cdot (u_\varepsilon c_\varepsilon)
------------------------------------------

\nabla \cdot (D\nabla c_\varepsilon)
=0,
\quad x\in \Omega_\varepsilon(t)
]

这是一个典型的 **对流-扩散方程**。

固-液界面上的通量条件为：

[
(-u_\varepsilon c_\varepsilon + D\nabla c_\varepsilon)\cdot n_\varepsilon
-------------------------------------------------------------------------

\varepsilon \rho_s^{-1} r_\varepsilon(c_\varepsilon-\rho_s)
=0,
\quad x\in \Gamma_\varepsilon(t)
]

初始条件为：

[
c_\varepsilon=c^0,
\quad x\in \Omega_\varepsilon(0)
]

其中：

| 符号              | 含义         |
| --------------- | ---------- |
| (D)             | 分子扩散系数     |
| (u_\varepsilon) | 孔隙尺度流速     |
| (n_\varepsilon) | 固-液界面单位法向量 |
| (\rho_s)        | 固体密度       |

这个界面条件本质上是：**固-液界面不能储存质量，所以流体侧的通量变化必须与矿物反应消耗/释放的质量一致**。

---

## 2.3 孔隙尺度流动方程：Stokes 问题

孔隙尺度流速 (u_\varepsilon) 和压力 (p_\varepsilon) 由 Stokes 型方程得到：

[
\varepsilon^2 \Delta u_\varepsilon
==================================

\frac{1}{\mu}\nabla p_\varepsilon,
\quad x\in \Omega_\varepsilon(t)
]

[
\nabla \cdot u_\varepsilon =0,
\quad x\in \Omega_\varepsilon(t)
]

固-液界面上的边界条件为：

[
u_\varepsilon=v_{n,\varepsilon} n_\varepsilon,
\quad x\in \Gamma_\varepsilon(t)
]

其中 (\mu) 是动力黏度。论文指出，这对应于随时间演化边界上的无滑移条件。

---

# 3. 上尺度：从孔隙尺度到宏观尺度

论文用均质化/双尺度渐近展开把孔隙尺度模型推导成微-宏耦合模型。

## 3.1 双尺度变量

引入宏观变量 (x) 和微观变量 (y)：

[
y=\frac{x}{\varepsilon}
]

因此空间微分算子展开为：

[
\nabla = \nabla_x+\frac{1}{\varepsilon}\nabla_y
]

[
\nabla\cdot = \nabla_x\cdot+\frac{1}{\varepsilon}\nabla_y\cdot
]

[
\Delta
======

\Delta_x
+
\frac{2}{\varepsilon}\nabla_x\cdot\nabla_y
+
\frac{1}{\varepsilon^2}\Delta_y
]

---

## 3.2 渐近展开假设

任意变量 (\phi_\varepsilon) 被展开为：

[
\phi_\varepsilon(t,x)
=====================

\phi_0(t,x,y)
+
\varepsilon \phi_1(t,x,y)
+
\varepsilon^2 \phi_2(t,x,y)
+\cdots
]

Level-set 函数也展开为：

[
L_\varepsilon(t,x)
==================

L_0(t,x,y)
+
\varepsilon L_1(t,x,y)
+
\varepsilon^2 L_2(t,x,y)
+\cdots
]

法向量展开为：

[
n_\varepsilon=n_0+\varepsilon n_1+O(\varepsilon^2)
]

其中：

[
n_0=
\frac{\nabla_y L_0}{|\nabla_y L_0|}
]

这部分假设很关键：**宏观尺度和微观尺度可以分离，并且微观结构的变化可以通过渐近展开传递到宏观有效参数中**。

---

# 4. 微-宏耦合模型

经过上尺度后，模型在宏观上求解平均浓度、平均流速和压力；在微观 REV 中更新几何结构和有效参数。

## 4.1 微观 level-set 方程

微观流体区域定义为：

[
Y_{l,0}(t,x)={y\in Y:L_0(t,x,y)<0}
]

固-液界面为：

[
\Gamma_0(t,x)={y\in Y:L_0(t,x,y)=0}
]

微观界面演化为：

[
\frac{\partial L_0}{\partial t}
+
v_{n,0}|\nabla L_0|
=0,
\quad \text{in } \Omega\times Y
]

---

## 4.2 宏观反应-运移方程

宏观组分浓度 (c_0) 满足：

[
\frac{\partial \phi c_0}{\partial t}
+
\nabla_x\cdot(u_0 c_0)
----------------------

# \nabla_x\cdot(D\nabla c_0)

-s r_0,
\quad \text{in } \Omega
]

其中：

| 符号     | 含义       |
| ------ | -------- |
| (\phi) | 孔隙率      |
| (u_0)  | 宏观平均流速   |
| (D)    | 有效扩散张量   |
| (s)    | 反应比表面积   |
| (r_0)  | 宏观等效反应速率 |

孔隙率和反应比表面积由 REV 几何计算：

[
\phi=\frac{|Y_{l,0}|}{|Y|}
]

[
s=\frac{|\Gamma_0|}{|Y|}
]

这就是这个模型的核心：**矿物反应改变 REV 几何，REV 几何改变 (\phi,s,K,D)，这些参数再反馈到宏观反应-运移方程。**

---

## 4.3 宏观流动方程：Darcy 定律

宏观平均流速由 Darcy 定律给出：

[
u_0
===

-\frac{K}{\mu}\nabla p_0,
\quad \text{in } \Omega
]

连续性方程为：

[
\nabla\cdot u_0=0,
\quad \text{in } \Omega
]

其中 (K) 是有效渗透率张量。

---

## 4.4 有效扩散张量

有效扩散张量定义为：

[
D_{ij}(t,x)
===========

\frac{1}{|Y|}
\int_{Y_{l,0}(t,x)}
D
\left(
\frac{\partial \zeta_j}{\partial y_i}
+
\delta_{ij}
\right)
dy
]

其中 (\zeta_j) 是浓度扰动场，满足辅助胞元问题：

[
-\nabla_y\cdot(\nabla_y \zeta_j)=0,
\quad \text{in } Y_{l,0}(t,x)
]

[
\nabla_y\zeta_j\cdot n_0=-e_j\cdot n_0,
\quad \text{on } \Gamma_0(t,x)
]

[
\zeta_j \text{ is periodic in } y,
\quad
\frac{1}{|Y|}\int_{Y_{l,0}}\zeta_j dy=0
]

---

## 4.5 有效渗透率张量

有效渗透率张量定义为：

[
K_{ij}(t,x)
===========

\frac{1}{|Y|}
\int_{Y_{l,0}(t,x)}
\omega_j^i dy
]

其中 (\omega_j) 是速度扰动场，(\pi_j) 是压力扰动场，二者满足 Stokes 型辅助胞元问题：

[
-\Delta_y \omega_j+\nabla_y \pi_j=e_j,
\quad \text{in } Y_{l,0}(t,x)
]

[
\nabla_y\cdot \omega_j=0,
\quad \text{in } Y_{l,0}(t,x)
]

[
\omega_j=0,
\quad \text{on } \Gamma_0(t,x)
]

[
\omega_j,\pi_j \text{ periodic in } y,
\quad
\frac{1}{|Y|}\int_{Y_{l,0}}\pi_jdy=0
]

这些胞元问题每次几何更新后都要重新求解，用于更新 (K) 和 (D)。

---

# 5. 数值耦合方法

论文的 RTSPHEM-P 模型用的是 **operator splitting / SNIA** 思路：

1. 宏观尺度求 Darcy 流速；
2. 宏观尺度求溶质运移；
3. 把各网格点浓度传给 PHREEQC；
4. PHREEQC 计算水化学平衡、矿物溶解/沉淀反应；
5. 根据反应速率更新 REV 中固-液界面；
6. 重新计算孔隙率 (\phi)、比表面积 (s)、渗透率 (K)、有效扩散张量 (D)；
7. 进入下一个时间步。

原始 RTSPHEM 直接用主组分总浓度参与反应，并采用 GIA；本文改进后的 RTSPHEM-P 把 IPhreeqc 接入 RTSPHEM，用 SNIA 处理反应-运移，并屏蔽了原 RTSPHEM 内部简化的地球化学反应模块。

---

# 6. 算例一：SPE10 派生的两矿物体系

这个算例用于比较原始 RTSPHEM 和改进 RTSPHEM-P。

## 6.1 几何与网格设置

宏观区域：

[
\Omega=22\times 6 \ \mathrm{mm}^2
]

SPE10 原始数据分辨率为 (220\times 60)，论文将每 (5\times 5) 个数据块粗化为一个块，得到：

[
44\times 12 = 528
]

个 REV 单元。

宏观网格使用：

[
11162
]

个三角形单元。

每个 REV 代表：

[
0.25\ \mathrm{mm}^2
]

的宏观区域。

微观几何假设为矩形孔隙结构，用矩形长宽两个自由度调节初始渗透率分布。

固体矿物组成假设为：

[
50% \ \text{calcite}
+
50% \ \text{dolomite}
]



---

## 6.2 原 RTSPHEM 中的矿物反应路径

方解石：

[
\mathrm{Calcite:}
\quad
\mathrm{CaCO_3}
\rightleftharpoons
\mathrm{Ca^{2+}}
+
\mathrm{CO_3^{2-}}
]

白云石：

[
\mathrm{Dolomite:}
\quad
\mathrm{CaMg(CO_3)_2}
\rightleftharpoons
\mathrm{Ca^{2+}}
+
\mathrm{Mg^{2+}}
+
2\mathrm{CO_3^{2-}}
]

论文指出，这种写法没有显式考虑 (\mathrm{H^+}) 消耗，会导致 pH 和矿物溶解速率模拟失真。

---

## 6.3 RTSPHEM-P 中补充的水相平衡反应

改进模型中进一步考虑：

[
\mathrm{CO_3^{2-}+H^+}
\rightleftharpoons
\mathrm{HCO_3^-}
]

[
\mathrm{CO_3^{2-}+2H^+}
\rightleftharpoons
\mathrm{CO_2+H_2O}
]

[
\mathrm{H_2O}
\rightleftharpoons
\mathrm{OH^-+H^+}
]

[
\mathrm{Ca^{2+}+CO_3^{2-}+H^+}
\rightleftharpoons
\mathrm{CaHCO_3^+}
]

[
\mathrm{Mg^{2+}+CO_3^{2-}+H^+}
\rightleftharpoons
\mathrm{MgHCO_3^+}
]

[
\mathrm{Ca^{2+}+H_2O}
\rightleftharpoons
\mathrm{Ca(OH)^+ + H^+}
]

[
\mathrm{Mg^{2+}+H_2O}
\rightleftharpoons
\mathrm{Mg(OH)^+ + H^+}
]

[
\mathrm{Mg^{2+}+CO_3^{2-}}
\rightleftharpoons
\mathrm{MgCO_3}
]

这些反应由 PHREEQC 数据库处理，因此 RTSPHEM-P 能够考虑更真实的水化学形态分配和活度效应。

---

## 6.4 两矿物体系的动力学方程

方解石反应速率：

[
r_{\mathrm{Calcite}}
====================

(k_1 c_{\mathrm{H^+}}+k_2)
\left(
1-
\frac{
c_{\mathrm{Ca^{2+}}}c_{\mathrm{CO_3^{2-}}}
}{
K_{\mathrm{Calcite}}
}
\right)
]

白云石反应速率：

[
r_{\mathrm{Dolomite}}
=====================

k c_{\mathrm{H^+}}^{0.5}
\left(
1-
\frac{
c_{\mathrm{Ca^{2+}}}
c_{\mathrm{Mg^{2+}}}
c_{\mathrm{CO_3^{2-}}}^2
}{
K_{\mathrm{Dolomite}}
}
\right)
]

参数为：

[
k_1=0.89\ \mathrm{mol/m^2/s}
]

[
k_2=6.7\times 10^{-7}\ \mathrm{mol/m^2/s}
]

[
k=4.5\times 10^{-4}\ \mathrm{mol/m^2/s}
]

[
K_{\mathrm{Calcite}}=10^{-8.234}
]

[
K_{\mathrm{Dolomite}}=10^{-16.5}
]

---

## 6.5 两矿物算例的溶质边界条件

外边界分为：

[
\Gamma_{\mathrm{left}},
\Gamma_{\mathrm{right}},
\Gamma_{\mathrm{up}},
\Gamma_{\mathrm{down}}
]

### 上下边界：无通量边界

[
(c^*u-D\nabla c^*)\cdot n=0,
\quad \text{on } \Gamma_{\mathrm{up}},\Gamma_{\mathrm{down}}
]

其中：

原 RTSPHEM：

[

* \in {\mathrm{Ca^{2+}},\mathrm{Mg^{2+}},\mathrm{CO_3^{2-}},\mathrm{H^+}}
  ]

RTSPHEM-P：

[

* \in {\mathrm{Ca},\mathrm{Mg},\mathrm{C},\mathrm{H^+}}
  ]

### 左边界：除氢离子外，其他组分无扩散-对流净通量

[
(c^#u-D\nabla c^#)\cdot n=0,
\quad \text{on } \Gamma_{\mathrm{left}}
]

原 RTSPHEM：

[
# \in {\mathrm{Ca^{2+}},\mathrm{Mg^{2+}},\mathrm{CO_3^{2-}}}
]

RTSPHEM-P：

[
# \in {\mathrm{Ca},\mathrm{Mg},\mathrm{C}}
]

### 左边界：氢离子输入

[
(c_{\mathrm{H^+}}u-D\nabla c_{\mathrm{H^+}})\cdot n
===================================================

c_{\mathrm{H^+}}u\cdot n,
\quad \text{on } \Gamma_{\mathrm{left}}
]

### 右边界：氢离子零扩散通量

[
D\nabla c_{\mathrm{H^+}}\cdot n=0,
\quad \text{on } \Gamma_{\mathrm{right}}
]

### 初始条件

其他组分初始浓度：

[
c^#=0,
\quad \text{in } \Omega
]

氢离子初始浓度：

[
c_{\mathrm{H^+}}=10^{-2}\ \mathrm{mol/m^3},
\quad \text{in } \Omega
]

这个浓度对应初始酸性流体环境，论文后面也把它解释为约 pH = 5 的背景条件。

---

## 6.6 两矿物算例的 Darcy 边界条件

左边界给定流入通量：

[
u\cdot n
========

-5\times 10^{-5}\ \mathrm{m/s},
\quad \text{on } \Gamma_{\mathrm{left}}
]

上下边界无流：

[
u\cdot n=0,
\quad \text{on } \Gamma_{\mathrm{up}},\Gamma_{\mathrm{down}}
]

右边界定压：

[
p=0,
\quad \text{on } \Gamma_{\mathrm{right}}
]

---

## 6.7 时间步长设置

RTSPHEM-P 使用 SNIA，所以需要满足 Courant 条件：

[
Cr=\frac{\Delta t \cdot u}{\Delta x}\leq 1
]

因此最大时间步长为：

[
\Delta t_{\max}=4\ \mathrm{s}
]

总模拟时间：

[
T=4.32\times 10^3\ \mathrm{s}
]

原 RTSPHEM 使用 GIA，最大时间步长设为：

[
\Delta t_{\max}=10000\ \mathrm{s}
]

总模拟时间为：

[
T=4.32\times 10^7\ \mathrm{s}
]

为了比较优先流通道形成，论文还把 RTSPHEM-P 的反应速率放大 1000 倍，使 (4.32\times 10^3\ \mathrm{s}) 等效为 (4.32\times 10^6\ \mathrm{s})。不过作者也承认，这种线性缩放如果不同时调整流速等参数，会违反质量守恒，因此主要用于定性展示模型捕捉优先流通道的能力。

---

# 7. 算例二：真实 SEM 图像多矿物体系

第二个算例是模型的扩展应用，用真实二维 SEM 图像构造多矿物系统。

## 7.1 SEM 图像与几何设置

原始 SEM 图像分辨率为：

[
9474\times 6947
]

像素尺寸为：

[
1.2\ \mu m
]

如果直接按原始分辨率模拟，网格数会达到：

[
65,815,878
]

计算成本太高，所以论文利用微-宏模型进行上尺度处理。

宏观区域为：

[
\Omega=11.28\times 8.16\ \mathrm{mm^2}
]

粗化后形成：

[
31\times 23=713
]

个典型 REV 单元。

每个典型单元代表：

[
0.0576\ \mathrm{mm^2}
]

的宏观区域。

宏观网格使用：

[
13778
]

个三角形单元。

微观几何演化网格：

[
100\times 100
]

有效参数计算网格：

[
50\times 50
]

固体区域渗透率设为：

[
10^{-16}\ \mathrm{m^2}
]

孔隙区域渗透率设为：

[
10^{-12}\ \mathrm{m^2}
]

矿物体积分数近似为：

| 矿物         | 固体骨架体积分数 |
| ---------- | -------: |
| Quartz     |    37.5% |
| Dolomite   |      25% |
| K-feldspar |      15% |
| Muscovite  |   11.25% |
| Kaolinite  |   11.25% |

Ilmenite 因含量很低被忽略。

---

## 7.2 多矿物反应路径

石英：

[
\mathrm{Quartz:}
\quad
\mathrm{SiO_2+2H_2O}
\rightleftharpoons
\mathrm{H_4SiO_4}
]

钾长石：

[
\mathrm{K\ feldspar:}
\quad
\mathrm{KAlSi_3O_8+8H_2O}
\rightleftharpoons
\mathrm{K^+ + Al(OH)_4^- + 3H_4SiO_4}
]

白云母：

[
\mathrm{Muscovite:}
\quad
\mathrm{KAl_3Si_3O_{10}(OH)_2+10H^+}
\rightleftharpoons
\mathrm{K^+ + 3Al^{3+}+3H_4SiO_4}
]

高岭石：

[
\mathrm{Kaolinite:}
\quad
\mathrm{Al_2Si_2O_5(OH)_4+6H^+}
\rightleftharpoons
\mathrm{H_2O+2Al^{3+}+2H_4SiO_4}
]

白云石的反应路径沿用前面的碳酸盐矿物溶解反应。

---

## 7.3 通用矿物溶解/沉淀动力学方程

论文采用一般形式：

[
r_m
===

\pm k_m
\left|
1-
\left(
\frac{Q_m}{K_m}
\right)^p
\right|^q
]

其中：

| 符号      | 含义        |
| ------- | --------- |
| (r_m>0) | 溶解        |
| (r_m<0) | 沉淀        |
| (k_m)   | 动力学反应速率常数 |
| (Q_m)   | 离子活度积     |
| (K_m)   | 平衡常数      |
| (p,q)   | 经验参数      |

矿物反应速率常数考虑中性、酸性和碱性三种机制：

[
k_m
===

k_{25}^{nu}
\exp\left[
-\frac{E_a^{nu}}{R}
\left(
\frac{1}{T}
-----------

\frac{1}{298.15}
\right)
\right]
]

[
+
k_{25}^{H}
\exp\left[
-\frac{E_a^{H}}{R}
\left(
\frac{1}{T}
-----------

\frac{1}{298.15}
\right)
\right]
a_H^{n_H}
]

[
+
k_{25}^{OH}
\exp\left[
-\frac{E_a^{OH}}{R}
\left(
\frac{1}{T}
-----------

\frac{1}{298.15}
\right)
\right]
a_{OH}^{n_{OH}}
]

由于该模拟假设恒温：

[
T=25^\circ C=298.15\ K
]

所以简化为：

[
k_m
===

k_{25}^{nu}
+
k_{25}^{H}a_H^{n_H}
+
k_{25}^{OH}a_{OH}^{n_{OH}}
]

论文表 1 给出的参数为：

| 矿物         |          (k_{25}^{nu}) |           (k_{25}^{H}) | (n_H) |          (k_{25}^{OH}) | (n_{OH}) |
| ---------- | ---------------------: | ---------------------: | ----: | ---------------------: | -------: |
| Quartz     | (1.023\times 10^{-14}) |                      — |     — |                      — |        — |
| Dolomite   |  (2.951\times 10^{-8}) |  (6.457\times 10^{-4}) |   0.5 |                      — |        — |
| K-feldspar | (3.890\times 10^{-13}) | (8.710\times 10^{-11}) |   0.5 | (6.310\times 10^{-22}) |   -0.823 |
| Kaolinite  | (6.918\times 10^{-14}) | (4.898\times 10^{-12}) | 0.777 | (8.913\times 10^{-18}) |   -0.472 |
| Muscovite  | (3.020\times 10^{-13}) | (7.762\times 10^{-12}) |   0.5 |                      — |        — |

单位主要为 (\mathrm{mol/m^2/s})。

---

## 7.4 SEM 多矿物算例的溶质边界条件

上下边界无通量：

[
(c^*u-D\nabla c^*)\cdot n=0,
\quad
\text{on } \Gamma_{\mathrm{up}},\Gamma_{\mathrm{down}}
]

其中：

[

* \in {\mathrm{K,Ca,Mg,Al,C,Si,H^+}}
  ]

左边界除 (\mathrm{H^+}) 外其他组分无通量：

[
(c^#u-D\nabla c^#)\cdot n=0,
\quad
\text{on } \Gamma_{\mathrm{left}}
]

其中：

[
#\in {\mathrm{K,Ca,Mg,Al,C,Si}}
]

左边界输入氢离子：

[
(c_{\mathrm{H^+}}u-D\nabla c_{\mathrm{H^+}})\cdot n
===================================================

c_{\mathrm{H^+}}u\cdot n,
\quad
\text{on } \Gamma_{\mathrm{left}}
]

右边界氢离子零扩散通量：

[
D\nabla c_{\mathrm{H^+}}\cdot n=0,
\quad
\text{on } \Gamma_{\mathrm{right}}
]

初始条件：

[
c^#=0,
\quad \text{in } \Omega
]

[
c_{\mathrm{H^+}}=10^{-2}\ \mathrm{mol/m^3},
\quad \text{in } \Omega
]

---

## 7.5 SEM 多矿物算例的 Darcy 边界条件

左边界流入：

[
u\cdot n
========

-1.5\times 10^{-7}\ \mathrm{m/s},
\quad
\text{on } \Gamma_{\mathrm{left}}
]

上下边界无流：

[
u\cdot n=0,
\quad
\text{on } \Gamma_{\mathrm{up}},\Gamma_{\mathrm{down}}
]

右边界定压：

[
p=0,
\quad
\text{on } \Gamma_{\mathrm{right}}
]

时间步长条件：

[
Cr\leq 1
]

最大时间步长：

[
\Delta t_{\max}=1000\ \mathrm{s}
]

总模拟时间：

[
T=1.296\times 10^6\ \mathrm{s}
]



---

# 8. 几何演化误差与质量守恒检验

因为 VIIM 界面追踪方法只有一阶精度，论文专门用质量守恒思想来估计几何演化误差。

## 8.1 白云石误差

论文认为在 SEM 算例中，硅酸盐矿物溶解速率远小于白云石，孔隙体积变化主要由白云石溶解决定。因此孔隙率变化应接近白云石体积分数变化。

白云石误差定义为：

[
Error_{\mathrm{Dolomite}}
=========================

\frac{
|\Delta \phi-\Delta V_{\mathrm{Dolomite}}|
}{
\Delta \phi
}
\cdot
P_{\mathrm{Dolomite}}
\times 100%
]

其中：

| 符号                             | 含义             |
| ------------------------------ | -------------- |
| (\Delta \phi)                  | 模拟前后孔隙率变化      |
| (\Delta V_{\mathrm{Dolomite}}) | 白云石体积分数变化      |
| (P_{\mathrm{Dolomite}})        | 白云石在固体骨架中的体积分数 |

## 8.2 其他硅酸盐矿物误差

对于 K-feldspar、Kaolinite、Muscovite、Quartz：

[
Error_*
=======

\frac{
|\Delta M_*|
}{
M_0^*
}
\times 100%
]

其中：

[

* \in
  {
  \mathrm{K\ feldspar, Kaolinite, Muscovite, Quartz}
  }
  ]

最终论文得到：

| 矿物         | 平均相对误差 |
| ---------- | -----: |
| Dolomite   |  1.80% |
| K-feldspar |  3.03% |
| Kaolinite  |  4.96% |
| Muscovite  |  4.96% |
| Quartz     |  0.99% |

加权平均总误差为：

[
2.39%
]

作者认为这可以近似视为满足质量守恒，说明 VIIM 对复杂多矿物典型单元几何演化是可用的。 

---

# 9. 关键假设汇总

这篇论文的模拟建立在下面这些核心假设上：

| 类别    | 假设                                                |
| ----- | ------------------------------------------------- |
| 几何假设  | 孔隙介质为二维矩形区域；每个宏观网格点附带一个 REV                       |
| 尺度假设  | 微观尺度与宏观尺度可分离，(\varepsilon\ll 1)                   |
| 周期性假设 | 微观固体结构具有周期性                                       |
| 饱和假设  | 模型讨论的是饱和孔隙介质中的反应运移                                |
| 界面假设  | 化学反应发生在固-液界面，界面本身不储存质量                            |
| 流动假设  | 孔隙尺度为 Stokes 流，宏观尺度为 Darcy 流                      |
| 边界假设  | 固-液界面满足无滑移/界面随反应移动条件                              |
| 化学假设  | RTSPHEM-P 用 PHREEQC 处理水相平衡、活度、矿物动力学               |
| 温度假设  | SEM 多矿物算例为恒温 (25^\circ C)                         |
| 数值假设  | 微-宏数据交换使用 operator splitting / SNIA               |
| 误差假设  | VIIM 为一阶精度，需要用质量守恒估计几何演化误差                        |
| 缩放假设  | 两矿物长时间模拟中，RTSPHEM-P 通过放大反应速率实现等效时间缩放，但该缩放主要用于定性展示 |
| 沉淀限制  | 模型可处理一定矿物沉淀，但多矿物同时沉淀会出现界面沉淀位置不确定问题，因此仍有限制         |

论文结论中也明确说，RTSPHEM-P 的优势在于：相比原 RTSPHEM，它能通过 PHREEQC 更真实地处理复杂水化学反应；同时能捕捉多矿物体系中因拓扑变化造成的渗透率、扩散张量和反应比表面积突变。不过作者也承认，VIIM 的一阶精度和多矿物同时沉淀问题仍是后续需要改进的地方。 

---

# 10. 用一句话概括模拟全过程

这篇文章的模拟逻辑可以压缩成：

[
\text{Darcy流}
\rightarrow
\text{溶质对流-扩散}
\rightarrow
\text{PHREEQC地球化学反应}
\rightarrow
\text{矿物溶解/沉淀速率}
\rightarrow
\text{固-液界面移动}
\rightarrow
\text{REV几何更新}
\rightarrow
(\phi,s,K,D)\text{更新}
\rightarrow
\text{进入下一时间步}
]

也就是说，它不是单纯做反应运移，而是把 **化学反应—孔隙结构演化—水力参数变化—宏观运移反馈** 连成了一个闭环。
