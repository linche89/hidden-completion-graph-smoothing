# 有限时线性网络计算与 rounds–bits–energy–compute 权衡：定理级尽调

> 截止日期：2026-09-13。范围：有限轮线性变换、网络线性方程、图滤波、通信下界、有限比特与率失真、in-network computation、无线传感器能耗、localized Kalman/SLS，以及 2026 年最接近的碰撞工作。本文把论文原结论与本文自行推出的基本引理严格分开；“本地 PDF p.x”均指 PDF 文件页序，而非期刊印刷页码。

## 0. 结论先行

在没有固定协议模型之前，不存在一个有意义的“最 relaxed 充要条件”。一旦模型固定，答案分为三层：

| 问题层次 | 截至 2026 年的状态 | 精确判断 |
|---|---|---|
| 同步、可靠、消息为任意维无限精度实数，节点 (i) 只需输出给定线性映射的第 (i) 个块 | **SOLVED** | (r) 轮可算当且仅当目标输出不依赖于 (r)-跳因果锥之外的数据；最少轮数就是所有非零依赖块的最大有向距离。环完全允许；无环不是必要条件。 |
| 每节点只有一个标量/固定维状态，每轮仅乘一个图稀疏矩阵 | **形式上 SOLVED，结构上 PARTIAL** | (r) 轮可算当且仅当目标矩阵有 (r) 个图稀疏因子的乘积分解。一般图上可检验、图论化的闭式刻画并未被找到；Hamilton 圈等只给 generic 充分结果。 |
| 有限比特且要求对任意实数输入零误差精确输出 | **不可能，除非跨割目标为常数** | 有限 transcript 只能区分有限个情形，不能精确表示连续值域。必须改成有限字母表、有限精度、概率误差或 MSE/失真。 |
| 任意图、任意统计源、任意目标、有限比特、轮数、物理能耗、计算与内存、随机/恶意失效的完整 Pareto 充要区域 | **OPEN-LOOKING，而且当前表述过宽** | 它严格包含仍开放的多终端 remote rate–distortion，也包含一般网络函数计算/网络编码与无线调度问题；因此不可能由一个仅依赖“拓扑 + 秩”的简单条件终结。 |
| localized Kalman/SLS | **PARTIAL** | 对固定线性估计器/控制器类，系统响应的 affine feasibility 是充要的；但 locality/FIR 可行性本身仍是实例依赖的优化问题，而且论文通常不同时建模 bits、joules、失效与计算。 |
| Zhao et al. 2026 的 design trilemma | **最接近但未闭环** | Theorem 1 的 NMI 只是充分条件；Remark 1 明说 sufficiency，结尾把 necessary conditions 列为未来工作。它不是资源 Pareto converse。 |

因此，2011 年论文中的“无环 + 每个 (A_i) 满列秩”确实过强。对静态线性估计，唯一全局状态通常只需堆叠后的全局信息矩阵满秩；对某个局部输出则只需该线性泛函可识别。对通信拓扑，最弱精确条件是**目标依赖对应的可达性/因果锥条件**，而非无环。可是，一旦同时要求有限比特、容错、MSE、RF 能量和计算，问题就不再只有一个拓扑充要条件。

## 1. 先固定协议，否则“局部”没有数学含义

令有向通信图为 (G=(V,E))，节点 (j) 初始持有 (z_j\in\mathbb R^{d_j})，只要求输出节点集合 (O\subseteq V) 中的节点计算

\[
y_i=\sum_{j\in V}T_{ij}z_j,\qquad i\in O.
\]

至少要同时固定以下六项：

1. **正确性**：零误差精确、ε-范数误差、均方误差、稳态协方差，还是失败概率不超过 δ；
2. **消息**：任意精度实数、每边每轮 (q_e(t)) 个实符号、有限域符号，还是 (b_e(t)) bit；
3. **时序**：同步轮、连续时间 ODE、每个采样时刻一次通信、异步事件触发，还是允许一个采样周期内无限次融合；
4. **输出范围**：单个 sink、任意指定子集，还是每个节点都要得到整个全局状态；
5. **失效**：独立丢包、相关衰落、crash、固定数 Byzantine，还是时变图；
6. **能耗**：发送次数、bit-hop、(E_{\rm elec}b+\epsilon b d^\nu)、MAC/干扰调度、空闲监听、重传和本地 FLOP 是否计入。

下文使用资源向量

\[
\mathsf p=(R,B_{\rm hop},E,C,M,D,p_{\rm fail}),
\]

分别表示轮数、bit-hop、物理能量、计算、峰值内存、失真和失败概率。不同论文通常只刻画其中一到两个坐标；把这些结果拼在一起不能自动得到一个共同 Pareto 面。

## 2. 四个不依赖文献的基准引理

以下是用于审计论文声称的“局部/有限时/低通信”的自包含结论，不主张它们是本文检索到的某一篇论文中的新定理。

### 引理 A：无限精度 LOCAL 模型的最弱轮数充要条件

**协议。** 同步可靠有向图；第 (t) 轮每个节点可向每条出边发送任意维、任意精度的消息，并保留无限内存；输出在第 (r) 轮结束后产生。协议可预先知道 (G,T)。

记 (d_G(j,i)) 为从 (j) 到 (i) 的有向距离。则存在 (r) 轮协议精确计算所有 (i\in O) 的 (y_i)，当且仅当

\[
T_{ij}=0\quad\text{对所有 }i\in O,\ d_G(j,i)>r.
\tag{A1}
\]

从而最少轮数为

\[
R_\star(T,G,O)=\max\{d_G(j,i):i\in O,\ T_{ij}\neq0\},
\tag{A2}
\]

约定不存在所需路径时 (R_\star=\infty)。

**必要性。** 归纳可知，(r) 轮后节点 (i) 的全部内部状态只可能依赖其 (r)-跳入邻域的数据；改变因果锥外的 (z_j) 不会改变输出。若 (T_{ij}\neq0)，就产生矛盾。

**充分性。** 沿边携带来源标签泛洪 (z_j)，第 (r) 轮后 (i) 已收到所有非零依赖块，再本地相乘求和。无需全图无环，也无需所有节点获得所有数据。

这个定理是“最 relaxed”的同时也揭示了模型过松：泛洪消息与内存可随 (n) 线性增长，且完全忽略 bit 数、数值精度、能量与离线设计成本。更聪明的树上聚合可以不传播原始数据，但**所需充分统计量的影响仍必须穿越相应割**。

对线性最小二乘，若

\[
\widehat x=(A^\top R^{-1}A)^{-1}A^\top R^{-1}b,
\]

则唯一全局估计只要求全局信息矩阵 (A^\top R^{-1}A\succ0)，并不要求每个局部 (A_i) 满列秩。若节点 (i) 只要 (L_i\widehat x)，把目标块写成 (T_{ij}=L_i(A^\top R^{-1}A)^{-1}A_j^\top R_j^{-1})，再用 (A1) 即可得到最少传播半径。目标行若稠密，精确估计不可避免地需要全网信息的影响；这不等于必须发送全部原始观测。

### 引理 B：固定维线性迭代就是稀疏因子分解

若每节点只有一个标量状态，且每轮

\[
x(t)=A_t x(t-1),\qquad (A_t)_{ij}=0\text{ whenever }j\notin N_i^-\cup\{i\},
\]

则 (T) 可在 (r) 轮实现，当且仅当

\[
T=A_rA_{r-1}\cdots A_1,\qquad A_t\in\mathcal A_G.
\tag{B1}
\]

这里 “当且仅当” 是精确的，却只是把协议设计转写成受稀疏约束的矩阵乘积分解。距离支撑条件仍然必要，但在这个受限状态模型里一般不充分。允许每节点隐藏状态/向量消息时，应把节点状态提升后再做 block-sparse factorization；可实现集合随隐藏维数而变。

### 引理 C：线性协议的割秩下界

对任意割 (S\subset V)，冻结 (S) 内输入，仅考察 (S^c) 输入对 (S\cap O) 输出的线性作用。若进入 (S) 的所有消息总实数维数为

\[
Q_S=\sum_{t=1}^R\sum_{e\in\delta^-(S)}q_e(t),
\]

则任何**线性**交互协议必有

\[
Q_S\ge \operatorname{rank} T_{S\cap O,S^c}.
\tag{C1}
\]

因为远端输入到本地输出的映射必须因子化为“远端输入 → 入割 transcript → 本地输出”。这是 message dimension/real-symbol-hop 的 converse，不是有限 bit、非线性编码或物理能耗的完整 converse；一般网络编码里，所有 cut 条件同时成立也未必充分。

### 引理 D：任意实输入的零误差精确计算不能用统一有限 bit

固定割内输入。若跨割 transcript 最多 (B<\infty) bit，接收侧至多看到 (2^B) 种 transcript；若 (T_{S,S^c}\neq0)，连续变化的远端实输入会产生不可数个正确输出。因此，不存在对所有实输入都零误差的统一有限-bit 协议。研究 bits–accuracy 时必须选有限字母表、量化范围/先验分布和失真；否则“最少 bit”不是有限数。

## 3. 有限轮任意线性变换、稀疏分解与图滤波

### 3.1 真正的有限轮线性变换

| 工作 | 协议与定理 | 是否全网传播/轮数 | 消息、内存、计算与能耗假设 | 碰撞判断 |
|---|---|---|---|---|
| Costello & Egerstedt, 2015 | 连续时间、固定连通无向图（含 self-loop）、每节点一个标量，ẋ=(W(t)x)，(W(t)) 图稀疏。Theorem 1：完整状态转移 (T) 可实现 iff (det T>0)（本地 PDF p.2）。[DOI](https://doi.org/10.1109/TAC.2014.2380643) | 所有节点最终形成 (Tx)；不是离散 LOCAL 轮。ODE 对任意正时间一般已有沿任意长路径的非零影响，没有有限传播锥。p.5 还说明任意非零标量线性泛函可在一个指定节点计算、其余输出自由。 | 实数、无限精度；全局离线设计 (W(t))。p.5 明说 shooting 只是可行性演示，并非 efficient/scalable。无 bits、内存、RF energy 或 failure。 | 是连续时间可控性结果，不是“有限通信轮 + 资源 Pareto”。奇异 consensus 在完整状态转移模型中反而不可能。 |
| Kar, Püschel & Moura, 2021 | 有向图、每节点一个复标量；Definition 1 正是 (B1)。Theorem 4.1：有 self-loop 的有向环上几乎所有 (T\in\mathbb C^{N\times N}) 可分成 (N) 个因子；Theorem 4.2 推广到含 Hamilton 有向圈的图（本地 PDF pp.2,4）。[arXiv](https://arxiv.org/abs/2104.01502) | 每节点得到自身目标行；至多 (N) 轮，随网络规模增长并允许全局信息绕圈传播。Hamilton 圈比强连通严格更强。 | 每边每轮一个标量；当前运行状态为一个标量，但若本地存全部时变权重，节点约存 (O(N(\deg i+1))) 个系数（由更新式推得）。求权重仍是全局、非凸的稀疏分解；p.5 将高效设计、少于 (N) 因子和放宽 Hamilton 圈列为后续问题。无 bits/energy/failure。 | 给出“几乎所有矩阵 + 特殊图”的 generic 充分构造，不是任意 (T)、任意强连通图的图论充要条件。 |
| Segarra, Marques & Ribeiro, 2017 | 以图 shift (S) 反复邻居交换；node-invariant filter (\sum_{l=0}^{L-1}c_lS^l)，node-variant filter (\sum_l\mathrm{diag}(c^{(l)})S^l)。Proposition 1 给 node-invariant 精确可实现的共同特征向量/重特征值条件；Proposition 4 给逐行 node-variant 条件；Corollary 2 在 (S) 特征值互异且特征向量无零元时用 (L=N) 实现任意 (B)（本地 PDF pp.3–7）。[DOI](https://doi.org/10.1109/TSP.2017.2703660) | (L-1) 次邻居交换；任意矩阵的 generic 情形需 (N-1) 轮，因此是全局尺度。 | 每次交换一个标量；node-variant 每节点存 (O(L)) 系数，在线 workspace 可常数级（由递推式推得）。系数设计需 (S) 的全局谱信息。无 finite bits/energy/failure。 | 是给定 shift 的代数刻画；不是任意局部协议，也没有资源 converse。 |
| Jiang et al., 2025 | 对 scaled all-ones/consensus 矩阵构造稀疏双随机因子（RHB/DSHB/SDS）。[DOI](https://doi.org/10.1137/24M1633790) | 目标是特殊 rank-one 变换与特定层次结构，不是任意 (T)/任意给定图。 | 实数稀疏因子；不刻画 bits、RF energy 或 failure。 | 是很近的 sparse factorization 文献，但没有终结一般问题。 |

本地全文：

- `literature/04_distributed_linear_computation/2015_costello_egerstedt_global_linear_to_local_edge_rules.pdf`
- `literature/04_distributed_linear_computation/2021_kar_pueschel_moura_finite_time_linear_transforms.pdf`
- `literature/03_sparse_inverse_graph_filters/2015_segarra_marques_ribeiro_distributed_linear_network_operators.pdf`

### 3.2 多项式图滤波与逆算子

Shuman et al. 用 Chebyshev 多项式近似稀疏对称 PSD 算子的谱乘子。一次 (K) 阶乘子需要 (K) 轮，并精确计为 (2K|E|) 个标量消息（本地 PDF p.6）；Propositions 4–5 给算子误差上界及在 (M+1) 阶光滑时 (O(K^{-M})) 的逼近率（pp.8–9）。算法只需多项式系数和谱上界，不需特征向量。[DOI](https://doi.org/10.1109/TSIPN.2018.2824239)；本地 `literature/03_sparse_inverse_graph_filters/2011_shuman_et_al_distributed_chebyshev_approximation.pdf`。

Emirov et al. 用多个可交换 graph shifts 构造多变量多项式逆滤波。Algorithm 2.2 每节点运算量

\[
O\!\left((\deg G+1)\prod_{k=1}^d(L_k+1)\right),
\]

内存

\[
O\!\left((\deg G+L_d+1)\prod_{k=1}^{d-1}(L_k+1)\right)
\]

（本地 PDF p.5），邻居通信同阶（p.6）；这些量仅在度和多项式阶数一致有界时才与 (N) 无关。Theorem 3.1 对其固定迭代给出收敛 iff (\rho(I-HG)<1)（p.7），是**算法特定**充要条件。迭代次数随精度增加，信息最终可扩散到全网。[DOI](https://doi.org/10.1007/s43670-021-00019-x)；本地 `literature/03_sparse_inverse_graph_filters/2020_emirov_cheng_jiang_sun_polynomial_inverse_graph_filter.pdf`。

结论：图滤波提供“每一阶 = 一跳”的透明实现及误差–轮数上界，但一般逆矩阵是稠密的；固定半径通常只能近似。它没有给任意图上 rounds–bits–energy 的匹配 converse。

## 4. 网络线性方程：有限时、有限 data rate 与资源账本

| 工作 | 精确结果与输出范围 | 传播与轮数 | 消息/内存/计算 | 缺失资源 |
|---|---|---|---|---|
| Yang et al., Automatica 2020 | (N) 个方程、未知量 (y\in\mathbb R^m)，堆叠 (H) 满列秩；无向连通或有向强连通。Theorem 3：任意选定节点可从自身状态历史与 Hankel 秩亏精确抽取 LS 解（本地 PDF p.5）。[DOI](https://doi.org/10.1016/j.automatica.2019.108798) | 第 (j) 个分量用 (2(D_{r,j}+1)) 个连续样本；(D_{r,j}+1) 是相应全局线性系统的可观测矩阵秩。并非由固定局部半径界定；最坏可随 (N,m) 增长。若每节点都运行抽取器，则都可输出。 | 每节点保有 (x_i,v_i\in\mathbb R^m)，每轮向邻居交换两个 (m)-向量；需存 (O(mD_r)) 历史，稠密 Hankel 求核朴素可到 (O(D_r^3))（后两项为算法结构推算）。 | 可靠无限精度实数；没有 bits、energy、failure。其步长 iff 仅针对所提迭代，不是所有协议的 iff。 |
| Gade, Liu & Vaidya, 2020 | TITAN 在可靠同步强连通有向图上用 top-(k) flooding；Theorem 2：(T\ge\mathrm{diam}(G)) 时在 (T\lceil m/k\rceil) 轮精确平均，每个节点恢复全部扰动输入（本地 PDF p.5）。LS 先聚合 (A_i^\top A_i,A_i^\top b_i)，再每节点求逆。[arXiv](https://arxiv.org/abs/2004.04680) | 明确全网传播且每节点重构全部充分统计；轮数在环/链上为 (O(m^2/k)) 量级。 | 输入维 (d) 时每节点内存 ((2k+m)d)，节点 (i) 通信量 (\deg_i^{out}(2kT\lceil m/k\rceil+1)d) 个实数单位（p.5）。按论文未压缩向量化，LS 的 (d=n^2+n)；每节点最终逆矩阵 (O(n^3))。 | 模运算仍按精确实数计，不是有限 bit；无 RF energy/丢包。隐私需要额外 weak vertex connectivity (\tau+1)，正确性只需强连通。 |
| Lei et al., SIAM J. Optim. 2020 | 无向连通图、每节点一条方程、每节点最终得到完整 (m)-维解。Theorems 3.4/3.6：任意 (K\ge1) 的有限级 zooming quantizer 均可选参数使一致线性方程指数收敛；Remark 3.7 强调三电平已足够（本地 PDF pp.8–9）。[DOI](https://doi.org/10.1137/19M1258864) | 精确等式只在 (t\to\infty)；ε-解所需轮数 (\tau(\epsilon)=O(\log 1/\epsilon))，常数含图谱和方程条件。 | 每边每轮传 (m) 个量化符号；论文利用“零由静默表达”，按每坐标 (\lceil\log_2(2K)\rceil) bit 计。Corollary 3.9：节点 (i) 的总量不超过 (m|N_i|\lceil\log_2(2K)\rceil\tau(\epsilon))，全网不超过 (2m|E|\lceil\log_2(2K)\rceil\tau(\epsilon))（p.10）。解码器需保存各邻居估计，内存 (O(m\deg i))。 | 无 joule/failure；(K,α,h) 的总 bit 最优化非凸且“global optimum generally unavailable”，论文明确列为 open（p.10）。不构成最优 Pareto。 |

本地全文：

- `literature/04_distributed_linear_computation/2018_yang_et_al_distributed_least_squares_solver.pdf`
- `literature/04_distributed_linear_computation/2020_gade_liu_vaidya_finite_time_linear_equations.pdf`
- `literature/04_distributed_linear_computation/2020_lei_peng_shi_anderson_finite_data_rates.pdf`

这些论文已经证明“无环 + 每节点局部满列秩”远非必要；但有限时精确方案要么使用无限精度状态历史，要么本质上把全网充分统计量传播到输出节点。有限 data-rate 方案则改为渐近/ε-精确，并没有匹配所有协议的 converse。

## 5. 分布式优化通信下界与有限 bit 统计估计

### 5.1 下界只对明确 oracle 模型成立

Scaman et al. 考虑连通无向图，每节点有 α-强凸、β-光滑的 (f_i:\mathbb R^d\to\mathbb R)；一次本地梯度/Fenchel 操作为 1 时间单位，传一个 (d)-向量到邻居耗时 τ，所有节点都要输出全局最优解。黑盒算法的内存只允许由既往局部一阶信息与邻居向量作有限线性张成（本地 PDF pp.2–3）。

- Corollary 1：存在难例使时间至少
  \[
  \Omega\!\left(\sqrt{\kappa_g}(1+\Delta\tau)\log\frac1\epsilon\right),
  \]
  其中 Δ 为直径。
- Corollary 2：固定 gossip eigengap (\gamma) 时至少
  \[
  \Omega\!\left(\sqrt{\kappa_l}\left(1+\frac{\tau}{\sqrt\gamma}\right)\log\frac1\epsilon\right).
  \]
- Theorem 4 的 MSDA 达到后一阶数（本地 PDF pp.4–5）。

这是匹配的 computation–communication time 下界，但只覆盖光滑强凸、无限精度 (d)-向量、一阶线性张成/gossip oracle；不计 bits、内存峰值、无线 energy 或失效。[PMLR 主页面](https://proceedings.mlr.press/v70/scaman17a.html)；本地 `literature/06_resource_tradeoffs/2017_scaman_et_al_optimal_distributed_optimization.pdf`。

Woodworth et al. 在 (M) 台机器、(R) 次间歇通信、每轮每机 (K) 次串行随机梯度的固定 oracle 下给 minimax 上下界（至多差对数因子）；这再次说明所谓“Pareto 最优”依赖 oracle 和网络结构，而不是一个无模型的普适结论。[PMLR 主页面](https://proceedings.mlr.press/v134/woodworth21a.html)。

Braverman et al. 允许 (m) 台机器与 coordinator 之间任意多轮交互，代价是 transcript 总 bit；Theorem 4.5/Corollary 4.8 对稀疏 Gaussian mean / Gaussian design 给近紧 minimax risk–communication 下界（本地 PDF pp.11,13）。它没有一般图、bit-hop 或 latency，输出也只在 coordinator。[DOI](https://doi.org/10.1145/2897518.2897582)；本地 `literature/06_resource_tradeoffs/2016_braverman_et_al_communication_lower_bounds_estimation.pdf`。

### 5.2 Gaussian CEO / 多终端率失真给出不可绕过的开放边界

Oohama 已精确求出**标量 Gaussian CEO** 的二次失真率失真函数，[DOI 10.1109/18.669162](https://doi.org/10.1109/18.669162)。这是假定长块编码、多个独立噪声观测编码器、无噪声 rate-limited 链路、单个中央 decoder 的特殊星型模型，不包含图轮数或 RF 能量。

Tavildar, Viswanath & Wagner 考虑 (N) 个 memoryless jointly Gaussian sources、分离编码、无损 bit pipes 到一个 decoder，并只重构 (x_1)：

- Theorem 1：binary Gauss–Markov tree 的完整 rate–distortion region 由“逐点 Gaussian vector quantization + Slepian–Wolf binning”达到（本地 PDF p.13）；
- Theorem 3：若观测可嵌入一般 Gauss–Markov tree，同一 separation 架构达到 many-help-one 的完整区域（p.14）；
- Section 4 仅给树条件必要性的 partial converse，并展示 separation 在一般相关结构上可任意差。

[DOI](https://doi.org/10.1109/TIT.2009.2034791)；本地 `literature/06_resource_tradeoffs/2010_tavildar_viswanath_wagner_gaussian_many_help_one.pdf`。这里的“tree”是**源统计的 Gauss–Markov 树**，不是通信拓扑无环。

更关键的是，Tang & Yang 2025 在引言中仍明确写道：一般 multi-terminal remote rate–distortion region remains open；其 Theorem 1 只是一般 achievable region，只有源在给定 decoder side information 后条件独立时内外界才在 Theorem 3 重合。[全文与 DOI](https://pmc.ncbi.nlm.nih.gov/articles/PMC12385534/)。向量 Gaussian CEO 的一般二次失真区域也长期未知；已有工作明确只给外界，[arXiv:1202.0536](https://arxiv.org/abs/1202.0536)。

**开放性归约。** 假设有人给出“任意图、任意联合统计、有限 bits、任意失真、任意局部编码”的完整资源 Pareto 充要刻画。把图选成所有传感器直接连到一个输出节点的星型，把 rounds 固定为一次长块传输，忽略 energy/compute 坐标，目标设为恢复潜在状态，就直接得到上述一般 multi-terminal remote RD region。由于这个严格子问题截至 2025/2026 仍开放，原始超宽目标至少同样难。这不是说任何受限版本都不可解，而是说明不能合理宣称用一个纯图论条件“彻底终结”全部版本。

## 6. In-network computation 与能量–时延：已经解决的是哪些窄模型

### 6.1 有限字母表、单 sink 的网络函数计算

Appuswamy et al. 考虑有限有向无环 multigraph、独立有限字母表 block sources、无噪声单位容量边、一个 receiver 与零误差 computing capacity。论文对 multi-edge tree 的任意目标以及一般网络的有限域线性目标给 tight 结果。[DOI](https://doi.org/10.1109/TIT.2010.2095070)。但 Huang et al. 2015 指出原论文面向一般函数的两个 cut-set bounds 并非普遍有效，并用目标函数诱导的等价关系修正；这说明即使在有限字母表 DAG + 单 sink 中，简单 min-cut 叙事也很微妙。[修正论文](https://arxiv.org/abs/1501.01084)。该模型没有 all-node state estimate、MSE、轮时延、内存或 RF energy。

### 6.2 树上有损线性计算

Yang, Grover & Kar 考虑无向树、每节点一个独立 Gaussian block、无噪声双向 bit pipes、指定根节点计算加权和；时间分槽且每槽仅一个节点沿一条边发送，协议调度 oblivious，可重复发边（本地 PDF pp.8–9）。

- Theorem 1：总 MMSE 等于沿树路径 incremental distortions 的累积；
- Theorem 2：由此得到比普通 cut-set 更强的 sum-rate 外界；
- Theorem 4：Gaussian random codebook 给内界；内外界仅在 (D\to0) 时相合，有限 (D) 的 gap 为 (O(\sqrt D))；
- Section V 才扩展为所有节点共识，代价相当于在每个根方向重复树计算。

因此它提供很有价值的 rate–distortion converse，却没有一般图上的精确区域，也不计物理能量、本地编码复杂度和有限 blocklength。[DOI](https://doi.org/10.1109/TIT.2017.2710059)；本地 `literature/06_resource_tradeoffs/2016_yang_grover_kar_lossy_in_network_computation.pdf`。

### 6.3 “time–energy 联合最优”确有先例，但 energy 的含义很窄

Karamchandani, Appuswamy & Franceschetti 2011 考虑 $\sqrt n\times\sqrt n$ grid geometric network，每节点一个 bit，指定 sink 计算 identity 或 symmetric function。无噪声模型允许边并行；noisy broadcast 模型中，有共同邻居的发送者不能并发，每个邻居收到独立 BSC 噪声副本。Theorems III.1–III.4、IV.1–IV.6 分别给 latency 和 transmission count 的下界及常数因子匹配方案；Theorems V.1–V.2 才把 symmetric function 的 transmission 界部分推广到一般连通图。[DOI](https://doi.org/10.1109/TIT.2011.2168902)。

这是与本项目最容易“撞车”的早期结果之一，但必须保留四个限定：

1. 输出仅为一个 sink；
2. 输入是 bit，目标是 identity/symmetric Boolean functions，不是任意实线性估计；
3. “energy complexity”被定义为**传输次数**，不是 joules；
4. 不计本地计算/内存，也没有一般拓扑的完整 latency–transmission Pareto。

Balister et al. 2011 则真正使用几何路径损耗能量：节点均匀落在增长的 (d)-维 Euclidean 区域，求和并送到指定 root，发射代价随距离的 ν 次方。Theorem 1 对附加 latency δ 给最优期望能量标度；特别地 ν<d 时为 (\Theta(n))，ν>d 时为

\[
\Theta\!\left(\max\{n,n^{\nu/d}(1+\delta)^{1-\nu}\}\right),
\]

ν=d 的临界情形带 logarithmic 项/间隙。其策略在随机网络上 order-optimal；Theorem 2 对可按 proximity-graph maximal cliques 分解的函数只在有限 latency 区间推广。[DOI](https://doi.org/10.1109/INFCOM.2011.5934949)，[arXiv](https://arxiv.org/abs/1101.0858)。它并非任意给定图、任意线性 (T)、有限 bit/MSE/compute/failure 的统一前沿。

Giridhar & Kumar 2005 在随机平面无线模型中给 divisible/type-sensitive/type-threshold 函数的 transport-capacity scaling，说明函数的可聚合结构会改变吞吐标度；同样是指定 sink、特定干扰与随机几何模型。[DOI](https://doi.org/10.1109/JSAC.2005.843543)。

### 6.4 “邻居通信必然更省电”是错误命题

LEACH 的一阶 radio model 为（本地 PDF p.2）

\[
E_{TX}(k,d)=E_{\rm elec}k+\epsilon_{\rm amp}kd^2,
\qquad E_{RX}(k)=E_{\rm elec}k,
\]

其中实验参数 $E_{\rm elec}=50$ nJ/bit、$\epsilon_{\rm amp}=100$ pJ/bit/m²。Eq. (7) 明确给出 direct transmission 比 nearest-neighbor multihop 更省能的参数区间，因为每一跳都重复支付发射电子学与接收成本（本地 PDF p.3）。所以“本地通信天然低功耗”并不成立；要看距离、包长、接收/监听、重传、MAC 和聚合压缩。[DOI](https://doi.org/10.1109/HICSS.2000.926982)；本地 `literature/06_resource_tradeoffs/2000_heinzelman_chandrakasan_balakrishnan_leach.pdf`。

Shi, Johansson & Murray 2008 的树上 Kalman 研究恰好展示资源假设反转结论：

- 固定 rooted sensor tree，逐 hop 一采样延迟；只建模 packet drop，不允许 node failure（Remark 2.3，本地 PDF p.2）；
- scheme 1 转发原始测量、几乎不本地计算；scheme 2 每个节点运行全状态 KF 并上送估计/协方差；
- 完美通信、无限计算时两者在 root 等价（Theorems 4.1–4.2，p.4）；
- 有丢包、无限计算时 scheme 2 的期望 covariance 不劣（Theorem 5.1）；
- 完美通信但本地计算受限时 scheme 1 反而不劣，因为 scheme 2 增加一步估计延迟（Theorem 5.4，p.6）；Theorem 5.5 给 packet-arrival threshold 的存在性。

[DOI](https://doi.org/10.3182/20080706-5-KR-1001.2609)；本地 `literature/06_resource_tradeoffs/2008_estimation_wireless_sensor_networks_tradeoffs.pdf`。论文没有 joule、bit 数或完整 Pareto，只是固定树与两种 architecture 的比较。

## 7. Localized Kalman / SLS：局部可实现性的真正数学核心

### 7.1 Wang–You–Matni 2015：系统响应 feasibility 是受限类中的充要条件

对固定 delayed estimator 结构，Proposition 1 证明有限 MSE 的闭环响应 ((M_w,M_v)) 可由某个稳定线性估计器诱导 iff

\[
M_w(zI-A)-M_vC=I,
\qquad M_w,M_v\in z^{-1}\mathcal{RH}_\infty
\]

（本地 PDF p.3）。加入空间 sparsity、通信 delay 和 FIR 集合 $\mathcal S$ 后，论文把“$\mathcal S$-localizable”**定义**为该 affine 约束与 $\mathcal S$ 的交非空（problem (19)，p.4）。这给出了线性稳定估计器类中的精确 feasibility test，但没有把 feasibility 化成仅由图与局部 rank 构成的闭式判据。

行可分性允许每个局部子问题只用局部 plant model；Table 1 的数值例在 2-hop、FIR (T=15) 时给 normalized MSE 1.01，并把传统 KF synthesis 记为 (O(n^3))、LDKF 在固定局部规模下记为 (O(n))（本地 PDF p.6）。这是特定稀疏网络族的 scaling，不是所有系统的定理；协议仍传实数局部信号，没有 bit/energy/failure converse。[DOI](https://doi.org/10.1016/j.ifacol.2015.10.306)；本地 `literature/05_localized_control_observability/2015_wang_you_matni_localized_distributed_kalman_filters.pdf`。

### 7.2 SLS：可达系统响应的 affine N&S，不等于 locality 总能可行

Anderson et al. 的 Theorem 2.1（有限时域）和 Theorem 4.1（无限时域）给 state-feedback achievable system responses 的 affine 充要参数化；Lemma 4.2 说明无限时域 affine 方程可行 iff ((A,B)) stabilizable（本地 PDF pp.6–7,17–18）。Definitions 2–4 再施加 (d)-localized 与 FIR-(T) 约束，并把 ((d,T))-localizable 定义为相应 SLS feasibility（pp.21–22）。若 locality region 大小对全局 (n) 一致有界，partially separable synthesis 的局部子问题规模可与 (n) 无关。

这是一套很强的数学抽象：本质是“仿射行为子空间 ∩ 时空支撑锥是否为空”。但它不自动保证任意 estimator/plant 在指定 (d,T) 下可行，也不刻画量化 bit、RF energy 或故障。[DOI](https://doi.org/10.1016/j.arcontrol.2019.03.006)；本地 `literature/05_localized_control_observability/2019_anderson_et_al_system_level_synthesis.pdf`。

Kjellqvist & Yu 2022 的 Definition 2 令支撑先随 (A^k) 扩散到 (d) 步后冻结，并**先假设**该 SLC feasible；Theorem 3 将最优无限时域 state-feedback SLS 分成局部有限时变 LQR 加 terminal Riccati，Kalman 由对偶得到。论文强调用 separation 构造的 output-feedback/LQG 一般 suboptimal；低内存 state-space 实现优于随 FIR horizon 线性增长的实现。[DOI](https://doi.org/10.1109/CDC51059.2022.9992443)；本地 `literature/05_localized_control_observability/2022_kjellqvist_yu_infinite_horizon_sls.pdf`。仍无 bits/energy/failure。

Arbelaiz et al. 2025 对 continuous-space spatially invariant Kalman gain 给了重要负面结构事实：compact support 才是真正有限通信半径；一般 Fourier square-root 不是 entire function，因此精确集中式 gain 没有 compact support，只是空间衰减（本地 PDF p.5）。Proposition IV.1 给完全去中心化 gain 的一个充分匹配条件；Theorem V.1 对高阶 diffusion 给精确指数衰减率（p.6）。因此 truncation 是可控近似，不是一般 exact locality。[DOI](https://doi.org/10.1109/TAC.2024.3504257)；本地 `literature/05_localized_control_observability/2025_arbelaiz_et_al_how_far_to_share_measurements.pdf`。

### 7.3 2026 Zhao et al. design trilemma：必须纠正的“充要”误读

[Zhao et al., EJC 2026, Article 101559](https://doi.org/10.1016/j.ejcon.2026.101559) 的模型是离散 LTI、静态有向图、每节点渐近重构**整个** (x)。三项假设是（本地 PDF p.2）：

1. **Marginally joint detectability**：堆叠 (C) 可检测，但删除任意一个 (C_j) 后均不可检测；这选择了“每个传感器都不可缺”的边界情形，并不是所有实例中最弱的必要 sensing 条件；
2. 每个 sampling instant 仅通信一次，但每节点发送当前**完整 (n_x)-维状态估计**；
3. 图为 minimal strong digraph，即强连通且删除任意边即失去强连通。

Theorem 1 的逻辑形式是：**若**存在 (Q\succ0,W,L,M) 满足给定 NMI，**则**所有节点渐近重构（本地 PDF p.4）。Remark 1 原文括注“sufficiency”（p.5）；Algorithm 1 是 cone-complementarity 式迭代 SDP，其收敛结论还要求 admissible set 非空及 complementary feasible point。结论段明确把“examining the necessary conditions of the design”列为未来工作（p.6）。

因此碰撞结论是：

- 它很好地揭示 sensing、每采样通信量/频率和 connectivity 的设计 tension；
- 它没有必要条件、没有最少轮数 converse、没有 bits/joules/compute/MSE Pareto；
- 它是 asymptotic observer stability，不是有限轮静态线性映射；
- 所以不能据此声称 2026 年已经得到“最 relaxed 充要条件”。

本地：`literature/05_localized_control_observability/2026_zhao_et_al_discrete_lti_design_trilemma.pdf`；开放版：[arXiv:2603.20144](https://arxiv.org/abs/2603.20144)，[UCL repository](https://discovery.ucl.ac.uk/id/eprint/10228689/)。

### 7.4 三篇 2026 近邻工作的 converse 审计

| 工作 | 协议与结果 | 为什么不是完整资源前沿 |
|---|---|---|
| Pérez-Salesa, Aldana-López & Sagüés, NAHS 62:101768 | 连续时间 stochastic system；传感器 measurement 和 estimator-to-estimator estimate 都异步事件触发。论文证明通过调小 event thresholds、增大 consensus gain，可把性能任意逼近 full-measurement centralized Kalman–Bucy。[DOI/摘要](https://doi.org/10.1016/j.nahs.2026.101768) | “任意逼近”对应越来越密/高增益的通信极限；仿真通信率只是占用时隙比例。没有 packet bit、bit-hop、joule、计算、匹配下界或 Pareto 必要性。 |
| Gao, Cheng, Liu & Shen, Neurocomputing 683:133440 | packet length 同时决定 quantization error 与 packet success probability；给 mean-square ultimate boundedness 的**充分条件**与误差上界，把包长/增益写成 nonlinear mixed-integer problem，用 PSO 求解。[DOI/摘要](https://doi.org/10.1016/j.neucom.2026.133440) | PSO 不证明全局最优；无必要条件、无所有协议 converse，也未同时刻画 rounds–energy–compute。 |
| Gao, Yang & Chai, Automatica 186:112837 | 时变网络、Byzantine agents、quality evaluation + multi-hop relaying，regular agents 估计完整状态并丢弃可疑消息；建立该专门 adversary/protocol 下的图与收敛条件。[DOI/摘要](https://doi.org/10.1016/j.automatica.2026.112837) | 研究重点是 resilience 与任意/有限时间收敛，不是 packet bits、RF energy、本地复杂度的 Pareto converse；结论不能外推到随机丢包或一般编码协议。 |

三篇都说明“通信更少/包更短/多跳容错”已成为显式设计变量，但没有任何一篇给任意图上的匹配 resource converse。

## 8. 为什么完整 rounds–bits–energy–compute Pareto 不能由图拓扑单独决定

至少有五个彼此独立的障碍：

1. **值域障碍**：任意实数精确估计与有限 bit 不兼容（引理 D）。
2. **网络编码障碍**：给定 edge capacities 的线性/函数计算可行性是 time-expanded network coding/稀疏多线性分解问题；cut-rank 必要但一般不自动充分。
3. **率失真障碍**：一般多终端 remote RD 本身仍开放；星型单 decoder 已是普适问题的子例。
4. **物理层障碍**：同一抽象图可有完全不同的节点间距离、path loss、干扰、接收与空闲功耗；因此 topology 不足以决定 joules。LEACH 的 Eq. (7) 甚至让 direct 比 nearest-neighbor multihop 更省电。
5. **容错障碍**：独立 erasure 的期望 MSE、worst-case crash 的可达性、Byzantine 的 graph robustness 是三个不同问题。允许“某种失败”而不说明 failure set/distribution，不存在共同充要条件。

还有一个经常被混淆的量纲问题：

- (R) 是 sequential latency；
- (∑_e b_e) 是总 bit；
- (∑_e b_e\,d_e^\nu) 才接近某类 RF transmit energy；
- 峰值节点能耗决定 network lifetime，未必等于总能耗最小；
- 本地计算节省 bit 可能增加 FLOP/内存和等待；
- 一次长 block code 的 rate 很低，但 blocklength latency 和编码复杂度可能极高。

故“轮越少、能量越低、拓扑要求越弱、计算越小”不是单调同向目标，而是一个真正的多目标设计问题。必须先定义 feasible region，再谈 Pareto；把所有坐标用一个加权和合并还可能漏掉非凸 Pareto 点。

## 9. 可发表、可证明的研究范围建议

### 建议 A：先做 deterministic linear exact 的“骨架定理”

这是最接近“原创充要条件”且不冒充解决信息论开放问题的范围：

- 固定任意有向图 (G)、输出子集 (O)、线性目标 (T)、轮数 (R)；
- 每边每轮允许 (q_e(t)) 个有限域符号，或另开一版允许实符号；
- 允许节点有明确上限的 hidden memory；
- energy 先定义为带权 symbol-hop (\sum_{t,e}w_eq_e(t))，不要称为通用 RF joule；
- failure 固定为一个明确残余图族 $\mathfrak F$，例如最多 $f$ 条 erasure/crash，暂不混入 Byzantine；
- 目标是 time-expanded graph 上 local encoder matrices 的代数 feasibility N&S，加上 causal-cone 与 cut-rank converse；在 tree/polytree/bounded-treewidth 图上争取把必要条件证明为充分，并用动态规划求完整 ((R,\text{symbol-hop},M)) Pareto。

任意图上“存在局部矩阵使端到端 transfer 等于 (T)”可以精确写成多项式方程组，逻辑上是 N&S，但算法可能困难；真正有价值的创新是找可解图类、FPT 参数或新的 sufficient cut/branch decomposition，而不是把分解定义本身包装成终结性定理。

### 建议 B：真实-valued state estimation 单独做 Gaussian approximate 层

在 A 的骨架之后再固定：线性 Gaussian plant、MSE (D)、量化范围/分布、有限 blocklength、独立 packet erasure、具体 first-order radio model。优先从 tree/polytree、单 sink 或指定输出子集开始；此时 belief propagation/Schur complement 可给局部 sufficient statistics，动态规划有希望联合优化轮数、bit allocation 与第一阶能量。不要一开始宣称任意相关源/任意图的 exact RD region。

### 建议 C：动态 localized estimation 用 SLS feasibility 作接口

固定 (d)-hop locality、FIR horizon (T) 与每采样一次消息，先研究

\[
\text{SLS affine constraints}\cap\text{locality cone}\cap\text{FIR cone}\neq\varnothing
\]

的图/系统结构；输出可只覆盖指定 state blocks。把 bits 与 RF energy 作为第二层近似/鲁棒扰动，而不是在首个 theorem 同时求全资源 converse。这样既承接 LDKF/SLS，又能明确超越 Zhao 2026 的纯充分 NMI：例如在一类 chordal/bounded-treewidth systems 上得到可验证 N&S 和最小 (d,T)。

### 最终状态标记

- **SOLVED**：固定同步、可靠、无限精度且消息/内存不限时，线性目标的最少轮数就是目标依赖的最大有向距离；固定 scalar-state linear iteration 时，可实现性等价于图稀疏因子分解。
- **PARTIAL**：Hamilton 圈上的 generic 任意线性变换、给定 shift 的图滤波、特定迭代的收敛 iff、树上/特殊 Gaussian RD、受限 oracle 的优化下界、localized SLS/Kalman feasibility、特定无线模型的能量–时延标度。
- **OPEN-LOOKING**：任意图 + 任意统计源/实数目标 + finite bits + finite rounds + physical energy + compute/memory + failures 的完整 Pareto 充要刻画。这个范围包含截至 2026 年仍开放的多终端 remote RD，故当前不能诚实地称为已解决，也不宜作为一篇论文的单一定理目标。

## Sources

- [Kar, Püschel & Moura, *Finite-time, in-network computation of linear transforms*, arXiv:2104.01502](https://arxiv.org/abs/2104.01502)
- [Anderson et al., *System Level Synthesis*, Annual Reviews in Control 47 (2019)](https://doi.org/10.1016/j.arcontrol.2019.03.006)
- [Scaman et al., *Optimal Algorithms for Smooth and Strongly Convex Distributed Optimization in Networks*, PMLR 70 (2017)](https://proceedings.mlr.press/v70/scaman17a.html)
- [Yang, Grover & Kar, *Rate Distortion for Lossy In-Network Linear Function Computation*, IEEE TIT 63(8)](https://doi.org/10.1109/TIT.2017.2710059)
- [Tang & Yang, *Rate-Distortion Analysis of Distributed Indirect Source Coding*, Entropy 27:844 (2025)](https://doi.org/10.3390/e27080844)
