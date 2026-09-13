# 网络化静态 WLS / nodewise state estimation 的直系谱系尽调（2011/2013--2026）

> 检索与复核截止：2026-09-13（Asia/Shanghai）  
> 任务边界：以用户给出的“图无环 + 每个本地 (A_i) 满列秩”的稿件为种子，追踪其静态、部分重构、邻居消息传递 WLS 直系谱系；同时检查容易被误认为“已经给出最宽充要条件”的相邻动态观测器、Kalman、Gaussian BP 工作。  
> 证据分级：A = 已逐页复核全文/定理；B = 出版商摘要或可见正文片段，全文未取得；C = 只用于划边界的相邻文献。凡属 B 级，本文不猜定理编号、页码或未公开假设。

## 0. 先给结论

1. **[证据 A：全文定理] “无环”不是任何分布式状态估计方法存在的必要条件。** 它只是这条谱系中 cavity / Gaussian-BP 型小消息精确消元在有限步内不重复计数的一个强而便利的充分条件。2015 年的 Richardson 型方法已能在一般含环图上渐近得到精确全局 WLS；2018 年 Wu 等还展示了在一般连通通信图上先分布式抽取生成树、再用全局状态维信息消息做有限步精确融合。[^4][^7]

2. **[证据 A：全文定理] “每个 (A_i) 满列秩”也不是必要条件。** 2015 年 Marelli--Fu 已把它明确放宽到聚合测量矩阵 (H)（其记号为 (A)）满列秩；其电力系统短文甚至明说本地区域可以不可单独辨识。[^4][^5] 对“只要节点 (i) 的目标分量”而言，最宽的信息论条件还可放宽为
   
   \[
   P_i\ker H=\{0\}
   \quad\Longleftrightarrow\quad
   \operatorname{row}(P_i)\subseteq\operatorname{row}(H),
   \]
   
   即所有与全体测量不可区分的状态，只要在目标 (P_i x) 上也不可区分即可；根本不要求整个 (x) 唯一。

3. **[证据 A + 检索判定] 直系谱系目前没有“所有局部算法存在”的最宽充要条件。** 找到的“iff”均有较窄量词：2020/2021 TAC 的 iff 是“该特定 GaBP 递推是否收敛”；2026 Fattore 等的 iff 是“在 Assumptions 1--2 下，是否存在式 (12)--(13) 这一结构的动态观测器”；2026 Jaiswal 等的 iff 是“在强连通/平衡或无向连通图假设下，是否存在式 (13) 的连续时间功能观测器”。它们都不是静态 nodewise WLS、有限半径、有限消息预算下的任意算法 iff。[^9][^14][^15] 这里“没有找到”是截至截止日的系统检索结论，不是对全世界文献不存在性的形式证明。

4. **[逐项证据见括注；2022/2023/2026 JFI 为 B，其余为 A] 直系工作实际到达的边界是：**

   - 向量节点、无环信息图：全局满列秩即可有限步精确（2015）。
   - 标量节点、pairwise edge measurements、任意连通含环图：在文中非退化边测量与至少一个 anchor 的特殊模型下，特定 BP/WLS 算法渐近精确（2020）。
   - 向量节点、含环图：2022 年 generalized block diagonal dominance 只给特定算法渐近精确的**充分条件**，不是必要条件。
   - 含环图上的有限步小消息算法：2018/2021 的结果给误差随 loop-free depth 指数变小，通常不是固定含环实例上的精确中央 WLS。
   - 2023 年 TCAS-II 的摘要仍明确把有限步结论限定为无环图；2026 年 JFI 工作改成动态公共状态及不完整测量，且只声称与中央 Kalman “comparable”，不是同一问题的精确 iff。[^6][^8][^10][^11][^13]

5. **[本文数学推导，不冒充既有文献定理]** 在**固定线性实例、固定目标、任意精度实数消息、允许泛洪**的明确计算模型下，可以写出一个干净而真正最宽的 (r)-轮 iff。令

   \[
   K_R=(H^\top R^{-1}H)^{-1}H^\top R^{-1},\qquad T_i=P_iK_R,
   \]

   并令 (S_{i,r}) 选择在 (r) 轮内能到达节点 (i) 的测量坐标，则

   \[
   \boxed{\exists\ r\text{-轮精确局部算法}}
   \quad\Longleftrightarrow\quad
   \boxed{\exists L_{i,r}:\ T_i=L_{i,r}S_{i,r}}.
   \]

   等价地，(T_i) 对所有 (r) 轮内不可见的测量扰动为零。必要性来自 LOCAL 模型的不可区分性；充分性由在 (r) 轮内把所需数据沿路由收集到 (i)，再应用 (L_{i,r}) 得到。它允许不做全局覆盖：中央答案确实不依赖的数据无须到达。**但一旦再要求消息维度与网络规模无关、内存 (O(\deg i))、容错概率、或者模型矩阵也只局部可知，问题就变成另一套资源受限 iff；直系文献没有把这些约束统一解决。**

6. **[明确的否定性尽调判定]** 截至 2026-09-13，未发现一个已发表定理同时覆盖：fixed-radius LOCAL 通信、模型矩阵仅局部持有、受限消息/内存，以及允许失败的 ((epsilon,delta)) 近似，并对“任意算法是否存在”给问题级 iff。A 级全文只分别解决其中若干切片；B 级摘要也未声称这一统一结论。因此本项目若要研究这个联合边界，不能把 2021 的谱半径 iff 或两篇 2026 动态观测器 iff 当作撞车。

## 1. 种子论文究竟是哪一年、解决什么问题

本地种子文件是 `E:\Code\MAS\[ACC_11]A Distributed Algorithm for State Estimation with Application to Smart Grid.pdf`。PDF 元数据显示创建于 2012-09-10，题名/作者元数据仍是模板占位符；没有证据表明它在 2011 年正式发表。可核实的正式出版物是 Tai、Lin、Fu、Sun 的 2013 American Control Conference 论文 *A New Distributed State Estimation Technique for Power Networks*，页 3338--3343，DOI `10.1109/ACC.2013.6580347`。[^1][^2]

两份文本的核心模型一致：

\[
z_i=A_i x_i+w_i,
\qquad
z_{ij}=B_{ij}x_i+B_{ji}x_j+w_{ij},
\]

噪声块独立、协方差已知。节点 (i) 只需输出自己分块 (x_i) 的估计，但拼接后的结果要**等于中央 WLS**。运行时使用本地节点测量、相邻边测量、局部拓扑和邻居发来的边界消息。原假设为：

- Assumption 1：物理/测量图 (G) 无环；
- Assumption 2：每个节点的自测矩阵 (A_i) 满列秩。

正式 ACC 版 Theorem 3.1（本地 PDF 第 4 页，证明跨第 4--5 页）说：节点 (i) 的协方差在其 eccentricity (ε_i) 轮后等于中央 WLS 对应块；借 Lemma 3.1，把协方差相等转成估计相等；全网在直径 (D) 轮后完成。消息是 ((\beta,\Phi))：对边 ((i,j))，(β) 是边测量维数 (r_{ij}) 的向量，(Phi) 是 (r_{ij}\times r_{ij}) 矩阵，因而不携带整个全局状态。[^2]

未发表的 8 页稿在 Theorem 2（PDF 第 4 页）给同一理想通信结论，并在仿真/结论中声称丢包、延迟、异步下仍几乎必然有限时间收敛；但正文只定义了 Assumptions 1--2，仿真段却写“Assumptions 1 to 3”，不存在 Assumption 3，也没有对应容错定理。因此该稿的丢包结论只能记为**主张 + 仿真，不是定理级保证**。真正的容错定理出现在 2015 SYSID 短文的 Assumption 3 / Corollary 6。[^1][^5]

## 2. 四种不能混用的“条件量词”

后续论文看似都在谈 convergence / optimality，但实际有四个不同问题：

1. **中央目标是否定义良好**：全状态唯一 WLS iff (H^\top R^{-1}H\succ0)，在 (R\succ0) 时等价于 (ker H=0)。
2. **某个目标是否可辨识**：只需 (P_i\ker H=0)，不必全状态唯一。
3. **某个已经指定的消息递推是否收敛**：例如谱半径小于 1；算法发散不代表别的算法不存在。
4. **在明确通信/消息/内存/轮数模型下，是否存在任意算法**：这是用户问的最宽存在性 iff，必须先固定资源模型。

本报告凡写“充要条件”，都会标明它属于哪一层。

## 3. 直系谱系：定理级总表

### 3.1 全文已复核的核心工作（A 级）

| 年份与论文 | 模型、每节点目标 | 拓扑与秩/可观条件 | 保证、轮数与定理量词 |
|---|---|---|---|
| 2013 Tai--Lin--Fu--Sun[^2] | 静态 node/edge 向量测量；节点 (i) 只输出 (x_i) | (G) 无环；每个自测 (A_i) 满列秩 | Theorem 3.1，PDF p.4：精确中央 WLS；节点 eccentricity 轮，全网 diameter 轮。是该算法的充分条件，不是存在性 iff。 |
| 2015 Marelli--Fu，Algorithm 1[^4] | 一般块测量 (y_i=\sum_j A_{ij}x_j+v_i)；节点只输出本块 | 聚合 (A) 满列秩、(R) 非奇异；允许一般含环信息图；能和测量邻居通信 | Richardson 递推，式 (6)--(8)，PDF pp.3--4；(0<\gamma<2/\|\Upsilon\|) 时**渐近精确**，无限轮。Theorem 7（p.5）给可安全选的缩放范围。树不是必要条件。 |
| 2015 Marelli--Fu，Algorithm 4[^4] | 同上；Schur-complement/cavity 消息 | 全局满列秩；Assumption 11；信息图 ((i,j)\iff\Psi_{ij}\ne0) 无环 | Theorem 15（p.6，证明 p.12）：节点 (i) 在半径 (ρ_i) 轮后精确，全网直径轮。无环是该小消息有限消元算法的充分条件。 |
| 2015 Marelli--Ninness--Fu[^5] | 电网 local/boundary measurement；每区只输出本区状态 | 区域图无环；聚合 (H) 满列秩、(R_\nu) 可逆；**不要求本区单独可观** | Theorem 3（刊页 565）：半径/直径轮精确中央 WLS。Assumption 3 + Corollary 6（刊页 565）：只要任意时刻之后都有同样新或更新的数据包最终到达且时钟错位有界，则几乎必然存在有限 (t^*) 后精确。 |
| 2018 Sui--Marelli--Fu--Lu[^6] | 2013 式 node/edge 向量测量；节点只输出本块 | 任意无向含环图；但 Assumption 2 仍要求 (C_i^TR_i^{-1}C_i\succ0)，即每个自测满列秩；另需稀疏度/SNR 收缩常数 (\rho<1) 或 (\kappa<1) | Theorem 11（p.6）协方差误差 (le c\rho^{\ell_i})；Theorem 17（p.7）状态误差 (le c\kappa^{\ell_i+1})，其中 (ell_i) 是 loop-free depth。是跨实例/节点的近似界；不等于固定含环图随轮数趋于中央解。 |
| 2018 Wu--Fu--Xu--Lu[^7] | 动态**单一公共全局状态**，每节点都估整个 (x(t))；WLS 是无先验特例 | 原通信图无向连通；主融合先用分布式 DFS 去掉通信环成为生成树；聚合系统可观 | Theorem 1（p.3）：生成树上 diameter 轮精确中央 KF/ML；Sec.3.3（pp.4--5）任意连通含环通信图可先抽树；Sec.4.4（p.5）给 WLS 特例。Lemma 4（p.4）在其“每邻居每轮一次”约束下给 diameter 下界。消息为整个公共状态维的矩阵/向量，故不解决 nodewise 稀疏消息问题。 |
| 2020 Yang--Zhang--Fu[^8] | 静态 self/edge measurements；向量或标量节点，各输出本块 | Assumptions 1--2（p.2）：连通、至少某节点有自测；标量循环结果还隐含每条边确实依赖两端（非零系数） | Theorem 1（p.4）与 GaBP 等价；Lemma 1（p.4）仅标量证明 comparison matrix (\bar\Psi\succ0)；Theorem 2（p.5）：向量树 diameter 轮精确，标量含环图渐近精确并有指数率。**标量结论很强但模型特殊；不是一般向量 iff。** DOI 的文章号是 **109091**，不是 109090。 |
| 2020 预印本 / 2021 TAC，Marelli--Sui--Fu--Sun[^9] | 静态 pairwise 向量 Gaussian 模型；各节点求本块 ML/WLS 边缘 | 一般含环图；Assumption 1（p.4）为局部严格信息支配。非叶节点可没有满秩自测，但叶节点仍须自测满秩 | Theorem 1（p.5）信息矩阵指数收敛。Theorem 2（p.6）：在 Assumption 1 下，**该算法**的均值估计收敛 iff 约化非回溯消息矩阵 (̅A(\infty)) 稳定；不稳定时几乎所有测量发散。Theorems 4--5（pp.8,10）给 loop-free depth 准确度界。这个 iff 既不保证固定含环图收敛值等于中央 WLS，也不排除其他算法。 |
| 2026 Fattore--Valcher--Gao--Yang，arXiv:2603.10656[^14] | 离散时间 LTI **动态观测器**；(y_i=C_ix) 无测量噪声；每节点估计整个公共 (x(t)) | Assumption 1（p.3）联合可检测；Assumption 2（p.4）为局部 Jordan 输出向量线性无关；有向图可含环。每个不稳定 Jordan miniblock 对应检测节点集和 grounded Laplacian | Theorem 4（p.6）：在 Assumptions 1--2 下，存在**式 (12)--(13)** 的观测器 iff 对每个不稳定 ((\ell,h)) 有实增益 (k^{\ell,h}) 使 (lvert1-k^{\ell,h}\mu\rvert<1/lvert\lambda_\ell\rvert) 对 grounded Laplacian 的所有特征值成立。渐近时域估计，不是静态 WLS、有限轮或任意算法 iff。 |
| 2026 Jaiswal--Berger--Tomar，arXiv:2604.00680[^15] | 连续时间 LTI 公共状态；节点 (i) 测 (y_i=C_ix+D_iu)；所有节点估同一个功能 (z=Kx) | Assumption 6（p.3）：有向图须平衡且强连通，或无向连通；聚合 (\tilde C) 对 (K) 联合部分可检测 | Theorem 10（p.4）给中央 partial detectability 的秩 iff；Lemma 15（p.6）换成聚合 (\tilde C)；Theorem 17（p.6，证明 p.7 起）：在 Assumption 6 下，存在**式 (13)** 的分布式功能观测器并渐近 omniscience iff 联合部分可检测。它不是噪声 WLS、不是离散通信轮、也不是每节点不同目标的定理。 |

### 3.2 2022、2023、2026 三个易误报节点（B 级）

| 论文 | 公开材料真正说了什么 | 是否放宽种子假设、是否“撞车” | 公开全文状态（截止检索日） |
|---|---|---|---|
| Yang--Zhang--Fu--Cai, SCL 2022[^10] | 研究 2015 Algorithm 4 在**向量节点、含环图**上的收敛；在 generalized block diagonal dominance 下，输出渐近趋于全局最优 WLS。 | 放宽“无环”，但代价是换成矩阵支配性；这是**特定算法正确收敛的充分条件**，摘要没有声称必要性。没有最宽 iff。 | 出版商正式页可读摘要/部分正文；OpenAlex 与 Semantic Scholar 未给合法 OA PDF。故主定理编号、页码、完整假设编号不在本报告臆测。 |
| Zhu--Wang--Sun, TCAS-II 2023[^11] | 摘要：基于 coupled measurements 和邻居信息，每节点所得等同中央 WLS；**当通信网络为无环图时**有限次收敛；abstract graph 的大小决定迭代数；另做虚假数据检测。 | 没有从摘要看到解除无环。是否放宽局部满秩也无法由公开摘要确认。标题中的 “optimal / finite steps” 不能解释成一般含环图定理。 | IEEE 正式落地页、DOI、摘要可核；未找到合法 OA PDF。定理号、页码、消息维度和秩假设必须等取得全文后再填。 |
| Zhu--Wang--Sui--Wu, JFI 2026[^13] | 离散时间、时变公共状态 (x_{k+1}=A_kx_k+w_k)，传感器测量含 Bernoulli/incomplete 因子；每节点估**整个公共状态**。无环图上做有限次消息迭代，含环图用 loop elimination；摘要只说精度与中央 Kalman **comparable**。 | 允许含环的方式仍是删通信环；模型、目标和“精确”标准均已改变。不是静态 nodewise WLS 的放宽，更没有该问题的必要充分条件。 | ScienceDirect 正式摘要和正文片段可读，但 PDF 需机构/购买；未找到合法 OA 版。不能从摘要推出“等于中央 KF”。 |

特别提醒：2023 JAI 的 Cai--Zhang--Fu 是综述/统一说明，不是新 maximal iff。它把低复杂度解释为每轮一次邻居交换、每节点计算和存储随度数线性，并重述树上 diameter 轮及某些含环情形的渐近结果。[^12]

## 4. 每篇的通信对象、消息维数与本地知识

| 工作 | 运行时可用信息 / 本地知识 | 每条有向边每轮消息量级 | 备注 |
|---|---|---|---|
| 2013 Tai 等 | 本地 (A_i,R_i,z_i)，相邻 (B_{ij},B_{ji},R_{ij},z_{ij})，邻居列表 | (β\in\mathbb R^{r_{ij}})、(Phi\in\mathbb R^{r_{ij}\times r_{ij}}) | separator/边测量维；树上不重复计数。 |
| 2015 Marelli--Fu Algorithm 1 | 测量节点知道自身整行块 (A_{i,j})、(R_i)；节点间按 (N_i) 转发状态分块；步长/预条件器需估计或界定全局谱量 | 主要传 (x_j(t)\in\mathbb C^{d_j})，以及初始化信息 | 任意图但无限轮；Richardson 稳定性，而非 BP 有限消元。 |
| 2015 Marelli--Fu Algorithm 4 | 节点 (i) 有 (alpha_i) 及相邻 (Psi_{j,i})；Assumption 11 明确邻居通信 | 发给 (j)：(gamma_{i,j}\in\mathbb C^{d_j})、(Gamma_{i,j}\in\mathbb C^{d_j\times d_j}) | 消息按接收方状态分块维数，不随总节点数直接增长。 |
| 2015 SYSID | 本地 reduced measurement 及相邻边界测量模型/协方差 | 边界维向量 (gamma) + 同维方阵 (Upsilon) | 有丢包时保留最新到达消息；Corollary 6 依赖最终新鲜度。 |
| 2018 Sui 等 | 与 2013 类似；每节点要能反演自测信息 | 发给 (j)：(alpha_{i\to j}\in\mathbb R^{n_j})、(Q_{i\to j}\in\mathbb R^{n_j\times n_j}) | 含环时构造 computation tree；小消息但一般近似。 |
| 2020 Yang 等 | 本地 self/edge measurements 及其协方差 | ((x_{i\to j},\Sigma_{i\to j}))，向量情形为 (n_i) 向量 + (n_i\times n_i) 矩阵；标量为 2 个标量 | 含环精确渐近只证明标量。 |
| 2021 TAC GaBP | 同类 self/edge Gaussian factors；运行时只用邻居信息参数 | factor-to-variable 信息向量 + 信息矩阵，约 (n_j+n_j^2) | 检查 Theorem 2 的全局非回溯谱条件本身不一定是纯本地一次性知识；Theorem 3 提供较易分布检查的充分条件。 |
| 2012 Pasqualetti 等[^3] | 每 monitor 知 (H_i,\Sigma_i,z_i)；一般连通通信图 | 传整个增广估计 (hat\xi_i\in\mathbb R^{n+p}) 及其 nullspace basis (K_i)，最坏可达矩阵级 | diameter 轮得到 (hat x(\epsilon))，但固定 (\epsilon>0) 不是精确 WLS。它证明“含环通信可有限传播”，却不是可扩展 nodewise 小消息解。 |
| 2018 Wu 等 | 所有节点最终都保有公共全状态估计；建树需唯一数值 ID、最大共识和 DFS | (Q_{i\to j}\in\mathbb R^{n\times n})、(alpha_{i\to j}\in\mathbb R^n) | (n) 是**全局目标状态维**；节点数增大而状态维固定时可称低复杂度，但不满足原问题按局部状态维扩展的含义。 |
| 2026 Fattore 等 | 运行时交换邻居 Jordan miniblock 估计；设计需系统 (A) 的统一 Jordan 变换、各 (C_i) 的可检测分解、grounded Laplacian 谱和可行增益 | 对每个未本地完整可观的 Jordan miniblock 交换其估计；总量可接近全状态维 | 运行时局部通信，离线设计并非仅本地知识。 |
| 2026 Jaiswal 等 | 运行时交换内部状态 (w_j)；Algorithm 1 设计要聚合 (C̃)、构造 (T,\bar L,\bar P)、计算 (L+L^T) 的谱隙并选全局 (gamma) | 邻居内部观测器状态 (w_j\in\mathbb R^q) | 连续时间共识，没有“第 (r) 轮精确”概念；离线需要全局系统/图谱量。 |

## 5. 两个原假设分别被放宽到了哪里

### 5.1 从“每个 (A_i) 满列秩”到目标函数可辨识

放宽链条是：

\[
\boxed{\forall i,\ \operatorname{rank}A_i=s_i}
\Longrightarrow
\boxed{\operatorname{rank}H=\dim x}
\Longrightarrow
\boxed{P_i\ker H=0}.
\]

- 第一箭头严格：每个本地自测已经锁定本地状态，当然保证聚合 (H) 满列秩，但远非必要。
- 2015 Marelli--Fu / SYSID 已到第二层：允许某个或多个节点单独不可辨识，只要全网联合测量矩阵满列秩。[^4][^5]
- 2020 Yang 等在**标量 pairwise** 特例进一步给出非常易检条件：连通 + 至少一个非零自测 anchor，再结合每条边对两端的非零相对信息，使比较矩阵正定。这只是保证全局满秩的一种特殊充分结构，不是一般矩阵的必要条件。[^8]
- 若只要求 (P_i x)，第三层才是信息论最宽条件。它等价于存在矩阵 (L_i) 使 (P_i=L_iH)，即目标可由测量确定。2026 Jaiswal 的“partial detectability”在**动态连续时间**语境表达了相似思想：不可检测的不稳定模态可以存在，只要 (Kx) 不依赖它；但它不是静态 WLS 定理。[^15]

若题目坚持输出某个指定的 Moore--Penrose 最小范数解 (H^\dagger z)，则即使 (H) 秩亏也可定义一个单值算法目标；但这已改变了“物理状态/所有 WLS 极小点共享的节点分量”这一估计语义，必须在论文中明说。

### 5.2 从“无环”到什么

无环的地位取决于资源模型：

- **不限制轮数**：2015 Richardson 方法在一般含环图上渐近精确；因此树绝非必要。[^4]
- **不限制消息维数、允许泛洪/抽树**：任意有限连通图都能先建立生成树并传播全部充分统计量；2018 Wu 等给了这一思路的正式实例。[^7]
- **坚持 cavity/GaBP 小消息**：树上有限步精确；标量 pairwise 情形 2020 年可在一般环图渐近精确；向量含环目前直系结果主要是谱/支配性充分条件或算法特定 iff。[^8][^9][^10]
- **坚持含环图、少轮、允许误差**：2018 和 2021 结果表明误差随节点周围的 loop-free depth 指数衰减，但这描述空间局部性/近似性，不是无条件精确。[^6][^9]
- **坚持有限轮且消息宽度仅与局部块有关**：自然参数不是“是否含环”这个二值量，而是消元 separator / treewidth。树宽 1 时原算法最漂亮；一般图可用 junction tree 精确化，但消息矩阵维数随 separator 状态维增长。直系 WLS 论文没有给出这一资源--拓扑联合 iff。

### 5.3 2020 Yang 向量树陈述的一个技术警告

该文 Assumption 2 字面只说图中“至少有某些节点具有一个 self measurement”。Lemma 1 用标量 (x_i\in\mathbb R)、非零边系数和连通性证明 (Ψ̄\succ0)，从而推出 (H) 满列秩；这个证明不能自动推广到任意向量块。Theorem 2 的树分支却写成向量也成立。若没有另加全局 (H) 满列秩或足够的矩阵非退化条件，简单维数反例即可使 WLS 不唯一。因此其向量树结论应按 2015 Theorem 15 的 well-posedness 条件阅读；不能把“连通 + 一个自测”误报为一般向量 WLS 唯一性的充分条件。[^8]

## 6. “最 relaxed”充要条件：必须把信息条件与通信条件拆开

### 6.1 中央问题层：全状态与局部功能

设 (R\succ0)。中央 WLS 目标

\[
J(x)=(z-Hx)^TR^{-1}(z-Hx)
\]

有唯一全状态极小点 iff

\[
H^TR^{-1}H\succ0
\iff \ker H=0.
\]

若只关心节点功能 (P_i x)，则所有 WLS 极小点给出同一答案 iff

\[
P_i\ker H=0.
\]

这是秩条件的最宽形式，和图是否有环无关。

### 6.2 固定实例、(r) 轮 LOCAL 层：线性目标的精确 factorization iff

先假设 (H,R) 固定且算法设计时已知；测量值按节点/边初始分布。令 ({\cal A}_{i,r}) 是在有向通信图上 (r) 轮内可到达 (i) 的测量坐标，(S_{i,r}) 是相应选择矩阵。唯一 WLS 时，中央 nodewise 输出是线性映射

\[
y_i^*(z)=T_i z,
\qquad
T_i=P_i(H^TR^{-1}H)^{-1}H^TR^{-1}.
\]

那么

\[
\exists\text{任意确定性 }r\text{-轮局部算法精确输出 }y_i^*(z)\ \forall z
\]

当且仅当

\[
\ker S_{i,r}\subseteq\ker T_i
\quad\Longleftrightarrow\quad
\operatorname{row}(T_i)\subseteq\operatorname{row}(S_{i,r})
\quad\Longleftrightarrow\quad
T_i=L_{i,r}S_{i,r}\text{ for some }L_{i,r}.
\]

**必要性**：两组全局测量若在 (i) 的 (r)-hop view 中完全相同，任何 (r)-轮算法的输出必相同；故它们之差必须落在 (ker T_i)。  
**充分性**：所有被 (S_{i,r}) 选中的坐标都可在 (r) 轮内沿路径送达 (i)，节点直接计算 (L_{i,r}S_{i,r}z)。

这就是允许“不全局覆盖”时的精确答案：只有 (T_i) 非零依赖的测量需要到达。若仅做无噪声可辨识，类似条件为

\[
\ker H_{{\cal A}_{i,r}}\subseteq\ker P_i,
\]

但它和“对每个噪声实现精确复制中央 WLS”不是同一个条件；后者必须检查完整的 (T_i) 支撑。

### 6.3 丢包/失败层

令 ({\cal A}_{i,t}(\omega)) 是失败样本路径 (omega) 下，截止时间 (t) 实际到达节点 (i) 的测量集合。则：

- 在给定失败样本路径上，时刻 (t) 精确 iff (T_i) 可通过 (S_{{\cal A}_{i,t}(\omega)}) 因子化；
- 几乎必然最终精确 iff 对几乎所有 (omega)，存在有限随机 (t_i(\omega)) 使上述因子化成立；
- 若每个“必要测量”都有一条存活路径且独立重传最终成功概率为 1，则可满足；
- 若某次失败永久切断一个 (T_i) 非零依赖的数据，则没有任何局部算法能保持对所有 (z) 精确。

2015 SYSID Assumption 3 / Corollary 6 正是“新鲜或更新数据最终到达”这一充分机制的算法化版本，但没有把它提升成一般 (T_i)-支撑 iff。[^5]

### 6.4 为什么还不能说“彻底终结”所有资源版本

上述 factorization iff 允许消息携带任意多实数、节点知道计算系数，并且不计建路由/拓扑发现成本。一旦限制：

- 每条边每轮最多 (b) 比特或固定维数；
- 节点存储不随网络规模增长；
- (H,R) 也只局部持有，算法要对一族未知实例统一工作；
- 只能 (r) 轮、容许 (epsilon) 误差或失败概率 (delta)；
- 能量按发送次数、距离、比特数或矩阵运算共同计价；

就需要通信复杂度、treewidth/separator、逼近误差及容错可达性的联合条件。不同预算下的“最 relaxed”条件不同，不存在脱离资源模型的单一图性质。直系论文分别优化其中一个切片，没有给出统一 Pareto 边界。

## 7. 容易混淆但不真正撞车的相邻文献

### 7.1 Pasqualetti--Carli--Bullo 2012

该文考虑一般 (z=Hx+v)、任意连通 monitor 通信图，每个 monitor 最终估计**整个全局状态**。Theorem 3.1（PDF p.6）用迭代投影得到一致线性方程的最小范数解；对噪声 WLS，Theorem 3.2（p.8）只证明 (hat x(\epsilon)\to x_{WLS}) 当 (epsilon\downarrow0)。Theorem 3.3（p.9）在图直径 (d) 轮得到固定 (epsilon) 的 (hat x(\epsilon))，Corollary 3.3（p.10）给异步 (dT)；因此固定 (epsilon>0) 仍是近似 WLS。Theorem 4.1（p.12）给局部误差指数衰减。消息含全局增广变量和 nullspace basis，故它证明环图通信并非障碍，却不满足原问题的局部块小消息可扩展性。[^3]

### 7.2 Wu--Fu--Xu--Lu 2018

这篇是最有价值的“反例控制”：它确实在任意连通含环通信图上有限步复现中央 ML/KF/WLS，但每节点目标是同一个全局 (n)-维状态，消息是 (n\times n) 信息矩阵加 (n)-向量，并先删除通信环。它终结的是“环图绝对不可能有限步精确”这一错误命题，不是“原 nodewise、separator-sized 消息在原含环图上何时存在”的问题。[^7]

### 7.3 Gaussian BP 的 iff

Du 等 2018 JMLR 对 vector-valued Gaussian BP 给出均值/方差收敛条件；2021 TAC 也给本谱系算法的约化消息矩阵稳定性 iff。二者都属于**指定递推的收敛性**。一个算法的迭代矩阵谱半径 (ge1)，只说明该递推不行，不说明可以泛洪、共轭梯度、嵌套分解或 junction tree 的其他局部算法不存在。[^9][^16]

### 7.4 2021 clustered WLS + GBP 与 2023 动态 accuracy

Živojević 等把 cluster 内 WLS 与 cluster 间 GBP 集成，主要是仿真上更快/更高成功率，结论把一般收敛条件列为未来工作，并未给 maximal iff。[^17] Zhu、Wang、Miao、Sui 2023 研究动态 MAP 在含环图上的精度界；其基准在无环情形本身也未必等于中央最优，故不是本静态精确 WLS 问题的回答。[^18]

### 7.5 两篇 2026 arXiv 新文献的严格量词

**Fattore et al., arXiv:2603.10656。** Theorem 4 的量词是：固定系统 (4)--(5)、固定有向通信图，先假设联合可检测和技术性 Assumption 2，再问对每个不稳定 Jordan miniblock 是否能选一个实 coupling gain，使所提出的 Luenberger + consensus 结构 Schur 稳定。其条件同时编码了：检测该 miniblock 的节点到其他节点必须有 rooted spanning forest，以及 grounded Laplacian 特征值能被同一个实增益映到半径 (1/|\lambda_\ell|) 的圆盘内。它比 Gao--Yang 2025 的统一 gain 条件宽，但只对该观测器结构是 iff；估计随物理时间渐近，非有限通信轮。公开 PDF：`https://arxiv.org/pdf/2603.10656`；本地归档：`E:\Code\MAS\literature\05_localized_control_observability\2026_fattore_et_al_discrete_lti_jordan_observer.pdf`。[^14]

**Jaiswal et al., arXiv:2604.00680。** 其秩条件使用

\[
D_{[A,C],n,\lambda}
=\operatorname{col}\big(C,C(\lambda I-A),\ldots,C(\lambda I-A)^{n-1},(\lambda I-A)^n\big).
\]

Theorem 10 表明对所有 (lambda\in\mathbb C_+)，

\[
\operatorname{rank}
\begin{bmatrix}D_{[A,C],n,\lambda}\\K\end{bmatrix}
=\operatorname{rank}D_{[A,C],n,\lambda}
\]

iff 中央系统对 (K) 部分可检测。Lemma 15 把 (C) 换成所有节点测量矩阵堆叠的 (C̃)。Theorem 17 再说：**在 Assumption 6 的图类内**，式 (13) 的分布式观测器存在并使所有节点渐近估计同一 (Kx) iff 联合部分可检测。论文引言所谓“无需额外系统/图假设”应读成“除已经规定的 Assumption 6 之外不再加条件”，不能读成任意弱连通/故障图。公开 PDF：`https://arxiv.org/pdf/2604.00680`；本地归档：`E:\Code\MAS\literature\05_localized_control_observability\2026_jaiswal_berger_tomar_partial_state_estimation.pdf`。[^15]

两篇论文都提供很好的“功能可观测性/模态可达性”数学语言，但均不与“一个静态 noisy WLS snapshot、节点只求自己的分块、有限 (r)-hop 通信、复制中央 WLS”撞车。

## 8. 年表式判定

| 命题 | 截止 2026 的证据判定 |
|---|---|
| 每个 (A_i) 满列秩是必要的 | **SOLVED：否。** 2015 已用全局满列秩替代；目标功能只需 (P_i\ker H=0)。置信度 0.99。 |
| 无环是估计存在的必要条件 | **SOLVED：否。** 一般环图可 Richardson 渐近解，也可生成树泛洪有限解。置信度 0.99。 |
| 标量 pairwise 含环图上，谱系小消息算法能否精确 | **SOLVED（受限模型）：能渐近精确。** 2020 Theorem 2。置信度 0.97。 |
| 一般向量含环图上，同一 BP/cavity 算法的简单 N&S | **PARTIAL。** 2021 有算法特定全局谱 iff；2022 GBDD 只是易检充分条件；准确度与中央正确性仍需分开。置信度 0.94。 |
| 一般含环图、固定局部消息宽度、有限轮、精确 nodewise WLS 的存在 iff | **OPEN-LOOKING。** 应以目标算子局部支撑、separator/treewidth 和消息预算共同刻画；直系论文未给。置信度 0.92。 |
| 同时允许失败、((\epsilon,\delta)) 近似、能量/轮数/内存多目标的最宽 iff/Pareto 边界 | **OPEN-LOOKING。** 找到的是分散的充分条件和个别下界，不是统一结论。置信度 0.90。 |
| 2022 / 2023 / 2026 已经彻底解决本题 | **否。** 2022 充分条件；2023 仍树；2026 三篇换了动态/common-state 目标或只给 comparable 精度。置信度 0.95（其中付费全文细节置信度较低）。 |

## 9. 最终判决

**总体判决：OPEN-LOOKING（置信度 0.92）。**

“无环 + 每个本地满秩”这组强假设早已分别被否定为必要条件，因此若论文只主张放宽这两条会撞车。尚未在直系谱系中发现的是：在清楚固定的 LOCAL/CONGEST 型资源模型下，对 nodewise 静态 WLS 的**目标依赖支撑 + 有限半径 + 消息宽度/内存 + 容错**给出统一必要充分条件，并进一步求通信能量、轮数和近似误差的 Pareto 最优算法。最稳妥的原创切入点不是再证明一种“含环也收敛”的充分条件，而是先以 (T_i=P_iK_R) 的局部 factorization/不可区分性给信息论 iff，再叠加 separator/treewidth 或比特通信约束；这样才能把“任何算法存在”与“某个 BP 更新收敛”彻底分开。

## 10. 主来源与公开链接

[^1]: 本地种子稿：`E:\Code\MAS\[ACC_11]A Distributed Algorithm for State Estimation with Application to Smart Grid.pdf`。PDF 创建日期 2012-09-10；未找到该题名的 2011 正式出版记录。

[^2]: X. Tai, Z. Lin, M. Fu, Y. Sun, “A New Distributed State Estimation Technique for Power Networks,” *2013 American Control Conference*, 3338--3343. DOI: [10.1109/ACC.2013.6580347](https://doi.org/10.1109/ACC.2013.6580347). 作者公开 PDF: [ACC_2013_2.pdf](https://www.eng.newcastle.edu.au/~mf140/home/Papers/ACC_2013_2.pdf). 本地：`E:\Code\MAS\literature\01_direct_lineage\2013_tai_lin_fu_sun_distributed_state_estimation_acc.pdf`。

[^3]: F. Pasqualetti, R. Carli, F. Bullo, “Distributed Estimation via Iterative Projections with Application to Power Network Monitoring,” *Automatica* 48(5), 747--758, 2012. DOI: [10.1016/j.automatica.2012.02.025](https://doi.org/10.1016/j.automatica.2012.02.025). [arXiv:1103.0579](https://arxiv.org/abs/1103.0579). 本地：`E:\Code\MAS\literature\01_direct_lineage\2012_pasqualetti_carli_bullo_iterative_projections.pdf`。

[^4]: D. E. Marelli, M. Fu, “Distributed Weighted Least-Squares Estimation with Fast Convergence for Large-Scale Systems,” *Automatica* 51, 27--39, 2015. DOI: [10.1016/j.automatica.2014.10.077](https://doi.org/10.1016/j.automatica.2014.10.077). [PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4308017/); [作者 PDF](https://www.eng.newcastle.edu.au/~mf140/home/Papers/Automatica2014_2.pdf). 本地：`E:\Code\MAS\literature\01_direct_lineage\2015_marelli_fu_distributed_wls_automatica.pdf`。

[^5]: D. E. Marelli, B. Ninness, M. Fu, “Distributed Weighted Least-Squares Estimation for Power Networks,” *IFAC-PapersOnLine* 48(28), 562--567, 2015. DOI: [10.1016/j.ifacol.2015.12.188](https://doi.org/10.1016/j.ifacol.2015.12.188). [作者 PDF](https://www.eng.newcastle.edu.au/~mf140/home/Papers/SYSID_2015.pdf)。

[^6]: T. Sui, D. E. Marelli, M. Fu, R. Lu, “Accuracy Analysis for Distributed Weighted Least-Squares Estimation in Finite Steps and Loopy Networks,” *Automatica* 97, 82--91, 2018. DOI: [10.1016/j.automatica.2018.07.016](https://doi.org/10.1016/j.automatica.2018.07.016). [arXiv:1806.09104](https://arxiv.org/abs/1806.09104); [作者 PDF](https://www.eng.newcastle.edu.au/~mf140/home/Papers/Automatica_2018.pdf). 本地：`E:\Code\MAS\literature\01_direct_lineage\2018_sui_marelli_fu_lu_loopy_wls_accuracy.pdf`。

[^7]: Z. Wu, M. Fu, Y. Xu, R. Lu, “A Distributed Kalman Filtering Algorithm with Fast Finite-Time Convergence for Sensor Networks,” *Automatica* 95, 63--72, 2018. DOI: [10.1016/j.automatica.2018.05.012](https://doi.org/10.1016/j.automatica.2018.05.012). [作者 PDF](https://www.eng.newcastle.edu.au/~mf140/home/Papers/Automatica_2018_3.pdf). 本地：`E:\Code\MAS\literature\01_direct_lineage\2018_wu_fu_xu_lu_distributed_kalman_loop_removal.pdf`。

[^8]: Q. Yang, Z. Zhang, M. Fu, “Distributed Weighted Least-Squares Estimation for Networked Systems with Edge Measurements,” *Automatica* 120, 109091, 2020. DOI: [10.1016/j.automatica.2020.109091](https://doi.org/10.1016/j.automatica.2020.109091). [arXiv:2002.11221](https://arxiv.org/abs/2002.11221). 本地：`E:\Code\MAS\literature\01_direct_lineage\2020_yang_zhang_fu_wls_edge_measurements.pdf`。

[^9]: D. Marelli, T. Sui, M. Fu, X. Sun, “Convergence and Accuracy Analysis for a Distributed Static State Estimator Based on Gaussian Belief Propagation,” *IEEE Transactions on Automatic Control* 66(10), 4785--4791, 2021. DOI: [10.1109/TAC.2020.3037454](https://doi.org/10.1109/TAC.2020.3037454). [arXiv:2004.01969](https://arxiv.org/abs/2004.01969). 本地：`E:\Code\MAS\literature\01_direct_lineage\2020_marelli_et_al_gaussian_bp_static_estimator.pdf`。

[^10]: Q. Yang, Z. Zhang, M. Fu, Q. Cai, “Asymptotic Convergence of a Distributed Weighted Least Squares Algorithm for Networked Systems with Vector Node Variables,” *Systems & Control Letters* 165, 105265, 2022. DOI/出版商正文页: [10.1016/j.sysconle.2022.105265](https://doi.org/10.1016/j.sysconle.2022.105265). 截止 2026-09-13 未定位到合法 OA PDF。

[^11]: M. Zhu, R. Wang, X.-M. Sun, “The Optimal Distributed Weighted Least-Squares Estimation in Finite Steps for Networked Systems,” *IEEE Transactions on Circuits and Systems II: Express Briefs* 70(3), 1069--1073, 2023. DOI: [10.1109/TCSII.2022.3215514](https://doi.org/10.1109/TCSII.2022.3215514). [IEEE 记录 9925114](https://ieeexplore.ieee.org/document/9925114). 截止 2026-09-13 未定位到合法 OA PDF。

[^12]: Q. Cai, Z. Zhang, M. Fu, “Distributed Computations for Large-Scale Networked Systems Using Belief Propagation,” *Journal of Automation and Intelligence* 2(2), 61--69, 2023. DOI/开放出版商页面: [10.1016/j.jai.2023.06.003](https://doi.org/10.1016/j.jai.2023.06.003).

[^13]: M. Zhu, Y. Wang, T. Sui, Y. Wu, “Finite-Time Distributed State Estimation for Networked Systems with Incomplete Measurements,” *Journal of the Franklin Institute* 363(8), 108674, 2026. DOI: [10.1016/j.jfranklin.2026.108674](https://doi.org/10.1016/j.jfranklin.2026.108674). [ScienceDirect 正式页面](https://www.sciencedirect.com/science/article/pii/S0016003226002747). 截止 2026-09-13 PDF 为机构/购买访问，未定位到合法 OA 全文。

[^14]: G. Fattore, M. E. Valcher, R. Gao, G.-H. Yang, “Distributed State Estimation of Discrete-Time LTI Systems via Jordan Canonical Representation,” arXiv:2603.10656v1, 2026（ECC 2026 会议稿的扩展版）. [摘要](https://arxiv.org/abs/2603.10656); [公开 PDF](https://arxiv.org/pdf/2603.10656); [arXiv DOI](https://doi.org/10.48550/arXiv.2603.10656). 本地：`E:\Code\MAS\literature\05_localized_control_observability\2026_fattore_et_al_discrete_lti_jordan_observer.pdf`。

[^15]: J. Jaiswal, T. Berger, N. K. Tomar, “Distributed Partial State Estimation for Linear State-Space Systems,” arXiv:2604.00680v1, 2026（稿件页眉为 “Preprint submitted to Automatica”）. [摘要](https://arxiv.org/abs/2604.00680); [公开 PDF](https://arxiv.org/pdf/2604.00680); [arXiv DOI](https://doi.org/10.48550/arXiv.2604.00680). 本地：`E:\Code\MAS\literature\05_localized_control_observability\2026_jaiswal_berger_tomar_partial_state_estimation.pdf`。

[^16]: J. Du, S. Ma, Y.-C. Wu, S. Kar, J. M. F. Moura, “Convergence Analysis of Distributed Inference with Vector-Valued Gaussian Belief Propagation,” *Journal of Machine Learning Research* 18(172), 1--38, 2018. [JMLR 正式开放页面/PDF](https://www.jmlr.org/papers/v18/16-556.html). 本地：`E:\Code\MAS\literature\02_graphical_models_local_inference\2018_du_et_al_vector_gaussian_bp_convergence.pdf`。

[^17]: D. Živojević, M. Delalić, D. Raca, D. Vukobratović, V. Švenda, “Distributed Weighted Least-Squares and Gaussian Belief Propagation: An Integrated Approach,” *IEEE SmartGridComm 2021*. DOI: [10.1109/SmartGridComm51999.2021.9632296](https://doi.org/10.1109/SmartGridComm51999.2021.9632296). [TechRxiv DOI](https://doi.org/10.36227/techrxiv.14781489).

[^18]: M. Zhu, R. Wang, X. Miao, T. Sui, “Accuracy Analysis for Distributed Dynamic State Estimation in Large-Scale Systems with a Cyclic Network Graph,” *Science China Information Sciences* 66, 190206, 2023. DOI: [10.1007/s11432-022-3846-1](https://doi.org/10.1007/s11432-022-3846-1).
