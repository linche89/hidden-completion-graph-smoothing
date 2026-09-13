# 局部 Gaussian/WLS 状态估计：定理级原创性与撞车审计

**核验截止：2026-09-13。**　本报告是独立的数学原创性审计，不是文献数量统计。对象限定为：稀疏线性/SPD 正规方程、有限轮局部通信、节点输出，以及固定模型或只知局部模型时的 minimax 误差。动态 observer 的 detectability、某一迭代算法的收敛条件、Bayes 平均 MSE 和这里的“任何算法是否存在”不是同一个量词问题。

> **后续专项结论。** 本报告收口时把一般块 exact SDP 标为 PENDING；随后完成的逐式审计已将其精确改写为单个 norm-bounded full block，并确认被非严格 Petersen lemma 直接包含。因此该候选作为原创主定理的最终判决也是 NO-GO。详见 [sdp_exact_novelty_audit.md](./sdp_exact_novelty_audit.md)。下文保留 PENDING 字样，用于记录审计在该时间点的证据链，而不是当前最终状态。

## 0. 结论先行

截至本次检索，没有找到一篇论文同时给出

> common local rule + arbitrary global completions + 对单位球输入的鲁棒误差 + 至少 \(1-\delta\) 节点 + 有限半径

这一整套 Gaussian/WLS 定理。然而，“没有整篇同款论文”并不等于当前候选定理有足够原创性。逐项拆开后：

1. fixed-model row-tail iff 是经典 information-based complexity / optimal recovery 的直接特例，**已撞车**；
2. completion 下的逐局部视图 Chebyshev 半径是把 completion 并入 problem element 后的同一经典定理，**直接组合**；
3. Schur/Dirichlet/oracle 精确分解是 block inverse、离散 DtN/transparent boundary 和 geometric resolvent identity 的改写；其范数 sandwich 是 restriction 加 submultiplicativity，**直接组合**；
4. 谱窗常数可锐化到 \((\kappa+1)/(2\sqrt\kappa)\)，但这是 Kantorovich–antieigenvalue 不等式加一个投影，**直接组合，不是新定理**；
5. \((1-\delta)\)-分位数只是逐点序关系取 quantile；若没有新的 pasting theorem，它是**定义级重述**；
6. grounded Laplacian 的 occupation、首次出球后的 remaining lifetime、Neumann survival tail 都是 Green 函数与强 Markov 性的直接结果；两独立游走的 intersection identity 也已有同式先例，**已知/直接组合**；
7. 任意 completion 的 \(\eta\)-Schur 与 pendant-tentacle no-go 似乎未见以 LOCAL-WLS 量词单独发表，但证明是一行 block inverse/KCL；它最多是**可能新的问题表述**，不是足以单独承载论文的深结论；
8. “rounds–messages–flops–memory–energy”目前只有成本核算，没有同一模型下的 matching lower/upper bound，因而还不是 tradeoff theorem；
9. “只知 exterior Schur 谱窗时的 robust boundary correction”在应用形式上仍有空隙，且 robust approximate inverse、robust linear estimator 和 optimized Schwarz 已覆盖大框架；但本项目新推得了一个重要升级：对**实矩阵、单个自伴 Loewner 谱窗、operator norm 和 locality-constrained \([H,0]\)**，一般 block 问题可由单约束 lossless S-lemma **exact 化为有限规模 SDP**，并非只能写成 semi-infinite convex program。这个结果在形式上已达到主定理门槛，但与 full-block S-procedure/robust control 的逐式原创性仍待专项查重。

**修订总判决：旧候选包若原样打包为“原创数学定理”，仍是 NO-GO（高置信度）；新得到的 block robust-boundary exact SDP 则单列为“形式上达到主定理门槛、原创性 PENDING”。** 它若通过 full-block S-procedure/robust-control 专项查重，可成为真正主定理；若被证明只是现成 lossless full-block theorem 的换元特例，则仍回到 NO-GO。其他可继续的创新核是 completion-pasting 的非平凡结构刻画、绝大多数根的 capacity/uniform-integrability iff、或同一通信模型下 matching 的资源 Pareto 定理。不能用整体 scope 很大来代替其中至少一个强定理。

## 1. 审计口径与证据等级

### 1.1 四道量词防火墙

| 轴 | 弱问题 | 强问题 | 不可偷换之处 |
|---|---|---|---|
| 模型知识 | 固定全局 \(\theta\)，系数可离线预装 | 节点只见 radius-\(r\) marked view，同一规则要跨所有 completion | 每个模型各有局部系数不推出一个 common rule |
| 输入风险 | 给定分布下 MSE | \(\sup_{\|z\|\le1}\) adversarial operator error | Gaussian screening/MSE 小不等于 operator row-tail 小 |
| 失败语义 | 每个输入可换一批坏节点 | 存在固定 \(1-\delta\) 好节点集、每个节点对所有输入正确 | row-norm quantile 对应后者，不对应前者 |
| 算法类 | 任意选择函数 | 可测/可计算、线性、有限精度或有限消息 | Chebyshev center 的集合论存在性不自动给低复杂度协议 |

尤其要严格区分：

- **算法收敛 iff**：如 walk-summability 保证某个 GaBP/Neumann 迭代；
- **问题可解 iff**：允许所有合法局部算法后仍能否达到误差 \(\varepsilon\)。

本报告只在证据确实跨过第二道量词时才称“问题级 iff”。

### 1.2 证据等级

- **A（全文）**：核对了公开全文中的定理、公式与页码；
- **B（正式摘要/出版页）**：核对官方摘要与书目信息，但未把它用作细粒度公式证据；
- **C（检索阴性）**：多术语检索未发现直接先例；只能支持“本次未找到”，不能证明不存在。

## 2. 总表：候选结论与最近先例的精确包含关系

| 候选核 | 最近的已知结果 | 精确包含关系 | 判定 | 置信度 |
|---|---|---|---|---|
| fixed-model row-tail iff | Werschulz 的 radius of information | 候选是线性解算子、单位球、坐标投影信息的特例 | **已撞车** | 高，A |
| common local rule across completions | radius of information + LOCAL indistinguishability | 把 \((\theta,z)\) 当作 problem element、把 marked view 与本地数据当作 information map 即得 | **直接组合** | 高，A/B |
| common **linear** rule \(\inf_\ell\sup_\theta\|T_\theta-\ell S\|\) | robust linear approximation / robust estimator | 候选是结构不确定线性算子族上的共同近似 | **框架已撞车；专门闭式未必** | 高/中，A |
| Schur exact factorization | block inverse；geometric resolvent；DtN/transparent BC | 候选等式就是大域 resolvent 通过边界耦合因子化 | **已知** | 高，A |
| \(g\le d\le\sqrt{1+\|F\|^2}g\) | 上述等式 + norm restriction/submultiplicativity | 两个不等式不需新工具 | **直接组合** | 高，A |
| 谱窗 sharp \(\gamma\) | Kantorovich/antieigenvalue | \(\|F\|\le(\kappa-1)/(2\sqrt\kappa)\) 是最大 turning angle 的投影推论 | **直接组合** | 高，A/B |
| 节点 \(1-\delta\) quantile | order statistic 定义 | pointwise inequality 取 quantile；固定模型 success fraction 是计数恒等式 | **定义级/直接组合** | 高 |
| killed-walk occupation/remaining life | Green potential + strong Markov | 三个误差量分别是标准 additive reward 在不同 stopping time 的期望 | **已知/直接组合** | 高，A |
| \(\ell_2\) Green tail = two-walk intersections | Ben-Hamou–Oliveira–Peres Eq. (3.1) | 正则/对称、截断或 killed 版本为同一展开 | **已撞车** | 高，A |
| arbitrary-completion \(\eta\)-Schur no-go | block inverse instability | 一行构造；未见相同 LOCAL-WLS 叙述 | **可能新表述，数学弱** | 中，A/C |
| pendant tentacle no-go | KCL/Green reciprocity/Kron reduction | pendant 零注入子网对端口等效、路径上转移 Green 系数恒定 | **folklore 组合** | 中高，A/C |
| robust optimal boundary correction | robust approximate inverse；optimized Schwarz/DtN；lossless S-lemma/full-block machinery | 思想与技术框架已知；新两列损失在实单 Loewner 谱窗下有 exact finite SDP，尚未在先例中逐式定位 | **主定理形态已达标；原创性 PENDING** | 数学正确性高；新颖性中低，A/C |
| local compute/communication tradeoff | CONGEST Laplacian solver matching rounds；distributed optimization oracle bounds | 已知 theorem 的误差、输出与带宽模型不同；当前候选仅 accounting | **证据不足/open-looking** | 高（“尚未证明”） |

## 3. Fixed-model row-tail iff：经典 optimal recovery 的特例

令 \(T_i:\mathcal Z\to\mathcal Y_i\) 是固定模型下节点 \(i\) 的目标算子，\(S\) 提取 radius-\(r\) 内可见输入。若 \(S\) 是正交坐标抽取，\(P=S^*S\)、\(Q=I-P\)，候选结论为

\[
\inf_f\sup_{\|z\|\le1}\|T_i z-f(Sz)\|=\|T_iQ\|. \tag{3.1}
\]

Werschulz 在 information-based complexity 的标准表述中定义

\[
r(N)=\sup_{y\in N(F)}\operatorname{rad} S\bigl(N^{-1}(y)\cap F\bigr)
=\inf_{\phi\text{ using }N}e(\phi,N), \tag{3.2}
\]

即所有共享同一信息值的 problem elements 的输出集合之 Chebyshev 半径，恰等于使用该信息的任何算法的最优最坏误差；公开技术报告 p.3 明列两式及 central algorithm。[^1] 取 problem class 为 Hilbert 单位球、solution operator 为 \(T_i\)、information 为坐标投影 \(S\)，零信息纤维是 \(\ker S\cap B\)，中心对称性给出半径 \(\|T_iQ\|\)，而 \(f(Sz)=T_iPz\) 达到它。

因此包含关系不是“相似”，而是

\[
\boxed{\text{fixed-model row-tail theorem}}
\subset
\boxed{\text{radius-of-information theorem}}.
\]

有限维下最优值取到；若模型类/纤维非闭或非紧，泛化为 completion 时应保留 \(\inf\)，或只声称任意 \(\eta>0\) 可达 \(r+\eta\)。可测选择、数值可计算性及通信编码也不由 (3.2) 自动给出。

**判定：已撞车（高置信度，A）。** 可以作为本文的基础 lemma 或把 WLS 译成 optimal recovery 语言，但不能声称新的“最宽 iff”。

## 4. Common local rule across completions：量词重要，但数学对象仍是 radius of information

### 4.1 任意非线性规则

把实例写成 \(w=(\theta,z)\)，令信息映射

\[
N_r(w)=\bigl(\text{rooted radius-}r\text{ marked model view},\;\text{visible data}\bigr),
\]

solution map 为 \(S_0(w)=T_i(\theta)z\)。任何 deterministic \(r\)-round LOCAL 算法对同一 \(N_r(w)\) 必须给相同输出；这是 LOCAL 模型的 view indistinguishability 原理，而不是 WLS 特性。Naor–Stockmeyer 奠定了常数半径局部计算模型；同半径同 view 的输出相同是模型定义的直接后果。[^2]

于是最优任意规则仍是

\[
\sup_v \operatorname{rad}\{T_i(\theta)z:N_r(\theta,z)=v,\ \|z\|\le1\}. \tag{4.1}
\]

这就是 (3.2) 对 compound problem element \(w=(\theta,z)\) 的应用。因此“逐 view 的 Chebyshev radius 是 common local algorithm 存在的 iff”在量词上是正确的，但原创性上是 **IBC 的直接实例化**。若每个 fibre 非空、闭、有界且有限维，则中心存在；否则“值 \(\le\varepsilon\)”与“存在恰好 \(\varepsilon\)-规则”的边界仍需 attainment/measurable-selection 条件。

### 4.2 限制为共同线性系数

共同线性规则的精确对象

\[
a_r^{\rm lin}(v)=\inf_\ell\sup_{\theta:\,\mathcal V_r(\theta)=v}
\|T_i(\theta)-\ell S\| \tag{4.2}
\]

是 robust linear approximation。若按可见/不可见输入分块 \(T_i(\theta)=[A_\theta,B_\theta]\)，则

\[
a_r^{\rm lin}(v)\ge
\max\left\{
\sup_\theta\|B_\theta\|,
\frac12\sup_{\theta,\theta'}\|A_\theta-A_{\theta'}\|
\right\}. \tag{4.3}
\]

第一项来自把输入限制在不可见块，第二项只是 triangle inequality。若 \(A_\theta\equiv A\)，取 \(\ell=A\) 又给上界，故

\[
a_r^{\rm lin}(v)=\sup_\theta\|B_\theta\|. \tag{4.4}
\]

这组式子很干净，但没有超出 robust approximation 的基本几何。El Ghaoui 已直接研究“对结构不确定矩阵族选择共同 approximate inverse，使最大 spectral-norm inverse error 最小”，Theorem 6.2（文章 p.184）给 structured 情形的 SDP 上界；unstructured 情形 Eq. (6.5) 给解析最优。[^3] Eldar–Ben-Tal–Nemirovski 也把不确定线性模型下的共同 minimax linear estimator 化为 SDP。[^4] 它们没有加入 graph view，但已经覆盖“跨模型选择一个共同线性算子”的核心优化范式。

**判定：任意规则的精确半径是直接组合；线性版的大框架已撞车（高置信度）。** 可能的贡献只能来自 completion class 的图结构使 (4.2) 可被局部、显式、低复杂度地计算，而不能来自重新定义 (4.2)。

## 5. Schur/Dirichlet 与 oracle row-tail：等式已知，sandwich 是一步推论

### 5.1 精确因子分解

对球 \(B\ni i\) 与外部 \(O\)，写

\[
J=\begin{bmatrix}A&E\\E^*&D\end{bmatrix}\succ0,
\qquad T=D-E^*A^{-1}E,
\]

并定义

\[
Y=P_iJ^{-1}P_O=-P_iA^{-1}ET^{-1},
\qquad F=E^*A^{-1}.
\]

直接 block inversion 给

\[
P_iJ^{-1}-P_iA^{-1}P_B=Y[-F,I]. \tag{5.1}
\]

这不是新的 Schur identity。离散 PDE/网络中，消去外部自由度后的 Schur 补就是 DtN/Steklov–Poincaré 边界算子；已知精确透明边界即可在截断域复现大域解。Gander–Zhang 的综述明确回顾：transparent boundary 产生 finite-step optimal Schwarz，而 approximate DtN 导致 optimized Schwarz。[^5] Kirsch 的 Theorem 5.20（p.40，Eq. (5.53)）把大域 resolvent 的内外矩阵元精确写成小域 resolvent、边界耦合与大域 resolvent的乘积，正是 (5.1) 的 operator/resolvent 近亲。[^6]

Gander–Jakabčin–Outrata 更直接把 domain truncation、Schur complement 与 artificial boundary 放在同一框架：Definition 2.1（文章 p.323）定义 Schur；Theorem 3.11（pp.330–331）证明有限 Dirichlet truncation 的 Schur symbol 是无限透明 Schur symbol 在无穷远的 ([i,i])-Padé 逼近。[^7]

### 5.2 oracle–Dirichlet sandwich

对 primitive RHS \(b\) 的单位 \(\ell_2\) 球，固定模型只见 \(b_B\) 的 oracle 最优误差是

\[
g_i(B)=\|P_iJ^{-1}P_O\|=\|Y\|.
\]

principal/zero-Dirichlet 局部解的误差是

\[
d_i(B)=\|Y[-F,I]\|.
\]

故

\[
g_i(B)\le d_i(B)
\le \sqrt{1+\|F\|^2}\,g_i(B). \tag{5.2}
\]

下界仅把 residual operator 限制到 outside-input 子空间；上界仅用

\[
\|Y[-F,I]\|\le\|Y\|\,\|[-F,I]\|,
\qquad \|[-F,I]\|=\sqrt{1+\|F\|^2}.
\]

所以即使没有找到文献逐字写出“oracle row-tail vs Dirichlet”这句话，(5.2) 仍是**标准等式的一步范数推论**，不是独立原创定理。对标量 \(B/O\) 分块可达到乘法因子，因此不能靠把常数称作 sharp 来抬高原创性。

### 5.3 sharp \(\gamma\) 的专项核查

若

\[
mI\preceq J\preceq MI,
\qquad \kappa=M/m,
\]

则最紧的只依赖谱窗的 block bound 是

\[
\|E^*A^{-1}\|
\le \frac{M-m}{2\sqrt{mM}}
=\frac{\kappa-1}{2\sqrt\kappa}, \tag{5.3}
\]

从而

\[
\gamma\le
\sqrt{1+\|F\|^2}
\le \frac{M+m}{2\sqrt{mM}}
=\frac{\kappa+1}{2\sqrt\kappa}. \tag{5.4}
\]

证明是已知 Kantorovich/antieigenvalue 几何的一个投影：对 SPD \(J\)，向量与其像的最小 cosine 为

\[
\inf_{x\ne0}\frac{\langle x,Jx\rangle}{\|x\|\,\|Jx\|}
=\frac{2\sqrt{mM}}{m+M};
\]

把 \(J[u;0]=[Au,E^*u]\) 投影到内/外子空间，即得 tangent bound (5.3)。二维仅含端点特征值的 rotation 例达到等号。Gustafson 1968 引入 operator angle；1999 论文专门给出 Kantorovich–Wielandt 的几何/antieigenvalue 解释。[^8]

因此此前较粗的 \((\kappa-1)/2\) 可被修正，但 **sharp \(\gamma\) 不是原创核**：

\[
\boxed{\text{Kantorovich sharp angle}}
+\boxed{\text{subspace projection}}
+\boxed{\text{submultiplicativity}}
\Longrightarrow (5.3)\text{--}(5.4).
\]

### 5.4 适用边界

(5.1)–(5.2) 精确对应的是节点直接获得 primitive normal-equation RHS \(b\)。原始 WLS 若写 \(x=J^{-1}Gz\)，则从 observation \(z\) 回到该定理还需：\(Gz\) 能否在同一球内形成、\(G\) 是否跨边界、输入范数经 \(G\) 如何变化。否则不能把 primitive-RHS sandwich 无条件称为 raw-measurement WLS theorem。

**判定：Schur 等式已知；sandwich 与 sharp \(\gamma\) 均为直接组合（高置信度，A）。**

## 6. 新补查：只知 exterior Schur 谱窗时的 robust boundary correction

这是本轮最值得保留、也最需要降调的候选核。

### 6.1 问题的准确形式

固定本地 \(A,E\)，令

\[
C=P_iA^{-1}E,\qquad F=E^*A^{-1},\qquad X=T^{-1},
\]

只知道

\[
aI\preceq X\preceq bI,
\qquad a=1/M,\quad b=1/m.
\]

允许在 zero-Dirichlet 本地系数 \(P_iA^{-1}\) 上加一个只依赖 \(A,E,a,b\) 的共同 correction \(H\)。对 local/exterior primitive RHS 的共同误差为

\[
\mathcal R(H)=
\sup_{aI\preceq X\preceq bI}
\bigl\|[CXF-H,\,-CX]\bigr\|_2. \tag{6.1}
\]

这里第二列是完全看不见的 exterior RHS，不能从目标中删掉。\(H=0\) 才是 zero-Dirichlet；精确 \(X\) 已知时可取 \(H=CXF\)，但这违反 common-completion 量词。

### 6.2 已有工作覆盖到哪里

存在三条高度相关的成熟谱系：

1. **uncertain approximate inverse**：El Ghaoui 2002 已把“选一个共同矩阵近似所有不确定 inverse、最小化最大 operator-norm error”定义并作 SDP/解析处理；Theorem 6.2 和 Eq. (6.5) 是直接撞车警报。[^3]
2. **robust linear estimator**：不确定 observation/covariance 下共同 minimax estimator 的 SDP 已知。[^4]
3. **optimized Schwarz / artificial boundary**：Bennequin–Gander–Halpern 研究 homographic minimax approximation；Theorems 2.2、2.6（文章 pp.186、189）给紧集上最优传输条件的存在与唯一性。[^9] Gander–Jakabčin–Outrata §4（pp.333–335）数值最小化 Robin truncation 的 \(\ell_\infty\) symbol error，并展示 equioscillation。[^7]

但 (6.1) 与标准 optimized Schwarz 目标并不相等：后者常最小化迭代反射/收敛因子

\[
\sup_\lambda\left|\frac{p-f(\lambda)}{p+f(\lambda)}\right|
\]

或直接近似 DtN symbol；(6.1) 则是一次性 nodewise solution-map error，并同时惩罚 unseen exterior RHS。另一方面，任意矩阵 \(H\) 也未必可实现为对称正定 Robin/DtN 边界算子；若论文使用“optimal boundary condition”措辞，还必须证明 realization。

因此最准确的包含判断是：

\[
\text{(6.1)}\ \text{属于 robust approximate-operator 范式，}
\quad
\text{但不是已检到 Schwarz 损失的逐字特例。}
\]

### 6.3 标量 cut 的闭式（本轮推导，尚未见同式先例）

若 cut dimension 为一，则 \(X=u\in[a,b]\)。令 \(q=\|F\|_2^2\)，\(C\ne0\)。把 \(H\) 正交投影到 \(C\) 的列空间、再把其行方向投影到 \(F\)，都不会增大 residual；故最优解可写 \(H=hCF\)。于是

\[
\mathcal R(hCF)
=\|C\|\max_{u\in[a,b]}
\sqrt{q(u-h)^2+u^2}. \tag{6.2}
\]

被最大化的平方关于 \(u\) 凸，所以只需比较两个端点。令

\[
q_0=\frac{a+b}{b-a}.
\]

当 \(q=0\) 时 \(H_*=0\)、\(\mathcal R_*=\|C\|b\)。当 \(q>0\) 时，一个最优参数为

\[
h_*=
\begin{cases}
b, & q\le q_0,\\[2mm]
\dfrac{(a+b)(q+1)}{2q}, & q>q_0,
\end{cases} \tag{6.3}
\]

且

\[
\mathcal R_*=
\begin{cases}
\|C\|b, & q\le q_0,\\[2mm]
\|C\|\sqrt{
b^2+\dfrac{(q(b-a)-(a+b))^2}{4q}}, & q>q_0.
\end{cases} \tag{6.4}
\]

zero-Dirichlet 的风险为

\[
\mathcal R(0)=\|C\|b\sqrt{q+1}, \tag{6.5}
\]

故 \(C,F\ne0\) 时共同 correction 严格改善它；第一种区间甚至达到由 unseen outside column 强制的下界 \(\|C\|b\)。

式 (6.3)–(6.4) 经代数和数值网格复核，但尚不是经同行评审的结果。专项检索了 *robust optimal recovery with uncertain operator*、*approximate inverse of uncertain matrices*、*optimized Schwarz/Robin minimax*、*transparent/artificial boundary approximation*、*optimal domain truncation*，未找到同一两列损失下的同式闭解。这个阴性结论仅为 C 级证据。

### 6.4 更新：一般实 cut 可 exact 化为有限 SDP

进一步攻击推翻了“这里只能得到 semi-infinite convex program”的旧判断。令 \(C\in\mathbb R^{p\times n}\)、\(F\in\mathbb R^{n\times k}\)、\(H\in\mathbb R^{p\times k}\)，并置

\[
h=\frac{a+b}{2},\qquad d=\frac{b-a}{2}.
\]

关键几何恒等式是：对任意 \(v\in\mathbb R^n\)，

\[
\{Xv:aI\preceq X\preceq bI,\ X=X^T\}
=\{z:\|z-hv\|_2\le d\|v\|_2\}. \tag{6.6}
\]

“\(\subseteq\)”来自 \(X=hI+dK\)、\(K=K^T\)、\(\|K\|\le1\)；反向可取一个把 \(v/\|v\|\) 映到目标方向的对称 Householder 矩阵，再作尺度收缩。它是 self-adjoint contraction completion 的专门几何；Davis–Kahan–Weinberger 给出更一般的 norm-preserving block completion 理论。[^28]

定义

\[
D_0=
\begin{bmatrix}
-I_n&hC^T\\
hC&-ab\,CC^T
\end{bmatrix},\qquad
E_0=
\begin{bmatrix}
0&0\\0&I_p
\end{bmatrix},\qquad
L_H=
\begin{bmatrix}
F^T&-H^T\\
I_n&0
\end{bmatrix}.
\]

则在 \(a<b\)、\(C\ne0\) 时，

\[
\boxed{
\inf_H\mathcal R(H)^2
=
\min_{\tau\ge0,\ H,\ \lambda\ge0}\ \tau
\quad\text{s.t.}\quad
\begin{bmatrix}
-\tau E_0+\lambda D_0&L_H^T\\
L_H&-I_{k+n}
\end{bmatrix}\preceq0 .
} \tag{6.7}
\]

证明量词是精确的：谱范数先写成左测试向量 \(u\) 的 supremum，令 \(v=C^Tu\)、\(z=Xv\)，式 (6.6) 把全部 Loewner completion 精确压成一个欧氏球，即单个齐次二次约束 \(w^TD_0w\ge0\)，其中 \(w=[z;u]\)。因 \(a<b\)、\(C\ne0\) 给严格可行点，单约束 S-lemma 是 lossless；再作一次 Schur complement 即得 (6.7)。S-lemma 的严格可行单约束版本是经典结果。[^29] 退化情形 \(C=0\) 或 \(a=b\) 应单独处理。

因此：

- 一般实 block、单 Loewner 谱窗的值已经是多项式规模 exact SDP，不应再称“只能 semi-infinite”或无依据称 NP-hard；
- 对固定 \(H\)，最坏 \(X\) 可取 \(aI+(b-a)P\)，但不能只查 \(aI,bI\) 两个各向同性端点；优化 \(H\) 后仍存在严格端点 gap 的 \(2\times2\) 解析反例，详见本地攻击报告；[^30]
- 任意额外的仿射 \(H\)-约束可直接加入 (6.7)，而不破坏这一 exactness；
- 当前证明限于实对称情形。复 Hermitian 情形不能不经证明直接沿用 Householder 球像恒等式。

与 El Ghaoui 2002 的差别必须准确表述：其 Eq. (1.4) 已提出共同 approximate inverse 的 min–max operator-norm 问题；Theorem 6.2/Eq. (6.3) 对一般 structured LFR 给 SDP bound，而论文明确把 exact analytic claim 留给 unstructured full perturbation 的 Eq. (6.5)。[^3] 式 (6.7) 的额外结构是 row-selected \(CX\)、locality-constrained approximant \([H,0]\)、耦合损失 \([CXF-H,-CX]\)，以及 self-adjoint Loewner interval 对单测试向量的球像。故它**不由已核到的 Eq. (6.5) 直接推出**；但证明仍属于 full-block geometry + lossless S-procedure 的成熟机器，不能在专项逐式查重前宣称原创。

**修订判定：数学上 exact，且形式上达到主定理门槛；原创性 PENDING（数学正确性高，文献新颖性中低置信度）。** 标量 cut 现在只是 (6.7) 的解析特例。完整推导、端点反例与 SDP 尺寸见 [block robust boundary attack](./block_robust_boundary_attack.md)。

## 7. \(1-\delta\) 节点：quantile 本身不是原创 theorem

### 7.1 固定模型、固定好节点集

令

\[
e_{i,r}=\inf_{f_i}\sup_{\|z\|\le1}
\|T_i z-f_i(S_{i,r}z)\|.
\]

对节点概率测度 \(\mu\)，存在至少 \(1-\delta\) 的**固定节点集**逐节点满足 worst-case error \(\le\varepsilon\)，当且仅当

\[
Q_{1-\delta}^{\mu}(e_{i,r})\le\varepsilon. \tag{7.1}
\]

这是 order statistic 的定义。若有 pointwise \(d_i\le\Gamma g_i\)，则

\[
Q_{1-\delta}(d)\le\Gamma Q_{1-\delta}(g) \tag{7.2}
\]

也只是 quantile 对逐点序的单调性。因此“Schur sandwich 再取 quantile”是直接组合，不能称新的多数节点定理。

### 7.2 两种常被混淆的失败语义

(7.1) 对应

\[
\mu\{i:\sup_{\|z\|\le1}\mathrm{err}_i(z)\le\varepsilon\}\ge1-\delta,
\]

即好节点集在输入之前固定。它一般强于

\[
\sup_{\|z\|\le1}\mu\{i:\mathrm{err}_i(z)>\varepsilon\}\le\delta, \tag{7.3}
\]

后者允许坏节点随 \(z\) 改变。误差算子 \(I_n\) 已能区分二者：每行 operator norm 是 1，但单位 \(\ell_2\) 输入不可能在很多坐标同时大。因而任何“允许一定失败”的 theorem 都必须先固定采用 (7.1) 还是 (7.3)。

### 7.3 completion 与 quantile 不可随意交换

真正的 common-rule 风险是

\[
\inf_{\{\ell_v\}}\sup_\theta
Q_{1-\delta}^{\mu_\theta}
\bigl(e_i(\theta;\ell_{v_i})\bigr). \tag{7.4}
\]

逐 view 先取最坏 completion、再在节点上取 quantile，通常只是 (7.4) 的充分上界；最坏 completion 可能随 view 改变，而这些局部最坏片段未必能同时拼成一个合法全局 \(\theta\)。只有 completion class 对 disjoint union/patching 有足够闭包时，才能交换“每个 view 的 sup”与“每个模型的 quantile”。

measured asymptotic expanders、measured Roe algebras确实已有“测度/大多数点”的成熟语言；例如 Li–Špakula–Zhang 的 Proposition 4.8、Corollary 4.21、Theorem 6.1 处理 measured expansion 与 Roe rigidity。[^10] 但它们不是 nodewise WLS row-quantile theorem。反过来，Špakula–Zhang Theorem 3.3 在 Property A 等条件下连接 quasi-local 与 finite-propagation/Roe approximation；Ozawa 的 Theorems A–B、Corollary C 又表明没有合适几何条件时 quasi-local 不必属于 uniform Roe algebra。[^11][^12] 这些结果提醒我们：不能从一个分位数定义无条件升级成图拓扑 iff。

**判定：quantile sandwich 是直接组合；固定模型 success-fraction 是 tautological characterization（高置信度）。** 真正可能新的只有 completion-pasting 的必要充分结构，或可检验的 measured spectral/capacity 条件。

## 8. Grounded Laplacian / M-matrix：概率表示几乎全部已知

### 8.1 三个误差量的 stopping-time 表示

令

\[
J=S(I-P),\qquad P\ge0,\quad \rho(P)<1,
\qquad N=(I-P)^{-1}=\sum_{t\ge0}P^t,
\]

把 substochastic \(P\) 增广一个 cemetery state，并令 \(\zeta\) 为首次被杀死的时刻。对 \(\|c\|_\infty\le1\)，正性使 row operator norm 等于行和。若 \(B=B_r(i)\)，令 \(\tau_B=\inf\{t:X_t\notin B\}\)，则：

\[
o_i(B)=\sum_{j\notin B}N_{ij}
=\mathbb E_i\sum_{t<\zeta}\mathbf1\{X_t\notin B\}, \tag{8.1}
\]

\[
d_i(B)=
\mathbb E_i[(\zeta-\tau_B)\mathbf1\{\tau_B<\zeta\}], \tag{8.2}
\]

\[
n_i(r)=\mathbb E_i[(\zeta-r-1)_+], \tag{8.3}
\]

分别对应 oracle outside tail、zero-Dirichlet 首次出球后的剩余 reward、以及 \(r\)-step Neumann truncation 的 survival tail。因此

\[
o_i(B)\le d_i(B)\le n_i(r). \tag{8.4}
\]

Lawler–Limic 的 Theorem 1.6.1 是 strong Markov property；§4.6 的 Lemma 4.6.1（book p.96）把 killed-domain Green 定义为出域前期望访问次数，Proposition 4.6.2（p.97）正是“总 Green = 出域前 Green + 出域后条件 Green”；§9.10（p.233）又明确写出 \(G=(I-Q)^{-1}\) 是严格 sub-Markov 链的 expected visits。[^13] 因而 (8.1)–(8.4) 是标准 Markov reward 分解在三个截断规则上的应用。

令 \(h_i(B)=\Pr_i(\tau_B<\zeta)\)。若所有可能 exit state \(x\) 的剩余寿命满足

\[
\sup_x\mathbb E_x\zeta\le L,
\]

则强 Markov 立即给

\[
h_i(B)\le o_i(B)\le d_i(B)\le Lh_i(B). \tag{8.5}
\]

这里 \(L\) 是 domination/uniform-integrability 型 envelope；它是合理且必要的稳定性语言，但 (8.5) 本身只是 conditional expectation，不是新“余寿命 criterion”。

### 8.2 \(\ell_2\) 行尾与两游走相交

两条从 \(i\) 出发的独立 killed walks \(X,Y\) 满足（对称情形；一般 reversible 情形需 stationary/degree 权重）

\[
\sum_jN_{ij}^2
=\mathbb E_i\sum_{s,t<\zeta}\mathbf1\{X_s=Y_t\}. \tag{8.6}
\]

Ben-Hamou–Oliveira–Peres 的 Eq. (3.1)（文章 p.380）已写出 truncated Green 平方和等于两独立随机游走 path intersections 的期望；正则图时就是无权版本。[^14] 把时间截断换成 killing、把求和限制到球外，不改变其乘法展开本质。

### 8.3 什么还可能新

Green/Poisson kernel、Dirichlet problem 与 finite-network potential theory 已有系统理论；Bendito–Carmona–Encinas 给出有限网络 Schrödinger operator 的 Green/Poisson 框架，Carmona 等 2024 又专门连接 symmetric M-matrix 与随机游走。[^15][^16]

因此可发表的新核不能只是给 (8.1)–(8.6) 换 WLS 名字。需要的是例如：

\[
\boxed{
\text{随机根可局部验证的 capacity/tightness 条件}
\iff
Q_{1-\delta}(o_i(r))\to0
}
\]

并证明 UI 假设的 sharpness、rare bottleneck 反例以及 completion 下的必要方向。Bordenave–Lelarge 的 local weak limit/resolvent 工作只控制谱测度/非实谱参量的 resolvent；靠近零点的 inverse 可被消失谱隙放大。[^17] Anantharaman–Sabri 也显式在 local weak tree limit 外加入 Green moment control。[^18] 这说明“local weak convergence alone”不够，却没有替我们给出上述 majority iff。

**判定：随机游走表示与 \(h\)–\(L\) sandwich 已知/直接组合；可检验的 majority capacity/UI iff 仍 open-looking。**

## 9. Arbitrary-completion tentacle/no-go：结论可信，原创深度有限

### 9.1 generic SPD 的 \(\eta\)-Schur 构造

固定可见 \(A,E\)，若 completion class 只要求 SPD，而不给 exterior stability，则可取

\[
D_\eta=E^*A^{-1}E+\eta I,
\qquad T_\eta=\eta I,
\qquad \eta\downarrow0. \tag{9.1}
\]

只要 \(P_iA^{-1}E\ne0\)，block inverse 中的 \(T_\eta^{-1}\) 使不可见系数和本地系数 ambiguity 同时发散。由此得到：固定局部 view 下，对所有 SPD completion 的 finite minimax error 通常为无穷。

这是正确且很有用的量词反例，但数学证明就是 Schur complement criterion 加 (5.1)。本次没有找到以“only-local-model WLS common rule no-go”原样陈述的正式 theorem；这只能支持“可能是新的表述”，不能使一行构造变成深原创。

### 9.2 grounded network 的 pendant tentacle

在可见边界端口接一条长度 \(N\)、内部无 grounding 的 pendant path。对从 tentacle 任一点注入的单位电流，全部电流必须经同一端口进入内部 grounded network；由 KCL 与 Green reciprocity，根 \(i\) 对每个 tentacle 节点的 transfer Green coefficient 相同。因此球外

\[
\ell_1\text{ tail}\asymp N,
\qquad
\ell_2\text{ tail}\asymp\sqrt N. \tag{9.2}
\]

Kron reduction 本来就是 Laplacian 的 Schur complement；Dörfler–Bullo Theorem III.8（文章 pp.157–159 附近）证明 boundary effective resistance 在 Kron reduction 下保持。[^19] pendant 零注入树对端口的等效行为则是基础电路/KCL。故 (9.2) 的 LOCAL-completion 叙述可能未被专门写过，但机制是 folklore。

该 no-go 依赖 completion 允许：规模任意增长、无统一 killing/grounding、无统一 \(mI\preceq J\)、并保持根的有限 view 不变。一旦模型类承诺统一 Schur 下界、统一余寿命或限制 exterior 体积，结论不再适用。

**判定：可能新的 LOCAL-WLS formulation（中置信度），但现有版本不足以当主数学贡献。** 真正的新定理应分类：哪些 completion closure 导致 minimax 无穷；哪些最小 stability promise 恰好使其有限。

## 10. 逆衰减、finite section、Gaussian screening 与本问题的边界

这些文献不是“漏掉的答案”，但会限制可声称的原创范围。

### 10.1 inverse decay / localized resolvent

- Demko–Moss–Smith Proposition 2.1、2.2 与 Theorem 2.4（pp.492–494）已从谱区间和 Chebyshev 逼近推出 banded SPD inverse 的指数衰减，且正定情形速率可达。[^20]
- Jaffard 1990 证明带指数/多项式 off-diagonal decay 的矩阵类在 \(\ell_2\) 逆下封闭，并明确把 inversion 称为“local” transformation。[^21]
- Gröchenig–Rzeszotnik–Strohmer Theorem 12（manuscript p.15）给 weighted decay algebra 中 finite-section 解的定量 tail error。[^22]
- Cheng–Jiang–Sun Theorem 5.3（manuscript p.14）是 inverse-closed Wiener lemma；Theorems 6.1–6.2（pp.15–16）连接全局 stability 与足够大局部 restrictions；Proposition 7.1（pp.17–18）给 localized finite solve 误差。[^23]

所以“uniform spectral gap/condition number + sparse finite propagation \(\Rightarrow\) inverse locality”已高度撞车。它们主要给结构性充分条件或条件化双向 stability，不给 arbitrary completion + \(1-\delta\) roots 的最宽 iff。

### 10.2 Gaussian screening

Stein 2011 的 Conjecture 1、Theorems 1–2 研究连续空间 Gaussian prediction 的 screening effect；这是相对 kriging risk 的渐近结论，不是图上 fixed-radius adversarial row norm。[^24] 因而可以作为统计动机，不能用来证明本项目的 operator iff，也不能声称“screening 从未研究过”。

### 10.3 operator quasi-locality

Property A / Roe 理论回答 whole-operator norm 下 quasi-local 是否可由 finite propagation 逼近；本项目的 nodewise row quantile 更弱。Ozawa 的反例又表明一般空间上 quasi-local 与 uniform Roe 不等价。故若想要无条件“远集相互作用趋零 iff finite-round approximation”，在 whole-operator 层面甚至是假的；多数节点版本必须另加 measure/capacity 结构。

## 11. Local computation / communication tradeoff：现在只有账本，没有 Pareto 定理

在 unlimited-message LOCAL 中，\(r\) 轮可泛洪整个 radius-\(r\) ball，但 ball size 可指数增长；把稠密局部系统显式收集后求解，朴素时间/空间约为 \(O(m_r^3)/O(m_r^2)\)。Neumann/多项式法则以一次 sparse matvec 对应一轮，存储小，但误差依赖谱隙/条件数。这些是正确的 resource accounting。

已有真正 matching 的通信下界必须在明确协议模型内陈述。例如 Anagnostides 等 Theorem 1 证明：在 CONGEST 与 Supported-CONGEST 上，固定图 \(G\) 的 Laplacian solver 即使只要求 \(\varepsilon\le1/2\)，仍需 \(\widetilde\Omega(\mathrm{SQ}(G))\) 轮；Theorem 2 给 Supported-CONGEST 中 \(n^{o(1)}\mathrm{SQ}(G)\log(1/\varepsilon)\) 上界。[^25] Andoni–Krauthgamer–Pogrow Theorems 1.1–1.2 则在 sublinear-query 模型给单坐标线性系统的上/下界。[^26] Scaman 等在 first-order distributed optimization oracle 中给 computation/communication time 的匹配界。[^27]

这些都不直接包含当前语义：nodewise output、只要求 \(1-\delta\) 节点、local-model completions、以及 row-wise adversarial error。不同模型的 lower bound 不能拼成同一 Pareto frontier。

若要把“能耗低、估计快、拓扑要求弱、计算低”变成数学贡献，至少要先固定

\[
(r,\ \text{bits/edge},\ \text{total bit-hop},\ \text{flops},\ \text{memory},\ \varepsilon,\ \delta)
\]

并在**同一协议与输入类**证明 matching lower/upper bound。没有 converse 的误差曲线不是“最优 tradeoff theorem”。

**判定：当前证据不足；目标 open-looking，但范围必须收窄。**

## 12. 剩余创新核：还必须多证明什么才不只是拼接

| 候选创新核 | 当前已有 | 必须新增的非平凡 lemma/theorem | 若做不到的判决 |
|---|---|---|---|
| robust common boundary correction | 一般实 cut 的 exact finite SDP (6.7)；标量闭式；优化后端点反例；robust approximate inverse/full-block 大框架已知 | 首要缺口已从“数学定理”转为“逐式原创性”：证明 (6.7) 不被既有 full-block S-procedure/DKW/robust performance theorem 直接包含；再给网络 completion 的物理 realization 或一个严格分离 corollary 可显著增强 | **形式上已过主定理门槛；若专项查重撞车则 NO-GO，若未撞车则主定理候选** |
| completion-pasting + \(1-\delta\) | Chebyshev fibre 与 quantile 定义；sup/quantile 不交换的观察 | 对具体 completion class 给 necessary-and-sufficient pasting/compactness 条件，使 per-view ambiguity 与 per-model quantile 恰好等价；或构造 sharp gap family 并定量刻画 | 没有结构 theorem：纯定义，NO-GO |
| majority capacity/UI | killed-walk identities；\(h\le d\le Lh\) | 用 radius-\(r\) 可见量、random-root law 或 Dirichlet capacity 给可检验 iff；证明 UI/tightness 恰好必要；rare bottleneck/tentacle 反例达到边界 | 只写 remaining lifetime：经典 Markov 推论，NO-GO |
| generic block-SPD completion | Schur transfer factorization与谱窗上界 | 允许 cancellation 的 transfer singular-value/capacity invariant；对 majority roots 给双向而非 block-norm majorant 的单向界 | 只有 submultiplicativity：NO-GO |
| arbitrary-completion no-go 分类 | \(T_\eta=\eta I\) 与 pendant tentacle | “minimax 有限 iff completion family 满足何种最小 Schur coercivity/tightness”完整分类；最好给必要性、充分性与可局部验证代理 | 两个反例只能作警示 lemma |
| resource Pareto | 轮数/ball-size/flops 账本；其他模型已有 lower bounds | 在同一 nodewise LOCAL/CONGEST 模型、同一 \((\varepsilon,\delta)\) 风险下给 matching rounds–bits–memory lower/upper；再接硬件能量模型 | 无 matching converse：不是优化定理 |

### 对两个“看似尖锐”点的最终裁决

1. **sharp \(\gamma\)**：确实比粗界好且可达，但来源是 antieigenvalue sharp angle + 投影 + submultiplicativity；它是漂亮的 sharp corollary，**不是独立原创核**。
2. **\(\delta\)-quantile**：若逐点已有界，取 quantile 是序保持；固定模型“至少 \(1-\delta\) 节点成功 iff quantile 小”是定义。除非加入 nontrivial pasting/capacity theorem，**不是原创核**。

### 强主定理的硬门槛

**一句判据：此前 package 明确 NO-GO；新式 (6.7) 已在同一实单谱窗模型中给出一般非交换 cut 的 exact robust minimax、有限可计算上界和 lossless converse，因而形式上跨过主定理门槛；现在的硬门槛变成“专项查重证明它不只是现成 full-block S-procedure 的直接实例”，否则仍按 NO-GO 处理。** 另外三条可独立跨门槛的路线仍是 completion-pasting 的 sharp iff、majority capacity/UI 的 sharp iff，或同一模型下 matching 的 rounds–bits–memory Pareto。

把这些方向并列成大 scope、把若干已知 lemma 串联、给标量特例、改善一个次乘法常数、或在不同通信模型间拼接上下界，都不跨过这道门槛。

## 13. 可安全声称与不应声称

### 可安全声称

- “我们把 optimal recovery 的 radius-of-information 专门化为 nodewise WLS row-tail，并显式区分 fixed model 与 common-completion rule。”
- “我们给出一个 WLS-specific Schur factorization，由此得到 oracle 与 principal solve 的常数因子比较。”
- “我们把三种局部截断误差统一翻译成 killed-walk stopping-time rewards。”
- “本次检索未发现同时覆盖 common local rule、arbitrary completions、\(1-\delta\) nodes 与 finite radius 的现成论文。”
- “标量 cut robust correction 的特定两列 minimax 闭式尚未找到同式先例；原创性仍待同行/更系统数据库核验。”
- “在实对称单 Loewner 谱窗下，两列 robust boundary 损失可由单约束 lossless S-lemma exact 化为有限 SDP；该数学结论已经内部逐式核验，但文献原创性仍在专项查重。”

### 不应声称

- “首次发现局部估计的充要条件是 inverse row tail / Chebyshev radius”；
- “Schur sandwich 或 sharp \(\gamma\) 是全新矩阵不等式”；
- “首次用随机游走 occupation/remaining lifetime 表示局部线性求解误差”；
- “取 \(1-\delta\) 分位数得到新的概率图论 iff”；
- “robust optimal boundary correction 这一思想首次提出”；
- “一个纯拓扑条件彻底终结任意 Gaussian/WLS 局部估计问题”。

## 14. 检索边界与最终判决

本次沿以下等价术语检索并核对主来源：optimal recovery/radius of information；LOCAL indistinguishability；robust approximate inverse/robust linear estimation；domain decomposition、DtN/Steklov–Poincaré、transparent/artificial boundary、optimized Schwarz/Robin；geometric resolvent、finite section、localized inverse、Jaffard algebra；Green/Poisson kernel、killed walk occupation、intersection local time、uniform integrability；Gaussian screening；Property A、quasi-local/uniform Roe；distributed Laplacian solver 与 communication lower bound。2025–2026 的相关新作也未出现与上述完整量词同款的 theorem。

阴性检索不是数学上的不存在证明；尤其是式 (6.7)，仍可能藏在 full-block S-procedure、robust performance、interval analysis、\(H_\infty\) model matching 或 matrix Chebyshev-center 文献的不同符号下。因此报告只给“本次未找到逐式先例”，不写“全球首次”。

最终分层判决：

- **SOLVED / 已知**：fixed-model optimal error 的 information-radius 表述；Schur/DtN/resolvent 因子化；killed Green occupation；给定 pointwise error 后的 quantile 计数。
- **PARTIAL / NOVELTY-PENDING**：completion no-go、WLS-specific oracle–Dirichlet comparison、标量 cut 闭式，以及数学上 exact 的一般实 block SDP (6.7)。前三者增量偏薄；(6.7) 的数学形态足够强，但尚未通过外部原创性审计。
- **OPEN-LOOKING**：completion-pasting 的可检验 iff；多数根 capacity/UI iff；同一 nodewise 协议下 matching 资源 Pareto；以及复 Hermitian/多独立不确定块是否仍能 lossless finite 化。
- **原创性总判决**：**旧 theorem package = NO-GO；新增 (6.7) = 主定理形态已达标、原创性 PENDING。专项查重若证明它被既有 full-block theorem 直接包含，则整体回到 NO-GO；若无直接包含，它是当前唯一足以支撑强主定理的候选。**

## Sources

[^1]: Arthur G. Werschulz, *An Overview of Information-Based Complexity*, Columbia Technical Report CUCS-022-02 (2002), p.3, radius formula and optimal error identity. [Official PDF](https://mice.cs.columbia.edu/getTechreport.php?format=pdf&techreportID=152). 另见其专著 chapter DOI [10.1093/oso/9780198535898.003.0004](https://doi.org/10.1093/oso/9780198535898.003.0004).

[^2]: Moni Naor, Larry Stockmeyer, *What Can Be Computed Locally?*, SIAM J. Comput. 24(6) (1995), 1259–1277. DOI [10.1137/S0097539793254571](https://doi.org/10.1137/S0097539793254571).

[^3]: Laurent El Ghaoui, *Inversion Error, Condition Number, and Approximate Inverses of Uncertain Matrices*, Linear Algebra Appl. 343–344 (2002), 171–193; Theorem 6.2 and Eq. (6.5), p.184. DOI [10.1016/S0024-3795(01)00273-7](https://doi.org/10.1016/S0024-3795(01)00273-7); [author PDF](https://people.eecs.berkeley.edu/~elghaoui/Pubs/InvErr_LAA02.pdf).

[^4]: Yonina C. Eldar, Aharon Ben-Tal, Arkadi Nemirovski, *Robust Mean-Squared Error Estimation in the Presence of Model Uncertainties*, IEEE Trans. Signal Process. 53(1) (2005), 168–181. DOI [10.1109/TSP.2004.838933](https://doi.org/10.1109/TSP.2004.838933); [author PDF](https://www.weizmann.ac.il/math/yonina/sites/math.yonina/files/Robust%20Mean-Squared%20Error%20Estimation.pdf).

[^5]: Martin J. Gander, Hui Zhang, *Schwarz Methods by Domain Truncation*, Acta Numerica 31 (2022), 1–134. DOI [10.1017/S0962492922000034](https://doi.org/10.1017/S0962492922000034); [arXiv:2207.09791](https://arxiv.org/abs/2207.09791).

[^6]: Werner Kirsch, *An Invitation to Random Schrödinger Operators*, Panoramas et Synthèses 25 (2008), 1–119; Theorem 5.20, p.40, Eqs. (5.53)–(5.55). [arXiv:0709.3707](https://arxiv.org/abs/0709.3707); [author PDF](https://www.fernuni-hagen.de/mi/fakultaet/emeriti/docs/kirsch/invitation.pdf).

[^7]: Martin J. Gander, Lukáš Jakabčin, Michal Outrata, *Domain Truncation, Absorbing Boundary Conditions, Schur Complements, and Padé Approximation*, ETNA 59 (2023), 319–341; Definition 2.1 p.323, Theorem 3.11 pp.330–331, Robin optimization §4 pp.333–335. [Official PDF](https://etna.ricam.oeaw.ac.at/vol.59.2023/pp319-341.dir/pp319-341.pdf).

[^8]: Karl Gustafson, *The Angle of an Operator and Positive Operator Products*, Bull. Amer. Math. Soc. 74 (1968), 488–492. DOI [10.1090/S0002-9904-1968-11974-3](https://doi.org/10.1090/S0002-9904-1968-11974-3). 另见 *The Geometrical Meaning of the Kantorovich–Wielandt Inequalities*, Linear Algebra Appl. 296 (1999), 143–151, DOI [10.1016/S0024-3795(99)00106-8](https://doi.org/10.1016/S0024-3795(99)00106-8).

[^9]: Daniel Bennequin, Martin J. Gander, Laurence Halpern, *A Homographic Best Approximation Problem with Application to Optimized Schwarz Waveform Relaxation*, Math. Comp. 78 (2009), 185–223; Theorem 2.2 p.186, Theorem 2.6 p.189. DOI [10.1090/S0025-5718-08-02145-5](https://doi.org/10.1090/S0025-5718-08-02145-5); [author PDF](https://www.math.univ-paris13.fr/~halpern/Publis/MathComp08.pdf).

[^10]: Kang Li, Ján Špakula, Jiawen Zhang, *Measured Asymptotic Expanders and Rigidity for Roe Algebras*, IMRN 2023(17), 15102–15154; Proposition 4.8, Corollary 4.21, Theorem 6.1. DOI [10.1093/imrn/rnac242](https://doi.org/10.1093/imrn/rnac242). 本地：`literature/03_sparse_inverse_graph_filters/2020_li_spakula_zhang_measured_asymptotic_expanders_roe.pdf`。

[^11]: Ján Špakula, Jiawen Zhang, *Quasi-Locality and Property A*, J. Funct. Anal. 278 (2020), 108299; Theorem 3.3, manuscript pp.9–13. DOI [10.1016/j.jfa.2019.108299](https://doi.org/10.1016/j.jfa.2019.108299). 本地：`literature/03_sparse_inverse_graph_filters/2019_spakula_zhang_quasi_locality_property_a.pdf`。

[^12]: Narutaka Ozawa, *Embeddings of Matrix Algebras into Uniform Roe Algebras and Quasi-Local Algebras*, Theorems A–B and Corollary C, pp.1–2. DOI [10.4171/JEMS/1672](https://doi.org/10.4171/JEMS/1672); [arXiv:2310.03677](https://arxiv.org/abs/2310.03677). 本地：`literature/03_sparse_inverse_graph_filters/2023_ozawa_matrix_embeddings_quasi_local_not_uniform_roe.pdf`。

[^13]: Gregory F. Lawler, Vlada Limic, *Random Walk: A Modern Introduction*, Cambridge UP (2010); Theorem 1.6.1, Lemma 4.6.1 p.96, Proposition 4.6.2 p.97, §9.10 p.233. Book DOI [10.1017/CBO9780511750854](https://doi.org/10.1017/CBO9780511750854); [author PDF](https://www.math.uchicago.edu/~lawler/srwbook.pdf).

[^14]: Anna Ben-Hamou, Roberto I. Oliveira, Yuval Peres, *Estimating Graph Parameters with Random Walks*, Math. Stat. Learn. 1 (2018), 375–399; Eq. (3.1), p.380. DOI [10.4171/MSL/9](https://doi.org/10.4171/MSL/9); [official PDF](https://ems.press/content/serial-article-files/28922).

[^15]: Enrique Bendito, Ángeles Carmona, Andrés M. Encinas, *Potential Theory for Schrödinger Operators on Finite Networks*, Rev. Mat. Iberoam. 21 (2005), 771–818. DOI [10.4171/RMI/435](https://doi.org/10.4171/RMI/435); [official PDF](https://ems.press/content/serial-article-files/38136).

[^16]: Ángeles Carmona, Andrés M. Encinas, María José Jiménez, Àngela Martín, *Random Walks Associated with Symmetric M-matrices*, Linear Algebra Appl. 693 (2024), 324–338. DOI [10.1016/j.laa.2023.10.009](https://doi.org/10.1016/j.laa.2023.10.009); [repository PDF](https://upcommons.upc.edu/server/api/core/bitstreams/8857a4d1-6399-42af-a6c9-9c631b772bc1/content).

[^17]: Charles Bordenave, Marc Lelarge, *Resolvent of Large Random Graphs*, Random Structures Algorithms 37 (2010), 332–352. DOI [10.1002/rsa.20313](https://doi.org/10.1002/rsa.20313); [arXiv:0801.0155](https://arxiv.org/abs/0801.0155).

[^18]: Nalini Anantharaman, Mostafa Sabri, *Quantum Ergodicity on Graphs: From Spectral to Spatial Delocalization*, Ann. Math. 189 (2019), 753–835. DOI [10.4007/annals.2019.189.3.3](https://doi.org/10.4007/annals.2019.189.3.3); [official PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v189-n3-p03-s.pdf).

[^19]: Florian Dörfler, Francesco Bullo, *Kron Reduction of Graphs with Applications to Electrical Networks*, IEEE TCAS-I 60(1) (2013), 150–163; Theorem III.8. DOI [10.1109/TCSI.2012.2215780](https://doi.org/10.1109/TCSI.2012.2215780); [author PDF](https://motion.me.ucsb.edu/pdf/2011d-db.pdf).

[^20]: Stephen Demko, William F. Moss, Philip W. Smith, *Decay Rates for Inverses of Band Matrices*, Math. Comp. 43 (1984), 491–499; Propositions 2.1–2.2, Theorem 2.4, pp.492–494. DOI [10.1090/S0025-5718-1984-0758197-9](https://doi.org/10.1090/S0025-5718-1984-0758197-9). 本地：`literature/03_sparse_inverse_graph_filters/1984_demko_moss_smith_decay_inverse_band_matrices.pdf`。

[^21]: Stéphane Jaffard, *Propriétés des matrices « bien localisées » près de leur diagonale et quelques applications*, Ann. Inst. H. Poincaré C 7(5) (1990), 461–476. DOI [10.1016/S0294-1449(16)30287-6](https://doi.org/10.1016/S0294-1449(16)30287-6); [official record/PDF](https://ems.press/journals/aihpc/articles/4077063). 本地：`literature/03_sparse_inverse_graph_filters/1990_jaffard_matrices_bien_localisees.pdf`。

[^22]: Karlheinz Gröchenig, Ziemowit Rzeszotnik, Thomas Strohmer, *Quantitative Estimates for the Finite Section Method*, Theorem 12, manuscript p.15. [arXiv:math/0610588](https://arxiv.org/abs/math/0610588); [author PDF](https://www.math.ucdavis.edu/~strohmer/papers/2006/fsm.pdf).

[^23]: Cheng Cheng, Min-Hsiu Jiang, Qiyu Sun, *Spatially Distributed Sampling and Reconstruction*, Appl. Comput. Harmon. Anal. 47 (2019), 109–148; Theorem 5.3, Theorems 6.1–6.2, Proposition 7.1. DOI [10.1016/j.acha.2017.07.007](https://doi.org/10.1016/j.acha.2017.07.007). 本地：`literature/03_sparse_inverse_graph_filters/2015_cheng_jiang_sun_spatially_distributed_sampling_reconstruction.pdf`。

[^24]: Michael L. Stein, *2010 Rietz Lecture: When Does the Screening Effect Hold?*, Ann. Statist. 39(6) (2011), 2795–2819; Conjecture 1, Theorems 1–2. DOI [10.1214/11-AOS909](https://doi.org/10.1214/11-AOS909); [arXiv:1203.1801](https://arxiv.org/abs/1203.1801).

[^25]: Ioannis Anagnostides, Christoph Lenzen, Bernhard Haeupler, Goran Zuzic, Themis Gouleakis, *Almost Universally Optimal Distributed Laplacian Solvers via Low-Congestion Shortcuts*, DISC 2022; Theorems 1–2, pp.6:5–6:6. DOI [10.4230/LIPIcs.DISC.2022.6](https://doi.org/10.4230/LIPIcs.DISC.2022.6); [official PDF](https://drops.dagstuhl.de/storage/00lipics/lipics-vol246-disc2022/LIPIcs.DISC.2022.6/LIPIcs.DISC.2022.6.pdf).

[^26]: Alexandr Andoni, Robert Krauthgamer, Yosef Pogrow, *On Solving Linear Systems in Sublinear Time*, ITCS 2019; Theorems 1.1–1.2 and Proposition 3.3. DOI [10.4230/LIPIcs.ITCS.2019.3](https://doi.org/10.4230/LIPIcs.ITCS.2019.3); [arXiv:1809.02995](https://arxiv.org/abs/1809.02995).

[^27]: Kevin Scaman et al., *Optimal Algorithms for Smooth and Strongly Convex Distributed Optimization in Networks*, ICML 2017. [PMLR paper](https://proceedings.mlr.press/v70/scaman17a.html). 本地：`literature/06_resource_tradeoffs/2017_scaman_et_al_optimal_distributed_optimization.pdf`。

[^28]: Chandler Davis, W. M. Kahan, Harry F. Weinberger, *Norm-Preserving Dilations and Their Applications to Optimal Error Bounds*, SIAM J. Numer. Anal. 19(3) (1982), 445–469. DOI [10.1137/0719029](https://doi.org/10.1137/0719029).

[^29]: Imre Pólik, Tamás Terlaky, *A Survey of the S-Lemma*, SIAM Review 49(3) (2007), 371–418; 单二次约束及严格可行条件下的 lossless 形式。DOI [10.1137/S003614450444614X](https://doi.org/10.1137/S003614450444614X).

[^30]: 本项目内部交叉核验，*Block robust boundary：端点反例、精确 SDP 与已知框架边界*，Theorem 4、§4、§8。路径：[block_robust_boundary_attack.md](./block_robust_boundary_attack.md)。这是数学推导证据，不是外部原创性证据。
