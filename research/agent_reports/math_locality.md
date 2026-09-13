# “局部可推断性”的数学本质：跨领域定理级尽调

**尽调日期：2026-09-13。** 本报告只讨论“中心估计能否由有限轮局部通信实现/逼近”的数学条件。重点区分：（i）固定实例还是图族；（ii）全节点同时正确还是允许按节点测度有 \(\delta\) 比例失败；（iii）连续 Gaussian/WLS 还是有限字母 Gibbs；（iv）算法可预先知道全局系数，还是只能看到半径 \(r\) 内的模型。页码以下优先指论文印刷页；没有印刷页时注明 PDF 页。检索结论中的“未发现”只表示截至本次所列主来源和交叉引文检索，没有发现相应定理，不是绝对不存在声明。

## 0. 结论先行

1. **固定线性估计器、允许离线知道全局系数时，最宽的有限轮充要条件非常简单，而且与无环、单点满秩都无关：中心估计矩阵相应的行必须支撑在该节点的 \(r\)-邻域内；近似版则由邻域外行核的对偶范数精确给出。** 这是信息集/线性代数结论，不是需要新论文才能成立的深定理。

2. **若要求一个图族上全网向量误差的统一 \(\ell_2\to\ell_2\) 控制，正确抽象是有限传播算子的算子范数闭包，即 uniform Roe algebra。** “可由有界轮局部线性算子统一逼近 iff 属于 uniform Roe algebra”本身是定义的展开，不应包装成结构定理。真正的结构定理是：在有界几何且有 Property A 时，quasi-locality（远隔集合间作用趋零）等价于这种有限传播逼近 [S14, Theorem 3.3, PDF pp. 9–13]。

3. **不能在任意有界度图族上宣称“相关性/远程影响衰减 iff 有限轮全网逼近”。** Ozawa 已证明：含 expander 序列的 uniformly locally finite 空间上存在 quasi-local、但不属于 uniform Roe algebra 的算子，因而它不能被任何有限传播算子在算子范数下一致逼近 [S15, Corollary C, pp. 1–2]。这直接卡死无额外几何假设的通用 operator-level 版本。不过 Ozawa 的算子并未被证明是稀疏 Gaussian/WLS 的逆或后验均值，所以它限制“一般算子定理”，不自动否定某个特殊 WLS 类定理。

4. **有限字母局部 Gibbs 推断领域已有最接近题目原话的真正结构性 iff。** Feng–Yin 的 Theorem 5.1 证明：任意 LOCAL 近似边缘推断算法蕴含 strong spatial mixing；反向在“locally admissible local Gibbs distributions”上成立，轮数就是达到目标 TV 衰减所需的距离加 \(O(1)\) [S1, Definition 5.1 and Theorem 5.1, pp. 21–23]。这是直接撞车，但它是有限字母、TV、任意可行边界、每个节点正确的版本，并不是连续 Gaussian/WLS，也不是允许 \(\delta\) 比例输出节点失败的版本。

5. **Gaussian/GaBP 文献解决的是特定消息迭代何时收敛，而不是“是否存在任意 \(r\)-轮局部估计器”。** 例如 Du et al. 对其 vector Gaussian BP 得到 \(\rho(Q)<1\) 当且仅当均值消息从任意初始化收敛到中心 MMSE [S5, Theorem 13 and Corollary 17, pp. 18, 24]；这是真 iff，但算法特定、无限迭代，并且模型仍假设每个相连局部系数块满列秩。Walk-summability \(\rho(|R|)<1\) 在一般 loopy 图上只是 GaBP 收敛的充分条件，不是必要条件 [S4, Proposition 1, p. 2040; Proposition 21 and discussion, pp. 2049–2052]。

6. **“允许 \(\delta\) 比例节点失败”尚未发现已有的通用 Gaussian/WLS 结构 iff。** measured coarse geometry 确实有字面上非常接近的 \(1-\delta\) 测度穷尽定理，但其对象是 block-rank-one averaging/projection 的准局部性和 Roe 代数成员关系；高测度子集是 expander core，不是“估计正确的节点集合”。因此除 rank-one 投影特例外，它只能借用语言，不能当作本问题已解决。

7. **最可能形成原创且可证明的主结果**不是再声称一个无条件“相关性衰减 iff 局部计算”，而是一个分层的 **measured localizability theorem**：用随机根节点的邻域外行核误差之“依测度收敛”精确刻画 \(1-\delta\) 节点的局部可估计性；在一致可积/统一行有界条件下，再等价于归一化 Hilbert–Schmidt 尾；若升级为全网谱范数，则必须另加 Property A、固定块秩、Jaffard/逆闭代数等结构，Ozawa 说明这些桥梁不能省。

## 1. 先把问题分成四个不同量词版本

令 \(G=(V,E)\)，每个测量块 \(z_j\) 归属节点 \(j\)，\(B_r(i)\) 是节点 \(i\) 的闭 \(r\)-球，\(\Pi_{B_r(i)}\) 是保留这些测量坐标的投影。中心线性估计写成

\[
\widehat x=Kz,\qquad
K=(H^*R^{-1}H+Q)^{-1}H^*R^{-1},
\]

其中 \(Q\succeq0\) 可表示 Gaussian 先验精度；无先验 WLS 取 \(Q=0\)，但需相应可识别条件。节点 \(i\) 只关心 \(C_i\widehat x\)，记 \(T_i=C_iK\)。

必须分别回答：

| 版本 | 算法知道什么 | 成功量词 | 自然数学对象 |
|---|---|---|---|
| A. 固定实例、全局离线设计 | 可预先知道整个 \(G,H,R,Q,K\) | 每个指定节点/所有数据 | 行支撑、有限传播算子 |
| B. 图族、统一全网误差 | 系数可按实例离线设计，但半径需与 \(n\) 无关 | \(\sup_n\|K_n-L_n\|_{2\to2}\) | uniform Roe algebra / decay algebra |
| C. 允许 \(\delta\) 节点失败 | 同 B | 对至少 \(1-\delta\) 的节点逐点达标 | 行尾误差依节点测度的分位数/依测度收敛 |
| D. 只允许局部模型 | 一个统一规则只读 rooted \(r\)-ball 的模型标签和数据 | A/B/C 任一语义 | LOCAL indistinguishability、rooted completions、local weak/local-global limits |

版本 A–C 都可能暗含“全局设计器已把正确系数下发给各节点”。这不等同于版本 D。许多所谓 distributed filtering 论文只证明在线乘法局部，却把特征值界、滤波系数、增益或预条件器的全局计算排除在通信成本之外。

## 2. 固定实例的精确充要条件：这是基线，不是结构定理

### 2.1 任意解码器也逃不过的行支撑定理

节点 \(i\) 在 \(r\) 轮后只看到 \(\Pi_Bz\)，其中 \(B=B_r(i)\)。即使允许任意非线性解码器 \(\phi_i\)，对所有 \(z\) 精确复现线性目标 \(T_i z\) 的充要条件都是

\[
\ker \Pi_B\subseteq\ker T_i
\iff T_i(I-\Pi_B)=0
\iff T_i=T_i\Pi_B. \tag{2.1}
\]

必要性来自两组具有相同局部观测、但邻域外数据相反的输入；充分性取 \(\phi_i(\Pi_Bz)=T_i\Pi_Bz\)。所以“无环 + 每个 \(A_i\) 满列秩”显然不是最宽条件：真正条件是**中心映射这一行有没有用到半径外信息**。

若输入以 \(\ell_2\) 单位球作最坏情形，任意局部解码器的精确 minimax 误差为

\[
\inf_{\phi_i}\sup_{\|z\|_2\le1}
\|T_i z-\phi_i(\Pi_Bz)\|_2
=\|T_i(I-\Pi_B)\|_{2\to2}. \tag{2.2}
\]

下界仍用 \(z_\perp\) 与 \(-z_\perp\) 的不可区分性，截断线性解码器达到上界。标量输出、输入 \(\ell_p\) 球时，右端就是邻域外系数行的 \(\ell_{p^*}\) 范数。因此必须先声明信号范数：\(\ell_2\)、每节点幅度有界的 \(\ell_\infty\)、Bayes MSE 会给出不同“最宽条件”。

全网向量若要求单个谱范数保证，令 \({\cal P}_r\) 为每个输出行只支撑于相应 \(r\)-球的线性算子，则最佳误差是

\[
d_r(T)=\inf_{L\in{\cal P}_r}\|T-L\|_{2\to2}. \tag{2.3}
\]

式 (2.1)–(2.3) 是可行信息集的线性代数展开；若把“\(r\)-轮局部可估计 iff \(d_r(T)\le\varepsilon\)”称为主定理，会是 tautology。真正需要研究的是哪些**可检查的图/谱/概率条件**等价于 \(d_r(T)\to0\)，以及半径对 \(\varepsilon,\delta\)、条件数和几何常数的定量依赖。

### 2.2 无噪声线性可识别性的真正最宽条件

若局部可用观测为 \(y_B=H_Bx\)，目标为 \(Fx\)，则从 \(y_B\) 对所有 \(x\) 精确恢复 \(Fx\) 的充要条件是

\[
\ker H_B\subseteq\ker F
\iff \operatorname{row}(F)\subseteq\operatorname{row}(H_B)
\iff \exists L:\;F=LH_B. \tag{2.4}
\]

只有在 \(F=I\)（恢复全部状态）时，它才退化为 \(H_B\) 满列秩。故“每个本地 \(A_i\) 满列秩”远强于必要条件；最多需要某个通信球内的**联合测量矩阵**对目标功能满秩。动态系统只需把 \(H_B\) 换成该传感器集合生成的 observability rows \([C_B;C_BA;\ldots]\)。这正是 functional observability 的行空间本质，见 §8。

但 (2.4) 只解决 noiseless identifiability。要求局部估计与中心 WLS/MMSE **完全相同**时，仍需 (2.1) 的中心增益行支撑；一个目标局部可识别，不代表忽略远端噪声观测后仍能取得中心最优方差。

### 2.3 Bayes/Gaussian 版本

对任意平方可积随机变量，局部信息 \({\cal F}_{B_r(i)}=\sigma(Z_{B_r(i)})\) 下的最佳均方估计是

\[
\widehat X_i^{(r)}=\mathbb E[X_i\mid Z_{B_r(i)}].
\]

由条件期望的正交投影恒等式，

\[
\mathbb E\|X_i-\widehat X_i^{(r)}\|^2-
\mathbb E\|X_i-\mathbb E[X_i\mid Z]\|^2
=\mathbb E\|\mathbb E[X_i\mid Z]-\widehat X_i^{(r)}\|^2. \tag{2.5}
\]

所以局部 Bayes risk 与中心 risk 的差恰好等于两个估计器之间的 MSE。精确相等 iff 中心后验均值对局部 \(\sigma\)-代数可测。联合 Gaussian 且相关条件协方差非退化时，这又等价于

\[
X_i\perp Z_{V\setminus B_r(i)}\mid Z_{B_r(i)}
\iff \operatorname{Cov}(X_i,Z_{V\setminus B}\mid Z_B)=0
\iff I(X_i;Z_{V\setminus B}\mid Z_B)=0. \tag{2.6}
\]

Gaussian 的 excess MSE 是

\[
\Delta_i(r)=\operatorname{tr}\!left[
\operatorname{Cov}(X_i\mid Z_B)-\operatorname{Cov}(X_i\mid Z)
\right]\ge0. \tag{2.7}
\]

这给出了概率论上的精确充要条件，但仍是 conditional sufficiency 的基本恒等式，而非只看原始图拓扑就能判断的结构定理。连续 Gaussian 若允许边界值任意大，有限字母 Gibbs 文献那种“对任意两组边界配置的统一 TV 衰减”通常不会原样成立；应改用单位扰动的 operator/Lipschitz/Wasserstein 影响或 (2.7)。

## 3. 真正直接撞车：strong spatial mixing 与 LOCAL inference

Feng–Yin [S1] 在标准 LOCAL 模型中允许每轮无界消息和无界本地计算；\(t\) 轮节点只能区分半径 \(t\) 的带标签实例。其 approximate inference 要求每个节点输出有限字母边缘分布，且对任意实例、任意可行部分边界条件都满足 TV 误差 \(\le\delta\) [Definition 2.2 之后，p. 6]。

其 Definition 5.1 定义 strong spatial mixing：对类中每个 \(n\)-节点分布、每个节点 \(v\)、任意边界集合 \(\Lambda\) 和两组可行配置 \(\sigma,\tau\)，若差异集合为 \(D\)，则

\[
d_{\rm TV}(\mu_v^\sigma,\mu_v^\tau)
\le \delta_n(d_G(v,D)). \tag{3.1}
\]

**Theorem 5.1（pp. 21–23）的精确量词：**

- 必要性：对任何联合分布类，若存在误差任意 \(\delta>0\)、轮数 \(t(n,\delta)\) 的 LOCAL 边缘推断算法，则该类有 SSM rate
  \[
  \delta_n(t)=2\min\{\delta:t(n,\delta)\le t-1\}.
  \]
  证明就是两个半径内不可区分实例加三角不等式。
- 充分性：若该类是 **locally admissible, local Gibbs distributions**，且有 rate \(\delta_n(t)\)，则存在 LOCAL 推断算法，轮数
  \[
  t(n,\delta)=\min\{t:\delta_n(t)\le\delta\}+O(1).
  \]

这是本次检索中最接近“局部推断存在 iff 相关性衰减”的非平凡、算法级定理。限制必须原样保留：有限字母；TV；边界条件意义的 SSM；必要方向一般，充分方向需要局部 Gibbs 与 locally admissible；正确性要求每个节点，不是允许 \(\delta\) 比例节点失败；没有通信字节/能耗复杂度。Corollary 5.3（p. 24）进一步在指数 SSM 下给出 \(O((1-\alpha)^{-1}\log^3n)\) 轮 exact sampling，但这也不是 Gaussian state estimation。

Chen–Peng–Liu [S2, Theorem 1, p. 4988; Theorem 2 and Corollary 2, pp. 4989–4990] 用 Dobrushin comparison 给出另一条严格充分路线：若 \(c=\max_i\sum_j C_{ij}<1\)，局部截断模型与全局模型在查询节点的边缘误差由 \(D=(I-C_{\alpha\alpha})^{-1}\) 和边界影响显式控制。它允许有环，但依赖全局收缩/有限离散状态，只是充分条件，不是所求的最宽 iff。

**撞车判断：** 若论文主张“对有限字母 local Gibbs，LOCAL 可推断性由 SSM 精确刻画”，基本已被 [S1] 做掉。若主张“连续 Gaussian/WLS、按节点测度允许失败、并同时处理轮数—误差—拓扑”，[S1] 只能作为离散先例和不可区分性证明模板。

## 4. Gaussian graphical models 与 GaBP：解决的是算法收敛，不是任意局部可实现性

### 4.1 computation tree 与相关性衰减

Weiss–Freeman [S3, NIPS version pp. 674–680] 证明有限次 loopy BP 等价于在相应深度的 unwrapped computation tree 上做精确推断；其式 (2)–(3) 将原图与展开树的均值/方差误差关联到根—叶条件协方差。足够快的根叶相关衰减保证 BP 收敛；一旦 Gaussian BP 收敛，均值正确，但方差一般不正确。这是重要充分机制，不是一般有限半径 iff。

Malioutov–Johnson–Willsky [S4] 对标准化精度矩阵 \(J=I-R\) 给出：

- Proposition 1（p. 2040）：walk-summability 等价于所有绝对 walk sums 收敛，等价于 \(\sum_{k\ge0}|R|^k\) 收敛，等价于 \(\rho(|R|)<1\)，等价于 \(I-|R|\succ0\)。这是 WS 类本身的 iff。
- Proposition 21（p. 2049）：WS \(\Rightarrow\) LBP well-posed，均值收敛到真值，方差收敛到 backtracking self-return walks 的和。
- p. 2051 明说：WS 对树和单环是相应收敛条件的必要充分条件，一般图只充分；有非-WS 但 LBP 收敛的例子。
- Proposition 25（pp. 2051–2052）：computation-tree 极限谱半径 \(\rho_\infty<1\) 保证方差收敛；\(\rho_\infty>1\) 最终失效；边界 \(=1\) 未由该命题解决。

所以不能把 \(\rho(|R|)<1\) 宣称为“一般 Gaussian 局部估计存在”的必要条件。

### 4.2 特定 vector Gaussian BP 的真正 iff

Du–Ma–Wu–Kar–Moura [S5] 的模型为每个 agent 的局部线性 Gaussian 观测
\(y_n=\sum_{i\in\{n\}\cup I(n)}A_{n,i}x_i+z_n\)，假设每个出现的 \(A_{n,i}\) 满列秩、先验和噪声协方差正定，且全局 \(A\) 满列秩（pp. 4–5）。

- Theorem 6（p. 14）：任意半正定信息消息初始化下，信息矩阵消息收敛到唯一正定不动点。
- 信息矩阵收敛后，均值消息成为 \(v^{(\ell)}=-Qv^{(\ell-1)}+b\)。Theorem 13（p. 18）证明从任意初始化收敛到唯一值 **iff** \(\rho(Q)<1\)。
- Corollary 17（p. 24）：belief mean 从任意允许初始化收敛到中心 MMSE **iff** \(\rho(Q)<1\)；factor graph 为 forest 加一个 single loop 是充分拓扑条件。

这是需要正面承认的直接相关 iff，但它只终结“这套 BP 无限迭代是否收敛”，不终结“是否存在另一套有限轮/近似局部算法”；\(Q\) 又由收敛后的消息矩阵构成，不是单纯图论性质；而且原问题想放松的局部满列秩仍在假设中。

Giscard et al. [S6, Theorem 2, pp. 4–5] 则对任意有限图上的正定 \(J\) 给出 \(J^{-1}\) 每个 block entry 的有限 path-sum/branched continued-fraction 精确表达，包含所有 simple paths/cycles。它表明有环不妨碍精确代数公式，但枚举使用全局图，最坏复杂度很高，不能据此得到固定半径 LOCAL 算法。

## 5. 稀疏 SPD 逆、Green 函数衰减与局部重构

### 5.1 有环完全允许；条件数给出充分的指数局部化

Demko–Moss–Smith [S7] 的经典结果已经排除“必须无环”：

- Proposition 2.1（p. 492）给出 \(1/x\) 在 \([a,b]\) 上的最佳多项式逼近，其几何率由
  \(q=(\sqrt\kappa-1)/(\sqrt\kappa+1)<1\)、\(\kappa=b/a\) 控制。
- Proposition 2.2 和 Theorem 2.4（pp. 492–494）对正定 banded、可逆算子给出逆矩阵条目的指数离对角衰减。
- Proposition 5.1（p. 498）用
  \(S_m(A)=\bigcup_{k=0}^m\{(i,j):(A^k)_{ij}\ne0\}\)
  推广到一般稀疏模式：在该集合外，\(|(A^{-1})_{ij}|\) 以 \(O(q^{m+1})\) 衰减。

若 \(A\) 的非零模式沿图一跳，则次数 \(m\) 的多项式 \(p_m(A)\) 有传播 \(m\)，且证明直接给出
\(\|A^{-1}-p_m(A)\|_2\le O(q^{m+1}/a)\)。这是全网谱范数、有限轮构造性的充分结果；关键是图族上有统一谱隙/条件数，而非无环。它不是必要条件：条件数变坏也可能由于特殊右端、目标行或坏节点很少而仍局部可估计。

Benzi–Razouk [S8, Theorem 3.3, p. 21; Theorem 3.4, pp. 24–25] 将其推广到稀疏可对角化矩阵的解析矩阵函数。若 \(f\) 在包含谱的更大复区域解析，存在 \(K>0,\lambda<1\) 使
\(|f(A)_{ij}|\le\kappa(X)K\lambda^{d(i,j)}\)，并有相应 polynomial/finite-propagation 的算子范数逼近。Theorem 3.3 中 Bernstein 的“几何多项式逼近 iff 可解析延拓”是**函数 \(f\)** 层面的 iff，不是“某个图族的逆局部化 iff 条件数有界”。

Motee–Sun [S10, Theorems 6.1–6.2, PDF p. 13; Theorems 8.1, 9.1–9.2, pp. 17, 22–24; Theorem 10.1, p. 25] 用 weighted Gröchenig–Schur/Jaffard 型衰减代数处理无限网络：适当 \(q\)-Banach 子代数逆闭；Lyapunov/Riccati 解与反馈增益留在同一局部化代数；截断误差满足 \(\|K-K^{(T)}\|\le Cw(T)^{-1}\)。这是强大的充分闭包框架，但起点就是矩阵属于指定 decay algebra，并非普适必要条件。

### 5.2 最接近“局部稳定性 iff 全局稳定性”的结果

Cheng–Jiang–Sun [S9] 研究 spatially distributed sampling/reconstruction，在 agent graph 计数测度 doubling、测量矩阵属于 Jaffard class \(J_\alpha\)、\(\alpha>d\) 等假设下：

- Theorem 5.3（manuscript p. 15）为 Wiener lemma：\(A\in J_\alpha\) 且在 \(\ell_2\) 上 boundedly invertible \(\Rightarrow A^{-1}\in J_\alpha\)。
- Theorem 6.1（p. 16）证明全局 \(\ell_2\) lower stability 推出所有足够大 quasi-main finite restrictions 的统一 lower bound。
- Theorem 6.2（p. 17）在显式 tail inequality 下给出反向：一个固定足够大半径的所有 quasi-main restrictions 有统一 lower bound \(\Rightarrow\) 全局稳定。
- Proposition 7.1 / Theorem 7.2（pp. 18–19）给出局部 finite-section least squares 对中心 LS 的 \(\ell_\infty\) 多项式/指数收敛率。

这是“局部有限块稳定性与全局稳定性”等价方向的高撞车结果，且允许有环；但它依赖 doubling、polynomial growth、Jaffard decay 和显式尾界，目标是 sampling matrix 稳定及特定 finite-section 算法，不是任意 WLS 中心算子的无假设 iff。

## 6. 图滤波：有限轮构造与受限表示类

次数 \(r\) 的图移位多项式 \(p_r(S)\) 天然只有传播 \(r\)，因此是有限轮算法，但“所有有限传播算子”远大于“同一个 shift 的 node-invariant polynomial”。

- Shuman et al. [S11, Section IV-B, pp. 6–7; Proposition 4, p. 8] 的 Chebyshev 递推用 \(K\) 轮、约 \(2K|E|\) 次标量通信实现次数 \(K\) 逼近；若 \(\eta\) 个谱乘子各自的 uniform scalar approximation error 最大为 \(B(K)\)，则整体分析算子误差 \(\le\sqrt\eta B(K)\)。这是充分构造，且系数、\(\lambda_{\max}\) 上界通常是全局设计输入。
- Segarra–Marques–Ribeiro [S12, Proposition 1, pp. 4–5] 对 node-invariant polynomial graph filter 的精确可表示性给出 simultaneous diagonalizability、重复特征值一致性和次数条件；Proposition 2（pp. 5–6）在该表示类内给出 Frobenius/MSE 最优系数和 operator-norm SDP；Corollary 3（p. 9）说明 connected graph 上 consensus averaging 可选择 graph-supported shift 后在至多 \(N-1\) 次交互精确实现。
- Emirov et al. [S13, Theorem 3.1, p. 7] 对其迭代多项式逆滤波证明：对每个右端指数收敛到 \(H^{-1}b\) **iff** \(\rho(I-HG)<1\)。这是指定迭代的无限步 iff，不是任意有限传播逼近存在性的结构 iff。

因此 graph-filter 文献给出了很成熟的“轮数—多项式次数—谱逼近误差”上界和受限类的 representability 条件，但没有替代 §2 的一般传播算子条件。

## 7. Uniform Roe / quasi-locality：全网图族版本的正确抽象及其边界

把图族 \(G_n\) 做 coarse disjoint union \(X=\bigsqcup_nV_n\)，并令 \(K=\bigoplus_nK_n\) 作用在 \(\ell_2(X)\)。传播 \(\le r\) 正是矩阵块满足 \(d(i,j)>r\Rightarrow L_{ij}=0\)。于是

\[
K\in C_u^*(X)
\iff
\forall\varepsilon>0\ \exists r<\infty,\ L:\operatorname{prop}(L)\le r,
\ \sup_n\|K_n-L_n\|_{2\to2}<\varepsilon. \tag{7.1}
\]

这是 uniform Roe algebra 的定义 [S14, Definitions 2.4, 2.7, pp. 5–6]，与“图族上、全网谱范数、轮数不随 \(n\) 增长的线性局部实现”完全同义，但没有提供易检验的图/模型条件。

Quasi-locality 要求：对每个 \(\varepsilon>0\)，存在 \(R\)，使任意相距大于 \(R\) 的集合 \(A,B\) 都有
\(\|\chi_AK\chi_B\|<\varepsilon\)。Špakula–Zhang [S14, Theorem 3.3, p. 9] 证明，在有界几何空间上：对 \(p\in(1,\infty)\)（包括 Hilbert 情形 \(p=2\)），quasi-local 与多种 commutator 条件等价；若空间还有 Property A，它们再等价于 Roe membership/有限传播范数逼近。对 \(p=0,1,\infty\) 不需 Property A。Corollary/Proposition 6.4（pp. 18–19）还给出相应 inverse-closed 性。

### 7.1 Ozawa 反例精确限制什么

Ozawa [S15] 对 uniformly locally finite 空间证明：

- Theorem A：\(\prod_nM_n\) 不嵌入任何 uniform Roe algebra；
- Theorem B：若 \(X\) 含一列 expanders，\(\prod_nM_n\) 可嵌入 quasi-local algebra；
- Corollary C（p. 2）：因此 \(C_u^*(X)\subsetneq C_{ql}^*(X)\)，存在 quasi-local operator 不能在算子范数下由有限传播算子逼近。

故以下命题在一般有界度图族上是假的：

> 只要任意远隔集合之间的影响趋零（quasi-local/correlation decay），就存在统一有限轮的全网谱范数近似。

可修复方式至少有四类：（a）图族有 Property A/适当分解性质；（b）算子属于 Jaffard/Gröchenig–Schur 等已知逆闭衰减代数；（c）限制 block rank；（d）把性能降为按节点测度而非全网谱范数。

反例也有边界：它使用可容纳增长维矩阵代数的算子，不说明该算子是局部稀疏 SPD 的逆、Gaussian posterior mean，亦不说明每行在随机根意义下不可局部。因此它不能被拿来否定带 WLS 结构或 \(1-\delta\) 逐节点语义的专门定理。

2026 年 Li–Zhang–Zhu 预印本 [S18, Main Theorem I, pp. 2–3; Theorems 4.1–4.2, pp. 14–16] 进一步显示“秩”是关键边界：对 uniformly bounded geometry sparse family 的 block-diagonal partial isometries/projections，若块秩统一 \(\le M\)，quasi-locality 与相应 Roe membership 等价；rank-one 投影还有定量有限传播逼近。其 Main Theorem II 给出 asymptotic expanders 上 quasi-local 但非-Roe 的 ghost projection/unitary，坏例需要非一致有界的秩。该文截至本报告仍是 2026-08 的新预印本，应标作未正式同行评议；它支持“固定低维目标/固定秩可能有更强结构 iff”，但不是一般 WLS 定理。

## 8. Measured coarse geometry 与“允许 \(\delta\) 节点失败”：字面接近，语义不同

这是本次额外重点审查的部分。

### 8.1 已有 measured theorem 的精确对象

Li–Špakula–Zhang [S16] 取 block-rank-one projection

\[
P=\bigoplus_nP_n,\qquad P_n\eta=\langle\eta,\xi_n\rangle\xi_n,
\quad m_n(x)=\|\xi_n(x)\|^2,
\]

其中 \(m_n\) 是 \(X_n\) 上概率测度。Lemma 4.3（PDF p. 15）给出特殊 rank-one 恒等式

\[
\|\chi_AP_n\chi_B\|=
\sqrt{m_n(A)m_n(B)}. \tag{8.1}
\]

由此 Proposition 4.4（p. 15）把 \(P\) 的 quasi-locality 精确化为远隔 \(A,B\) 的质量乘积趋零；Proposition 4.8（pp. 15–16）证明

\[
P\text{ quasi-local}
\iff \{(X_n,d_n,m_n)\}\text{ 是 measured asymptotic expanders}. \tag{8.2}
\]

Definition 1.2（pp. 3–4）的 measured asymptotic expansion 量词是：对每个 \(\alpha\in(0,1/2]\)，存在与 \(n\) 无关的 \(c_\alpha,R_\alpha\)，使每个 \(\alpha\le m_n(A)\le1/2\) 的集合满足
\(m_n(\partial_{R_\alpha}A)>c_\alpha m_n(A)\)。Theorem 6.1（pp. 26–27）再证明 bounded geometry 下，对 block-rank-one \(P\)，quasi-local iff \(P\in C_u^*(X)\)。论文首页的 Theorem B 汇总为：rank-one 投影属于 uniform Roe iff 对应测度空间是 measured asymptotic expanders。

最像 \(1-\delta\) 的是 **Corollary 4.21（pp. 22–23）**：measured asymptotic expansion iff 存在 \(\alpha_k\downarrow0\)，使每个 \(X_n\) 都含有通过统一 Lipschitz 嵌入得到的 bounded-valency measured expander graph core \(V_{n,k}\)，且

\[
m_n(V_{n,k})\ge(1-\alpha_k)m_n(X_n), \tag{8.3}
\]

另有固定 \(k\) 的 expansion 常数和相邻点测度比控制；反向不需要该测度比条件。Lu–Wang–Zhang [S17, Theorem B/Thms. 3.9, 3.18, pp. 3–4, 11, 14; Proposition 4.7, p. 17; Theorem D/5.3, pp. 4, 20] 在 groupoid/dynamical Roe 框架中也得到“asymptotic expansion in measure iff 由质量趋一的 expansion domains 穷尽”，以及 averaging/rank-one projection 的 dynamical quasi-local iff Roe iff expansion。

### 8.2 为什么它不等于 \(\delta\)-node estimation

二者有四个不可偷换的差别：

1. **测度含义不同。** \(m_n(x)=\|\xi_n(x)\|^2\) 是 rank-one 向量的振幅质量；只有 \(\xi_n\) 均匀时才是节点计数比例。估计问题中的 \(\mu_n\) 通常是均匀随机节点或业务权重。
2. **高质量 core 的语义不同。** (8.3) 的 complement 不是“允许输出错误的节点”，core 是证明整个 measured space 有 expansion 的结构证人。
3. **误差对象不同。** [S16] 最终得到整个 \(P\) 的全局算子范数有限传播逼近，而不是只压缩输出行 \(\chi_SP\) 后允许 \(S^c\) 任意坏。
4. **秩极特殊。** (8.1) 完全依赖 rank one；一般 WLS/Kalman 增益的秩随网络增长，不能把证明直接搬过去。2026 的 bounded-rank 结果 [S18] 也恰好说明这一障碍。

因此，**直接撞车只发生在中心算子本身是 block-rank-one averaging/projection，且所用节点权重就是 \(|\xi_n|^2\) 的特殊片段。** 对一般多状态 WLS，“measured asymptotic expander iff \(1-\delta\) 节点可局部估计”不是已有结论，也不是 [S16]/[S17] 的推论。

### 8.3 \(\delta\)-node 版本目前能无条件得到的精确公式

给节点测度 \(\mu_n\)（均匀计数只是特例），定义固定实例的逐点最佳 \(r\)-轮误差

\[
e_{n,i}(r)=\|T_{n,i}(I-\Pi_{B_r(i)})\|_{2\to2}. \tag{8.4}
\]

允许为每个实例离线下发系数时，存在一个 \(r\)-轮算法使至少 \(1-\delta\) 测度的节点在所有 \(\|z\|_2\le1\) 上误差 \(\le\varepsilon\)，**充要条件恰为**

\[
\mu_n\{i:e_{n,i}(r)\le\varepsilon\}\ge1-\delta. \tag{8.5}
\]

对固定 \(n\)，最优轮数就是误差分布的分位数；图族的统一轮数再对 \(n\) 取上确界：

\[
r^*_{\varepsilon,\delta}
=\inf\{r:Q_{1-\delta}(e_{n,I_n}(r))\le\varepsilon\},
\quad I_n\sim\mu_n. \tag{8.6}
\]

Bayes/Gaussian 版本把 \(e_{n,i}\) 换成 (2.7) 的 \(\Delta_{n,i}\)。若还要求一个 good set 上的**联合向量谱范数**，自然量是

\[
\inf_{S:\mu_n(S)\ge1-\delta}\;
\inf_{L\in{\cal P}_r}\|\chi_S(T_n-L)\|_{2\to2}, \tag{8.7}
\]

但 (8.7) 不能简单化成逐行截断误差的最大值；多行误差可相干叠加。

式 (8.5) 是严格 iff，却仍属显式化定义，不是图结构定理。它揭示最适合 \(\delta\)-node 的数学拓扑是**随机根节点下行尾核的依测度收敛**：

\[
\forall\varepsilon,\delta>0\ \exists r<\infty:\quad
\sup_n\Pr_{I_n\sim\mu_n}[e_{n,I_n}(r)>\varepsilon]\le\delta. \tag{8.8}
\]

若输出为标量（或输出块维数一致有界）且逐行算子范数平方一致可积，则 (8.8) 等价于邻域外行尾的归一化 Hilbert–Schmidt 均方趋零：

\[
\lim_{r\to\infty}\sup_n
\frac1{|V_n|}\|T_n-\operatorname{trunc}_r(T_n)\|_F^2=0 \tag{8.9}
\]

（一般权重用加权行平方和）。(8.9) \(\Rightarrow\) (8.8) 用 Markov；反向用好节点上 \(\varepsilon^2\) 加坏节点上统一上界，再令 \(\varepsilon,\delta\downarrow0\)。这一路径不与 Ozawa 冲突，因为它是 row-in-measure / normalized Schatten-2，而非全网 operator norm。

一个说明范数差异的例子是 \(P_n=\mathbf1\mathbf1^*/n\)：每行 \(\ell_2\) 范数为 \(n^{-1/2}\)，故对全局 \(\ell_2\) 单位输入，零轮零估计器的每节点误差都趋零；但 \(\|P_n\|_{2\to2}=1\)，全网谱范数误差丝毫不降。若输入改为 \(\ell_\infty\) 单位球，每行 \(\ell_1\) 范数为 1，结论又不同。

### 8.4 “失败”至少有四种语义

必须在新工作中固定：

- 每个实例可有一个确定的坏节点集合，\(\mu(S^c)\le\delta\)；
- 均匀随机根节点成功概率 \(\ge1-\delta\)；
- 随机算法对每节点或对全网的失败概率；
- Byzantine/传感器故障，即有 \(\delta\) 节点提供恶意数据。

前三者已不等价，第四者更是鲁棒估计问题。Measured coarse geometry 对应“节点/振幅测度”，不是随机算法失败；Byzantine distributed observer 文献也不能视为解决“允许一些输出节点不准确”。

## 9. Graph limits / local-global convergence：可用来表述，不是已有估计 iff

Hatami–Lovász–Szegedy [S19] 的 local-global convergence 比 Benjamini–Schramm local weak convergence 更强：它考察所有有限着色后 rooted neighborhood 分布的可能极限，可表达局部算法在均匀随机节点上的输出统计。Göös–Hirvonen–Suomela 等 LOCAL lower-bound 工作中的 \((\alpha,r)\)-homogeneous 图，也用“有 \(\alpha\) 比例节点看到同一半径视图”制造不可区分性。

这些理论可以自然承载 (8.8)：把模型参数、测量类型、可能的估计器系数作为 marks，在随机根极限上讨论误差事件。但现有 graph-limit 定理并不自带中心 WLS 算子 \(K_n\)，也没有给出“中心算子可按 \(1-\delta\) 节点局部逼近 iff 某个标准 graphon/graphing 性质”。目前判断：**是可借的紧致性、随机根和下界语言，不是直接撞车。**

若进一步要求“只知道局部模型”，必须增加一个 completion indistinguishability 条件：任意两个合法全局实例，只要在根的半径 \(r\) marked neighborhood（包括所见数据）同构，同一个 deterministic LOCAL 算法就输出相同值；因此两实例的中心目标若相差 \(>2\varepsilon\)，根不可能在两边都成功。对所有 completion 的目标集合取 Chebyshev radius 可得到 all-node 的必要充分选择条件；但若每个实例可丢弃不同的 \(\delta\) 坏节点，量词变成分布式选择/耦合问题。把“存在一个 neighborhood-type selector”重新命名为 iff 仍是 tautology，真正的新内容应是用 graph limit、局部谱数据或 completion 的稳定性给出可检验上界。

## 10. Functional observability / distributed observers：动态、无限时间版本已有强 iff

Park–Martins [S20, Theorem 3.2, PDF p. 6] 考察固定有向通信图上的 LTI plant。若图的 source strongly connected components 为 \(V_1,\dots,V_s\)，则其 LTI distributed observer 存在实现所有节点 asymptotic omniscience 的参数 **iff** 每个 \((A,H_{V_i})\) detectable。论文 §III-A（PDF p. 7）进一步说明若此条件失败，则在同一通信图下，连 nonlinear/time-varying scheme 也不能实现渐近 omniscience。这是真正的图论—线性系统 iff，并远比“图无环 + 每点满秩”宽。

但它的语义是：动态无噪声/有限二阶噪声下，时间趋于无穷，每步相邻 observer 交换内部状态，所有节点最终重构全 plant state。它不回答静态 WLS 在给定有限轮内等于中心估计、误差—轮数—能耗 tradeoff 或允许 \(\delta\) 输出节点失败。

Functional observability 本身有精确行空间条件。Montanari et al. [S21, Eq. (3), main text p. 3] 给出

\[
(A,C,F)\text{ functionally observable}
\iff \operatorname{rank}\!\begin{bmatrix}O\\F\end{bmatrix}
=\operatorname{rank}(O)
\iff \operatorname{row}(F)\subseteq\operatorname{row}(O). \tag{10.1}
\]

其 Supporting Information Theorem 2（SI pp. 19–20）在“\(F,C\) 每行各选择一个 target/sensor、自由参数 generic、目标与传感器独立”等假设下，证明 structural functional observability iff：（i）每个 target 有通向 sensor 的路径；（ii）target 集不与所有 minimal dilation sets 的并相交。这是 generic identifiability 的图论 iff，不是 noisy center-equivalence。

Mitra–Sundaram 的 distributed functional observer [S22] 以 centralized Darouach rank 条件和 feasible leader 为基础，在强连通图上给出构造性充分定理；论文结论仍把 minimum-order distributed functional observer 留作开放问题。因此不能引用它声称 distributed functional estimation 已有完整 N&S。

较新空间局部 Kalman 工作 Arbelaiz et al. [S23, Proposition IV.1, PDF p. 5; Theorem V.1, p. 6] 在 matching spectral/spatially invariant PDE 等特殊结构下给出 exact decentralization 或核的指数衰减；2026 Zhao et al. [S24] 明确研究 observability–communication–connectivity trilemma，但其主要网络测量条件仍是充分条件，论文也把必要条件列为未来工作。它们说明优化方向活跃，但尚未给出本报告所定义的通用 \(\delta\)-node finite-round iff。

## 11. “已解决到哪里”对照表

| 命题版本 | 已有结果 | 必要/充分 | 尚缺什么 |
|---|---|---|---|
| 固定 \(K\)，单节点，最坏 \(\ell_2\) | 邻域外行核范数 (2.2) | 精确 iff/minimax | 只是线性代数，不是结构判据 |
| 固定 noiseless \(H\)，目标 \(F x\) | \(\ker H_B\subseteq\ker F\) (2.4) | 精确 iff | 不等于 noisy central optimum |
| Bayesian/Gaussian 单节点 | conditional expectation / covariance (2.5)–(2.7) | 精确 iff/risk gap | 需从局部模型导出可检验 decay |
| finite-alphabet local Gibbs、每节点 | SSM \(\leftrightarrow\) LOCAL inference [S1] | 必要一般；充分需 locally admissible local Gibbs | 非连续、非 \(\delta\)-node |
| 指定 vector GaBP、无限轮 | \(\rho(Q)<1\) [S5] | 算法特定 iff | 非任意算法、非有限轮 |
| 稀疏 SPD 逆、统一谱隙 | 多项式/指数衰减 [S7,S8] | 强充分条件 | 非必要；局部坏区/\(\delta\) 未刻画 |
| Jaffard + doubling sampling | local/global stability [S9] | 条件化的双向定理 | 非一般 WLS/operator iff |
| 图族、全网谱范数 | uniform Roe membership (7.1) | 定义性 iff | 结构 iff 需 Property A 等；Ozawa 阻止无条件版 |
| quasi-local \(\Rightarrow\) finite propagation | Property A 下 [S14]；rank-one/有界秩特例 [S16,S18] | 条件化 iff | 一般增长秩为假 [S15] |
| \(1-\delta\) 节点逐点误差 | 行尾分位数 (8.5)–(8.8) | 固定系数模型精确 iff | 未发现已有可检验图/谱结构 theorem |
| measured expansion 的 \(1-\delta\) core | [S16,S17] | rank-one 投影结构 iff | core 不是“成功节点”，不可泛化偷换 |
| 动态 LTI、无限时间 omniscience | source-SCC detectability [S20] | 强 iff | 非静态有限轮/noisy center-equivalence |

## 12. 最可能原创、又能避开已有反例的数学主线

建议把拟议结果命名为 **Measured Localizability / Local Estimation in Measure**，而不是笼统“locality iff correlation decay”。可分三层，第一层应能完整证明，后两层才是结构创新。

### 层 I：精确的 measured row-kernel theorem（低撞车、可立即证明）

对 uniformly bounded-degree marked graph family、节点测度 \(\mu_n\)、一致有界中心线性估计器 \(T_n\)，证明以下等价：

1. 对每个 \(\varepsilon,\delta>0\)，存在与 \(n\) 无关的 \(r\)，使至少 \(1-\delta\) 测度节点可由 \(r\) 轮、离线已知系数的局部算子达到逐节点最坏 \(\ell_2\) 误差 \(\varepsilon\)；
2. 随机根 \(I_n\sim\mu_n\) 的邻域外行算子范数 \(e_{n,I_n}(r)\) 随 \(r\to\infty\) 一致依概率趋零；
3. 在平方一致可积条件下，截断误差在加权 normalized Hilbert–Schmidt 范数中趋零。

并给出最优半径的分位数公式 (8.6)。这把“允许失败”变成精确量词，不受 Ozawa 影响，也能直接比较 \(\ell_2\)、\(\ell_\infty\)、Bayes risk。其数学难度不高，原创点主要是把 distributed estimation 的问题定义干净；不能夸为深结构 theorem。

### 层 II：WLS/Gaussian 的可检查充分—必要桥梁（最值得做）

在 \(J_n=H_n^*R_n^{-1}H_n+Q_n\) 图局部、最大度有界的条件下，寻找仅允许 \(\delta\) 异常根节点的局部谱条件，例如：随机根附近的 Dirichlet restrictions 有高概率统一 lower bound、局部 condition number/tail resolvent 一致可积，且坏 bottleneck 集测度 \(\le\delta\)。目标是证明它与

\[
\|(C_iJ_n^{-1}H_n^*R_n^{-1})(I-\Pi_{B_r(i)})\|
\to0\quad\text{in }\mu_n\text{-probability}
\]

等价或两侧有匹配常数。技术可结合：Cheng–Jiang–Sun 的 quasi-main restriction local/global stability、Combes–Thomas/Demko resolvent decay、Schur complement、random-root/graphing 紧致性。这里最可能出现真正新定理，因为现有结果几乎都要求**每处**统一谱隙/增长条件，而题目只要求**绝大多数根**。

一个可检验 conjecture 方向是“\(\mu\)-local invertibility + uniform integrability of Green rows iff inverse kernels localize in measure”；必须认真构造 rare bottleneck、expander core、低秩长程模态反例，决定还需哪些 tightness 条件。仅写“局部 block 满秩”不够，因为 Schur complement 可通过很远的近零模态放大。

### 层 III：从离线全局系数升级为 only-local-model

再增加 **local coefficient reconstructibility**：截断中心行不仅尾小，而且半径 \(r\) 内的系数能由半径 \(R(r,\varepsilon)\) 的 marked model 稳定决定；等价的必要下界来自两种 globally different completions 的 rooted-neighborhood indistinguishability。Feng–Yin 的 SSM proof 可作为离散模板，graphings/local-global convergence 用于 \(1-\delta\) 量词。真正目标应是一个 Gaussian completion-stability theorem，而非把“存在 local rule iff 存在 local rule”写成定理。

### 全网谱范数作为单独 corollary，而非默认主语义

若要升级到所有节点同时、\(\ell_2\to\ell_2\) 算子范数，应明确加上以下之一：

- coarse space 有 Property A，先证 Gaussian/WLS 算子 quasi-local，再用 [S14]；
- 直接证算子进入某个 inverse-closed Jaffard/Schur algebra [S9,S10]；
- 估计目标的 block rank 一致有界，尝试利用 [S16,S18]；
- 直接构造 polynomial/rational local approximants [S7,S8]。

没有这些桥梁，Ozawa [S15] 使一般命题为假。

## 13. 撞车风险与可证伪缺口

### 撞车风险

- **高：** “有限字母 Gibbs 的 LOCAL inference iff strong spatial mixing”——[S1] 已有精确量词和轮数。
- **高：** “指定 GaBP 收敛 iff 一个迭代矩阵谱半径小于 1”——[S5] 已有；一般 stationary inverse iteration 也有 [S13]。
- **高：** “稀疏 SPD、条件数有界 \(\Rightarrow\) 逆指数衰减/多项式局部逼近”——[S7,S8] 是经典结论。
- **中高：** “局部 finite-section stability 与全局 sampling stability”——[S9] 在 doubling + Jaffard 条件下已有双向结果。
- **高但只是定义：** “全网有限轮逼近 iff 中心算子在有限传播算子闭包”——就是 uniform Roe 定义。
- **中：** measured expansion 的 \(1-\delta\) core——[S16,S17] 已有，但仅 rank-one projection；若新稿研究 rank-one averaging 会直接相撞。
- **低至中：** 一般 Gaussian/WLS 的 row-in-measure、\(1-\delta\) 输出节点、并区分离线系数与 only-local-model 的完整 iff——本次未找到直接定理。

### 最关键的可证伪测试

1. **稀有近零模态：** 只有 \(o(n)\) 节点附近存在小谱隙，是否仍能让正比例节点的 Green rows 长程化？若能，纯“局部坏点比例小”不充分。
2. **expander / non-Property-A 图族：** 构造由稀疏 SPD inverse 或 posterior mean 实现的 Ozawa 型 quasi-local-not-Roe 算子；若成功，全网 Gaussian 版本也必须降范数或加几何假设。
3. **rank-one global mode：** \(\mathbf1\mathbf1^*/n\) 区分逐节点 \(\ell_2\)、全网谱范数和 \(\ell_\infty\) 输入，防止性能范数偷换。
4. **同半径局部视图、不同 global completion：** 在 cycles/covers/lifts 上构造中心 WLS 行相差常数的成对实例，给 only-local-model 轮数下界。
5. **好集边界：** 即使 \(1-\delta\) 节点局部谱良好，坏集的 Schur complement 是否通过大边界影响好集？需要容量/capacity 而不只是节点计数。
6. **Bayes 与 adversarial 输入分离：** 平均 MSE 小可能由远端模态概率小造成，最坏算子误差仍大；两种 theorem 不应混写。

最终判断：**截至 2026-09，已有文献分别解决了离散 Gibbs 的 SSM iff、特定 GaBP 的谱半径 iff、稀疏逆的充分衰减、全网 operator-locality 的 Roe/Property-A 理论、以及动态 observer 的 detectability iff；但没有发现把“一般连续 Gaussian/WLS + 有限轮 + 至少 \(1-\delta\) 节点 + only-local-model”四者合在一起的最宽结构充要条件。** 这正是最可信的研究空白，但必须以上述分层方式推进，不能用一个未经限定的“correlation decay iff locality”口号宣称彻底终结。

## 主来源

- **[S1]** Weiming Feng, Yitong Yin, *On Local Distributed Sampling and Counting*, arXiv:1802.06686. https://arxiv.org/abs/1802.06686 。本地：`literature/02_graphical_models_local_inference/2018_feng_yin_local_distributed_sampling_counting.pdf`。
- **[S2]** Yuxin Chen, Jian Peng, Qiang Liu, *Efficient Localized Inference for Large Graphical Models*, IJCAI 2018, DOI: [10.24963/ijcai.2018/692](https://doi.org/10.24963/ijcai.2018/692)。本地：`literature/02_graphical_models_local_inference/2018_chen_peng_liu_efficient_localized_inference.pdf`。
- **[S3]** Yair Weiss, William T. Freeman, *Correctness of Belief Propagation in Gaussian Graphical Models of Arbitrary Topology*, NeurIPS proceedings version. https://proceedings.neurips.cc/paper/1999/hash/10c272d06794d3e5785d5e7c5356e9ff-Abstract.html 。本地：`literature/02_graphical_models_local_inference/2001_weiss_freeman_gaussian_bp_arbitrary_topology.pdf`（归档文件名按 2001 扩展版命名，所核验 PDF 为 proceedings pp. 674–680）。
- **[S4]** Dmitry M. Malioutov, Jason K. Johnson, Alan S. Willsky, *Walk-Sums and Belief Propagation in Gaussian Graphical Models*, JMLR 7 (2006), 2031–2064. https://jmlr.org/papers/v7/malioutov06a.html 。本地：`literature/02_graphical_models_local_inference/2006_malioutov_johnson_willsky_walk_sums_gaussian_bp.pdf`。
- **[S5]** Jian Du, Shengli Ma, Yih-Fang Wu, Soummya Kar, José M. F. Moura, *Convergence Analysis of Distributed Inference with Vector-Valued Gaussian Belief Propagation*, JMLR 18(172), 2018. https://jmlr.org/papers/v18/16-556.html 。本地：`literature/02_graphical_models_local_inference/2018_du_et_al_vector_gaussian_bp_convergence.pdf`。
- **[S6]** Pierre-Louis Giscard et al., *Exact Inference on Gaussian Graphical Models of Arbitrary Topology using Path-Sums*, JMLR 17 (2016). https://jmlr.org/papers/v17/14-445.html 。本地：`literature/02_graphical_models_local_inference/2016_giscard_et_al_path_sums_exact_gaussian_inference.pdf`。
- **[S7]** Stephen Demko, William F. Moss, Philip W. Smith, *Decay Rates for Inverses of Band Matrices*, Mathematics of Computation 43 (1984), 491–499, DOI: [10.1090/S0025-5718-1984-0758197-9](https://doi.org/10.1090/S0025-5718-1984-0758197-9)。本地：`literature/03_sparse_inverse_graph_filters/1984_demko_moss_smith_decay_inverse_band_matrices.pdf`。
- **[S8]** Michele Benzi, Nader Razouk, *Decay Bounds and O(n) Algorithms for Approximating Functions of Sparse Matrices*, ETNA 28 (2007), 16–39. https://etna.ricam.oeaw.ac.at/vol.28.2007-2008/pp16-39.dir/pp16-39.pdf 。本地：`literature/03_sparse_inverse_graph_filters/2007_benzi_razouk_decay_sparse_matrix_functions.pdf`。
- **[S9]** Cheng Cheng, Min-Hsiu Jiang, Qiyu Sun, *Spatially Distributed Sampling and Reconstruction*, Applied and Computational Harmonic Analysis 47 (2019), 109–148, DOI: [10.1016/j.acha.2017.07.007](https://doi.org/10.1016/j.acha.2017.07.007), arXiv:1511.08541. 本地规范目录暂未归档。
- **[S10]** Nader Motee, Qiyu Sun, *Sparsity and Spatial Localization Measures for Spatially Distributed Systems*, SIAM J. Control Optim. 55(1) (2017), 200–235, DOI: [10.1137/15M1049294](https://doi.org/10.1137/15M1049294), arXiv:1402.4148. 本地规范目录暂未归档。
- **[S11]** David I. Shuman, Pierre Vandergheynst, Daniel Kressner, Pascal Frossard, *Distributed Signal Processing via Chebyshev Polynomial Approximation*, arXiv:1111.5239. https://arxiv.org/abs/1111.5239 。本地：`literature/03_sparse_inverse_graph_filters/2011_shuman_et_al_distributed_chebyshev_approximation.pdf`。
- **[S12]** Santiago Segarra, Antonio G. Marques, Alejandro Ribeiro, *Optimal Graph-Filter Design and Applications to Distributed Linear Network Operators*, IEEE TSP 65(15) (2017), DOI: [10.1109/TSP.2017.2703660](https://doi.org/10.1109/TSP.2017.2703660)。本地：`literature/03_sparse_inverse_graph_filters/2015_segarra_marques_ribeiro_distributed_linear_network_operators.pdf`。
- **[S13]** Emil Emirov, Cheng Cheng, Min-Hsiu Jiang, Qiyu Sun, *Polynomial Approximation of Inverse Graph Filters*, DOI: [10.1007/s43670-021-00019-x](https://doi.org/10.1007/s43670-021-00019-x), arXiv:2003.11152. 本地：`literature/03_sparse_inverse_graph_filters/2020_emirov_cheng_jiang_sun_polynomial_inverse_graph_filter.pdf`。
- **[S14]** Ján Špakula, Jiawen Zhang, *Quasi-Locality and Property A*, Journal of Functional Analysis 278 (2020), 108299, DOI: [10.1016/j.jfa.2019.108299](https://doi.org/10.1016/j.jfa.2019.108299), arXiv:1809.00532. 本地：`literature/03_sparse_inverse_graph_filters/2019_spakula_zhang_quasi_locality_property_a.pdf`。
- **[S15]** Narutaka Ozawa, *Embeddings of Matrix Algebras into Uniform Roe Algebras and Quasi-Local Algebras*, JEMS online 2025, DOI: [10.4171/JEMS/1672](https://doi.org/10.4171/JEMS/1672), arXiv:2310.03677. 本地规范目录暂未归档。
- **[S16]** Kang Li, Ján Špakula, Jiawen Zhang, *Measured Asymptotic Expanders and Rigidity for Roe Algebras*, IMRN 2023, 15102–15154, DOI: [10.1093/imrn/rnac242](https://doi.org/10.1093/imrn/rnac242), arXiv:2010.10749. 本地规范目录暂未归档。
- **[S17]** Xinlu Lu, Qin Wang, Jiawen Zhang, *Asymptotic Expansion for Groupoids and Roe-Type Algebras*, JNCG online first 2026, DOI: [10.4171/JNCG/666](https://doi.org/10.4171/JNCG/666)。本地规范目录暂未归档。
- **[S18]** Kang Li, Jiawen Zhang, Jingming Zhu, *K-Theoretic Comparison of Roe and Quasi-Local Algebras via Projections*, arXiv:2608.22439. https://arxiv.org/abs/2608.22439 。截至 2026-09 为预印本；本地规范目录暂未归档。
- **[S19]** Hamed Hatami, László Lovász, Balázs Szegedy, *Limits of Locally-Globally Convergent Graph Sequences*, arXiv:1205.4356. https://arxiv.org/abs/1205.4356 。本地规范目录暂未归档。
- **[S20]** Se Young Park, Nuno C. Martins, *Design of Distributed LTI Observers for State Omniscience*, IEEE TAC 62(2) (2017), DOI: [10.1109/TAC.2016.2560766](https://doi.org/10.1109/TAC.2016.2560766), arXiv:1401.0926. 本地：`literature/05_localized_control_observability/2014_park_martins_lti_distributed_observers_iff.pdf`。
- **[S21]** A. N. Montanari et al., *Functional Observability and Target State Estimation in Large-Scale Networks*, PNAS 119 (2022), DOI: [10.1073/pnas.2113750119](https://doi.org/10.1073/pnas.2113750119)。本地：`literature/05_localized_control_observability/2022_functional_observability_target_state_estimation.pdf`。
- **[S22]** Aritra Mitra, Shreyas Sundaram, *Distributed Functional Observers for LTI Systems*, arXiv:1705.10891. https://arxiv.org/abs/1705.10891 。本地：`literature/05_localized_control_observability/2017_mitra_sundaram_distributed_functional_observers.pdf`。
- **[S23]** Julen Arbelaiz et al., *How Far to Share Measurements?*, IEEE TAC, DOI: [10.1109/TAC.2024.3504257](https://doi.org/10.1109/TAC.2024.3504257), arXiv:2406.14781. 本地：`literature/05_localized_control_observability/2025_arbelaiz_et_al_how_far_to_share_measurements.pdf`。
- **[S24]** Zhao et al., 2026 work on the discrete-LTI observability–communication–connectivity design trilemma, European Journal of Control, DOI: [10.1016/j.ejcon.2026.101559](https://doi.org/10.1016/j.ejcon.2026.101559), arXiv:2603.20144. 本地：`literature/05_localized_control_observability/2026_zhao_et_al_discrete_lti_design_trilemma.pdf`。
