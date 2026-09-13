# 主定理决策稿（scope frozen）：未知隐藏网络的精确局部逼近

> 状态日期：2026-09-14  
> 当前判决：**GO，研究边界已冻结**。两项独立数学审计、专项原始文献查重与 44 项
> deterministic checks 均已完成。本文不声称终结一般分布式状态估计；它精确解决
> scalar、undirected、nonnegative、grounded graph completion 上的 common linear
> local-rule 问题。dense-SPD、directed、signed、block-valued 与 nonlinear 版本不进入主文。

## 1. 论文真正要解决的问题

固定 \(q\) 个可见端口 \(\Gamma\)，端口以外允许出现任意有限数量的隐藏节点、
任意无向含环拓扑和任意非负边权。令

\[
D_G=\operatorname{diag}(D_\Gamma,D_Z)\succ0,
\qquad
T_G=D_G^{1/2}(D_G+L_G)^{-1}D_G^{1/2}.
\]

\(T_G\) 是 grounding-normalized graph-Tikhonov smoother。局部规则只能读取端口输入，
但性能必须对全部隐藏输入和所有 admissible completions 成立。给定输出矩阵
\(C\in\mathbb R^{p\times q}\) 和受通信/协议约束的
\(Q\in\mathcal Q_{\rm loc}\)，完整误差算子是

\[
\mathcal A_G(Q)=C[I_\Gamma\ 0]T_G-Q[I_\Gamma\ 0].
\]

在 unit-grounding 主模型中，\(T_Gy\) 正是 quadratic state-smoothing problem

\[
\arg\min_x\left\{\frac12\|x-y\|_2^2+
\frac12\sum_{(i,j)}w_{ij}(x_i-x_j)^2\right\}.
\]

因此这里的 state estimation 含义是明确的 graph-regularized linear estimation，
而不是把结论外推到任意动态 Kalman/WLS estimator。

问题是：是否存在一个共同的局部 \(Q\)，使

\[
\sup_G\|\mathcal A_G(Q)\|_2\le \varepsilon,
\]

以及这个无限拓扑 supremum 能否精确计算，而不是用一般 full-block uncertainty 做保守外包。

论文只使用两种嵌套 completion 类：

- **有限预算模型**：所有节点 unit grounding，hidden vertices 至多为 \(h\)；
- **无界模型**：hidden 数量、total grounding 与 edge conductance 均无统一上界。

有限预算模型是物理主模型；无界模型给出拓扑要求最弱的 sharp limit。其余模型不再加入
本论文，以免继续移动边界。

## 2. 有限隐藏预算主定理

对端口 full partition \(\pi\) 及整数分配
\(k_B\ge0,\sum_{B\in\pi}k_B\le h\)，定义

\[
X_{\pi,k}
=\sum_{B\in\pi}\frac{\mathbf1_B\mathbf1_B^T}{|B|+k_B},
\qquad
\mathcal X_{q,h}=\{X_{\pi,k}\}.
\]

### 定理 A：finite-\(h\) completion 的精确闭凸包

\[
\boxed{
\overline{\operatorname{conv}}
\{\big[(I+L_G)^{-1}\big]_{\Gamma,\Gamma}:|Z|\le h\}
=\operatorname{conv}\mathcal X_{q,h}.}
\]

“至多 \(h\)”与“恰好 \(h\)”给出同一闭凸包；要求每个有限图 connected 也不改变
结果。每个场景均可由 connected、maximum-degree-2 paths 逼近。因此 acyclicity 不再是
问题级必要条件，cycles 和任意无向拓扑已被精确覆盖。

候选场景数恰为

\[
N(q,h)=\sum_{s=1}^q S(q,s){h+s\choose s}.
\]

若 \(h\ge1\)，其 exposed vertices 恰为 \(k=0\) 或
\(\sum_Bk_B=h\) 的场景，数量为

\[
N_{\rm vert}(q,h)
=\sum_{s=1}^qS(q,s)
\left[1+{h+s-1\choose s-1}\right].
\]

### 定理 B：完整隐藏输入、所有 Schatten norms 的精确归约

令 \(T_G=(I+L_G)^{-1}\)、\(E_\Gamma=[I_q\ 0]\)，并定义

\[
\mathcal A_G(Q)=CE_\Gamma T_G-QE_\Gamma.
\]

则对每个 \(1\le p\le\infty\)，

\[
\boxed{
\sup_{G:|Z|\le h}\|\mathcal A_G(Q)\|_{S_p}
=\max_{X\in\mathcal X_{q,h}}
\left\|[CX-Q,\ C(X-X^2)^{1/2}]\right\|_{S_p}.}
\]

这里 \(C(X-X^2)^{1/2}\) 正是不能丢弃的 hidden-input 能量。对 \(p\ge2\)，风险是
\(X\) 的凸函数，可只查 \(N_{\rm vert}\) 个 vertices；对 \(1\le p<2\)，精确等式仍成立，
但应保留全部 \(N(q,h)\) 个 candidates。

谱范数设计可写成无保守有限 SDP。令 \(\rho\) 表示平方半径；对每个 retained scenario
\(X\) 施加

\[
\boxed{
\begin{bmatrix}
\rho I_p-CXC^T+CXQ^T+QXC^T&Q\\
Q^T&I_q
\end{bmatrix}\succeq0.}
\]

因此任何由 zero pattern、\(r\)-hop 可访问集、共享系数、无偏约束或量化字典凸松弛定义的
\(Q\in\mathcal Q_{\rm loc}\)，都有一个 finite necessary-and-sufficient certificate。

### 定理 C：预算--保守性的精确速率

记 finite hull 为 \(K_h\)，无界 partial-partition hull 为 \(K_\infty\)。则

\[
\boxed{d_H^{\rm op}(K_h,K_\infty)=\frac q{q+h}.}
\]

并且对 fixed \(Q\) 的 exact spectral radius，

\[
0\le R_\infty(Q)^2-R_h(Q)^2
\le
\bigl(\|C\|_2^2+2\|C\|_2\|Q\|_2\bigr)\frac q{q+h}.
\]

这不是一个只有渐近阶的说法：terminal uncertainty set 的 Hausdorff 距离对所有
\(q,h\) 都恰好等于 \(q/(q+h)\)。

## 3. 无界隐藏 completion 的 sharp limit

对非空 \(A\subseteq\Gamma\)，令

\[
v_A=D_\Gamma^{1/2}\mathbf1_A,
\qquad
P_A^{(d)}=\frac{v_Av_A^T}{v_A^Tv_A}.
\]

对 active subset \(U\subseteq\Gamma\) 的任意分割
\(\pi=\{A_1,\ldots,A_s\}\)，定义

\[
P_{\pi,U}^{(d)}=\sum_{A\in\pi}P_A^{(d)},
\]

并在 inactive ports 上补零。所有这类投影构成
\(\mathcal P_\Gamma^\partial(d)\)，数量恰为 Bell 数 \(B_{q+1}\)。

### 定理 D：无界 completion 的精确端口闭凸包

若 completion 类允许隐藏节点数和隐藏 total grounding 不受统一上界约束，则

\[
\boxed{
\overline{\operatorname{conv}}
\{T_G[\Gamma,\Gamma]:G\}
=\operatorname{conv}\mathcal P_\Gamma^\partial(d).}
\]

此外：

- 每个 \(P\in\mathcal P_\Gamma^\partial(d)\) 都是 exposed vertex；
- 每个 vertex 都能由连通、最大度数为 \(2\) 的加权路径 completion 渐近实现；
- 因而对任意只依赖端口压缩的连续凸损失 \(\Phi\)，

  \[
  \sup_G\Phi(T_G[\Gamma,\Gamma])
  =\max_{P\in\mathcal P_\Gamma^\partial(d)}\Phi(P).
  \]

这同时给出 upper bound 与 matching converse：环、任意隐藏规模和连续边权不会制造比这些
有限场景更坏的端口行为；反过来，每个场景确实都来自允许的低度网络极限。

### 定理 E：无界 completion 的完整隐藏输入归约

若 \(\psi:\mathbb S_+^p\to\mathbb R\) 连续、凸且对 Loewner 序单调不减，则

\[
\boxed{
\sup_G\psi\!\left(\mathcal A_G(Q)\mathcal A_G(Q)^T\right)
=\max_{P\in\mathcal P_\Gamma^\partial(d)}
\psi\!\left((CP-Q)(CP-Q)^T\right).}
\]

关键不是删除 hidden columns，而是对每个 full forest projection
\(\widehat P\) 使用

\[
\begin{aligned}
&\big(C[I\ 0]\widehat P-Q[I\ 0]\big)
 \big(C[I\ 0]\widehat P-Q[I\ 0]\big)^T\\
&\qquad=CXC^T-CXQ^T-QXC^T+QQ^T,
\qquad X=\widehat P[\Gamma,\Gamma],
\end{aligned}
\]

右侧对 \(X\) 是仿射的。于是特别得到

\[
\sup_G\|\mathcal A_G(Q)\|_{S_p}
=\max_{P\in\mathcal P_\Gamma^\partial(d)}\|CP-Q\|_{S_p},
\qquad 2\le p\le\infty.
\]

这一 Schatten 范数范围是 sharp 的：对每个 \(1\le p<2\)，一个 port、一个 hidden node
的两节点例子就使 full completion 的值等于 \(2^{1/p-1/2}>1\)，而两个 partial scenarios
都等于 \(1\)。因此 universal finite-scenario equality 在 Schatten family 中成立当且仅当
\(p\ge2\)。

### 推论 F：无界 completion 的 finite exact LMI

存在 \(Q\in\mathcal Q_{\rm loc}\) 使所有 completion 的 spectral error 不超过
\(\varepsilon\)，当且仅当

\[
\boxed{
\begin{bmatrix}
\varepsilon I_p&CP-Q\\
(CP-Q)^T&\varepsilon I_q
\end{bmatrix}\succeq0
\quad
\text{for every }P\in\mathcal P_\Gamma^\partial(d).}
\]

所以无限个隐藏规模、拓扑和连续参数被精确压缩成 \(B_{q+1}\) 个有限 LMI。变量规模和
每个 LMI 的阶数不依赖隐藏网络大小。

## 4. 为什么这不是把经典结果换名字

必须引用、不能宣称原创的部分包括：

- grounded Laplacian resolvent 的 matrix-forest 展开；
- graph-Tikhonov smoother 是随机森林分量平均算子的期望；
- DPP 的 mean-projection identity；
- 全分割 projection matrices 与 ratio-cut polytope；
- operator norm epigraph 的 Schur LMI。

直接先例是 Pilavci 等对
\((L+D)^{-1}D\) 的 random-spanning-forest 无偏表示；一般 mean-projection 定理可见
Dereziński--Khanna--Mahoney 与 Kassel--Lévy。[^1][^2][^3]

当前没有定位到直接先例、因此仍可作为候选原创贡献的组合是：

1. finite-\(h\) resolvent port hull 及其**完整 exposed-vertex characterization**；
2. finite-\(h\) full-input all-Schatten exact reduction；
3. \(d_H(K_h,K_\infty)=q/(q+h)\) 的 exact budget hierarchy；
4. 任意隐藏节点 completion 的 weighted partial-partition **精确闭凸包**；
5. 每个 vertex 的 connected degree-2/path converse；
6. common local rule 的 finite exact LMI iff 与无界 Schatten \(p=2\) sharp threshold；
7. 每个场景在 universal certificate 中的不可删除性。

其中第 7 点可在 Frobenius loss 下直接看出：对任意目标场景 \(P\)，取
\(C=I,Q=I-P\)，则

\[
\|P'-(I-P)\|_F^2=q-\|P-P'\|_F^2,
\]

所以 \(P\) 是唯一最坏场景。

截至 2026-09-14 的原始论文定向检索没有发现同时覆盖上述组合的工作；最危险的近邻是
Pilavcı 等的 random-forest smoother，以及 Kenyon--Wilson/Lam 的 grove/cactus electrical
compactification。[^4][^5] 这个结论只能写成 “we are not aware of”，不能当作全球不存在证明。
独立审计给出的综合原创置信约为 0.60：足以支持应用数学/GSP/control 论文，但不宜宣传为
新的 matrix-forest 基础理论。

## 5. 应用收益必须怎样证明

### 5.1 真正的拓扑放宽

定理允许任意 cycles、任意有限 hidden size、任意非负 weights；甚至 matching lower
witness 只需 paths。它不是“从树推广到某个谱半径充分条件”，而是把整个允许的隐藏拓扑类
精确消掉。对原始 2013 工作而言，这说明 acyclicity 不是问题级必要条件，而只是其有限步
cavity 算法的条件。

### 5.2 可测量的计算收益

无界模型的离线设计只依赖 boundary size \(q\)，而不依赖 hidden node/edge count：

| \(q\) | \(B_{q+1}\) |
|---:|---:|
| 2 | 5 |
| 3 | 15 |
| 4 | 52 |
| 5 | 203 |
| 6 | 877 |
| 7 | 4,140 |
| 8 | 21,147 |

这是一种 boundary-size fixed-parameter exact algorithm，不是关于 \(q\) 的 polynomial-time
算法。在线阶段只执行已经得到的稀疏/共享 \(Q\)，可直接计量 FLOPs、memory、rounds 与
bit-hops。有限预算模型提供 \(N(q,h)\) 个安全场景；spectral、Frobenius 与
Schatten-\(p\ge2\) 可精确裁剪到 \(N_{\rm vert}(q,h)\) 个。实验必须同时报告“更紧
bound”与“更多离线场景”这项真实 tradeoff。

### 5.3 严格而非口头的误差收益

有限预算已经给出一个不依赖人为约束的解析收益。取

\[
q=2,\quad h=1,\quad
C=P_{\rm av}=\frac12\mathbf1\mathbf1^T.
\]

精确有限预算设计为

\[
Q_* = \frac23P_{\rm av},
\qquad
R_1^*=\frac{\sqrt2}{3}=0.471404\ldots.
\]

无界模型的最优半径为 \(R_\infty^*=1/2\)。已知 \(h=1\) 因而使半径下降
\(5.72\%\)，平方风险从 \(1/4\) 降到 \(2/9\)，下降 \(11.11\%\)；\(h=0\) 时误差为零。
这验证 finite-\(h\) theorem 确实把结构先验转化为更低的最坏误差。

无界模型中，若 \(Q\) 完全自由，则
\(0,I\in\mathcal P_\Gamma^\partial(d)\) 导致一个必须公开的退化：

\[
Q_*=C/2,
\qquad
R_*=\|C\|_2/2.
\]

因此无界模型中的非平凡设计必须让 \(\mathcal Q_{\rm loc}\) 真实表达通信零模式、跨节点共享
opcode、保常/无偏条件、量化字典或 nominal-performance 约束。

一个完全解析的 shared-opcode 例子是

\[
C=\operatorname{diag}(1,2),
\qquad Q=\alpha I_2.
\]

若漏掉 connected partial scenario \(P=\mathbf1\mathbf1^T/2\)，会得到
\(\alpha=1,R=1\)。精确五场景问题则给出

\[
\alpha_*=\frac43-\frac{2\sqrt{10}}{15},
\qquad
R_*=\frac{10+2\sqrt{10}}{15}=1.088303\ldots.
\]

也就是说，粗糙场景集把 robust radius 低估 \(8.83\%\)，把平方风险低估约
\(18.44\%\)，并输出错误的共享系数。这是理论必须带来可观测收益的最小基准；大规模
实验还要比较 full-block SDP、随机拓扑 sampling 与 exact scenarios。

## 6. 不能越过的 claim 边界

- 不声称首次发现 matrix-forest 或 random-forest Tikhonov identity。
- 不声称无界模型的任意 full-matrix convex loss 都可归约；hidden-column directional loss
  与 Schatten \(p<2\) 有反例。finite-\(h\) 定理对所有 Schatten \(p\ge1\) 成立，但
  \(p<2\) 不能只检查 terminal-hull vertices。
- heterogeneous grounding 必须使用 grounding-normalized energy；未归一化 Euclidean
  hidden RHS 的 worst case 甚至可随 hidden count 发散。
- 不把 graph-smoother input gain 冒充一般 raw sensor-noise WLS MSE；本文应用对象是明确定义的
  quadratic graph state smoother。
- 不声称自由 \(Q\) 的 direct-smoother design 很复杂；其闭式退化必须写出。
- 不声称 Bell enumeration 对大 boundary 是多项式算法。
- 不声称解决 directed、signed、block-noncommutative 或 nonlinear decoder 的一般问题。
- 不声称“终结所有分布式状态估计”；精确终结的是本文明确定义的 scalar undirected
  hidden-completion + common linear local rule 子问题。

## 7. 论文结构与投稿定位

建议把论文写成一个主定理链，而不是堆很多小结论：

1. 端口 completion 模型与 local-rule constraints；
2. finite-\(h\) exact hull、vertex characterization 与 path converse；
3. full-input all-Schatten theorem、exact SDP 与预算收益例；
4. \(h\to\infty\) Hausdorff rate、partial-partition limit 与 sharp \(p=2\) boundary；
5. scenario solver 与 boundary-size scaling；
6. synthetic/PGLib/graph-signal experiments；
7. NO-GO boundaries。

以当前强度，优先定位 applied matrix analysis、graph signal processing 或 control/estimation；
不再为了追求纯数学期刊继续扩展模型。

## 8. 下一阶段实验的通过标准

1. 对 \(q=2,\ldots,8\) 记录 scenario generation、建模、solver、peak memory 与 active
   constraints；
2. 随机生成大量 hidden graphs、cycles、weights、groundings，确认 exact certificate
   从不被 sample 击穿；
3. 用 degree-2 construction 数值逼近 SDP 返回的 active scenario，验证 converse；
4. 比较 exact partial scenarios、只查 endpoints、随机 completion sampling 和
   \(0\preceq T\preceq I\) full-block relaxation 的误差与时间；
5. 在 PGLib 和图信号数据上画 rounds / cut size / error / byte-hop / FLOPs Pareto 曲线；
6. 将理论预言与实际 CPU/GPU crossover 分开报告，避免把 kernel throughput 当成统计收益。

## 9. 当前证据入口

- 完整边界证明与反例：`research/agent_reports/boundary_partial_partition_audit.md`
- finite-\(h\) 精确定理与应用例：`research/agent_reports/finite_hidden_budget_theorem.md`
- 原始文献专项查重：`research/agent_reports/direct_smoother_novelty_audit.md`
- 一般 factor / Schur / WLS 推导：`research/agent_reports/structured_completion_main_theorem_search.md`
- rank-one inverse 顶点碰撞审计：`research/agent_reports/rank_one_inverse_vertex_theorem.md`
- 多数节点 spectral-tail supporting theorem：`research/theory/spectral_tail_majority_theorem.md`
- 数值上界与 degree-2 converse：`experiments/theory_search/hidden_partial_partition_checks.py`
- norm 与输入几何反例：`experiments/theory_search/boundary_partial_partition_counterexamples.py`
- finite-\(h\) 六组复核：`experiments/theory_search/finite_hidden_budget_checks.py`

## References

[^1]: Y. Pilavci, P.-O. Amblard, S. Barthelmé, and N. Tremblay,
    “Graph Tikhonov Regularization and Interpolation via Random Spanning Forests,”
    *IEEE Transactions on Signal and Information Processing over Networks* 7 (2021),
    359--374. [DOI 10.1109/TSIPN.2021.3084879](https://doi.org/10.1109/TSIPN.2021.3084879).

[^2]: M. Dereziński, R. Khanna, and M. W. Mahoney,
    “Improved Guarantees and a Multiple-Descent Curve for Column Subset Selection and the
    Nyström Method,” NeurIPS 2020, Supplemental Lemma 5.
    [Primary paper page](https://proceedings.neurips.cc/paper/2020/hash/342c472b95d00421be10e9512b532866-Abstract.html).

[^3]: A. Kassel and T. Lévy, “On the Mean Projection Theorem for Determinantal Point
    Processes,” *ALEA* 20 (2023), 497--504.
    [DOI 10.30757/ALEA.v20-17](https://doi.org/10.30757/ALEA.v20-17).

[^4]: T. Lam, “Electroid Varieties and a Compactification of the Space of Electrical
    Networks,” *Advances in Mathematics* 338 (2018), 549--600.
    [DOI 10.1016/j.aim.2018.09.014](https://doi.org/10.1016/j.aim.2018.09.014).

[^5]: R. W. Kenyon and D. B. Wilson, “Boundary Partitions in Trees and Dimers,”
    *Transactions of the AMS* 363(3) (2011), 1325--1364.
    [DOI 10.1090/S0002-9947-2010-04964-5](https://doi.org/10.1090/S0002-9947-2010-04964-5).
