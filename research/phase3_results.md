# 第三阶段结果：强主定理门槛、数学碰撞与本地规模验证

> 状态日期：2026-09-13  
> 当前判决：旧理论包不能作为原创主定理；一般块鲁棒边界的 exact SDP 虽然数学正确，但已被 single-full-block Petersen lemma / lossless S-procedure 直接包含，同样不能独立承担主定理。集群和 MPI 不在本阶段声称范围内。

> **2026-09-14 收口更新：** 本文件记录被淘汰路线及既有 HPC 结果；论文主线已转为并冻结在
> grounded graph smoother 的 finite-hidden-budget exact hull 与 unbounded sharp limit。
> 最终定理、原创性边界和验证入口见 `main_theorem_decision_v2.md`，不再由本文件的旧候选结论决定。

## 0. 结论先行

本阶段把以下四点设为论文主定理的硬门槛：

1. 非平凡的必要充分条件或等价刻画；
2. 能实际构造估计器的有限算法；
3. 说明更弱条件为何失败的匹配反例或下界；
4. 与最近已知数学定理逐变量比较后，仍有实质原创余量。

按此门槛，information radius、Schur 分解、节点分位数、Green 函数和条件数衰减不能靠组合包装成强主定理。它们分别是 optimal recovery、block inverse/Dirichlet-to-Neumann、序统计、potential theory 和 Kantorovich 型不等式的特例或直接推论。[^1][^2][^3][^4][^5]

本阶段还推导出一个正确且实用的工具：对远端 Schur 补落在单个 Loewner 谱窗内的 completion 类，一般实 block 的最优共同局部线性误差可由有限 SDP 精确求出。该结果给出 iff、最优局部系数和最坏 completion 的恢复路线，并有二维解析反例排除“只检查两个各向同性端点”的错误简化。

但是专项查重已经确认，这个 exact SDP 可改写成单个 norm-bounded full block，继而由非严格 Petersen lemma 一行推出。[^6][^7] 局部系数只是 nominal block 中的仿射决策变量，并没有产生新的 lossless 数学机制。因此该结果保留为算法工具和应用命题，不作为原创主定理。

计算侧已经完成独立实跑：PGLib 14 至 10,000 节点的稀疏 DC-WLS 全部通过；约一百万状态的稀疏网格在本机 RTX 5080 上跑通批量 SpMV 和固定轮局部 Richardson。GPU 对 10k 单右端 CG 明显慢于 CPU，而对百万状态、常驻设备的 40 步局部迭代取得跨复跑约 52--62 倍的 kernel 吞吐优势。正负结果均保留。

## 1. Exact SDP 对应的网络问题

取根节点 $i$ 的半径 $r$ 区域 $B$，将全局 SPD 信息矩阵和正规方程右端分块为

$$
J=
\begin{bmatrix}
A&E\\
E^T&D
\end{bmatrix},
\qquad
b=
\begin{bmatrix}
b_B\\ b_O
\end{bmatrix},
\qquad A\succ0.
$$

令 $R_i$ 抽取根输出，并定义

$$
\Sigma=D-E^TA^{-1}E,
\qquad
C=R_iA^{-1}E,
\qquad
F=E^TA^{-1}.
$$

节点经过 $r$ 轮泛洪后知道 $A,E,b_B$，但不知道远端 $D,b_O$。所有只使用 $b_B$ 的线性估计器都可写成

$$
\widehat x_i=(R_iA^{-1}+Q)b_B.
$$

块逆公式给出它相对中心 WLS 输出的误差算子

$$
\bigl[C\Sigma^{-1}F-Q,\;-C\Sigma^{-1}\bigr].
\tag{1.1}
$$

第二列来自通信半径外的不可见右端项；局部算法不能用自由矩阵把它消掉。

## 2. 一般块谱窗的精确 SDP

考虑 completion 类

$$
mI\preceq\Sigma\preceq MI,
\qquad 0<m<M,
$$

以及共同局部误差

$$
R(C,F;m,M)=
\inf_Q\sup_{mI\preceq\Sigma\preceq MI}
\left\|[C\Sigma^{-1}F-Q,-C\Sigma^{-1}]\right\|_2.
\tag{2.1}
$$

置

$$
\alpha=M^{-1},\quad \beta=m^{-1},\quad
h=\frac{\alpha+\beta}{2},\quad
d=\frac{\beta-\alpha}{2},
$$

并定义

$$
D_0=
\begin{bmatrix}
-I&hC^T\\
hC&-\alpha\beta CC^T
\end{bmatrix},
\quad
E_0=
\begin{bmatrix}
0&0\\0&I
\end{bmatrix},
\quad
L_Q=
\begin{bmatrix}
F^T&-Q^T\\
I&0
\end{bmatrix}.
$$

已证明

$$
\boxed{
R(C,F;m,M)^2=\min_{\tau,Q,\lambda}\tau
}
\tag{2.2}
$$

满足

$$
\tau\ge0,\qquad\lambda\ge0,\qquad
\boxed{
\begin{bmatrix}
-\tau E_0+\lambda D_0&L_Q^T\\
L_Q&-I
\end{bmatrix}\preceq0.}
\tag{2.3}
$$

所以，存在同一个局部线性系数 $Q$，使所有合法 completion 和所有单位范数右端项的根误差不超过 $\varepsilon$，当且仅当 SDP (2.2)--(2.3) 的最优值不超过 $\varepsilon^2$；最优 $Q^\star$ 直接构造估计器。

若 $C=0$，根与割边精确解耦，最优值为零。若 $m=M$，completion 唯一，应直接求普通 operator-norm 近似，不应在退化谱窗上机械套用严格可行的 S-lemma。

### 2.1 正确性的短证明

令 $S=\Sigma^{-1}$。对任意向量 $v$，

$$
\{Sv:\alpha I\preceq S\preceq\beta I\}
=
\{z:\|z-hv\|_2\le d\|v\|_2\}.
\tag{2.4}
$$

正向包含来自自伴收缩；反向可用把 $v/\|v\|$ 映到目标方向的对称 Householder 算子构造。对式 (1.1) 的谱范数取左测试向量 $u$，再令 $v=C^Tu,z=Sv$，鲁棒条件变成一个齐次二次约束蕴含另一个二次不等式。因为 $m<M,C\ne0$ 时有严格可行点，单约束 S-lemma 无损；Schur complement 即给出式 (2.3)。[^8]

### 2.2 端点谱成立，只查两个端点矩阵却会失败

Loewner 区间的极点为

$$
S=\alpha I+(\beta-\alpha)P,\qquad P=P^T=P^2.
$$

因此可选择一个最坏 completion，使 $\Sigma$ 的特征值全在 $\{m,M\}$，但投影 $P$ 的方向不能省略。

取

$$
m=1,\quad M=2,\quad C=[1\;0],\quad F=[0\;4]^T.
$$

只检查 $S=\frac12I,I$ 得到的最优值是 $1$；完整 Loewner 区间上，即使先对 $Q$ 优化，真实值仍为

$$
R=\frac{5\sqrt{17}}{16}\approx1.288470508.
\tag{2.5}
$$

最坏点是非各向同性的 rank-one 投影极点，形成约 $28.8\%$ 的严格 gap。解析值、最优 $Q=0$ 和 SDP multiplier $\lambda=17$ 已由确定性脚本复核。

## 3. 为什么 exact SDP 仍不足以做主定理

表面上，式 (2.2)--(2.3) 同时有 exact iff、可计算构造和严格反例。但令

$$
B_Q=
\begin{bmatrix}
hF^TC^T-Q^T\\hC^T
\end{bmatrix},
\qquad
G=d
\begin{bmatrix}
F^T\\I
\end{bmatrix},
\qquad H=C^T,
$$

式 (2.4) 会把原问题精确改写为

$$
\inf_Q\sup_{\Delta^T\Delta\preceq I}
\|B_Q+G\Delta H\|_2.
\tag{3.1}
$$

这正是单 norm-bounded full block 的 robust-performance/model-matching 问题。van Waarde 等 2023 的 Proposition 4.16(b) 给出非严格 Petersen lemma 的必要充分 LMI；将 $B_Q,G,H$ 代入，再作 Schur complement，就得到与式 (2.3) 等价的证书。[^7] 共同 approximate inverse 的 min--max 量词则更早见于 El Ghaoui。[^6]

判决是：网络语义和不可见列使这个应用问题有意义，但没有产生新的 full-block elimination theorem。任意单一矩阵区间的椭球推广也落入既有 QMI image/contraction-completion 机制；标量谱窗下把 exterior 压到 cut 子空间再加至多一个方向，则是旋转对称带来的计算化简。二者均不足以挽救主定理原创性。

## 4. 其余数学结果的去留

| 结果 | 数学状态 | 安全位置 |
|---|---|---|
| 固定模型的最优球外行尾 | classical radius of information 特例 | 背景引理 |
| 同一局部 view 跨 completion 的 Chebyshev 半径 | optimal recovery + LOCAL indistinguishability | 问题定义 |
| Schur/oracle--Dirichlet 分解 | block inverse、DtN、resolvent | 推导工具 |
| sharp 条件数常数 $(\kappa+1)/(2\sqrt\kappa)$ | Kantorovich/antieigenvalue 推论 | baseline |
| 节点 $1-\delta$ 分位数 | 逐点误差的序统计恒等式 | 直接 corollary |
| killed-walk occupation/remaining lifetime | Green 函数 + strong Markov | 图论解释 |
| arbitrary-completion zero/infinity no-go | 一行 Schur instability | 假设边界 |
| scalar-cut robust 闭式 | 正确；一般 SDP 的一维特例 | 示例或附录 |
| 一般块 exact SDP | 正确；Petersen lemma 直接实例 | 算法工具，不作主定理 |

多数节点分位数若没有新的 completion-pasting 或低谱质量结构定理，只是计数恒等式，不能单独包装为贡献。

## 5. 下一主定理必须逃离 single-full-block 包含

下一轮只保留四类可能跨过硬门槛的方向：

1. **稀疏且可实现的 completion。** 强制远端由有界度图、局部测量因子或 M-matrix/Laplacian 网络实现；给出新的 exact iff，或证明复杂性并给可证最优近似。
2. **跨节点共享协议。** 重叠局部球的系数不能逐根独立选择，而必须由同一有限 opcode 或有限参数集生成；需要刻画联合约束。
3. **多数节点的低谱质量条件。** 对节点谱测度 $\mu_i$，研究小特征值对 $J^{-1}$ 行范数的贡献在除 $\delta$ 节点外消失，是否恰好等价于存在共同有限度局部多项式，并给显式迭代与轮数界。
4. **同一协议内的资源 converse。** 在同一输入类和失败语义下证明 rounds、bit-hop、memory 的下界，并与可执行算法匹配。

任意矩阵区间和 cut 维数压缩仍值得作为实现优化，但它们本身仍是 full-block 几何的直接推论。没有完成上述至少一项，不进入论文主张阶段。

## 6. 本地 CPU：PGLib DC-WLS 已跑通

PGLib-OPF v23.07 提供固定版本的电网拓扑与参数，但不是现成状态估计 benchmark；本项目用公开脚本冻结 DC 测量、噪声、持有者和参考节点。[^9] 从因子到 $H,J=H^TR^{-1}H$ 全程直接构造 COO/CSR，不形成 dense $H,J,J^{-1}$。

| case | buses / state | measurements | nnz(J) | case total | 中心相对残差 |
|---|---:|---:|---:|---:|---:|
| IEEE 14 | 14 / 13 | 26 | 73 | 0.0041 s | $4.56\times10^{-16}$ |
| IEEE 118 | 118 / 117 | 226 | 707 | 0.0126 s | $1.19\times10^{-14}$ |
| IEEE 300 | 300 / 299 | 499 | 1,741 | 0.0160 s | $1.15\times10^{-13}$ |
| PEGASE 1354 | 1,354 / 1,353 | 2,395 | 7,359 | 0.0685 s | $1.08\times10^{-13}$ |
| PEGASE 2869 | 2,869 / 2,868 | 5,442 | 16,300 | 0.1510 s | $1.55\times10^{-12}$ |
| GOC 10000 | 10,000 / 9,999 | 16,202 | 49,557 | 0.4745 s | $8.72\times10^{-14}$ |

六个 case 的抽样 whitened tail 均随半径非增；独立复跑复现了输入 SHA-256、矩阵维度、非零元、残差和该单调性。这里的时间只说明稀疏管线正确且规模可承受，不代表 AC block-WLS 或稳定 microbenchmark 性能。

另有一个 40 节点、块维 2、仅 $5\%$ 节点带全秩 anchor 的 synthetic case。全局 $J$ 可解，而 $95\%$ 节点没有本地 anchor，直接离开“每个 $A_i$ 满列秩”。

## 7. 本地 GPU：正结果和负结果均保留

项目隔离环境使用 CuPy 的 CUDA 12.x wheel及其 SPD sparse CG 接口。[^10][^11] 本机识别 RTX 5080，CPU 固定为单线程口径。

在同一个约 $9999\times9999$、约 49k 非零元的 PGLib GOC 矩阵上，Jacobi-CG 使用相同 RHS、相同容差：

| 指标 | CPU | GPU |
|---|---:|---:|
| 迭代次数 | 7,470 | 7,476 |
| 最终相对残差 | $9.28\times10^{-9}$ | $8.68\times10^{-9}$ |
| algorithm time | 0.241 s | 1.807 s |

GPU 明显更慢，因为每步工作量太小，Python 驱动和 kernel launch 占主导。

百万状态实验直接构造 $L_{\rm grid}+0.2I$ 和单锚点 grounded grid；约 4,996,000 个非零元，CSR 三数组约 60.99 MiB。跨最终运行与独立复跑：

- 标量 SpMV resident-kernel 相对单线程 CPU 约快 32--52 倍；
- 40 步 Richardson resident-kernel 约快 52--62 倍；
- CPU/GPU 数值相对差约 $10^{-16}$；
- $n=10,000$ 时 GPU 仍更慢，观测 crossover 位于 10k 与约 50k 状态之间；
- 首次 CUDA/cuSPARSE 初始化约 0.13--1.8 s，若计入冷启动，短任务结论会逆转。

单锚点压力例在 40 步后相对残差约 0.0028，但状态误差约 0.09；稳定谱隙例的残差约 0.015、状态误差约 0.013。病态系统中小残差不能替代状态误差，也不能单凭 GPU 吞吐证明局部可估计性。

## 8. 当前 go/no-go

- **旧理论包：NO-GO 作为主定理。**
- **一般块 exact SDP：NO-GO 作为主定理。** 它是 Petersen lemma 的直接实例；保留为计算边界规则。
- **新主线：尚未封口。** 只继续稀疏 completion 的结构 iff、跨节点共同 opcode、多数节点低谱质量 iff，或同一协议下 matching 资源下界。
- **实验路线：GO。** 稀疏 CPU/GPU 管线足以支持后续理论验证；近期不需要集群，也不作 MPI 声称。

## 9. 可复现入口

- 定理包：research/theory/theorem_package_v1.md
- 定理总审计：research/agent_reports/theorem_novelty_audit.md
- Exact SDP 与反例：research/agent_reports/block_robust_boundary_attack.md
- Petersen lemma 专项审计：research/agent_reports/sdp_exact_novelty_audit.md
- 矩阵区间与 cut 压缩证明：research/theory/candidate_main_theorem.md
- 反例脚本：experiments/theory_search/block_boundary_counterexample.py
- 数学回归：experiments/test_theory_identities.py
- 稀疏实现报告：research/agent_reports/scalable_experiments_implementation.md
- 百万状态 GPU 报告：research/agent_reports/local_gpu_scale.md

## Sources

[^1]: Arthur G. Werschulz, *An Overview of Information-Based Complexity*, Columbia Technical Report CUCS-022-02 (2002), p. 3. [Official PDF](https://mice.cs.columbia.edu/getTechreport.php?format=pdf&techreportID=152).

[^2]: Moni Naor and Larry Stockmeyer, “What Can Be Computed Locally?”, *SIAM Journal on Computing* 24(6) (1995), 1259--1277. [DOI 10.1137/S0097539793254571](https://doi.org/10.1137/S0097539793254571).

[^3]: Martin J. Gander and Hui Zhang, “Schwarz Methods by Domain Truncation,” *Acta Numerica* 31 (2022), 1--134. [DOI 10.1017/S0962492922000034](https://doi.org/10.1017/S0962492922000034).

[^4]: Gregory F. Lawler and Vlada Limic, *Random Walk: A Modern Introduction*, Cambridge University Press (2010). [DOI 10.1017/CBO9780511750854](https://doi.org/10.1017/CBO9780511750854).

[^5]: Karl Gustafson, “The Angle of an Operator and Positive Operator Products,” *Bulletin of the AMS* 74 (1968), 488--492. [DOI 10.1090/S0002-9904-1968-11974-3](https://doi.org/10.1090/S0002-9904-1968-11974-3).

[^6]: Laurent El Ghaoui, “Inversion Error, Condition Number, and Approximate Inverses of Uncertain Matrices,” *Linear Algebra and its Applications* 343--344 (2002), 171--193. [DOI 10.1016/S0024-3795(01)00273-7](https://doi.org/10.1016/S0024-3795(01)00273-7).

[^7]: Henk J. van Waarde, M. Kanat Camlibel, Jaap Eising, and Harry L. Trentelman, “Quadratic Matrix Inequalities with Applications to Data-Based Control,” *SIAM Journal on Control and Optimization* 61(4) (2023), 2251--2281; Proposition 4.16(b). [DOI 10.1137/22M1486807](https://doi.org/10.1137/22M1486807).

[^8]: Imre Pólik and Tamás Terlaky, “A Survey of the S-Lemma,” *SIAM Review* 49(3) (2007), 371--418. [DOI 10.1137/S003614450444614X](https://doi.org/10.1137/S003614450444614X).

[^9]: IEEE PES Power Grid Library, [PGLib-OPF v23.07](https://github.com/power-grid-lib/pglib-opf/tree/v23.07) and [license](https://github.com/power-grid-lib/pglib-opf/blob/v23.07/LICENSE).

[^10]: CuPy, “[Installation](https://docs.cupy.dev/en/stable/install.html).” Official documentation.

[^11]: CuPy, “[cupyx.scipy.sparse.linalg.cg](https://docs.cupy.dev/en/latest/reference/generated/cupyx.scipy.sparse.linalg.cg.html).” Official documentation.
