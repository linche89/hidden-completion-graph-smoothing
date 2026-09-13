# Block robust boundary：端点反例、精确 SDP 与已知框架边界

## 摘要结论

考虑实矩阵

\[
R(C,F;m,M)=\inf_Q\ \sup_{mI\preceq\Sigma\preceq MI}
\left\|[C\Sigma^{-1}F-Q,\,-C\Sigma^{-1}]\right\|_2,
\qquad 0<m\le M,
\]

其中 \(C\in\mathbb R^{p\times n}\)、\(F\in\mathbb R^{n\times k}\)、
\(Q\in\mathbb R^{p\times k}\)。严格攻击得到四点：

1. 现有 scalar-cut 分段闭式正确，包括阈值、最优 \(Q\) 和退化情形。
2. 对固定 \(Q\)，最坏逆 Schur 补可以限制到
   \(S=\Sigma^{-1}=\alpha I+(\beta-\alpha)P\)，其中 \(P\) 是任意正交投影、
   \(\alpha=1/M,\beta=1/m\)。因此每个特征值只需取谱端点，但不能只检查两个各向同性矩阵 \(\alpha I,\beta I\)；投影方向是连续且本质的。
3. 存在一个完全解析的 \(2\times2\) 反例：只查 \(mI,MI\) 得到 1，而真实 robust 值为
   \(5\sqrt{17}/16\approx1.28847\)。这个 gap 在优化过 \(Q\) 后仍存在。
4. 更重要的是，一般 block 问题在当前精确假设下并非只能写成 semi-infinite SDP。谱范数的左测试向量把矩阵区间压成一个精确欧氏球，再用单约束、lossless S-lemma，可得到一个尺寸 \(2n+p+k\) 的有限 SDP。因而没有依据对这个特例声称 NP-hard；matrix-cube relaxation 也没有必要。

第 4 点是一个有用的 exact reduction，但不能据此直接声称原创。El Ghaoui 已系统建立 uncertain inverse 的 common approximate inverse、LFR、elimination 和 S-procedure/SDP 框架，并明确区分 structured 上界与 unstructured exact 情形。[^1] 本文的精确 LMI 利用了更窄的“自伴 Loewner interval + operator norm + locality-constrained \([Q,0]\)”结构；它至少在写法上不等于照抄其 Eq. (6.5)，但很可能应被视作 lossless full-block S-procedure 的专门化，正式 novelty 判断仍需针对 robust-control 文献查重。

## 1. 变量替换与退化情形

令

\[
S=\Sigma^{-1},\qquad \alpha=\frac1M,\quad \beta=\frac1m,
\quad h=\frac{\alpha+\beta}{2},\quad d=\frac{\beta-\alpha}{2}.
\]

逆映射把不确定集双射为

\[
\mathcal S_{\alpha,\beta}
=\{S=S^T:\alpha I\preceq S\preceq\beta I\}.
\]

若 \(C=0\)，取 \(Q=0\) 即有 \(R=0\)。若 \(m=M\)，completion 唯一，问题就是一个普通 operator-norm best approximation；后文连续窗的 S-lemma 证明依赖 \(d>0\)，所以该情形应单独处理，而不是在无 Slater 点时硬套 S-lemma。

## 2. Scalar cut 闭式核验

设输出和 cut 状态均为一维，而 \(F=f u^T\)，其中 \(f=\|F\|_2\ge0\)、
\(\|u\|_2=1\)，\(C=c\in\mathbb R\)。问题变成

\[
\inf_q\sup_{s\in[\alpha,\beta]}
\sqrt{\|csF-q\|_2^2+c^2s^2}.
\]

### 命题 1（scalar cut 的精确值）

若 \(c=0\)，则 \(R=0\)。若 \(c\ne0\)，则

\[
R=|c|\begin{cases}
\beta,
&f^2\le \dfrac{\alpha+\beta}{\beta-\alpha},\\[2mm]
\sqrt{\displaystyle
\beta^2+\frac{[f^2(\beta-\alpha)-(\alpha+\beta)]^2}{4f^2}},
&f^2\ge \dfrac{\alpha+\beta}{\beta-\alpha}.
\end{cases}
\]

第一支可取

\[
q^*=c\beta F,
\]

第二支可取

\[
q^*=c\frac{(f^2+1)(\alpha+\beta)}{2f^2}F.
\]

当 \(f=0\) 时落在第一支，取 \(q^*=0\)。两支在阈值处连续且给出同一最优解。

### 证明

把 \(q\) 分解到 \(\operatorname{span}\{u^T\}\) 及其正交补。正交补只会给每个场景增加相同的平方误差，因此最优 \(q\) 必平行于 \(F\)。除去 \(|c|\) 并吸收 \(c\) 的符号，得到

\[
\inf_{\nu\in\mathbb R}\max_{s\in[\alpha,\beta]}
\{(fs-\nu)^2+s^2\}.
\]

固定 \(\nu\) 后，这是 \(s\) 的凸二次函数，所以最大值必在
\(\alpha,\beta\) 之一。若取 \(\nu=f\beta\)，\(\beta\) 端损失为 \(\beta^2\)，
\(\alpha\) 端损失为

\[
f^2(\beta-\alpha)^2+\alpha^2.
\]

后者不超过前者恰当且仅当
\(f^2\le(\alpha+\beta)/(\beta-\alpha)\)。另一方面，不可见块
\(-cs\) 在 \(s=\beta\) 已给出不可突破的下界 \(|c|\beta\)，故第一支最优。

超过阈值后，最优点必须平衡两个端点。解两端损失相等得到

\[
\nu^*=\frac{(f^2+1)(\alpha+\beta)}{2f}.
\]

此时两个 active quadratic 关于 \(\nu\) 的导数异号，故零属于 active subgradient 的凸包；该 equalizer 是全局 minimax，而非只是一处交点。代回 \(\beta\) 端即得第二支。证毕。

## 3. Loewner interval 的极点：谱端点够，两个标量端点不够

### 命题 2（极点归约）

若 \(\alpha<\beta\)，则

\[
\operatorname{ext}(\mathcal S_{\alpha,\beta})
=\{\alpha I+(\beta-\alpha)P:P=P^T=P^2\}.
\]

因此对每个固定 \(Q\)，

\[
\sup_{S\in\mathcal S_{\alpha,\beta}}
\|[CSF-Q,-CS]\|_2
=\sup_{P=P^T=P^2}
\|[C(\alpha I+(\beta-\alpha)P)F-Q,
-C(\alpha I+(\beta-\alpha)P)]\|_2.
\]

换回 \(\Sigma\)，可以选一个最坏 completion，使其全部特征值都属于
\(\{m,M\}\)。但投影 \(P\) 的秩和方向一般不能预先固定。

### 证明

仿射缩放后只需考虑 \(0\preceq K\preceq I\)。若 \(K\) 有一个位于
\((0,1)\) 的特征值，就可沿相应特征投影作正负小扰动，故它不是极点。反之，若 \(K=P\) 是投影且
\(P=(K_1+K_2)/2\)、\(0\preceq K_j\preceq I\)，分别在
\(\ker P\) 和 \(\operatorname{ran}P\) 上取二次型，正定序迫使两个
\(K_j\) 在这两个子空间及交叉块上都与 \(P\) 相同，所以 \(P\) 是极点。

固定 \(Q\) 时目标是 \(S\) 的仿射函数的范数，因此关于 \(S\) 连续且凸；有限维紧凸集上的凸函数至少有一个极点最大化者。证毕。

这条命题否定“需要搜索内部特征值”，但没有把不确定集变成有限集合：固定秩的投影构成 Grassmann 流形，方向仍连续。

## 4. 优化后仍成立的 2×2 端点反例

取

\[
m=1,\quad M=2,\quad
\alpha=\frac12,\quad\beta=1,
\quad C=[1\ 0],\quad F=[0\ 4]^T,
\]

且 \(Q=q\) 为标量。

若只检查 \(S=\alpha I,\beta I\)，则 \(CSF=0\)，所以

\[
\inf_q\max\{\sqrt{q^2+1/4},\sqrt{q^2+1}\}=1,
\]

唯一最优中心为 \(q=0\)。

对完整 Loewner interval，令 \(D=\operatorname{diag}(1,-1)\)。映射
\(S\mapsto DSD\) 保持不确定集和 \(\|CS\|\)，同时把 \(CSF\) 变号。因此最坏损失关于 \(q\) 是偶凸函数，仍在 \(q=0\) 处最小。

考虑 rank-one 极点

\[
S=\frac12I+\frac12uu^T,qquad
u=(\cos\theta,\sin\theta)^T,qquad y=\cos^2\theta.
\]

在 \(q=0\) 时，平方损失恰为

\[
16S_{12}^2+\|[S_{11},S_{12}]\|_2^2
=\frac14+\frac34y+4y(1-y)
=\frac14+\frac{19}{4}y-4y^2.
\]

它在 \(y=19/32\) 达到 \(425/256\)。故

\[
\boxed{R=\frac{\sqrt{425}}{16}=\frac{5\sqrt{17}}{16}
\approx1.288470508>1.}
\]

对应 \(S\) 的特征值确实是 \(\{1/2,1\}\)，故对应 \(\Sigma=S^{-1}\) 的特征值是 \(\{1,2\}=\{m,M\}\)。失败的不是“谱端点原则”，而是“只查两个各向同性端点”。复现脚本
`experiments/theory_search/block_boundary_counterexample.py` 同时检查了解析值和后文 SDP 的 \(\lambda=17\) 证书，最大 LMI 特征值约 \(7.3\times10^{-16}\)。

## 5. 有限 completion：epigraph SDP、dual 与 KKT

给定有限场景 \(S_1,\ldots,S_N\)。令

\[
r=k+n,\qquad
A_j^0=[CS_jF,-CS_j]\in\mathbb R^{p\times r},\qquad
P_0=[I_k\ 0]\in\mathbb R^{k\times r}.
\]

则 \(A_j(Q)=A_j^0-QP_0\)。精确 primal epigraph SDP 为

\[
\begin{array}{ll}
\text{minimize}_{t,Q}&t\\
\text{subject to}&
Z_j(t,Q):=
\begin{bmatrix}
tI_p&A_j^0-QP_0\\
(A_j^0-QP_0)^T&tI_r
\end{bmatrix}\succeq0,
\quad j=1,\ldots,N.
\end{array}
\tag{P-fin}
\]

该 LMI 与 \(\|A_j(Q)\|_2\le t\) 等价。取任意 \(Q\) 和充分大的 \(t\) 即有严格可行点，所以 Slater 条件成立。

采用 Frobenius 内积，令 dual 变量按如下符号分块：

\[
Y_j=\begin{bmatrix}U_j&-W_j\\-W_j^T&V_j\end{bmatrix}\succeq0,
\quad W_j\in\mathbb R^{p\times r}.
\]

精确 dual 是

\[
\begin{array}{ll}
\text{maximize}&2\displaystyle\sum_{j=1}^N\langle W_j,A_j^0\rangle\\
\text{subject to}&
Y_j\succeq0,\quad j=1,\ldots,N,\\
&\displaystyle\sum_{j=1}^N(\operatorname{tr}U_j+\operatorname{tr}V_j)=1,\\
&\displaystyle\sum_{j=1}^N W_jP_0^T=0.
\end{array}
\tag{D-fin}
\]

符号可直接由

\[
\langle Y_j,Z_j\rangle
=t(\operatorname{tr}U_j+\operatorname{tr}V_j)
-2\langle W_j,A_j^0-QP_0\rangle
\]

核验。Primal Slater 给出强对偶和 dual attainment。

一组变量最优当且仅当满足：

1. primal feasibility：\(Z_j(t,Q)\succeq0\)；
2. dual feasibility：上述 PSD、迹归一和 \(Q\)-stationarity 条件；
3. complementarity：
   \(\langle Y_j,Z_j(t,Q)\rangle=0\) 对每个 \(j\)。由于两者 PSD，这也等价于
   \(Y_jZ_j(t,Q)=0\)。

非零 \(Y_j\) 标识支撑最优值的 active completions；
\(\sum_jW_jP_0^T=0\) 是共同 \(Q\) 的平衡条件。若再给 \(Q\) 加线性局部性/稀疏约束，primal 仍是 SDP，但 dual stationarity 要加入该仿射可行空间的法锥，不能继续原样照抄最后一式。

## 6. 连续谱窗其实可以 exact 有限化

### 引理 3（矩阵区间作用于单向量的像是球）

对任意 \(v\in\mathbb R^n\)，

\[
\{Sv:\alpha I\preceq S\preceq\beta I\}
=\{z:\|z-hv\|_2\le d\|v\|_2\}.
\tag{1}
\]

### 证明

写 \(S=hI+dK\)，则 \(K=K^T\)、\(\|K\|_2\le1\)，故右侧包含关系立即成立。
反过来，写 \(z=hv+dw\)、\(\|w\|\le\|v\|\)。若 \(w=0\)，取 \(K=0\)。否则令
\(r=\|w\|/\|v\|\le1\)，选择一个对称正交 Householder 矩阵 \(H\)，使
\(H(v/\|v\|)=w/\|w\|\)，再取 \(K=rH\)。则 \(K\) 是自伴收缩且
\(Kv=w\)，从而构造出所需 \(S\)。证毕。

这与自伴 contraction completion 的经典几何一致；Davis–Kahan–Weinberger 给出了更一般的 norm-preserving block completion 公式。[^2] 但式 (1) 本身已有上述一行 Householder 构造，不需要调用更重的定理。

### 定理 4（连续 Loewner 谱窗的 exact SDP）

假设 \(m<M\) 且 \(C\ne0\)。定义

\[
D=\begin{bmatrix}
-I_n&hC^T\\
hC&-\alpha\beta CC^T
\end{bmatrix},\qquad
E=\begin{bmatrix}0_{n\times n}&0\\0&I_p\end{bmatrix},
\]

以及

\[
L_Q=
\begin{bmatrix}
F^T&-Q^T\\
I_n&0
\end{bmatrix}
\in\mathbb R^{(k+n)\times(n+p)}.
\]

则

\[
\boxed{
R(C,F;m,M)^2=
\min_{\tau,Q,\lambda}\ \tau
}
\]

满足

\[
\tau\ge0,\qquad\lambda\ge0,qquad
\boxed{
\begin{bmatrix}
-\tau E+\lambda D&L_Q^T\\
L_Q&-I_{k+n}
\end{bmatrix}\preceq0.}
\tag{SDP-cont}
\]

该 LMI 的阶数为 \(2n+p+k\)，所有变量均仿射出现。任何额外的仿射
\(Q\)-约束都可以原样加入而不破坏 exactness。

### 证明

固定 \(S,Q\)，谱范数可由左奇异向量写成

\[
\|[CSF-Q,-CS]\|_2^2
=\sup_{\|u\|=1}
\bigl(\|F^TSC^Tu-Q^Tu\|^2+\|SC^Tu\|^2\bigr).
\]

两个 supremum 可以交换。对固定 \(u\)，令 \(v=C^Tu\)、\(z=Sv\)。由引理 3，

\[
\|z-hC^Tu\|^2\le d^2\|C^Tu\|^2.
\]

令 \(w=[z;u]\)。上式等价于

\[
w^TDw\ge0,
\]

因为 \(d^2-h^2=-\alpha\beta\)。而待控制的平方输出正是
\(w^TL_Q^TL_Qw\)。由于约束齐次，除零向量外的 feasible \(w\) 均有
\(u\ne0\)，所以 robust 上界 \(\sqrt\tau\) 等价于二次型蕴含

\[
w^TDw\ge0
\quad\Longrightarrow\quad
w^T(L_Q^TL_Q-\tau E)w\le0.
\tag{2}
\]

当 \(d>0,C\ne0\) 时存在 \(w\) 使 \(w^TDw>0\)：取
\(C^Tu\ne0\) 且 \(z=hC^Tu\)。因此单二次约束的 S-lemma 在 Slater 条件下是 lossless，[^3]
式 (2) 当且仅当存在 \(\lambda\ge0\) 使

\[
L_Q^TL_Q-\tau E+\lambda D\preceq0.
\]

对 \(-I_{k+n}\) 作 Schur complement，恰得到 (SDP-cont)。最后对
\(Q,\tau\) 最小化。证毕。

这个证明还解释了为何没有 quantifier 偷换：
\(\sup_S\sup_u=\sup_{(S,u)}=\sup_u\sup_S\)，而每个 \(u\) 下允许最坏
\(S\) 随 \(u\) 变化，正是原联合 supremum 的含义。

## 7. 计算与最坏 completion 恢复

有限场景问题用 (P-fin) 即可。连续谱窗应直接解 (SDP-cont)，而不是枚举
Grassmann 流形或采用 matrix-cube relaxation。标准 SDP 求解器给出
\(Q^*,\tau^*,\lambda^*\)，值为 \(R=\sqrt{\tau^*}\)。

若还要恢复最坏 completion，可在固定 \(Q^*\) 后求齐次 generalized trust-region 问题

\[
\sup_{u,z}
\frac{\|F^Tz-Q^{*T}u\|^2+\|z\|^2}{\|u\|^2}
\quad\text{s.t.}\quad
\|z-hC^Tu\|\le d\|C^Tu\|.
\]

单二次约束使它可由同一 S-lemma/一维 multiplier 搜索求解。得到最坏
\((u,z)\) 后，令 \(v=C^Tu\)、\(w_0=(z-hv)/d\)，用引理 3 的 Householder 构造
\(K\)，再取 \(S=hI+dK\)、\(\Sigma=S^{-1}\)。最大化凸二次函数时可选球面点，因而可选对称正交 \(K\)，这与命题 2 的 endpoint-spectrum 极点一致。
若最坏方向恰有 \(C^Tu=0\)，则必有 \(z=0\)，此时该方向的损失与
completion 无关，任选一个 endpoint-spectrum completion 即可；不应在
\(v=0\) 上使用 Householder 除法。

目前没有得到一般 \(C,F\) 下 \(Q^*\) 的类似 scalar-cut 分段闭式；但这不再是“值不可计算”的缺口，而只是“能否解析消去 SDP”的问题。最小未解 lemma 可表述为：

> 刻画哪些 \((C,F)\) 允许 (SDP-cont) 的最优 \(Q\) 由 \(C,F,\alpha,\beta\) 的有限个谱投影闭式表示，并给出最坏投影的秩与方向。

在没有这条闭式前，(SDP-cont) 已是精确、多项式规模的可计算方案。不能再把本问题描述为一般 NP-hard，也不能把 structured robust SDP 的保守上界当作唯一方案。

## 8. 与 El Ghaoui 2002 的实质比较

El Ghaoui 的 Eq. (1.4) 已定义“选择一个共同矩阵去最小化不确定逆族的最大 operator-norm 误差”；Theorem 6.2/Eq. (6.3) 通过 LFR 和 multiplier SDP 优化 structured inversion-error 上界；只有 unstructured full perturbation 才明确声称条件必要，并在 Eq. (6.5) 给出解析 approximate inverse。[^1] 其摘要也明确区分“一般 structured SDP bounds”和“unstructured exact analytic expression”。

本问题与其重合之处是核心 min–max approximate-inverse 思想、operator norm、LFR/elimination/S-procedure 技术路线。以下内容因此不能作为新贡献：

- 提出 common robust inverse/center；
- 对有限场景写 operator-norm SDP；
- 对一般 structured uncertainty 使用 SDP multiplier 上界；
- 仅仅说“用 S-procedure 处理不确定逆”。

本问题的额外结构是：

- 只观察选定行 \(CS\)，并含耦合块 \([CSF,-CS]\)；
- approximant 被局部通信强制为 \([Q,0]\)，并非自由完整 inverse；
- uncertainty 不是一般 additive/LFR block，而是自伴正定矩阵的 Loewner 谱窗；
- 这个谱窗对单个左测试向量的像恰为一个球，从而只剩一个 lossless quadratic constraint。

因此 (SDP-cont) 比“套用 El Ghaoui structured upper bound”更精确，也不由其 unstructured Eq. (6.5) 直接给出。然而它的证明本质仍是 full-block 几何加 lossless S-lemma；El Ghaoui–Oustry–Lebret 对 uncertain SDP 的工作也明确指出 full uncertainty 下某些 robust 条件可以 necessary and sufficient。[^4] 在完成针对 full-block robust performance、Douglas/Parrott/Davis–Kahan–Weinberger completion 和 lossless S-procedure 文献的逐式查重前，最稳妥状态是：

> **数学上 exact，原创性未确认；很可能是已知 lossless full-block machinery 在 locality-constrained Schur window 上的一个干净专门化。**

## 9. 最终判定

- **Scalar cut：通过。** 原分段式、阈值和最优系数正确。
- **只查 \(mI,MI\)：否决。** 解析 \(2\times2\) 反例给出严格 28.8% gap。
- **只查 endpoint spectra：通过。** 最坏 \(\Sigma\) 可选全部特征值在
  \(\{m,M\}\)，但仍需连续优化投影方向。
- **有限 completion SDP/dual/KKT：通过。** 上述 primal、dual 和 stationarity/complementarity 量词一致。
- **连续谱窗只能 semi-infinite：否决。** 当前目标可 exact 化为 (SDP-cont)。
- **NP-hard：没有依据。** 对当前单 Loewner interval 和 spectral norm，已有多项式规模 exact SDP；只有改变 uncertainty 结构、引入多个独立二次约束或非凸离散限制后，才可能进入 matrix-cube/NP-hard 领域。
- **原创性：未通过。** exact reduction 值得保留，但在 robust-control/full-block S-procedure 谱系中高度可疑，不能因“未见同一公式”就宣称新定理。

## Sources

[^1]: Laurent El Ghaoui, “[Inversion Error, Condition Number, and Approximate Inverses of Uncertain Matrices](https://doi.org/10.1016/S0024-3795(01)00273-7),” *Linear Algebra and its Applications* 343–344 (2002), 171–193. [Author publication page](https://people.eecs.berkeley.edu/~elghaoui/pubs_inv_err.html).

[^2]: Chandler Davis, W. M. Kahan, and H. F. Weinberger, “[Norm-Preserving Dilations and Their Applications to Optimal Error Bounds](https://doi.org/10.1137/0719029),” *SIAM Journal on Numerical Analysis* 19(3) (1982), 445–469.

[^3]: Imre Pólik and Tamás Terlaky, “[A Survey of the S-Lemma](https://doi.org/10.1137/S003614450444614X),” *SIAM Review* 49(3) (2007), 371–418. The one-constraint result under a strict feasible point is lossless; the report uses precisely that case, not a multi-constraint relaxation.

[^4]: Laurent El Ghaoui, François Oustry, and Hervé Lebret, “[Robust Solutions to Uncertain Semidefinite Programs](https://doi.org/10.1137/S1052623496305717),” *SIAM Journal on Optimization* 9(1) (1998), 33–52. The paper distinguishes general sufficient robust SDP conditions from full-uncertainty cases where they become necessary and sufficient.
