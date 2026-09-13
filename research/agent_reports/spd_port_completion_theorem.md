# 任意 SPD 端口基准下的隐藏图 completion：定理、边界与原创性审计

## 0. 结论先行

本报告审计的推广在下述精确量词下是**正确的**：端口基准可以从正对角
grounding 放宽为任意

\[
\Lambda\succ0,
\]

同时允许任意有限隐藏节点、任意正对角隐藏基准、任意无向含环图和任意非负边权。
若

\[
H_Z=\operatorname{diag}(\Lambda,D_Z),\qquad
T_G=H_Z^{1/2}(H_Z+L_G)^{-1}H_Z^{1/2},
\]

则全部端口压缩 \([T_G]_{\Gamma\Gamma}\) 的闭凸包，恰好等于
\(B_{q+1}\) 个由 partial partitions 索引的正交投影的凸包。对完整隐藏输入的
direct-smoother 误差，Schatten--\(p\) 精确有限归约对所有实例成立**当且仅当**
\(p\ge2\)；这里 \(p=\infty\) 即谱范数。每个右侧场景均可由最大度数不超过
2 的路径森林逼近。

这个推广真正新增的代数困难是：\(\Lambda\) 非对角时，不同图分量在白化坐标中
不再正交。候选证明中提出的第二层 DPP/volume-sampling 分解恰好修复了这一点：

\[
X=V(V^TV+\operatorname{diag}h)^{-1}V^T
\]

仍是列子集空间投影的凸组合。\(h_j=0\) 时通过 \(h_j\downarrow0\) 和有限凸包
闭性处理，不需要假设每个分量都含隐藏节点。

不过有三个必须同时写进论文的限制。

1. 若 \(Q\) 完全自由，robust direct-smoother 设计退化为
   \(Q=C/2\)，最优值恒为 \(\lVert C\rVert_2/2\)。实际价值来自
   sparsity、通信半径、共享 opcode 参数、标称精确性等真实约束，而不是来自无约束
   \(Q\)。
2. 任意 SPD \(\Lambda\) 的自然语义是小端口集上的相关 fidelity/information
   geometry 或已经消元的 Schur base；它一般不是逐节点独立 grounding，也不自动
   保持传统 raw-coordinate graph-Tikhonov 的输入语义。
3. 两层 projection mixture 的数学工具已有明确先例。当前未检索到把“任意 SPD
   端口 base + 任意隐藏 graph completion + partial-partition convex hull + robust
   local smoother + sharp Schatten threshold”合在一起的原始论文，但不能把 DPP
   mean-projection、random-forest representation 或 infinite-conductance shorting
   本身宣传为原创。

我的判决是：

> **数学上 GO，且显著扩展 diagonal-grounding 主定理；最可信的论文定位是一个强的
> robust graph-smoother / structured inverse-completion 定理，而不是声称发现新的
> matrix-forest theorem。纯数学原创强度目前是中等，应用数学组合的价值更高。**

---

## 1. 精确模型

固定端口集合

\[
\Gamma=\{1,\ldots,q\},\qquad \Lambda\in\mathbb S_{++}^q.
\]

一次 completion 可以选择任意有限集合 \(Z\)（允许为空）、任意

\[
D_Z=\operatorname{diag}(d_z:z\in Z)\succ0,
\]

以及顶点集 \(V=\Gamma\mathbin{\dot\cup}Z\) 上的任意有限无向加权图 \(G\)。图可以
含环、不连通；边权有限且非负。令

\[
H=\begin{bmatrix}\Lambda&0\\0&D_Z\end{bmatrix},\qquad
E_\Gamma=[I_q\ 0],
\tag{1.1}
\]

并定义 symmetric normalized resolvent

\[
T_G=H^{1/2}(H+L_G)^{-1}H^{1/2}
\quad\text{及}\quad
X_G=E_\Gamma T_GE_\Gamma^T.
\tag{1.2}
\]

因为

\[
T_G=(I+H^{-1/2}L_GH^{-1/2})^{-1},
\]

所以 \(0\prec T_G\preceq I\)。这里的平方根始终取 symmetric principal square
root；由于 \(H\) 是 block diagonal，

\[
H^{1/2}=\operatorname{diag}(\Lambda^{1/2},D_Z^{1/2}).
\]

### 1.1 SPD partial-partition projections

取任意 active subset \(A\subseteq\Gamma\)，再取 \(A\) 的 set partition

\[
\pi=\{B_1,\ldots,B_k\}.
\]

令 \(U_\pi\in\{0,1\}^{q\times k}\) 的第 \(j\) 列为 \(\mathbf1_{B_j}\)，并定义

\[
P_{\pi,A}^{\Lambda}
=\Lambda^{1/2}U_\pi
(U_\pi^T\Lambda U_\pi)^{-1}
U_\pi^T\Lambda^{1/2}.
\tag{1.3}
\]

当 \(A=\varnothing\) 时约定 \(P_{\varnothing}^{\Lambda}=0\)。记所有这些矩阵的
集合为 \(\mathcal P_q^\partial(\Lambda)\)。

注意：当 \(\Lambda\) 非对角时，不能把 (1.3) 错写成各 block 的 rank-one
projector 之和。向量

\[
\Lambda^{1/2}\mathbf1_{B_1},\ldots,
\Lambda^{1/2}\mathbf1_{B_k}
\]

一般不正交；(1.3) 是它们**联合列空间**上的 Euclidean orthogonal projection。

每个 partial partition 对应一个不同子空间。因为 \(\Lambda^{1/2}\) 可逆，不同
indicator subspaces 在乘以 \(\Lambda^{1/2}\) 后仍不同。因此

\[
|\mathcal P_q^\partial(\Lambda)|
=\sum_{a=0}^q{q\choose a}B_a
=B_{q+1}.
\tag{1.4}
\]

把 inactive ports 与一个 cemetery symbol 放在同一 block 可得到 (1.4) 的双射。
这只是索引双射；并不是把普通 \((q+1)\)-point partition matrix 直接取主子块。

两个场景永远存在：

\[
0\in\mathcal P_q^\partial(\Lambda),\qquad
I_q\in\mathcal P_q^\partial(\Lambda).
\tag{1.5}
\]

第二个等式取全部端口 active 且每个端口为 singleton，此时 \(U_\pi=I_q\)，故

\[
\Lambda^{1/2}\Lambda^{-1}\Lambda^{1/2}=I_q.
\]

---

## 2. 主定理包

### 定理 2.1：端口压缩的精确闭凸包

令 \(\mathfrak G(\Lambda)\) 为第 1 节中全部有限 completion 的并。则

\[
\boxed{
\overline{\operatorname{conv}}
\{X_G:G\in\mathfrak G(\Lambda)\}
=\operatorname{conv}\mathcal P_q^\partial(\Lambda).}
\tag{2.1}
\]

而且：

- 每个 \(P\in\mathcal P_q^\partial(\Lambda)\) 都是右侧凸包的 exposed vertex；
- 每个 \(P\) 都是某个 maximum-degree-2 path-forest completion 序列的极限；
- 若坚持每个 completion 连通，可用趋于零的弱边把各条路径首尾连接，仍保留
  maximum degree 2 和同一极限；
- 对任意连续凸函数 \(\Phi:\mathbb S^q\to\mathbb R\)，

  \[
  \boxed{
  \sup_{G\in\mathfrak G(\Lambda)}\Phi(X_G)
  =\max_{P\in\mathcal P_q^\partial(\Lambda)}\Phi(P).}
  \tag{2.2}
  \]

这里必须写 supremum，因为非零 partial projections 通常需要边权趋于无穷，
而 inactive ports 还需要隐藏 grounding 总量趋于无穷；有限参数图未必达到极值。

### 定理 2.2：完整隐藏输入的 robust smoother

给定 \(C\in\mathbb R^{r\times q}\) 和
\(Q\in\mathbb R^{r\times q}\)，定义保留全部 hidden-input columns 的误差算子

\[
\mathcal A_G(Q)
=CE_\Gamma T_G-QE_\Gamma
\in\mathbb R^{r\times(q+|Z|)}.
\tag{2.3}
\]

若 \(\psi:\mathbb S_+^r\to\mathbb R\) 连续、凸且关于 Loewner 序单调不减，则

\[
\boxed{
\sup_{G\in\mathfrak G(\Lambda)}
\psi(\mathcal A_G(Q)\mathcal A_G(Q)^T)
=\max_{P\in\mathcal P_q^\partial(\Lambda)}
\psi((CP-Q)(CP-Q)^T).}
\tag{2.4}
\]

因此，对 \(2\le p<\infty\)，

\[
\boxed{
\sup_G\|\mathcal A_G(Q)\|_{S_p}
=\max_{P\in\mathcal P_q^\partial(\Lambda)}\|CP-Q\|_{S_p},}
\tag{2.5}
\]

并且 (2.5) 对 \(p=\infty\) 也成立。特别地：

- \(p=2\)：Frobenius norm；
- \(p=\infty\)：spectral/operator norm。

这个 Schatten 范围是 sharp 的。

### 定理 2.3：sharp Schatten 分类

在 Schatten norms 的通常范围 \(1\le p\le\infty\) 内，full-input universal finite
reduction (2.5) 对所有 \((q,\Lambda,C,Q)\) 和所有 completion 成立，当且仅当

\[
\boxed{p\ge2.}
\tag{2.6}
\]

“仅当”不是 proof-technique 的缺口：第 5 节给出同一个两节点例子，对**每一个**
\(1\le p<2\) 都严格违反 (2.5)。

---

## 3. 定理 2.1 的逐式证明

### 3.1 第一层：任意 SPD base 仍有 forest projection mixture

给图任取一个 oriented incidence matrix \(B\)，令

\[
L_G=BWB^T,
\qquad
A=H^{-1/2}BW^{1/2}.
\tag{3.1}
\]

于是

\[
T_G=(I+AA^T)^{-1}.
\tag{3.2}
\]

在 edge subsets 上取 L-ensemble

\[
\Pr(S)
=\frac{\det(A_S^TA_S)}{\det(I+A^TA)}.
\tag{3.3}
\]

一般 selected-column mean-projection identity 给出

\[
\mathbb E P_{\operatorname{col}(A_S)}
=A(I+A^TA)^{-1}A^T.
\tag{3.4}
\]

由 Woodbury identity，

\[
I-A(I+A^TA)^{-1}A^T=(I+AA^T)^{-1},
\]

故

\[
\boxed{
T_G=\mathbb E P_{\ker A_S^T}.}
\tag{3.5}
\]

\(H^{-1/2}\) 可逆，所以它不改变 incidence columns 的线性依赖关系。
\(\det(A_S^TA_S)>0\) 当且仅当所选 edges 线性独立，即 \(S\) 是 forest。因此
(3.5) 确实是 forest projections 的混合，尽管 \(H\) 不再是 diagonal matrix。

这一步是一般 DPP/volume-sampling 的已知结果，而不是本项目可认领的原创部分。

### 3.2 一个 forest projection 的端口块

固定一片 forest \(F\)。若它的 component indicators 构成矩阵 \(U_F\)，则

\[
\ker A_F^T
=\operatorname{col}(H^{1/2}U_F).
\tag{3.6}
\]

证明很直接：\(A_F^Tx=0\) 等价于
\(B_F^TH^{-1/2}x=0\)，即 \(H^{-1/2}x\) 在 forest 的每个 component 内为常数。
故对应 orthogonal projection 是

\[
\widehat P_F
=H^{1/2}U_F(U_F^THU_F)^{-1}U_F^TH^{1/2}.
\tag{3.7}
\]

只保留含端口的 components。设这些 components 在端口上的非空交集给出一个
full partition

\[
\pi=\{B_1,\ldots,B_k\}
\]

并令

\[
U=[\mathbf1_{B_1}\ \cdots\ \mathbf1_{B_k}],
\quad
V=\Lambda^{1/2}U,
\quad
h_j=\sum_{z\in C_j\cap Z}d_z.
\tag{3.8}
\]

hidden-only components 与端口块正交且没有端口贡献。对 port-containing
components，隐藏部分的 supports 互不相交，因此

\[
U_F^THU_F
=V^TV+\operatorname{diag}(h_1,\ldots,h_k).
\tag{3.9}
\]

虽然隐藏贡献仍是 diagonal，\(V^TV=U^T\Lambda U\) 一般不是 diagonal。
由 (3.7)，forest projection 的端口压缩恰为

\[
\boxed{
X_F
=V(V^TV+\operatorname{diag}h)^{-1}V^T.}
\tag{3.10}
\]

这一步明确保留了 \(\Lambda\) 造成的 component cross terms；若逐 block 分解，证明
会在这里出错。

### 3.3 第二层：非正交 components 的 column-subset DPP

先假设所有 \(h_j>0\)，并写

\[
D_h=\operatorname{diag}h,
\qquad M=VD_h^{-1/2}.
\]

则

\[
X_F=M(I+M^TM)^{-1}M^T.
\tag{3.11}
\]

再次在 \(M\) 的列子集上用 L-ensemble：

\[
\Pr(J)
=\frac{\det(M_J^TM_J)}{\det(I+M^TM)}.
\tag{3.12}
\]

由同一个 mean-projection identity，

\[
\boxed{
X_F=\mathbb E_JP_{\operatorname{col}(M_J)}
=\mathbb E_JP_{\operatorname{col}(V_J)}.}
\tag{3.13}
\]

第二个等号利用了正的列尺度不改变列空间。每个列子集 \(J\) 就是从 full partition
\(\pi\) 中选择一部分 active blocks；其 projection 正是

\[
P_{\pi_J,A_J}^{\Lambda}
=V_J(V_J^TV_J)^{-1}V_J^T.
\]

因此

\[
X_F\in\operatorname{conv}\mathcal P_q^\partial(\Lambda).
\tag{3.14}
\]

#### \(h_j=0\) 不构成漏洞

若某个 component 没有 hidden vertex，则相应 \(h_j=0\)。因为 \(U\) 的列线性
独立且 \(\Lambda\succ0\)，矩阵 \(V\) 满列秩，所以即使部分或全部 \(h_j=0\)，

\[
V^TV+\operatorname{diag}h\succ0.
\]

取

\[
h_j^{(\varepsilon)}=
\begin{cases}
h_j,&h_j>0,\\
\varepsilon,&h_j=0,
\end{cases}
\]

则 (3.13) 对每个 \(\varepsilon>0\) 成立，且

\[
V(V^TV+\operatorname{diag}h^{(\varepsilon)})^{-1}V^T
\longrightarrow X_F.
\]

右侧 scenario family 有限，故其 convex hull 是 compact/closed；所以极限仍在
同一个凸包中。这完整处理了 \(h_j=0\)，无需引入伪逆或假设每个端口 component
都有隐藏 grounding。

### 3.4 对 forest mixture 再取期望

由 (3.5)，

\[
X_G=\mathbb E_F X_F.
\]

结合 (3.14)，对每个有限 completion 都有更强的 pointwise inclusion：

\[
\boxed{X_G\in\operatorname{conv}\mathcal P_q^\partial(\Lambda).}
\tag{3.15}
\]

因此 (2.1) 的“\(\subseteq\)”成立；同时 (2.2) 的上界由 Jensen inequality 立即
得到。

### 3.5 sharp converse 与 varying dimension

固定一个 partial partition

\[
\pi=\{B_1,\ldots,B_k\}
\quad\text{of}\quad A\subseteq\Gamma.
\]

构造如下 forest。

- 对每个 active block \(B_j\)，用一条经过其中所有端口的 terminal path 把它们
  连起来。
- 对每个 inactive port \(i\in\Gamma\setminus A\)，加一个 hidden vertex \(z_i\)，
  并用一条边连接 \(i\) 与 \(z_i\)。
- 令每个 \(z_i\) 的 grounding 等于 \(s\)，先令该 forest 的全部 edge
  conductances 趋于无穷，再令 \(s\to\infty\)。

每个图都是 maximum-degree-2 forest；事实上 inactive components 只是单边路径。

令

\[
V_A=\Lambda^{1/2}U_\pi,
\qquad
V_I=\Lambda^{1/2}
[e_i:i\in\Gamma\setminus A].
\]

无限导通后的 full component projection 具有 Gram matrix

\[
G_s=
\begin{bmatrix}
V_A^TV_A&V_A^TV_I\\
V_I^TV_A&V_I^TV_I+sI
\end{bmatrix}.
\tag{3.16}
\]

标准 block inverse estimate 给出

\[
[G_s^{-1}]_{AA}\to(V_A^TV_A)^{-1},\qquad
[G_s^{-1}]_{AI}=O(s^{-1}),\qquad
[G_s^{-1}]_{II}=O(s^{-1}).
\tag{3.17}
\]

故其完整 port rows 满足

\[
E_\Gamma\widehat P_s
\longrightarrow
\left[
V_A(V_A^TV_A)^{-1}V_A^T,
0
\right]
=\left[P_{\pi,A}^{\Lambda},0\right]
\tag{3.18}
\]

in operator norm。隐藏列趋零是因为 (3.17) 的 inactive columns 还要乘
\(\sqrt s\)，结果为 \(O(s^{-1/2})\)。因此证明没有偷偷丢掉 hidden-input
columns。

对每个固定 \(s\)，若所有 forest edge weights 同乘 \(t\)，则

\[
T_{s,t}\to\widehat P_s\qquad(t\to\infty)
\]

in operator norm。用 diagonal sequence：先选 \(s_m\to\infty\)，再选足够大的
\(t_m\)，即可得到有限参数图序列满足 (3.18)。这也严谨处理了 completion dimension
或 base 随序列改变的问题。

若不允许隐藏 grounding 单点变大，但允许加入任意多个 unit-grounded hidden
vertices，可把每个 inactive edge 换成长 hidden path；总 grounding 仍趋于无穷，
相同 diagonal argument 成立。若要求每个图连通，把各个 path components 用权重
\(\varepsilon_m\downarrow0\) 的边首尾相连即可。

于是每个 scenario 属于 \(\{X_G\}\) 的 closure，得到 (2.1) 的反向 inclusion 和
(2.2) 的下界。

### 3.6 每个 scenario 都是 exposed vertex

固定 \(P\in\mathcal P_q^\partial(\Lambda)\)，定义线性泛函

\[
\ell_P(Y)=\langle2P-I,Y\rangle_F.
\]

对任意另一个 orthogonal projection \(R\)，

\[
\ell_P(P)-\ell_P(R)
=\operatorname{rank}P+\operatorname{rank}R-2\operatorname{tr}(PR)
=\|P-R\|_F^2.
\tag{3.19}
\]

不同 partial partitions 给出不同 projection，故 (3.19) 在 \(R\ne P\) 时严格为正。
所以没有 scenario 可以从 universal exact finite list 中删除。

---

## 4. 完整 smoother 误差的证明

由 (3.5)，写

\[
T_G=\mathbb E\widehat P_F,
\qquad
A_F=CE_\Gamma\widehat P_F-QE_\Gamma.
\]

于是

\[
\mathcal A_G(Q)=\mathbb EA_F
\]

且 matrix variance identity 给出

\[
\mathcal A_G(Q)\mathcal A_G(Q)^T
\preceq\mathbb E[A_FA_F^T].
\tag{4.1}
\]

对一个固定 forest projection，令

\[
X_F=E_\Gamma\widehat P_FE_\Gamma^T.
\]

利用 \(\widehat P_F^2=\widehat P_F\) 和 \(E_\Gamma E_\Gamma^T=I\)，完整 hidden
columns 的 Gram matrix 精确化为

\[
\begin{aligned}
A_FA_F^T
&=CX_FC^T-CX_FQ^T-QX_FC^T+QQ^T\\
&=:M_Q(X_F).
\end{aligned}
\tag{4.2}
\]

关键点是 \(M_Q(X)\) 对 \(X\) **affine**。第 3.3 节给出

\[
X_F=\sum_j\theta_jP_j,
\qquad P_j\in\mathcal P_q^\partial(\Lambda),
\]

所以

\[
M_Q(X_F)=\sum_j\theta_jM_Q(P_j).
\tag{4.3}
\]

对 (4.1) 使用 \(\psi\) 的 Loewner monotonicity 和 convexity，再对 (4.3) 使用
convexity，得到 (2.4) 的上界。对 projection \(P_j^2=P_j\)，

\[
M_Q(P_j)=(CP_j-Q)(CP_j-Q)^T.
\tag{4.4}
\]

第 3.5 节又给出

\[
E_\Gamma T_{G_m}\to[P_j,0],
\]

故每个右侧值都可由左侧逼近，(2.4) 成立。

对有限 \(p\ge2\)，取

\[
\psi_p(M)=\operatorname{tr}(M^{p/2}).
\]

因为 \(p/2\ge1\)，\(\psi_p\) 在 PSD cone 上 convex 且 Loewner monotone。
取 \(p\)-th root 得到 (2.5)。对 \(p=\infty\)，取

\[
\psi_\infty(M)=\lambda_{\max}(M)
\]

并取 square root。Frobenius case 则是 \(\psi_2(M)=\operatorname{tr}M\)。

这里第一层 forest averaging 也不可省略：\(T_G\) 通常不是 projection；(4.1) 正是
把 deterministic resolvent 拉回到 forest projections 的步骤。

---

## 5. \(p<2\) 的统一反例

只需一个 port 和一个 hidden vertex：

\[
q=1,\qquad \Lambda=1,\qquad D_Z=1.
\]

用 conductance \(t\to\infty\) 的单边把二者连接，则

\[
T_t\to\widehat P
=\frac12\begin{bmatrix}1&1\\1&1\end{bmatrix}.
\tag{5.1}
\]

取两维输出

\[
C=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
Q=\begin{bmatrix}1\\0\end{bmatrix}.
\tag{5.2}
\]

一个端口的 partial scenarios 只有 \(P=0\) 与 \(P=1\)，且

\[
\|CP-Q\|_{S_p}=1
\qquad(P=0,1)
\tag{5.3}
\]

对每个 \(p\ge1\) 都成立。另一方面，full hidden-input limiting error 是

\[
C[1/2,1/2]-[Q,0]
=\begin{bmatrix}-1/2&1/2\\1/2&1/2\end{bmatrix}.
\tag{5.4}
\]

它的两个 singular values 都是 \(1/\sqrt2\)，所以

\[
\|\text{(5.4)}\|_{S_p}
=2^{1/p-1/2}.
\tag{5.5}
\]

当且仅当 \(p<2\) 时，(5.5) 严格大于 1。有限 \(t\) 的值连续趋于该极限，故

\[
\sup_G\|\mathcal A_G(Q)\|_{S_p}
>\max_{P\in\{0,1\}}\|CP-Q\|_{S_p}
\qquad(1\le p<2).
\tag{5.6}
\]

这同时否定 nuclear norm 和每个 \(1<p<2\)，并证明 (2.6) 是 iff。注意：
这个反例只针对**完整 hidden-input error**。若损失只依赖端口 compression \(X_G\)，
定理 2.1 对任意连续凸 \(\Phi(X_G)\) 仍然成立。

---

## 6. 优化 \(Q\)：无约束时退化，有 locality 时才有设计价值

### 6.1 自由 \(Q\) 的闭式解

对 spectral norm 考虑

\[
\min_Q\max_{P\in\mathcal P_q^\partial(\Lambda)}\|CP-Q\|_2.
\tag{6.1}
\]

由 (1.5)，\(P=0,I\) 都是场景，因此

\[
\max_P\|CP-Q\|_2
\ge\max\{\|Q\|_2,\|C-Q\|_2\}
\ge\frac12\|C\|_2.
\tag{6.2}
\]

取 \(Q=C/2\)。每个 scenario 都是 orthogonal projection，所以

\[
\|CP-C/2\|_2
=\|C(P-I/2)\|_2
\le\|C\|_2\|P-I/2\|_2
=\frac12\|C\|_2.
\tag{6.3}
\]

故

\[
\boxed{
Q^*=C/2,
\qquad
\min_Q\sup_G\|\mathcal A_G(Q)\|_2=\|C\|_2/2.}
\tag{6.4}
\]

因此不能用“优化了一个自由 dense \(Q\)”来宣传应用收益；那只是两个 endpoint 的
midpoint。

### 6.2 合理的 constrained design 与 exact LMI

令 \(\mathcal Q_{\rm loc}\) 是非空 closed affine/convex feasible set，例如：

- 通信半径给出的固定 zero mask；
- 多个 agent 共用同一组 opcode coefficients；
- 对 nominal topology 的 reproduction/unbiasedness 约束；
- 固定 block locality 或 symmetry tying。

则逐点应用定理 2.2 得

\[
\boxed{
\min_{Q\in\mathcal Q_{\rm loc}}
\sup_G\|\mathcal A_G(Q)\|_2
=\min_{Q\in\mathcal Q_{\rm loc}}
\max_{P\in\mathcal P_q^\partial(\Lambda)}\|CP-Q\|_2.}
\tag{6.5}
\]

因为 \(P=0\) 使目标至少为 \(\|Q\|_2\)，closed feasible set 上存在 minimizer。
给定 \(\tau\ge0\)，(6.5) 的 feasibility 等价于有限组 LMI：

\[
\boxed{
\begin{bmatrix}
\tau I_r&CP-Q\\
(CP-Q)^T&\tau I_q
\end{bmatrix}\succeq0
\quad
\forall P\in\mathcal P_q^\partial(\Lambda).}
\tag{6.6}
\]

scenario 数是 \(B_{q+1}\)，与隐藏节点数、隐藏边数、环数和连续 edge weights
完全无关。这是 boundary-size fixed-parameter exact method；Bell number 增长很快，
不能写成对大 \(q\) 的 polynomial algorithm。

### 6.3 一个 intermediate partition 真正 active 的严格例子

取 \(q=r=2\)、\(\Lambda=I_2\)，令

\[
C=\operatorname{diag}(1,2),
\qquad
Q=\alpha I_2.
\tag{6.7}
\]

\(Q=\alpha I\) 是一个很自然的 shared-opcode constraint：两个端口使用同一标量
coefficient。若错误地删除 connected-block scenario

\[
P_{12}=\frac12\mathbf1\mathbf1^T,
\]

其余四个场景的 minimax 解为

\[
\alpha=1,\qquad R=1.
\tag{6.8}
\]

加入全部五个 partial scenarios 后，精确解变为

\[
\boxed{
\alpha^*=\frac43-\frac{2\sqrt{10}}{15}
=0.9116963119\ldots,}
\tag{6.9}
\]

\[
\boxed{
R^*=\frac{10+2\sqrt{10}}{15}
=1.0883036880\ldots.}
\tag{6.10}
\]

在最优点，\(P_{12}\) 与 singleton/full scenarios 共同 active。代入

\[
CP_{12}-\alpha I
=\begin{bmatrix}1/2-\alpha&1/2\\1&1-\alpha\end{bmatrix}
\]

并令其最大 singular value 等于 \(2-\alpha\)，得到

\[
15\alpha^2-40\alpha+24=0;
\]

取可行的小根即 (6.9)，随后得到 (6.10)。这个例子说明 constrained local design
不能只检查 \(0/I\) endpoints，也不能只检查 singleton partitions。

---

## 7. 应用解释与不能越过的边界

### 7.1 广义 graph-Tikhonov 的自然解释

考虑

\[
\widehat x
=\arg\min_x
\frac12(x-y)^TH(x-y)+\frac12x^TL_Gx.
\tag{7.1}
\]

其解为

\[
\widehat x=(H+L_G)^{-1}Hy.
\]

在 normalized coordinates

\[
u=H^{1/2}y,
\qquad
\widetilde x=H^{1/2}\widehat x,
\]

恰有

\[
\widetilde x=T_Gu.
\tag{7.2}
\]

因此 arbitrary SPD \(\Lambda\) 可以表示端口观测误差的 correlated precision，或
端口上的一般 quadratic fidelity。若 \(q\) 是一个小 boundary/interface，这个模型
是自然的。

但 \(\Lambda^{1/2}\) 通常 dense；计算 normalized port input 会混合全部 \(q\) 个
端口。若 \(q\) 横跨整个大网络，就不能把它称为逐节点 local operation。合理叙述是：
\(\Gamma\) 是已经可访问的小局部区域或 separator，\(\mathcal Q_{\rm loc}\) 再编码
真正的通信限制。

### 7.2 “固定 local core”可以编码到什么程度

若 fixed core 已经消元到端口上的 SPD information/energy matrix \(\Lambda\)，且未知
exterior 仅通过 scalar difference edges 增加 Laplacian \(L_G\)，则本文代数完全适用。
例如 grounded resistor/core graph 的 Kron reduction 给出的端口 Schur matrix 属于
SDDM/M-matrix 子类，是 arbitrary SPD theorem 的一个物理子类。

但是有一个输入端必须区分。标准 smoother 若 central matrix 是

\[
D_\Gamma+L_{\rm core}+L_G
\]

而 raw RHS 仍是 \(D_\Gamma y\)，把
\(\Lambda=D_\Gamma+L_{\rm core}\) 放进 (1.1) 会把 normalized direct map 的 RHS
改为 \(\Lambda y\)。要恢复原问题，需要保留固定 right factor，例如

\[
(\Lambda+L_G)^{-1}D_\Gamma
=\Lambda^{-1/2}T_G\Lambda^{-1/2}D_\Gamma.
\]

端口压缩定理 2.1 仍可用于只含 port-input 的固定 \(C,F,Q\) convex objective；
但完整 hidden-input direct-smoother theorem 不能不经推导直接改写。

### 7.3 Schur/WLS base 的范围

一般线性 WLS 在消去已知 local variables 后，端口 Schur complement 可以是任意
SPD matrix，这为 \(\Lambda\) 提供另一种解释。此时 \(T_G\) 描述以 \(\Lambda\)
geometry 归一化的 inverse/score response。

不过 partial-partition conclusion 依赖未知 factors 全部是 graph differences

\[
(e_u-e_v)(e_u-e_v)^T.
\]

若未知 measurement rows 是一般 rank-one \(a_ea_e^T\)，极限由一般 linear-matroid
flats 索引，而不是 set partitions。若 fixed base 在 port 与 hidden variables 之间已有
非零 cross block，(3.9) 也会改变。因此本定理不能声称覆盖任意 network WLS。

### 7.4 admissible family 的量词

精确下界使用：

- edge conductance 无统一上界；
- hidden grounding 总量可趋于无穷，或者 hidden node count 可增长；
- completion size 不预先固定。

若 edge weights 有上界、hidden count 固定、总 grounding 有上界，(3.15) 仍提供
universal upper envelope，但 scenario vertices 未必可达，(2.1)/(2.4) 的 equality
一般会变成保守上界。研究受限 budget 下更细的 reachable set，可能比继续抽象推广
更有应用价值。

最重要的一点是：原图 \(G\) 完全不要求 acyclic。forest 只是 DPP proof 中的随机
辅助对象，而不是网络拓扑假设。这确实消除了 2011 模型中“全图无环”的强限制，
但研究目标已明确收缩到 scalar Laplacian completion / normalized smoother class。

---

## 8. Primary-literature collision audit（截至 2026-09-14）

本次检索重点覆盖 matrix-forest/DPP、graph Tikhonov、Kron reduction、boundary
groves/electrical compactification、partition polytopes 和 uncertain-topology graph
filtering。下面是与候选定理最接近的原始文献及逐对象判断。

| 文献 | 已有结论 | 与本文候选结论的关系 |
|---|---|---|
| Chebotarev--Shamis, *Matrix-Forest Theorems* | \((I+L)^{-1}\) 的 entries/矩阵由 rooted spanning forests 表示。[^1] | forest representation 是经典起点；没有 arbitrary SPD port compression、hidden-completion hull 或 robust \(Q\)。 |
| Pilavcı--Amblard--Barthelmé--Tremblay 2021 | 对 diagonal \(Q\)，明确研究 \((L+Q)^{-1}Q\)，给出 component 内按 \(q_i\) 加权的 random-forest smoother 及其无偏 Monte Carlo 算法。[^2] | diagonal-grounding random-forest identity 直接相撞；其目标是已知图上的求解/抽样，不是未知任意 hidden completion 的 sharp convex envelope/minimax。没有 dense SPD \(\Lambda\)。 |
| Dereziński--Khanna--Mahoney 2020；Kassel--Lévy 2023 | selected-column/DPP mean projection 对任意 vectors 成立。[^3][^4] | (3.13) 的第二层 arbitrary-\(V\) identity 已知；这正是 SPD 推广证明的核心工具，所以不能把该 identity 当新定理。 |
| Dörfler--Bullo 2013 | loopy Laplacians 在 Kron/Schur reduction 下闭合，研究 boundary clique、effective resistance 和谱。[^5] | 解释端口 Schur base 的物理来源，但不描述 grounded resolvent compression 的 convex hull 或 robust smoother LMI。 |
| Kenyon--Wilson 2011 | 对 circular planar groves，boundary connection-partition probabilities 可由 response/DtN matrix 表示。[^6] | boundary partitions 是危险近邻，但研究箭头相反：从 response 求 grove probabilities；没有本报告的 normalized inverse convexification。 |
| Lam 2018, electroid/cactus compactification | compactify circular planar electrical networks；infinite-conductance shorts 由 noncrossing boundary partitions/cactus networks 描述，并用 grove coordinates/Grassmannian 分层。[^7] | 证明“边界节点无限导通后按 partition 识别”绝非新机制。但其对象是 planar ungrounded DtN response 的 projective compactification，不是 arbitrary hidden grounded resolvent 的 convex hull。详见下节。 |
| De Rosa--Khajavirad 2022 | 研究 \(\sum_B\mathbf1_B\mathbf1_B^T/|B|\) 型 full partition matrices 的 ratio-cut polytope/facets。[^8] | \(\Lambda=I\)、无 inactive block 的 scenario matrices 已知于 clustering polytope。arbitrary SPD joint projections、hidden electrical realization和 smoother minimax不在该文目标中。 |
| Moura--Yaman--Leus 2026 | fixed host graph 上 labelled connected subpartition incidence-vector polytope及其 facets/separation complexity。[^9] | “subpartition”组合对象相邻，但变量不是 projection matrix；宿主图固定，而本文允许任意 hidden completion。 |

### 8.1 与 Kenyon--Wilson / Lam compactification 的严格区分

这是最需要提前回应 reviewer 的近邻。

对一个传统 electrical network，boundary response/DtN matrix 是

\[
\mathcal R_G
=L_{\Gamma\Gamma}
-L_{\Gamma Z}L_{ZZ}^{-1}L_{Z\Gamma},
\tag{8.1}
\]

它把 boundary voltages 映到 boundary currents，通常 singular、row sums 为零。
Kenyon--Wilson 用 \(\mathcal R_G\) 表达 planar grove 的 boundary partition
probabilities。Lam compactify 的是 circular planar network 的 electrical-equivalence
class / projective grove-coordinate space；边界 shorting 只能给 noncrossing partitions，
compactification 本身仍有正维连续 cells。

本报告的对象是

\[
X_G
=\Lambda^{1/2}
\left[
\begin{pmatrix}\Lambda&0\\0&D_Z\end{pmatrix}+L_G
\right]^{-1}_{\Gamma\Gamma}
\Lambda^{1/2}.
\tag{8.2}
\]

若把括号内矩阵对 hidden variables 做 Schur complement，(8.2) 是一个**grounded
Kron operator 的 inverse**再归一化，而不是 DtN response 本身。本文还取的是全部
nonplanar completions 的**闭凸包**，并允许 inactive ports 被 hidden grounding
吸收；所以 scenarios 是全部 partial partitions，共 \(B_{q+1}\) 个，而不是 planar
noncrossing partitions。

此外，一般 SPD \(\Lambda\) 可能有任意符号的 off-diagonals，通常根本不是 resistor
network response/grounded Laplacian。因而 Lam 的 compactification 不实质覆盖定理
2.1。

反过来，也不能据此忽视 Lam：第 3.5 节中“高导通使端口融合为 partition blocks”
是 cactus/electrical compactification 早已系统研究的退化机制。可主张的新组合只能是
两层 DPP 后得到的 exact SPD partial-projection hull、inactive dilution、full-input
robust minimax 和 sharp norm classification。

2026 年新出现的 electroid-variety 工作继续研究 Lam compactification 的代数几何与
分层，而不是 (8.2) 的 robust convex envelope。[^10]

### 8.2 检索结论应如何措辞

截至报告日期，本次没有找到一篇原始论文逐式包含以下整组结论：

1. 任意 \(\Lambda\succ0\) 的 fixed port base；
2. 任意有限 hidden graph/size/cycles/weights/positive diagonal hidden base；
3. (8.2) 的 closed convex hull 等于 \(B_{q+1}\) 个 SPD partial-partition projections；
4. 每个场景由 degree-2 paths sharp realize；
5. full-input local smoother 的 exact finite minimax/LMI；
6. Schatten exactness iff \(p\ge2\)。

这支持“可能是新 theorem package”的保守判断，不构成不存在先例的数学证明。投稿前
仍应做一次 citation graph / MathSciNet / zbMATH 的人工交叉核验，特别沿着 DPP mean
projection、electrical compactification 和 partition-subspace arrangement 三条引用链。

---

## 9. 原创性与论文价值判决

### 9.1 哪些不能认领

- resolvent 是 random forest/component projections 的期望；
- arbitrary vectors 的 DPP mean-projection identity；
- infinite-conductance edges 把节点 short 成 partition blocks；
- \(\Lambda=I\) full partition matrices 与 clustering polytopes 的基本对象；
- Kron reduction 的网络解释。

### 9.2 哪些目前可以作为候选贡献

- 对 arbitrary SPD port geometry，任意 hidden completion 的端口像恰好闭合为
  SPD partial-partition projection polytope；
- 处理 nonorthogonal graph components 的两层 projection argument；
- 与上界匹配的 degree-2 inactive-dilution converse；
- 保留所有 hidden-input columns 后，robust direct smoother 精确归约到同一有限场景；
- Schatten--\(p\) 的 universal iff threshold \(p\ge2\)；
- locality-constrained \(Q\) 的 exact finite SDP，而复杂度只依赖 interface size \(q\)。

### 9.3 强度判断

arbitrary SPD 推广使定理不再依赖 component indicator 的坐标正交性，数学表达明显更
完整；它也统一了 correlated fidelity、preconditioned inverse 和 Schur-base 场景。
但证明的两个核心 identity 都能由已知 DPP mean-projection 很短地推出。因此：

- 作为 graph signal processing / robust estimation / applied matrix analysis 论文的
  主定理：**强且值得推进**；
- 作为只讲 convex hull 的纯矩阵理论论文：**原创强度仍偏中等**；
- 若实验能显示 constrained opcode SDP 相比 full-contraction/Petersen bounds 明显更小、
  并展示隐藏规模不影响 certificate cost，应用收益会很有说服力；
- 若想投更强的纯应用数学期刊，最好再加入受限 hidden budget 的精确 hull、有效
  separation/cut-generation，或一个非平凡 complexity theorem，而不是继续堆更抽象名词。

---

## 10. 数值审计

独立脚本：

```text
experiments/theory_search/spd_port_completion_checks.py
```

运行命令：

```powershell
python experiments\theory_search\spd_port_completion_checks.py
```

截至 2026-09-14，六项检查全部通过：

1. arbitrary SPD 下 scenario 数为 \(B_4=15\)，全部 symmetric/idempotent 且互异，
   并显式包含 \(0,I\)；
2. nonorthogonal component vectors 的第二层 DPP mixture 逐矩阵重构 (3.10)，并检查
   \(h_j=0\) 的正则化极限；
3. 40 个随机 SPD-base hidden graphs 的端口块均通过 convex-hull LP feasibility，
   full smoother 的 \(p=2,3,4,\infty\) 风险均不超过 finite scenarios；
4. 非对角 \(\Lambda\) 下的 active block + inactive port degree-2 lower witness 收敛；
5. 统一两节点反例对 \(p=1,1.2,1.5,1.9\) 全部严格失败，在 \(p=2\) 正好相等；
6. shared-opcode 例重现 (6.9)--(6.10)，确认 connected intermediate scenario active。

这些检查不是证明；它们专门锁定最容易出错的 component cross-Gram、hidden columns、
zero hidden mass、norm threshold 和 constrained-design 退化点。

---

## Sources

[^1]: Pavel Chebotarev and Elena Shamis, “[Matrix-Forest Theorems](https://arxiv.org/abs/math/0602575),” arXiv:math/0602575; earlier journal version, *Automation and Remote Control* 58(9) (1997), 1505--1514. Local archive: `literature/09_network_completion/2006_chebotarev_shamis_matrix_forest_theorems.pdf`.

[^2]: Yusuf Yiğit Pilavcı, Pierre-Olivier Amblard, Simon Barthelmé, and Nicolas Tremblay, “[Graph Tikhonov Regularization and Interpolation via Random Spanning Forests](https://arxiv.org/abs/2011.10450),” *IEEE Transactions on Signal and Information Processing over Networks* 7 (2021), 359--374, DOI [10.1109/TSIPN.2021.3084879](https://doi.org/10.1109/TSIPN.2021.3084879), especially equations (10)--(11), (19)--(21), and Proposition 2. Local archive: `literature/09_network_completion/2021_pilavci_et_al_graph_tikhonov_random_forests.pdf`.

[^3]: Michał Dereziński, Rajiv Khanna, and Michael W. Mahoney, “[Improved Guarantees and a Multiple-Descent Curve for Column Subset Selection and the Nyström Method](https://arxiv.org/abs/2002.09073),” NeurIPS 2020, especially Lemma 5. Local archive: `literature/09_network_completion/2020_derezinski_khanna_mahoney_dpp_projection.pdf`.

[^4]: Adrien Kassel and Thierry Lévy, “[On the Mean Projection Theorem for Determinantal Point Processes](https://doi.org/10.30757/ALEA.v20-17),” *ALEA* 20 (2023), 497--504; [arXiv:2203.04628](https://arxiv.org/abs/2203.04628). Local archive: `literature/09_network_completion/2023_kassel_levy_mean_projection_dpp.pdf`.

[^5]: Florian Dörfler and Francesco Bullo, “[Kron Reduction of Graphs with Applications to Electrical Networks](https://arxiv.org/abs/1102.2950),” *IEEE Transactions on Circuits and Systems I* 60(1) (2013), 150--163, DOI [10.1109/TCSI.2012.2215780](https://doi.org/10.1109/TCSI.2012.2215780). Local archive: `literature/09_network_completion/2013_dorfler_bullo_kron_reduction.pdf`.

[^6]: Richard W. Kenyon and David B. Wilson, “[Boundary Partitions in Trees and Dimers](https://arxiv.org/abs/math/0608422),” *Transactions of the American Mathematical Society* 363(3) (2011), 1325--1364, DOI [10.1090/S0002-9947-2010-04964-5](https://doi.org/10.1090/S0002-9947-2010-04964-5).

[^7]: Thomas Lam, “[Electroid Varieties and a Compactification of the Space of Electrical Networks](https://arxiv.org/abs/1402.6261),” *Advances in Mathematics* 338 (2018), 549--600, DOI [10.1016/j.aim.2018.09.014](https://doi.org/10.1016/j.aim.2018.09.014).

[^8]: Antonio De Rosa and Aida Khajavirad, “[The Ratio-Cut Polytope and K-Means Clustering](https://arxiv.org/abs/2006.15225),” *SIAM Journal on Optimization* 32(1) (2022), 173--203, DOI [10.1137/20M1348601](https://doi.org/10.1137/20M1348601). Local archive: `literature/09_network_completion/2022_de_rosa_khajavirad_ratio_cut_polytope.pdf`.

[^9]: Phablo F. S. Moura, Hande Yaman, and Roel Leus, “[On the Connected (Sub)partition Polytope](https://arxiv.org/abs/2401.01716),” *Mathematical Programming* (2026), DOI [10.1007/s10107-025-02321-1](https://doi.org/10.1007/s10107-025-02321-1). Local archive: `literature/09_network_completion/2026_moura_leus_yaman_connected_subpartition_polytope.pdf`.

[^10]: Dawei Shen, Mia Smith, and David E. Speyer, “[Algebraic Geometry of Electroid Varieties](https://arxiv.org/abs/2607.05576),” arXiv:2607.05576 (submitted July 2026).
