# 隐藏节点 grounded Laplacian 的端口极点与局部 smoother 审计

## 0. 结论先行

候选结论在修正量词和输入几何之后是**正确的**，但原先那句“任意 convex local-WLS loss 都只需检查 \(B_{q+1}\) 个场景”是错误的。最稳妥的最终分层如下。

1. **端口压缩层：成立且最干净。** 固定 \(q\) 个端口及其正 grounding，允许任意有限数量的隐藏节点、任意正隐藏 grounding、任意无向含环拓扑和非负边权。对
   
   \[
   T_G=D_G^{1/2}(D_G+L_G)^{-1}D_G^{1/2},
   \]
   
   所有端口主子块 \(T_G[\Gamma,\Gamma]\) 的闭凸包，恰好是 weighted partial-partition projections 的凸包。每个这样的投影都是 exposed vertex；每个 vertex 都可由最大度数 \(2\) 的连通路径网络逼近。因此，对任意连续凸函数 \(\Phi(T_G[\Gamma,\Gamma])\)，全拓扑最坏值恰等于 \(B_{q+1}\) 个有限场景中的最大值。

2. **完整隐藏输入层：只对合适的 Gram 风险成立。** 对直接 graph-Tikhonov smoother 的局部线性逼近，以及 common Schur/primitive-score 误差，谱范数、Frobenius 范数和 Schatten-\(p\) 范数 \(p\ge 2\) 都有同样的精确 \(B_{q+1}\)-场景归约。一个更统一的安全表述是：误差 Gram 矩阵上的连续、凸、Loewner 单调损失可以精确归约。

3. **一般 full-matrix convex loss：NO-GO。** 隐藏列含有端口主子块没有记录的方向信息。即使 \(q=1\)，取一个 port 与一个 hidden vertex 高导通融合，full projection 为 \(\mathbf1\mathbf1^T/2\)，对“隐藏坐标绝对值”这一凸损失，真实值为 \(1/2\)，而两个 padded partial 场景都给 \(0\)。Schatten-\(p\) 在所有 \(1\le p<2\) 都有严格反例。

4. **异质 grounding 必须归一化。** 上述 theorem 对 normalized inverse \(T_G\) 成立。若仍以未归一化 \((D_G+L_G)^{-1}\) 的全部 primitive RHS 列使用普通 Euclidean 范数，结论不仅失真，而且在隐藏 grounding 可趋于零、隐藏节点数不受限时 worst case 可以发散。

5. **自由 \(Q\) 的 smoother 设计会退化。** 若 \(Q\) 没有任何局部性或共享结构限制，则 \(Q=C/2\) 恒为最优，值恒为 \(\|C\|_2/2\)。精确 certificate 仍然成立，但要产生非平凡算法设计，必须加入通信零模式、shared opcode、量化或 nominal constraints。第 4 节给出一个 connected partial scenario 改变最优 shared gain、并产生 \(8.83\%\) 半径差的解析例子。

6. **原创性判决：数学正确，独立原创强度中等，应用组合值得继续。** Matrix-forest identity，以及 \((L+D)^{-1}D\) 是随机森林分量内加权平均算子的期望，已是经典结果；Pilavci 等 2021 年的论文已经逐式给出后者。[^1][^2] 全分割的 normalized partition matrices 及其 convex hull 也已出现在 clustering polytope 文献。[^7] 本次未找到直接覆盖“任意隐藏 completion 的端口 partial hull + degree-2 sharp converse + port-only robust smoother minimax/LMI”的原始论文。这个组合可以作为应用数学、graph signal processing 或 network estimation 论文的中心定理；若只把 forest mixture 换一种符号写成纯数学论文，则强度不足。

因此，本次审计给候选方向的判决是：

> **经修正后 GO，但不要宣称发现了新的 matrix-forest theorem，也不要宣称终结一般 WLS/state-estimation。最可信的 headline application 是：对未知隐藏网络 completion，端口局部 graph smoother 的最优鲁棒线性逼近可以精确化为 \(B_{q+1}\) 个有限 LMI。**

**证明状态。** 第 3--7 节的 diagonal-grounding 主结论、适用 loss 范围、三个 NO-GO、Schatten 阈值、自由-\(Q\) 退化和 shared-opcode 例子均已在本报告中给出完整代数证明，并由脚本作独立数值复核。第 11 节的 dense-\(\Lambda\) 结论和第 13 节的固定 hidden-budget 结论给出了可闭合的证明骨架，但建议在投稿主文中先列作 extension，待把 DPP/奇异极限的量词展开成正式引理。文献中的“未找到直接先例”只是截至本次原始论文检索的阴性证据，不是可证明的全球 novelty 命题。

## 1. 精确模型与不能省略的量词

固定一组有标签的端口

\[
\Gamma=\{1,\ldots,q\},\qquad
D_\Gamma=\operatorname{diag}(d_1,\ldots,d_q)\succ0.
\]

一个 completion 可以添加任意有限的隐藏集合 \(Z\)，并选择

\[
D_Z\succ0\ \text{diagonal},\qquad
L_G=L_G^T\succeq0
\]

其中 \(L_G\) 是 \(V=\Gamma\mathbin{\dot\cup}Z\) 上任意无向加权图的 Laplacian；边权非负，图可含环、可不连通。令

\[
D_G=\operatorname{diag}(D_\Gamma,D_Z),\qquad
T_G=D_G^{1/2}(D_G+L_G)^{-1}D_G^{1/2}.
\tag{1}
\]

若要求每个 admissible graph 连通，后面的反向构造仍成立：只需用趋于零的弱边把各条路径串成一条总路径。若要求固定 host graph、固定隐藏节点数、隐藏 grounding 总和有上界、边权有统一上界，反向等式一般不再成立；这时本文的多面体只是 sharp universal upper envelope，而不是该受限 family 的精确可达集。

令 \(E_\Gamma=[I_q\ 0]\) 为端口 selection matrix，并记

\[
X_G:=E_\Gamma T_GE_\Gamma^T=T_G[\Gamma,\Gamma].
\tag{2}
\]

对非空 \(A\subseteq\Gamma\)，定义

\[
v_A:=D_\Gamma^{1/2}\mathbf1_A,
\qquad d(A):=\mathbf1_A^TD_\Gamma\mathbf1_A,
\qquad P_A^{(d)}:=\frac{v_Av_A^T}{d(A)}.
\tag{3}
\]

取任意 active subset \(U\subseteq\Gamma\) 和 \(U\) 的 set partition 
\(\pi=\{A_1,\ldots,A_s\}\)，定义 weighted partial-partition projection

\[
P_{\pi,U}^{(d)}:=\sum_{A\in\pi}P_A^{(d)},
\tag{4}
\]

并在 
\(\Gamma\setminus U\) 上补零。各 \(v_A\) 支撑不交，所以 (4) 是 Euclidean orthogonal projection。记这些矩阵的集合为

\[
\mathcal P_\Gamma^\partial(d).
\]

其数量为

\[
|\mathcal P_\Gamma^\partial(d)|
=\sum_{u=0}^q{q\choose u}B_u=B_{q+1}.
\tag{5}
\]

把 inactive ports 与一个 cemetery symbol 放在同一 block，确实给出与 (q+1) 元素 set partitions 的双射；但这只是**组合索引双射**。普通 (q+1) 元素 partition matrix 的 cemetery block 在 inactive ports 上有非零项，而 (4) 把它整个删除，所以二者不是同一个矩阵 polytope，也不是简单的 principal-submatrix equality。

## 2. 经典起点：weighted random-forest projection identity

Chebotarev--Shamis 的 matrix-forest theorem 把 \((I+L)^{-1}\) 的条目写成 rooted spanning forests 的归一化权重。[^1] 对任意正 diagonal \(D\)，Pilavci 等明确使用

\[
K_G=(L_G+D_G)^{-1}D_G
\tag{6}
\]

并证明：给定随机 rooted spanning forest 的 component partition，条件 root distribution 在 component \(C\) 内与 \(d_j\) 成正比；相应的分量内加权平均矩阵 \(S_{\mathcal F}\) 满足

\[
(S_{\mathcal F})_{ij}
=\frac{d_j}{d(C)}\mathbf1\{i,j\in C\},
\qquad
\mathbb E S_{\mathcal F}=K_G.
\tag{7}
\]

这正是该文 equations (19)--(21) 与 Proposition 2 的矩阵内容。[^2] 作相似变换

\[
P_{\mathcal F}^{(d)}
:=D_G^{1/2}S_{\mathcal F}D_G^{-1/2}
=\sum_{C\in\pi(\mathcal F)}
\frac{v_Cv_C^T}{d(C)},
\qquad v_C=D_G^{1/2}\mathbf1_C,
\tag{8}
\]

便得到 orthogonal projection，且

\[
\boxed{T_G=\mathbb E P_{\mathcal F}^{(d)}.}
\tag{9}
\]

一般 DPP 的 mean-projection identity 也有更抽象的版本；Derezinski--Khanna--Mahoney 给出离散 selected-column projection 的期望公式，Kassel--Lévy 则系统研究 finite-rank DPP 的 mean projection theorem。[^3][^4] 这些工作说明“一个 contraction 是随机 projection 的均值”并非新的数学母题。

因此，论文里应把 (9) 标为引用的基础工具，不应作为原创 theorem。

## 3. 主数学结论：端口压缩的闭凸包

### 定理 3.1（任意隐藏 completion 的端口闭凸包）

令 \(\mathfrak G(D_\Gamma)\) 是第 1 节的全部有限 completion，并假定该类允许加入总 grounding 可趋于无穷的隐藏路径，例如任意多个 unit-grounded hidden vertices。则

\[
\boxed{
\overline{\operatorname{conv}}
\{X_G:G\in\mathfrak G(D_\Gamma)\}
=\operatorname{conv}\mathcal P_\Gamma^\partial(d).}
\tag{10}
\]

而且：

- 每个 \(P\in\mathcal P_\Gamma^\partial(d)\) 都是右侧多面体的 exposed vertex；
- 每个这样的 \(P\) 都是某个**连通、最大度数 \(2\)** 的 completion 序列的极限；事实上 witnesses 可以全取为路径；
- uniform grounding \(D_G=\mu I\) 时，
  
  \[
  P_A^{(d)}=\frac{\mathbf1_A\mathbf1_A^T}{|A|},
  \]
  
  因而得到未加权 partial-partition projections。

#### 证明

**上界。** 固定一个 full forest projection (8)，并令 \(C\) 遍历其 components。若

\[
A=C\cap\Gamma\ne\varnothing,
\]

则该 component 对端口主子块的贡献为

\[
\frac{v_Av_A^T}{d(C)}
=\alpha_A P_A^{(d)},
\qquad
\alpha_A:=\frac{d(A)}{d(C)}\in(0,1].
\tag{11}
\]

不同 \(A\) 两两不交。令 \(B_A\sim\operatorname{Bernoulli}(\alpha_A)\) 相互独立，则

\[
\sum_A\alpha_AP_A^{(d)}
=\mathbb E\!\left[\sum_A B_AP_A^{(d)}\right].
\tag{12}
\]

括号中的每个 realization 都是一个 partial-partition projection。因此每个 full forest projection 的 port compression 位于
\(\operatorname{conv}\mathcal P_\Gamma^\partial(d)\)。再对 (9) 取期望即得

\[
X_G\in\operatorname{conv}\mathcal P_\Gamma^\partial(d)
\quad\text{对每个有限 }G.
\tag{13}
\]

这个论证完整保留了 hidden components：它们不是被忽略，而是通过系数 \(\alpha_A=d(A)/d(C)\) 把 port block 连续衰减到 \(0\)。

**反向实现。** 固定 \(P_{\pi,U}^{(d)}\)。对每个 active block \(A\in\pi\)，用一条只经过 \(A\) 中 ports 的高导通路径连接它们。对每个 inactive port \(i\notin U\)，从 \(i\) 接一条含 \(h\) 个 unit-grounded hidden vertices 的路径。先把这些 path components 的内部 conductance 送到无穷。因为

\[
T=(I+D^{-1/2}LD^{-1/2})^{-1},
\]

高导通极限就是投影到每个 component 的
\(\operatorname{span}\{D^{1/2}\mathbf1_C\}\)，即 (8)。active component 不含 hidden vertex，所以端口块恰为 \(P_A^{(d)}\)。inactive component 的 port row 的平方范数为

\[
\frac{d_i}{d_i+h}\longrightarrow0.
\tag{14}
\]

于是再令 \(h\to\infty\)，端口行对可见列和所有隐藏列同时消失，而不是只让 principal block 消失。对每个 \(h\) 选择足够大的有限 conductance，取 diagonal sequence，便得到所需极限。

若要求 graph connected，把上述有限条 paths 依次以 conductance \(\varepsilon_h\to0\) 的弱边连接其 endpoints。每个内部 path endpoint 原度数至多 \(1\)，连接后总度数仍至多 \(2\)，整个 graph 是一条路径。由 inverse 对正定矩阵条目的连续性，可以把 \(\varepsilon_h\) 取得足够小而不改变极限。这证明每个 partial projection 都属于左侧 closure。

**vertex 与 exposedness。** 对任意 distinct orthogonal projections \(P,Q\in\mathcal P_\Gamma^\partial(d)\)，

\[
\langle2P-I,Q\rangle
=\operatorname{tr}P-\|P-Q\|_F^2
<\operatorname{tr}P
=\langle2P-I,P\rangle.
\tag{15}
\]

所以线性泛函 \(X\mapsto\langle2P-I,X\rangle\) 唯一暴露 \(P\)。证毕。

### 推论 3.2（任意端口 convex loss 的精确有限极值）

若 \(\Phi:\mathbb S^q\to\mathbb R\) 连续且凸，则

\[
\boxed{
\sup_{G\in\mathfrak G(D_\Gamma)}\Phi(X_G)
=\max_{P\in\mathcal P_\Gamma^\partial(d)}\Phi(P).}
\tag{16}
\]

上界来自 (13) 和 convexity，下界来自每个 vertex 的 degree-2 approximation。这里的 “sup” 往往不能改成 “max”，因为 inactive vertex \(0\) 需要隐藏 grounding 总和趋于无穷，或隐藏节点数趋于无穷。

这条推论是最宽、也最不容易被 reviewer 击穿的“任意凸损失”版本：损失必须只依赖 terminal compression \(X_G\)。

## 4. Headline application：端口局部 graph smoother 的 exact minimax

### 4.1 问题

在 normalized coordinates 中，global graph-Tikhonov smoother 是

\[
u\longmapsto T_Gu.
\]

给定 \(C\in\mathbb R^{p\times q}\)，我们关心 global smoother 的 port output

\[
C E_\Gamma T_Gu,
\]

但局部算法只能读取 port input \(E_\Gamma u\)，所以采用共同线性规则

\[
Q E_\Gamma u,
\qquad Q\in\mathcal Q_{\rm loc}\subseteq\mathbb R^{p\times q}.
\]

完整误差算子为

\[
\mathcal A_G(Q)
:=CE_\Gamma T_G-QE_\Gamma
\in\mathbb R^{p\times(|\Gamma|+|Z|)}.
\tag{17}
\]

注意 (17) 保留了所有 hidden-input columns。

### 定理 4.1（固定 \(Q\) 的精确场景归约）

令 \(\psi:\mathbb S_+^p\to\mathbb R\) 连续、凸且关于 Loewner 序单调不减。则对每个固定 \(Q\)，

\[
\boxed{
\sup_{G\in\mathfrak G(D_\Gamma)}
\psi\!\left(\mathcal A_G(Q)\mathcal A_G(Q)^T\right)
=
\max_{P\in\mathcal P_\Gamma^\partial(d)}
\psi\!\left((CP-Q)(CP-Q)^T\right).}
\tag{18}
\]

特别地，(18) 给出

\[
\sup_G\|\mathcal A_G(Q)\|_2
=\max_P\|CP-Q\|_2,
\tag{19a}
\]

\[
\sup_G\|\mathcal A_G(Q)\|_F
=\max_P\|CP-Q\|_F,
\tag{19b}
\]

以及对所有有限 \(2\le p_0<\infty\)，

\[
\sup_G\|\mathcal A_G(Q)\|_{S_{p_0}}
=\max_P\|CP-Q\|_{S_{p_0}}.
\tag{19c}
\]

这里 (19a) 取 \(\psi(M)=\lambda_{\max}(M)\)，(19b) 取 \(\psi(M)=\operatorname{tr}M\)，(19c) 取
\(\psi(M)=\operatorname{tr}(M^{p_0/2})\)。最后一个函数在 \(p_0/2\ge1\) 时是 convex 且 Loewner 单调。

#### 证明

由 (9)，写 \(T_G=\mathbb E\widehat P\)，其中 \(\widehat P\) 是 full forest projection，并令

\[
A_{\widehat P}:=CE_\Gamma\widehat P-QE_\Gamma.
\]

则 \(\mathcal A_G(Q)=\mathbb E A_{\widehat P}\)，且矩阵方差恒等式给出

\[
\mathcal A_G(Q)\mathcal A_G(Q)^T
\preceq\mathbb E[A_{\widehat P}A_{\widehat P}^T].
\tag{20}
\]

对一个固定的 full projection，令

\[
X=E_\Gamma\widehat P E_\Gamma^T.
\]

使用 \(\widehat P^2=\widehat P\) 和 \(E_\Gamma E_\Gamma^T=I_q\)，得到关键 Gram identity

\[
\begin{aligned}
A_{\widehat P}A_{\widehat P}^T
&=CXC^T-CXQ^T-QXC^T+QQ^T\\
&=:M_Q(X).
\end{aligned}
\tag{21}
\]

因此 full hidden columns 的全部能量只通过 \(X\) 出现，并且 \(M_Q(X)\) 对 \(X\) 是 affine。由第 3 节，\(X=\sum_j\theta_jP_j\) 是 partial projections 的凸组合，所以

\[
M_Q(X)=\sum_j\theta_jM_Q(P_j).
\tag{22}
\]

依次对 (20) 使用 \(\psi\) 的 Loewner monotonicity、convexity，再对 (22) 使用 convexity，得到左侧不超过右侧。对 partial projection \(P^2=P\)，

\[
M_Q(P)=(CP-Q)(CP-Q)^T.
\]

第 3 节的 degree-2 graph sequence 还使完整端口 rows 在 \(\ell_2\) 中收敛到
\([P\ 0]\)，故每个右侧场景都可从左侧逼近。等式成立。证毕。

### 4.2 优化 \(Q\)、attainment 与 exact LMI

若 \(\mathcal Q_{\rm loc}\) 是非空 closed affine set，例如由通信半径决定的固定 sparsity pattern，则逐个固定 \(Q\) 的 (19a) 立即给出

\[
\boxed{
\min_{Q\in\mathcal Q_{\rm loc}}
\sup_G\|\mathcal A_G(Q)\|_2
=
\min_{Q\in\mathcal Q_{\rm loc}}
\max_{P\in\mathcal P_\Gamma^\partial(d)}\|CP-Q\|_2.}
\tag{23}
\]

因为 \(P=0\) 是一个场景，右侧目标至少为 \(\|Q\|_2\)，故在 closed feasible set 上 coercive，minimum 达到。

这里有一个必须主动说明的退化情形。若 \(Q\) 在 
\(\mathbb R^{p\times q}\) 中完全自由，则设计问题恒有闭式解

\[
\boxed{
\min_Q\max_{P\in\mathcal P_\Gamma^\partial(d)}\|CP-Q\|_2
=\frac12\|C\|_2,
\qquad Q_*=\frac12C.}
\tag{23a}
\]

事实上，\(P=0,I\) 给出

\[
\max\{\|Q\|_2,\|C-Q\|_2\}\ge\frac12\|C\|_2,
\]

而每个 partial \(P\) 都是 orthogonal projection，所以

\[
\left\|CP-\frac12C\right\|_2
\le\|C\|_2\left\|P-\frac12I\right\|_2
=\frac12\|C\|_2.
\]

因此 fixed-\(Q\) exact certificate 仍有内容，但“自由选择任意 \(Q\)”本身没有复杂的设计 tradeoff。实际应用必须让
\(\mathcal Q_{\rm loc}\) 编码通信零模式、共享系数/共享 opcode、量化字典、constant-preserving constraint 或 nominal-model constraint，并使 \(C/2\) 不可行。

下面是一个最小、完全解析的非退化例子。取 \(q=p=2\)，

\[
C=\operatorname{diag}(1,2),
\qquad
\mathcal Q_{\rm loc}=\{Q=\alpha I_2:\alpha\in\mathbb R\}.
\tag{23b}
\]

这个约束表示两个输出必须使用同一个标量 opcode。若错误地漏掉 connected block

\[
P_{\rm av}=\frac12\mathbf1\mathbf1^T
\]

而只检查 \(0,I,e_1e_1^T,e_2e_2^T\)，则目标化为

\[
\min_\alpha\max\{|\alpha|,|2-\alpha|\}=1,
\qquad \alpha=1.
\]

完整五场景问题的精确解则为

\[
\boxed{
\alpha_*=\frac43-\frac{2\sqrt{10}}{15}
=0.911696\ldots,
\qquad
R_*=\frac{10+2\sqrt{10}}{15}
=1.088304\ldots.}
\tag{23c}
\]

证明如下。令

\[
g(\alpha)=\|CP_{\rm av}-\alpha I\|_2.
\]

若写 
\(\alpha=3/4+t\)，直接计算给出

\[
g(\alpha)^2
=\frac{11}{16}+t^2
+\frac12\sqrt{\frac58+10t^2},
\tag{23d}
\]

所以 \(g\) 在 
\([3/4,\infty)\) 单调不减。令 \(g(\alpha)=2-\alpha\)，其特征方程为

\[
\det\!\left((CP_{\rm av}-\alpha I)^T
(CP_{\rm av}-\alpha I)-(2-\alpha)^2I\right)
=\frac{15\alpha^2-40\alpha+24}{4}=0.
\]

较小根正是 (23c) 的 \(\alpha_*\)。当 \(\alpha\le\alpha_*\) 时，场景
\(e_2e_2^T\) 给下界 \(2-\alpha\ge R_*\)；当
\(\alpha\ge\alpha_*\) 时，(23d) 给 \(g(\alpha)\ge R_*\)。在
\(\alpha_*\) 处，\(e_2e_2^T\)、\(I\) 与 \(P_{\rm av}\) 同时达到 \(R_*\)，其余场景更小，故全局最优性得证。

这说明 connected partial scenario 不是装饰项：漏掉它会把 robust radius 低估
\(8.83\%\)，平方风险低估约 \(18.44\%\)，并给出错误的 shared coefficient。

给定 \(\tau\ge0\)，存在一个 local \(Q\in\mathcal Q_{\rm loc}\) 使全拓扑误差不超过 \(\tau\)，当且仅当

\[
\boxed{
\begin{bmatrix}
\tau I_p&CP-Q\\
(CP-Q)^T&\tau I_q
\end{bmatrix}\succeq0
\quad
\text{对所有 }P\in\mathcal P_\Gamma^\partial(d).}
\tag{24}
\]

因此最优 robust local smoother 是一个含 \(B_{q+1}\) 个、每个大小 \(p+q\) 的有限 SDP。这里不存在 scenario relaxation gap；唯一的数值误差来自 SDP solver。

此外，在 universal fixed-\(Q\) certificate 中，\(B_{q+1}\) 个场景一般一个也不能删除。固定任意目标
\(P\in\mathcal P_\Gamma^\partial(d)\)，在 Frobenius 实例中取

\[
C=I_q,\qquad Q=I_q-P.
\]

对任意另一个 partial projection \(P'\)，

\[
\|P'-(I-P)\|_F^2
=q-\|P-P'\|_F^2.
\tag{24a}
\]

所以目标场景 \(P'=P\) 的平方损失为 \(q\)，且它是唯一 worst scenario。换言之，虽然某个具体 \(C,Q\) 往往只有少数 active constraints，但不存在一个与实例无关、更小的 universal scenario 子集。

### 4.3 uniform 与 heterogeneous grounding 的应用解释

若 \(D_G=\mu I\)，则

\[
T_G=\mu(\mu I+L_G)^{-1},
\]

正是通常的 symmetric graph-Tikhonov smoother，(19)--(24) 使用普通 Euclidean input/output energy。

若 \(D_G\) 异质，Pilavci 等使用的物理 smoother 是

\[
K_G=(D_G+L_G)^{-1}D_G=D_G^{-1/2}T_GD_G^{1/2}.
\]

令 normalized signal \(u=D_G^{1/2}y\)、normalized estimate
\(\widetilde x=D_G^{1/2}x\)，则 \(x=K_Gy\) 等价于

\[
\widetilde x=T_Gu.
\]

所以 weighted theorem 对应 grounding-weighted signal energy，并非人为的代数换元。若应用坚持对原始 \(y\) 和 \(x\) 使用未加权 Euclidean norm，则必须把两个 \(D^{\pm1/2}\) 因子保留，不能直接引用 (19)。

## 5. 更一般的 Schur/primitive-score 误差：成立的范围

令 \(F\in\mathbb R^{q\times k}\)，并考虑

\[
\Delta_G(Q)=
\left[
CE_\Gamma T_GE_\Gamma^TF-Q,
\ -CE_\Gamma T_G
\right].
\tag{25}
\]

对 full forest projection \(\widehat P\)，仍令
\(X=E_\Gamma\widehat PE_\Gamma^T\)。由 \(\widehat P^2=\widehat P\)，

\[
\Delta_{\widehat P}(Q)\Delta_{\widehat P}(Q)^T
=(CXF-Q)(CXF-Q)^T+CXC^T
=:N_Q(X).
\tag{26}
\]

与 direct smoother 不同，\(N_Q(X)\) 不再 affine；但它是 matrix convex。若
\(X=\sum_j\theta_jP_j\)，则

\[
N_Q(X)\preceq\sum_j\theta_jN_Q(P_j),
\tag{27}
\]

因为 \(A\mapsto AA^T\) 满足矩阵 Jensen 不等式，而第二项对 \(X\) 线性。结合第一次 forest mixture 的 Gram Jensen，与定理 4.1 完全相同的论证给出：

\[
\sup_G\|\Delta_G(Q)\|_{S_{p_0}}
=\max_{P\in\mathcal P_\Gamma^\partial(d)}
\left\|[CPF-Q,-CP]\right\|_{S_{p_0}}
\tag{28}
\]

对 \(p_0=2\)、\(p_0=\infty\) 和所有 \(p_0\ge2\) 成立；也可用一般的 continuous convex Loewner-monotone Gram loss 表述。谱范数版同样有有限 Schur LMI，只是右侧矩阵宽度变成 (k+q)。

这证明 full hidden primitive-input 并不自动毁掉 theorem；毁掉的是“任意 full-matrix convex loss”这一过宽推广。

## 6. 三个必须写进论文的 NO-GO 边界

### 6.1 一般 full-matrix convex loss 不由端口主子块决定

取一个 port 和一个 hidden vertex，grounding 均为 \(1\)，并令两者之间 conductance 趋于无穷。full projection 极限为

\[
\widehat P=\frac12
\begin{bmatrix}1&1\\1&1\end{bmatrix}.
\]

考察端口行的负值 \(-e_1^T\widehat P=(-1/2,-1/2)\)，并定义 convex loss

\[
\ell(a_1,a_2)=|a_2|.
\]

真实 full loss 为 \(1/2\)。一个端口的 partial projections 只有 \(0\) 与 \(1\)；在同一两列空间中补零后，对应 rows 是 \((0,0)\) 与 \((-1,0)\)，两者 loss 都为 \(0\)。因此

\[
\sup_G\ell(\text{full error})
\ne\max_{P\in\{0,1\}}\ell(\text{padded partial error}).
\tag{29}
\]

这个 \(q=1\) 反例已经彻底否定“任意 convex full-matrix loss”。

### 6.2 Schatten threshold 是 sharp 的：所有 \(1\le p<2\) 都有反例

取 \(q=1\)，一个 port 与一个 hidden vertex 高导通融合，仍有

\[
\widehat P=\frac12\mathbf1\mathbf1^T.
\]

令 output dimension 为 \(2\)，并在 direct smoother (17) 中取

\[
C=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
Q=\begin{bmatrix}1\\0\end{bmatrix}.
\tag{30}
\]

两个 partial scenarios \(P=0,1\) 的 error 分别是 \(-Q=-e_1\) 与
\(C-Q=e_2\)，所以对每个 Schatten-\(p\) norm，其值都等于 \(1\)。full hidden-limit error 却是

\[
\mathcal A_{\widehat P}(Q)
=\begin{bmatrix}-1/2&1/2\\1/2&1/2\end{bmatrix},
\]

它的两个 singular values 均为 \(1/\sqrt2\)。因此

\[
\|\mathcal A_{\widehat P}(Q)\|_{S_p}
=2^{1/p-1/2}
\begin{cases}
>1,&1\le p<2,\\
=1,&p=2,\\
<1,&p>2.
\end{cases}
\tag{31}
\]

结合定理 4.1 的正面证明可得一个 sharp 判定：在 Schatten norms
\(1\le p\le\infty\) 中，“对所有 \(C,Q\) 和所有隐藏 completion 都能精确归约到 partial scenarios”当且仅当 \(p\ge2\)（把 \(p=\infty\) 解释为 spectral norm）。特别地 nuclear norm \(p=1\) 失败，而阈值 \(2\) 不是证明技术造成的空隙。

### 6.3 未归一化 heterogeneous RHS 失败并可发散

先看最小有限反例。一个 port 的 grounding 为 \(3\)，一个 hidden vertex 的 grounding 为 \(1\)，两者以 conductance \(t\to\infty\) 相连，则

\[
\lim_{t\to\infty}
e_1^T
\left(
\begin{bmatrix}3&0\\0&1\end{bmatrix}
+t\begin{bmatrix}1&-1\\-1&1\end{bmatrix}
\right)^{-1}
=\begin{bmatrix}1/4&1/4\end{bmatrix}.
\]

其 Euclidean row norm 为 \(\sqrt2/4>1/3\)，超过无隐藏节点的 raw endpoint \(1/3\)。

更强地，令 port grounding 为 \(1\)，加入 \(h\) 个 hidden vertices，每个 grounding 为 \(h^{-2}\)，并把整个 component 高导通融合。未归一化 inverse 的 port row 极限为每个坐标均等于

\[
\frac1{1+h\cdot h^{-2}}=\frac1{1+1/h},
\]

所以 full row norm 为

\[
\frac{\sqrt{h+1}}{1+1/h}\longrightarrow\infty.
\tag{32}
\]

这不是 normalized theorem 的漏洞，而是原始 Euclidean RHS 给大量低-grounding hidden coordinates 赋予了不合理的总能量。若物理模型确实采用 raw coordinates，就必须固定隐藏节点预算、grounding 下界或另行推导 bound。

## 7. WLS 风险解释：为什么 primitive RHS theorem 不是普通 measurement-noise MSE

考虑 WLS normal matrix 的块分解

\[
J=
\begin{bmatrix}A&E\\E^T&D\end{bmatrix},
\qquad
A\succ0,
\qquad
\Sigma=D-E^TA^{-1}E\succ0,
\tag{33}
\]

并定义

\[
F=E^TA^{-1},
\qquad
C=R_iA^{-1}E.
\]

centralized target row 与使用共同 local correction \(Q\) 的误差，对 primitive normal-equation RHS \(b\) 的系数是

\[
D_b(Q)=
\left[C\Sigma^{-1}F-Q, -C\Sigma^{-1}\right].
\tag{34}
\]

因此 
\(\|D_b(Q)\|_2\) 确实描述 Euclidean-bounded score/message corruption。但是若原始 whitened measurement noise 为 \(\xi\)，且

\[
b=H^T\xi,
\qquad J=H^TH,
\]

则真正的 noise transfer 是 \(D_b(Q)H^T\)，不是 \(D_b(Q)\)。利用 exact factorization

\[
J=
\begin{bmatrix}I&0\\F&I\end{bmatrix}
\begin{bmatrix}A&0\\0&\Sigma\end{bmatrix}
\begin{bmatrix}I&F^T\\0&I\end{bmatrix}
\tag{35}
\]

以及

\[
D_b(Q)
\begin{bmatrix}I&0\\F&I\end{bmatrix}
=\left[-Q,-C\Sigma^{-1}\right],
\]

得到

\[
\boxed{
D_b(Q)JD_b(Q)^T
=QAQ^T+C\Sigma^{-1}C^T.}
\tag{36}
\]

若 \(0\in\mathcal Q_{\rm loc}\)，则 \(Q=0\) 对每个 completion 都在 Loewner 序上最优；因 \(A\succ0\)，它还是唯一最优。这意味着：

- primitive-score operator norm 的非零最优 \(Q\) 可以是有意义的 robust solver/message-disturbance 设计；
- 它**不能**被直接宣传为普通 whitened sensor-noise 下的非平凡 state-estimation improvement；
- 要得到非平凡的 true-state risk，需要加入 signal prior/bias、模型失配或独立 score/message corruption。此时应逐式证明新的 Gram convexity，不能凭口头类比套用 (28)。

因此，论文的第一应用最好使用第 4 节的 direct smoother；WLS 放在第二应用，并把 primitive RHS、raw measurement noise 和 true-state MSE 三者明确区分。

## 8. 计算复杂度与实际收益

端口数 \(q\) 对应的 scenario 数为：

| \(q\) | \(B_{q+1}\) |
|---:|---:|
| 1 | 2 |
| 2 | 5 |
| 3 | 15 |
| 4 | 52 |
| 5 | 203 |
| 6 | 877 |
| 7 | 4,140 |
| 8 | 21,147 |
| 9 | 115,975 |
| 10 | 678,570 |

精确 SDP 的变量规模主要是 (pq+1)，每个 direct-smoother LMI 的大小是 (p+q)，并且完全不依赖 hidden node count、hidden edge count 或 topology 数量。场景可天然并行生成和检查。与“枚举所有隐藏拓扑/连续边权”相比，这是一项实质性、可测量的收益。

代价同样必须诚实写出：Bell number 超指数增长，所以该 formulation 是 boundary-size fixed-parameter exact method，而不是对大 \(q\) 的 polynomial algorithm。实践上 \(q\le6\) 或 \(7\) 很舒适，\(q=8\) 仍可批量并行，\(q\ge10\) 往往需要 symmetry reduction、active-scenario generation 或问题特定 separation oracle；本次审计没有证明一般的高效 separation theorem。

若 topology 被预先限制为一棵只含 \(q\) 个端口的 tree，没有 hidden vertices，则 full connected partitions 由删去 \(q-1\) 条边得到，场景数是 \(2^{q-1}\)，比 universal \(B_{q+1}\) 小。反过来，universal theorem 的价值正是用 \(B_{q+1}\) 个端口场景替代任意多 hidden vertices 与含环 topology；它不应被描述成所有受限 topology 下都最小的 scenario set。

## 9. 文献碰撞审计

| 文献方向 | 原始论文已经解决的内容 | 是否直接包含本文候选组合 |
|---|---|---|
| Matrix-forest | Chebotarev--Shamis 把 \((I+L)^{-1}\) 写成 rooted-forest 权重比；weighted multigraph 也在其范围内。[^1] | **包含基础 identity，不包含 hidden-port convex hull/minimax。** |
| Graph Tikhonov + RSF | Pilavci 等给出 \(K=(L+D)^{-1}D\)、component 内按 \(d_j\) 加权平均的随机矩阵 \(S\)，并证明 \(\mathbb ES=K\)；研究 Monte Carlo 无偏性、方差、SURE/LOOCV、运行时间和应用。[^2] | **最近的直接先例。** 本次全文术语核对未见 port-only minimax \(Q\)、unknown hidden completion、partial-partition hull、degree-2 converse 或 finite robust LMI。 |
| DPP mean projection | Derezinski 等的 Lemma 5 与 Kassel--Lévy 的 mean projection theorem 给出更一般的随机 projection 期望公式。[^3][^4] | **覆盖 projection-mixture 母题，不覆盖该 graph-completion image。** |
| Kron reduction | Dörfler--Bullo 证明 loopy Laplacian 对 Schur/Kron reduction 闭合，并研究 boundary clique、self-loop、effective resistance、谱和敏感性。[^5] | **说明模型自然，但没有给出 (10)、(18) 或 minimax LMI。** |
| Groves/electrical response | Kenyon--Wilson 研究 planar groves 的 boundary connection partitions，并证明连接概率由 response matrix/Dirichlet-to-Neumann map 决定。[^6] | **研究方向相邻但箭头相反：已知 response 求 partition probability；这里是未知 hidden completion 下 normalized inverse 的 universal convex envelope。** 其 planar/noncrossing 约束也不同。 |
| Electrical-network compactification | Lam 用 cactus networks 紧化 circular-planar electrical networks；无限电导会按 noncrossing partition 识别 boundary vertices，并以 grove measurements/electroid strata 描述边界。[^11] | **这是必须引用的相邻先例，但仍未直接覆盖本文组合。** 它研究 planar response/grove-coordinate 空间及其退化，而非 grounded normalized inverse 的任意非平面 hidden completion、inactive-port dilution、closed convex hull 或 port-only minimax。 |
| Partition-matrix polytopes | De Rosa--Khajavirad 研究 \(Z=\sum_A\mathbf1_A\mathbf1_A^T/|A|\) 型全分割矩阵及 ratio-cut polytope 的 facial structure。[^7] | **full unweighted projection polytope 已知；inactive/cemetery deletion、weighted hidden attenuation 与 electrical realization不在该文目标中。** |
| Connected subpartition polytope | Moura--Yaman--Leus 定义 fixed host graph 上 labeled connected \(k\)-subpartitions 的 incidence-vector convex hull，研究 facets 与 separation complexity。[^8] | **“subpartition”对象相邻，但变量表示、fixed-host connectivity 和 normalized Gram matrices均不同。** |
| Random-cluster / arboreal gas | Fortuin--Kasteleyn 建立以 component count 加权的 random-cluster model；后来的 arboreal-gas 工作研究其 \(q\downarrow0\) forest limit 与 percolation。[^9] | **提供随机 component partition 背景，不给 resolvent port hull 或 local approximation theorem。** |
| Random spanning forest probability | Avena--Gaudillière 研究 rooted forest 的 root process、谱、hitting time、coalescence/fragmentation。[^10] | **概率背景已知；没有本次 robust optimization 结论。** |

### 9.1 对 Pilavci 等的逐式比较

最容易被 reviewer 指出的 collision 是 Pilavci et al.：

- 他们的 \(S\) 在每个 sampled component 内把 observation 替换为按 \(d_j\) 加权的平均；
- 他们 equations (19)--(21) 和 Proposition 2 证明 \(\mathbb ES=(L+D)^{-1}D\)；
- 乘上 \(D^{1/2}\) 与 \(D^{-1/2}\) 后，本文所用 orthogonal forest projection identity (9) 立即得到。

因此以下话术不可用：

- “首次发现 graph Tikhonov inverse 是 partition projections 的均值”；
- “首次把 arbitrary diagonal \(D\) 推广到 weighted forests”；
- “首次用随机森林无偏计算 Tikhonov smoother”。

本次在该论文公开全文中检索并核对了 minimax、worst-case、robust、unknown topology、local approximation、partial partition、hidden node、convex hull、semidefinite/SDP 等术语，没有发现它研究如下问题：只允许读取固定 ports，在任意未知隐藏 completion 上选一个共同 \(Q\)，求其 sharp worst case，并证明所有 worst scenarios 由 \(B_{q+1}\) 个 partial projections 完备刻画。这个“未找到”属于全文阴性证据，不等于数学上的全球不存在证明。

## 10. 哪些是经典，哪些组合可能原创

### 明确经典

- grounded Laplacian inverse 的 rooted-forest 展开；
- graph-Tikhonov smoother 的 random-forest weighted averaging 与无偏 Monte Carlo；
- DPP/random projection 的 mean identity；
- full partition matrices 与相关 clustering polytopes；
- Kron reduction、response matrix、grove boundary partitions；
- convex function 在有限 polytope 上可在 vertices 达到最大这一一般原理；
- operator-norm epigraph 的 Schur LMI。

### 本次未定位到直接先例、因而可能形成论文贡献的组合

1. 任意数量 hidden vertices 下，normalized terminal inverse 的**精确 partial-partition closed hull**；
2. 每个 partial vertex 的**connected degree-2/path sharp converse**，尤其 inactive ports 通过总 grounding 发散而在完整 hidden-input row norm 中消失；
3. 任意 continuous convex terminal loss 的 topology-independent exact \(B_{q+1}\) reduction；
4. direct smoother 的 Gram identity (21) 与 port-only common-\(Q\) minimax equality (23)；
5. exact finite LMI iff (24)，以及 hidden-size-independent complexity；
6. spectral/Frobenius/Schatten-\(p\ge2\) 的正面结果与 nuclear/full-convex counterexamples 组成的清晰 norm frontier；
7. raw WLS noise identity (36) 对应用宣称的纠偏。

这些结果的证明主要由已知 forest identity、一个独立 Bernoulli activation、projection Gram identity 和 singular limit 组成，技术链条短而整齐。它们在 applied matrix analysis / graph signal processing 中足以成为强中心 theorem；要投偏纯数学的顶级期刊，还需要进一步给出该 polytope 的非平凡 facet/separation theory、受限 host graph 的 sharp characterization，或扩展到 block/matrix-valued edges 且保持 exactness。

## 11. 数学扩展：端口 base 可推广到任意 \(\Lambda\succ0\)

对角 grounding 不是凸包结论所需的最大代数范围。令

\[
H_G=\operatorname{diag}(\Lambda,D_Z),
\qquad \Lambda\in\mathbb S_{++}^q,
\qquad D_Z\succ0\ \text{diagonal},
\]

并定义

\[
T_G^{(H)}=H_G^{1/2}(H_G+L_G)^{-1}H_G^{1/2}.
\tag{37}
\]

对 partial partition \(\pi=\{A_1,\ldots,A_s\}\) of
\(U\subseteq\Gamma\)，令 \(U_\pi=[\mathbf1_{A_1}\ \cdots\ \mathbf1_{A_s}]\)，并定义

\[
P_{\pi,U}^{(\Lambda)}
=\Lambda^{1/2}U_\pi
(U_\pi^T\Lambda U_\pi)^{-1}
U_\pi^T\Lambda^{1/2}.
\tag{38}
\]

这是到
\(\operatorname{span}\{\Lambda^{1/2}\mathbf1_A:A\in\pi\}\) 的 orthogonal projection。即使 blocks 不交，dense
\(\Lambda\) 下这些生成向量通常不正交；所以不能再使用第 3 节逐 block 的独立 scalar activation。正确的替代是第二层 volume-sampling/DPP projection identity。

下面给出核验后的证明骨架。写 \(L_G=BWB^T\) 并令

\[
Y=H_G^{-1/2}BW^{1/2}.
\]

则

\[
T_G^{(H)}=(I+YY^T)^{-1}.
\]

对 columns of \(Y\) 采用 \(L\)-ensemble DPP，标准 mean-projection identity 给出

\[
I-T_G^{(H)}
=\mathbb E P_{\operatorname{span}(Y_S)},
\qquad
T_G^{(H)}
=\mathbb E P_{\ker(Y_S^T)}.
\tag{39}
\]

有正概率的 edge subsets 必须线性独立，因而是 forests。若其 component-indicator matrix 是
\(\widehat U\)，则 (39) 中的 projection 为

\[
H_G^{1/2}\widehat U
(\widehat U^TH_G\widehat U)^{-1}
\widehat U^TH_G^{1/2}.
\tag{40}
\]

令 \(U=\widehat U[\Gamma,:]\)、
\(V=\Lambda^{1/2}U\)，并令 \(h_C\) 是 component \(C\) 内 hidden groundings 的总和。由于 \(D_Z\) diagonal，(40) 的 terminal block 正是

\[
X=V\bigl(V^TV+\operatorname{diag}(h_C)\bigr)^{-1}V^T.
\tag{41}
\]

当所有 \(h_C>0\) 时，对 columns of
\(V\operatorname{diag}(h_C)^{-1/2}\) 再应用一次 mean-projection identity，得到 \(X\) 是
\(P_{\operatorname{span}(V_S)}\) 的凸组合；有 \(h_C=0\) 时用
\(h_C+\varepsilon\) 并令 \(\varepsilon\downarrow0\)。每个 selected subset \(S\) 恰好保留若干 forest components，所以这些 projections 正是 (38)。因此

\[
\overline{\operatorname{conv}}
\{T_G^{(H)}[\Gamma,\Gamma]\}
=\operatorname{conv}
\{P_{\pi,U}^{(\Lambda)}\}.
\tag{42}
\]

反向仍使用第 3 节的 degree-2 paths。active components 不放 hidden vertices；inactive components 的 hidden grounding 总和送到无穷。对 (41) 作 block inverse，可见 active terminal block 收敛到 (38)，而 terminal-to-hidden block 是
\(O(h_C^{-1/2})\)，所以完整 hidden-input rows 也消失。弱边再把 components 串成一条连通路径。所有 (38) 仍是 distinct orthogonal projections，故 exposedness、direct-smoother theorem、Schatten sharp threshold 和 scenario-minimality 原样成立。

这个扩展在数学上通过审计，但物理定位必须克制：dense \(\Lambda\) 本身已在 ports 之间引入 baseline coupling，且 normalized coordinate \(\Lambda^{1/2}y_\Gamma\) 混合多个端口；它不再是“每个节点独立 scalar grounding”。因此建议把 (42) 放作数学 extension，而 headline communication/state-estimation model 仍使用 diagonal \(D_\Gamma\)。两层 DPP identity 都是已知工具，潜在原创点仍是 hidden terminal image、path converse 与 robust local approximation 的组合，而不是 DPP identity 本身。[^3][^4]

## 12. 推荐的可发表主定理表述

建议正文主定理使用普通语言，不创造新名词：

> **Theorem (Exact robust local approximation of grounded graph smoothers).** Fix \(q\) observed vertices and their positive grounding weights. Over all finite undirected weighted graph completions with arbitrarily many hidden vertices, the closed convex hull of the observed principal blocks of \(D^{1/2}(D+L)^{-1}D^{1/2}\) is the convex hull of the \(B_{q+1}\) weighted partial-partition projections. Every vertex is exposed and is approached by connected path completions. Consequently, for any fixed local linear map \(Q\), the worst-case spectral, Frobenius, or Schatten-\(p\) error (\(p\ge2\)) in approximating the observed output of the global Tikhonov smoother from observed inputs is attained in the limiting sense at one of these \(B_{q+1}\) scenarios. The optimal common local map is therefore given exactly by a finite convex program; for the spectral norm it is the SDP (24).

论文应在 theorem 后立刻列出以下 assumptions：

- scalar、undirected、nonnegative edge weights；
- fixed positive port groundings；
- hidden count unbounded，admissible family 含 total hidden grounding 发散的 paths；
- heterogeneous case 使用 grounding-normalized coordinates；
- equality 是 over-all-completions supremum，不是每个 fixed topology 的 equality；
- arbitrary convex loss 只对 terminal compression；完整 error 只声称规定的 Gram/norm classes。

## 13. 下一步最有价值的理论与实验工作

一个已经可以严格写出的、比“任意多 hidden vertices”更物理的 refinement 是固定隐藏预算。先取 uniform/unit grounding，并限制 \(|Z|\le h\)。对 \(\Gamma\) 的一个 full partition \(\pi\)，以及每个 block 的非负整数分配
\[
k=(k_B)_{B\in\pi},\qquad \sum_{B\in\pi}k_B\le h,
\]
定义
\[
X_{\pi,k}
=\sum_{B\in\pi}\frac{|B|}{|B|+k_B}P_B
=\sum_{B\in\pi}\frac{\mathbf1_B\mathbf1_B^T}{|B|+k_B}.
\tag{43}
\]
则同一个 forest-mixture 证明给出精确 terminal equality
\[
\boxed{
\overline{\operatorname{conv}}
\{T_G[\Gamma,\Gamma]:|Z|\le h\}
=
\operatorname{conv}\{X_{\pi,k}\}.}
\tag{44}
\]
理由是：任一 full forest component 在 ports 上诱导 \(\pi\) 的一个 block，并含有恰好 \(k_B\) 个 hidden vertices；hidden-only components 只消耗未使用预算，对 terminal block 没有贡献。反过来，每个 \((\pi,k)\) 都由高导通 path forest 逼近，要求 connected 时仍可用弱边串联。有限边权只影响 closure，不影响 equality。场景总数至多且实际上按该参数化为
\[
\sum_{r=1}^q S(q,r){h+r\choose r}.
\tag{45}
\]
这些 scenarios 未必全是 vertices，所以不能照搬第 3 节的 exposedness 宣称。当 \(h\to\infty\) 时，令 active blocks 取 \(k_B=0\)、inactive blocks 取 \(k_B\to\infty\)，(44) 的极限恢复 partial-partition hull。

对 full-input direct smoother，(44) 仍给 exact finite reduction，但 \(X_{\pi,k}\) 一般不是 projection，正确场景损失是
\[
\lambda_{\max}^{1/2}\!\left(
CXC^T-CXQ^T-QXC^T+QQ^T
\right),
\]
等价地使用 canonical error
\[
\left[CX-Q,\ C(X-X^2)^{1/2}\right].
\]
不能把它错误简写为 \(\|CX-Q\|_2\)。这个 finite-budget 版本既保留隐藏列能量，又通常不再有自由 \(Q=C/2\) 的普遍退化，是值得优先发展的物理主模型。

1. **先把定理写成正式 LaTeX proof。** 主文采用定理 3.1 + 定理 4.1 + LMI (24)；三个 NO-GO 放 proposition/remark，避免 reviewer 误以为我们遗漏。
2. **做 exact SDP scaling。** 对 \(q=2,\ldots,8\)，记录 Bell scenario generation、建模时间、solver 时间、内存、active LMIs 数；与随机 graph/weight sampling 比较 certification gap。
3. **用真实图数据做隐藏 completion stress test。** 固定一个小 boundary，逐步加入大量 hidden nodes、cycles、heterogeneous grounding，验证 SDP certificate 对所有 samples 成立，并用 degree-2 constructions逼近预测的 active scenario。
4. **证明或否定受限 hidden budget 的精确集合。** 固定 hidden 数 \(h\) 或总 grounding \(M\) 后，系数 \(\alpha_A\) 不再能到任意 \(0\)；这可能产生比 universal partial polytope 更细、也更有物理意义的参数化 hull。
5. **不要急于推广到 block Laplacian。** 非交换 block weights 通常不再产生 component-wise rank-one orthogonal projections；必须先找最小反例或新的 dilation，不应把 scalar proof 逐符号复制。

## 14. 可复现实证据

本项目已有两组审计脚本：

- `experiments/theory_search/hidden_partial_partition_checks.py`：随机隐藏图 upper bound、degree-2 dilution converse、heterogeneous normalized theorem；
- `experiments/theory_search/boundary_partial_partition_counterexamples.py`：一般 full-matrix convex loss、Schatten sharp threshold、shared-opcode active scenario、direct/Schur nuclear norm、未归一化 heterogeneous RHS 的反例与解析值检查。

截至本报告写入时，两组脚本的全部检查均通过。数值检查不是证明；它们的作用是锁定最容易藏错的 hidden-column Gram 项与极限量词。

## Sources

[^1]: Pavel Chebotarev and Elena Shamis, “[Matrix-Forest Theorems](https://arxiv.org/abs/math/0602575),” arXiv:math/0602575 (manuscript posted 2006; weighted graph forest identities). See also their earlier *Automation and Remote Control* 58(9) (1997), 1505--1514, [official record](https://www.mathnet.ru/eng/at2672). Local archive: `literature/09_network_completion/2006_chebotarev_shamis_matrix_forest_theorems.pdf`.

[^2]: Yusuf Yiğit Pilavcı, Pierre-Olivier Amblard, Simon Barthelmé, and Nicolas Tremblay, “[Graph Tikhonov Regularization and Interpolation via Random Spanning Forests](https://arxiv.org/abs/2011.10450),” *IEEE Transactions on Signal and Information Processing over Networks* 7 (2021), DOI [10.1109/TSIPN.2021.3084879](https://doi.org/10.1109/TSIPN.2021.3084879). In particular equations (10)--(11), (19)--(21), and Proposition 2. Local archive: `literature/09_network_completion/2021_pilavci_et_al_graph_tikhonov_random_forests.pdf`.

[^3]: Michał Dereziński, Rajiv Khanna, and Michael W. Mahoney, “[Improved Guarantees and a Multiple-Descent Curve for Column Subset Selection and the Nyström Method](https://arxiv.org/abs/2002.09073),” NeurIPS 2020, especially Lemma 5 on the expected selected-column projection. Local archive: `literature/09_network_completion/2020_derezinski_khanna_mahoney_dpp_projection.pdf`.

[^4]: Adrien Kassel and Thierry Lévy, “[On the Mean Projection Theorem for Determinantal Point Processes](https://doi.org/10.30757/ALEA.v20-17),” *ALEA* 20 (2023), 497--504; [arXiv:2203.04628](https://arxiv.org/abs/2203.04628). Local archive: `literature/09_network_completion/2023_kassel_levy_mean_projection_dpp.pdf`.

[^5]: Florian Dörfler and Francesco Bullo, “[Kron Reduction of Graphs with Applications to Electrical Networks](https://arxiv.org/abs/1102.2950),” *IEEE Transactions on Circuits and Systems I* 60(1) (2013), 150--163, DOI [10.1109/TCSI.2012.2215780](https://doi.org/10.1109/TCSI.2012.2215780). Local archive: `literature/09_network_completion/2013_dorfler_bullo_kron_reduction.pdf`.

[^6]: Richard W. Kenyon and David B. Wilson, “[Boundary Partitions in Trees and Dimers](https://arxiv.org/abs/math/0608422),” *Transactions of the American Mathematical Society* 363(3) (2011), 1325--1364, DOI [10.1090/S0002-9947-2010-04964-5](https://doi.org/10.1090/S0002-9947-2010-04964-5). Local archive: `literature/09_network_completion/2006_kenyon_wilson_boundary_partitions.pdf`.

[^7]: Antonio De Rosa and Aida Khajavirad, “[The Ratio-Cut Polytope and K-Means Clustering](https://arxiv.org/abs/2006.15225),” *SIAM Journal on Optimization* 32(1) (2022), 173--203, DOI [10.1137/20M1348601](https://doi.org/10.1137/20M1348601). Local archive: `literature/09_network_completion/2022_de_rosa_khajavirad_ratio_cut_polytope.pdf`.

[^8]: Phablo F. S. Moura, Hande Yaman, and Roel Leus, “[On the Connected (Sub)partition Polytope](https://arxiv.org/abs/2401.01716),” *Mathematical Programming* (2026), DOI [10.1007/s10107-025-02321-1](https://doi.org/10.1007/s10107-025-02321-1). Local archive: `literature/09_network_completion/2026_moura_leus_yaman_connected_subpartition_polytope.pdf`.

[^9]: C. M. Fortuin and P. W. Kasteleyn, “[On the Random-Cluster Model: I. Introduction and Relation to Other Models](https://doi.org/10.1016/0031-8914(72)90045-6),” *Physica* 57(4) (1972), 536--564. For the forest/\(q\downarrow0\) connection in modern primary literature, see Roland Bauerschmidt, Nicholas Crawford, Tyler Helmuth, and Andrew Swan, “[Random Spanning Forests and Hyperbolic Symmetry](https://arxiv.org/abs/1912.04854),” *Communications in Mathematical Physics* 381 (2021), 1223--1261.

[^10]: Luca Avena and Alexandre Gaudillière, “[Two Applications of Random Spanning Forests](https://doi.org/10.1007/s10959-017-0771-3),” *Journal of Theoretical Probability* 31 (2018), 1975--2004.

[^11]: Thomas Lam, “[Electroid Varieties and a Compactification of the Space of Electrical Networks](https://arxiv.org/abs/1402.6261),” *Advances in Mathematics* 338 (2018), 549--600, DOI [10.1016/j.aim.2018.09.014](https://doi.org/10.1016/j.aim.2018.09.014). Local archive: `literature/09_network_completion/2018_lam_electroid_varieties.pdf`.
