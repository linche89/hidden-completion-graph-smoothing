# 稀疏 pairwise completion 的精确零空间与连通分割归约

## 结论先行

这条线得到了一条可以继续作为**候选主定理**推进的结果，而不是把
single-full-block Petersen lemma 换一套符号再写一次。

先考虑由统一 unary grounding 和任意标量 pairwise factors 组成的 exterior：

\[
\Sigma(w)=\mu I+\sum_{e\in E}w_e b_eb_e^T,
\qquad \mu>0,\quad w_e\ge0,
\]

其中每个 \(b_e\) 至多支撑两个变量；取 \(b_e=e_u-e_v\) 就回到
grounded Laplacian。对固定局部 transfer
\(C,F\) 和所有 completion 共用的局部系数 \(Q\)，Schur 行误差为

\[
\mathcal E_w(Q)
=\bigl[C\Sigma(w)^{-1}F-Q,\,-C\Sigma(w)^{-1}\bigr].
\tag{1}
\]

本报告证明：

1. 对任意固定 pairwise factor family，所有非负权重的 worst convex loss
   精确等于有限个公共零空间投影上的最大值；不同场景由所表示线性拟阵的 flats
   索引。该结论直接允许 signed/general pairwise factors，而不要求 M-matrix。
2. 若 exterior 只在一个固定端口宿主图 \(H\) 上改变非负边权，则
   \(\sup_G\|\mathcal E_G(Q)\|_2\) **精确等于**在 \(H\) 的所有连通块分割所对应的块平均投影上取最大；反向最坏 completion 由各块内高导通的 spanning trees 逼近。
3. 若允许任意多隐藏 exterior 节点，并允许任意稀疏 Laplacian 拓扑，则精确极端场景是端口集合的所有“部分分割投影”：一些端口块被融合，其余端口的影响被稀释到零。每个极端场景都能由**最大度 2 的路径森林**逼近。因此，仅要求稀疏或有界度并不会消除最坏 completion。
4. 上一条还能推广到 heterogeneous positive grounding：在 grounding-normalized input
   geometry 下，普通块平均投影被 \(\sqrt d\)-weighted partial projections 取代，exactness
   与 degree-2 converse 都保留；未归一化 Euclidean RHS 则有一端口严格反例。
5. 由此得到共同局部线性规则存在的有限、非保守 necessary-and-sufficient LMI；端口数为 \(q\) 时共有 \(B_{q+1}\) 个场景，所以问题对 cut size \(q\) 是 exact FPT。固定宿主树时只有 \(2^{q-1}\) 个场景。
6. 一个显式 two-sparse raw-WLS **normal-matrix** realization lemma 表明该 Schur
   completion 类不是抽象伪模型，并自动给出 \(C=R_iF^T\)；它不声称 score norm 等于
   raw-measurement norm。三端口 path 反例则说明：优化共同 \(Q\) 后，中间割仍可能是
   active worst completion。
7. 一个满足 \(C=R_iF^T\) 的二端口例子给出严格 Petersen gap：grounded-network
   的精确半径是 \(\sqrt2\)，而 enclosing full self-adjoint contraction 的 Petersen 值
   是 \(3/2\)。Petersen lemma 对外包球仍然 lossless，但该外包对网络 completion
   已经有严格保守性。

必须严厉划清已有结果：Pilavcı--Amblard--Barthelmé--Tremblay (2021)
Proposition 2 已明确证明 \((Q+L)^{-1}Q\) 是随机森林分区平均矩阵的期望；
Dereziński--Khanna--Mahoney (2020) Lemma 5 已证明一般 L-ensemble 下所选列空间
投影的期望公式。因此，本文的“凸组合/期望恒等式”**不是原创贡献**。
定向检索截至 2026-09-14 尚未找到的是把已知 mean-projection identity 推到：

- 对所有非负权重取 supremum 的 sharp converse；
- representable-matroid flats / graph connected partitions 的 exact worst-scenario reduction；
- 保留 hidden-input columns 的 partial-partition theorem 与 degree-2 converse；
- common Schur rule 的 finite exact LMI、Petersen strict gap 和 raw pairwise-WLS realization。

这些剩余部分仍是短而非平凡的组合定理；阴性检索不能证明首次性。因此最稳妥的判断是：

> **主定理候选：GO；正式原创性声明：仍待投稿级逐定理查重。**

它显著强于此前的 full-block exact SDP，但已知 mean-projection 近邻使原创性风险
高于初版判断，当前主定理原创性置信度下调到约 **0.45**，不能写成“已确认首次发现”。
另一个同等重要的 scope 限制是：以下 operator norm 针对 bounded normal-equation
RHS/score perturbations。对标准 whitened raw WLS measurement noise，(1e) 会让可调
\(Q\) 正交解耦并在 \(0\) 处达到最小；所以本报告尚不能直接宣称解决了非平凡的
raw-sensor MSE minimax。

## 1. 最小模型与 Schur 误差

固定 \(q\) 个 exterior ports \(\Gamma=\{1,\ldots,q\}\)。局部球消元后留下

\[
C\in\mathbb R^{p\times q},\qquad
F\in\mathbb R^{q\times k},\qquad
Q\in\mathcal Q_{\rm loc}\subseteq\mathbb R^{p\times k}.
\]

这里 \(Q\) 必须在看到远端 completion 以前选定；
\(\mathcal Q_{\rm loc}\) 可以编码线性局部性、固定零模式或其他 affine 约束。
对一个 exterior conditional precision \(\Sigma\succ0\)，代入标准块逆公式后误差
正是 (1)。若来自原始正规方程分块

\[
J=\begin{bmatrix}A&E\\E^T&D\end{bmatrix},
\qquad
F=E^TA^{-1},\qquad C=R_iA^{-1}E,
\]

则自动有 \(C=R_iF^T\)。后文所有结论针对 primitive normal-equation RHS；
若输入是原始传感器观测 \(z\)，仍须保留
\(b=\widehat H^TR_z^{-1}z\) 的输入映射，不能把
二者的 operator norm 无条件等同。

### 引理 1.1（显式 raw pairwise-WLS realization）

给定 \(F\in\mathbb R^{q\times k}\)、任意 \(R_i\in\mathbb R^{p\times k}\)，令
\(C=R_iF^T\)。给定任意 \(\mu>0\) 以及每列至多有两个非零元的
\(b_e\in\mathbb R^q\)。则存在 \(a>0\) 和只含 unary 或 two-coordinate
measurement rows 的矩阵 \(H(w)\)，使

\[
J(w)=H(w)^TH(w)
=\begin{bmatrix}A&E\\E^T&D(w)\end{bmatrix},
\quad
A=aI_k,\quad E=aF^T,
\tag{1a}
\]

并且

\[
D(w)-E^TA^{-1}E
=\mu I_q+\sum_e w_eb_eb_e^T,
\qquad
E^TA^{-1}=F,
\qquad
R_iA^{-1}E=C.
\tag{1b}
\]

所以本文的 completion 确实可由稀疏 pairwise WLS factors 生成，而不是先假定一个
自由 full block 再事后解释。

#### 构造与证明

对每个 \(F_{j\ell}\ne0\)，加入一个只支撑 local coordinate \(y_\ell\) 和
port \(x_j\) 的 row

\[
h_{j\ell}=s_{j\ell}e_{y_\ell}
+\frac{aF_{j\ell}}{s_{j\ell}}e_{x_j}.
\tag{1c}
\]

若第 \(\ell\) 列有 \(d_\ell\) 个非零元，取
\(s_{j\ell}^2=a/(2\max\{1,d_\ell\})\)，再用一个 local unary row 补足
\(A_{\ell\ell}=a\)。这样 cross block 恰为 \(E=aF^T\)，而这些 rows 在 exterior
diagonal 上消耗 \(aR_F\)，其中

\[
R_F=\operatorname{diag}(r_1,\ldots,r_q),\qquad
r_j=\sum_{\ell:F_{j\ell}\ne0}
2\max\{1,d_\ell\}F_{j\ell}^2.
\]

目标 raw exterior block 是

\[
D_0=\mu I_q+aFF^T.
\]

矩阵 \(D_0-aR_F\) 的非对角元均为 \(a(FF^T)_{ij}\)。先把每个非对角元
\(m_{ij}\) 写成 two-coordinate square

\[
|m_{ij}|(x_i+\operatorname{sgn}(m_{ij})x_j)^2,
\]

再用 unary squares 补对角。令

\[
M_F=\max_i\left[
r_i-(FF^T)_{ii}+\sum_{j\ne i}|(FF^T)_{ij}|
\right]_+.
\]

任取 \(0<a<\mu/(1+M_F)\)，就使 \(D_0-aR_F\) strictly diagonally dominant，
剩余 unary coefficients 非负。这就用 two-sparse rows 精确合成 \(D_0\)。最后加入
rows \(\sqrt{w_e}b_e^T\)；Schur complement 立即给出 (1b)。由于
\(A\succ0\) 且 Schur complement \(\succ0\)，\(J(w)\succ0\)。证毕。

这里写成 \(H(w)^TH(w)\) 只是采用了 whitened measurement rows。若 raw design 为
\(\widehat H\)、噪声协方差为 \(R_z\succ0\)，令
\(H=R_z^{-1/2}\widehat H\)，则
\(J=H^TH=\widehat H^TR_z^{-1}\widehat H\)；这与输出映射 \(R_i\) 是两个不同对象。

### 输入范数审计：normal-equation RHS 与 raw measurement 不能混同

令

\[
D_b(Q)=\left[C\Sigma^{-1}F-Q,-C\Sigma^{-1}\right].
\tag{1d}
\]

本文的 \(\|D_b(Q)\|_2\) 是对 Euclidean-bounded normal-equation RHS
\(b=[b_{\rm loc}^T,b_{\rm ext}^T]^T\) 的 gain。它**不是** whitened raw measurement
\(\xi\) 的 gain，因为 WLS 中 \(b=H^T\xi\)。事实上，令
\(M=C\Sigma^{-1}\)，利用
\(E=AF^T\) 与 \(D=\Sigma+FAF^T\)，可逐项消去交叉项得到

\[
\boxed{
D_b(Q)J D_b(Q)^T
=Q A Q^T+C\Sigma^{-1}C^T.}
\tag{1e}
\]

因此

\[
\|D_b(Q)H^T\|_2^2
=\lambda_{\max}\!\left(QAQ^T+C\Sigma^{-1}C^T\right).
\tag{1f}
\]

若 \(0\in\mathcal Q_{\rm loc}\)，则 \(Q=0\) 对**每个 completion** 都在 Loewner
序下最小化 (1f)，所以把当前 minimax 直接称为“raw sensor-noise WLS operator
optimization”会使 \(Q\)-设计退化。引理 1.1 只证明 normal matrix 的 pairwise
measurement realizability，不会消除这个输入范数变化。

这项退化只针对“纯 raw noise 下 local output 与 centralized output 的差”。若把
noiseless state variation 也纳入 central-WLS approximation error，则另有精确恒等式

\[
D_b(Q)J=\left[-QA,-QE-C\right],
\tag{1g}
\]

其右侧是相对 centralized noiseless solution 的 bias map，通常不会在 \(Q=0\) 最优。
第 4 节据此给出一个同时包含 state bias、raw measurement variance 和 score/message
disturbance 的非平凡 exact partition corollary。

当前非平凡 loss 的正确应用语义是：central normal-equation solution 对 bounded score
perturbations、量化误差或消息误差的近似。更一般地，若 score ellipsoid 是固定
\(G\succeq0\)，则 loss 可写为 \(\|D_b(Q)G^{1/2}\|_2\)；在 fixed-dimensional
factor model 中它仍是 \(T\) 的连续凸函数，所以定理 2.1 仍适用。

若论文要评价“对真实 state 的估计”而非 central-solution approximation，应另行采用
非平凡的 bias--variance loss。对 raw model
\(z=H(w)x+v\)、\(\operatorname{Cov}(x)=\Pi\)、
\(\operatorname{Cov}(v)=R_v\)，令目标选择矩阵
\(\widetilde R_i=[R_i\ 0]\)，并取 locality-constrained estimator \(L_Qz\)。一个正确的
worst-energy 指标是

\[
\left\|
\left[(L_QH(w)-\widetilde R_i)\Pi^{1/2},\ L_QR_v^{1/2}\right]
\right\|_2,
\tag{1h}
\]

对应 MSE 则是两个块乘自身转置后的 trace。该目标一般不会仅由
\(\Sigma^{-1}\) 的 affine image 决定，故**不能**无证明地套用后文的 flat/partition
exactness；它是下一阶段需要单独闭合的 state-estimation theorem，而不是本报告已经
解决的 claim。

另一个更直接的应用是 normalized graph-Tikhonov smoother 本身：比较
\(C[T_G]_{\Gamma,:}\) 与 local row \([Q,0]\)。该 objective 不含本节的
\((A,E,F)\) Schur blocks，也不是 (1) 的 raw-WLS 重解释；应作为独立 smoother theorem
陈述和审计。本报告以下仍只证明 Schur/score 模型及推论 4.3，不把两条应用线混写。

若进一步坚持所有 pairwise rows 都必须是纯差分 \(x_u-x_v\)，而不是任意
two-coordinate linear factor，则任意 signed \(F\) 不再可实现。一个完全差分的非空
子类是：每个 local \(y_j\) 仅用 weight \(a_j\) 的 matching edge 连到 port
\(x_j\)，外部再加 grounded graph。此时
\(A=\operatorname{diag}(a_j),E=-A\)，Schur block 正好是 grounded exterior，
且 \(F=-I\)。因此“任意 pairwise factor”与“纯 M-matrix difference network”在
scope 上必须分开写。

需要先明确一个量词：允许边权趋于 \(+\infty\)，故下面一般是 supremum，
而不保证某个有限权图达到 maximum。这个开放谱上界正是 sharp converse 能只用路径完成的原因；若每条边还有固定上界，结论会改变，见第 10 节。

## 2. 一般 pairwise factors：精确归约到拟阵 flats

固定 \(B=[b_1\ \cdots\ b_m]\in\mathbb R^{q\times m}\) 和任意 fixed base
\(\Lambda\succ0\)，令

\[
\Sigma(w)=\Lambda+B\operatorname{diag}(w)B^T,
\qquad
T(w)=\Lambda^{1/2}\Sigma(w)^{-1}\Lambda^{1/2},
\tag{2a}
\]

其中 \(w\in\mathbb R_+^m\)。定理本身不需要列稀疏性；本文的 network/WLS 应用
再把每列限制为至多 two-sparse。因此 dense columns 可视为一个纯线性代数扩展，
而不能宣传成局部 measurement model。写
\(\bar B=\Lambda^{-1/2}B\)。对 \(S\subseteq[m]\) 定义

\[
P_S=P_{\ker \bar B_S^T}.
\tag{2b}
\]

投影只依赖于 \(\operatorname{span}\{\bar b_e:e\in S\}\)，所以互不相同的场景可由
\(B\) 所表示线性拟阵的 flats 索引，而不必枚举所有 \(2^m\) 个 subsets。

### 定理 2.1（sharp flat/nullspace extremum）

对任意连续凸函数 \(\Phi:\mathbb S^q\to\mathbb R\)，

\[
\boxed{
\sup_{w\ge0}\Phi(T(w))
=\max_{S\subseteq[m]}\Phi(P_S)
=\max_{\mathcal F\in\operatorname{Flats}(B)}
\Phi(P_{\ker \bar B_{\mathcal F}^T}).}
\tag{2c}
\]

而且每个右侧投影都是可实现极限：令 \(S\) 中权重等于 \(t\)、其余为零，则
\(T(w(t))\to P_S\)。

#### 证明及已知边界

令 \(A=\bar B\operatorname{diag}(\sqrt w)\)，并在 factor index set 上取
L-ensemble

\[
\Pr(\mathsf S=S)
=\frac{\det(A_S^TA_S)}{\det(I+A^TA)}.
\tag{2d}
\]

线性相关的 subset 自动具有零概率。Dereziński--Khanna--Mahoney (2020), Lemma 5
给出已知 mean-projection identity

\[
\mathbb E\,P_{\operatorname{span}A_{\mathsf S}}
=A(I+A^TA)^{-1}A^T
=I-(I+AA^T)^{-1}.
\tag{2e}
\]

故

\[
T(w)=(I+AA^T)^{-1}
=\mathbb E\,P_{\ker A_{\mathsf S}^T}
=\sum_S\Pr(\mathsf S=S)P_S.
\tag{2f}
\]

Jensen 给出 (2c) 的“\(\le\)”方向。反过来，固定 \(S\)，
\((I+t\bar B_S\bar B_S^T)^{-1}\) 在
\(\operatorname{span}\bar B_S\) 上趋于零、在其正交补上恒等，因此趋于 \(P_S\)。
连续性给出反向不等式。证毕。

这里 (2e)--(2f) 是已知 DPP/volume-sampling 数学，不能作为本文原创性；候选新增
部分只是“对整个 weight orthant 的 convex supremum + every-flat converse + Schur
minimax”这一层。Kassel--Lévy (2023) 的 mean projection theorem 是更广的相邻结果，
也进一步说明不能把“随机投影的均值”包装成新概念。

### 异质 grounding 与共同 Schur 规则

令

\[
\bar C=C\Lambda^{-1/2},\qquad
\bar F=\Lambda^{-1/2}F.
\]

则对每个 \(P_S\) 定义

\[
A^{\Lambda}_{S}(Q)=
\left[
\bar C P_S\bar F-Q,
-\bar C P_S\Lambda^{-1/2}
\right].
\tag{2g}
\]

因为 Schur error 是 \(T\) 的 affine map，定理 2.1 立即给出

\[
\boxed{
\inf_{Q\in\mathcal Q_{\rm loc}}
\sup_{w\ge0}\|\mathcal E_w(Q)\|_2
=\inf_{Q\in\mathcal Q_{\rm loc}}
\max_{S\subseteq[m]}\|A^{\Lambda}_{S}(Q)\|_2.}
\tag{2h}
\]

所以固定 vertex/factor family 时，heterogeneous unary grounding 不破坏 exactness；
更一般的 fixed SPD affine shift \(K_0\) 也不破坏，只须取
\(\Lambda=K_0\)。它们只把 ordinary nullspace projections 变成 whitening 后的
nullspace projections。若 \(K_0\) 本身由已知 unary/pairwise factors 合成，模型仍完全
保持稀疏 WLS 解释。
当 \(\Lambda=\mu I\) 时，(2g) 退化为后文的
\([\mu^{-1}CP_SF-Q,-\mu^{-1}CP_S]\)。

## 3. Graph incidence 特例：精确归约到连通分割

给定无向宿主图 \(H=(\Gamma,E_H)\)。允许 completion 在每条宿主边上选择任意
\(w_e\ge0\)，并令

\[
L_H(w)=\sum_{e\in E_H}w_e b_eb_e^T,
\qquad
T_H(w)=\mu(\mu I+L_H(w))^{-1}.
\tag{2}
\]

对集合分割 \(\pi=\{B_1,\ldots,B_s\}\) 定义块平均投影

\[
P_\pi=\sum_{B\in\pi}\frac{\mathbf 1_B\mathbf 1_B^T}{|B|}.
\tag{3}
\]

它是投到“每个块上为常数”的子空间的正交投影。记

\[
\Pi_c(H)=\{\pi:\ H[B]\text{ 对每个 }B\in\pi\text{ 都连通}\}.
\]

### 推论 3.1（固定宿主图的 sharp convex extremum）

对任意连续凸函数
\(\Phi:\mathbb S^q\to\mathbb R\)，有

\[
\boxed{
\sup_{w\ge0}\Phi\bigl(T_H(w)\bigr)
=\max_{\pi\in\Pi_c(H)}\Phi(P_\pi).}
\tag{4}
\]

而且对每个 \(\pi\in\Pi_c(H)\)，存在只保留各块内一棵 spanning tree
的权重序列，使 \(T_H(w)\to P_\pi\)。因此等式右边没有加入不能由网络逼近的虚假 extreme points。

#### 证明

Chebotarev--Shamis 的 weighted matrix-forest theorem 说明

\[
(I+L)^{-1}_{ij}
=\frac{\text{使 }i,j\text{ 位于同一棵以 }j\text{ 为根之树的 rooted forests 总权重}}
{\text{全部 spanning rooted forests 总权重}}.
\tag{5}
\]

将 rooted forest 按其 underlying unrooted forest \(\mathcal F\) 分组。
给定 \(\mathcal F\) 的一个连通块 \(B\)，其根在 \(B\) 中均匀选择；所以该
unrooted forest 对矩阵的条件贡献恰为 \(P_{\pi(\mathcal F)}\)。因此存在概率
\(\lambda_{\mathcal F}\) 使

\[
T_H(w)=\sum_{\mathcal F}\lambda_{\mathcal F}
P_{\pi(\mathcal F)},
\qquad
\lambda_{\mathcal F}\ge0,\qquad
\sum_{\mathcal F}\lambda_{\mathcal F}=1.
\tag{6}
\]

每个 forest component 在 \(H\) 中连通，所以
\(\pi(\mathcal F)\in\Pi_c(H)\)。Jensen 不等式给出 (4) 的“\(\le\)”方向。

反过来，固定 \(\pi\in\Pi_c(H)\)，在每个 \(H[B]\) 中选一棵 spanning
tree；把这些 tree edge 的权重统一置为 \(t\)，其余边置零。随着
\(t\to\infty\)，在每个块的常数正交补上，
\(\mu(\mu I+tL)^{-1}\to0\)，在常数子空间上恒为 identity，故
\(T_H(w)\to P_\pi\)。连续性给出反向不等式。证毕。

### 对共同局部规则的直接推论

固定 \(Q\)，函数

\[
\Phi_Q(T)=
\left\|
\left[\mu^{-1}CTF-Q,\,-\mu^{-1}CT\right]
\right\|_2
\tag{7}
\]

是 \(T\) 的连续凸函数。因此

\[
\boxed{
\inf_{Q\in\mathcal Q_{\rm loc}}
\sup_{w\ge0}\|\mathcal E_{H,w}(Q)\|_2
=
\inf_{Q\in\mathcal Q_{\rm loc}}
\max_{\pi\in\Pi_c(H)}
\left\|
\left[\mu^{-1}CP_\pi F-Q,\,-\mu^{-1}CP_\pi\right]
\right\|_2.}
\tag{8}
\]

这里没有交换不合法的 \(\inf\) 与 \(\sup\)：推论 3.1 对**每个固定
\(Q\)** 都给出等式，最后才在两边取同一个 \(\inf_Q\)。

若 \(H\) 是一棵树，则每个连通分割唯一对应于删除一组 tree edges，所以
场景数恰为 \(2^{q-1}\)。这已经给出一个真正由 cut size 控制、而不是由整个
exterior 大小控制的 exact formulation。

## 4. 任意隐藏节点：精确归约到部分分割

现在允许 completion 添加任意有限隐藏节点 \(Z\)。令

\[
V=\Gamma\mathbin{\dot\cup}Z,\qquad
K_G=\mu I_V+L_G,\qquad
\widetilde C=[C\ 0],\qquad
\widetilde F=\begin{bmatrix}F\\0\end{bmatrix}.
\tag{9}
\]

隐藏节点的 primitive right-hand side 也属于 adversarial input，所以误差必须保留
全部 exterior 列：

\[
r_G(Q)=
\left\|
\left[
\widetilde C K_G^{-1}\widetilde F-Q,
-\widetilde C K_G^{-1}
\right]
\right\|_2.
\tag{10}
\]

这点很重要：只取 \((K_G^{-1})_{\Gamma\Gamma}\) 而删除 hidden-input columns
会换掉问题。下面的证明没有做这种偷换。

对任意 \(U\subseteq\Gamma\) 及 \(U\) 的一个分割 \(\pi\)，定义

\[
P_{\pi,U}=\sum_{B\in\pi}\frac{\mathbf1_B\mathbf1_B^T}{|B|},
\tag{11}
\]

并在 \(\Gamma\setminus U\) 上补零。记所有这类矩阵组成
\(\mathcal P_q^{\partial}\)。它们都是正交投影，共有

\[
|\mathcal P_q^{\partial}|
=\sum_{u=0}^q{q\choose u}B_u=B_{q+1}
\tag{12}
\]

个；最后一个等式可由添加一个 cemetery symbol 得到。

### 定理 4.1（任意隐藏 grounded completion 的精确最坏值）

令 \(\mathfrak G_q\) 包含所有以 \(\Gamma\) 为标记端口、允许任意有限隐藏节点和任意非负边权的有限图。则对每个固定 \(Q\)，

\[
\boxed{
\sup_{G\in\mathfrak G_q} r_G(Q)
=
\max_{P\in\mathcal P_q^{\partial}}
\left\|
\left[\mu^{-1}CPF-Q,\,-\mu^{-1}CP\right]
\right\|_2.}
\tag{13}
\]

更强地，(13) 的左边即使只允许**路径的 disjoint unions**，上确界也不变；
所以 lower witnesses 的最大图度为 2。

因此对任意非空 closed affine \(\mathcal Q_{\rm loc}\)，共同线性 minimax 是

\[
\boxed{
R_{\rm net}(C,F,\mu;\mathcal Q_{\rm loc})
=\min_{Q\in\mathcal Q_{\rm loc}}
\max_{P\in\mathcal P_q^{\partial}}
\left\|
\left[\mu^{-1}CPF-Q,\,-\mu^{-1}CP\right]
\right\|_2.}
\tag{14}
\]

由于 \(P=0\) 是一个场景，目标至少为 \(\|Q\|_2\)，故在 closed feasible
set 上不存在 \(Q\to\infty\) 的逃逸，minimum 达到。

#### 证明

**第一步：forest mixture。** 对固定 \(G\)，令
\(T_G=\mu K_G^{-1}\)。matrix-forest theorem 给出

\[
T_G=\sum_{\mathcal F}\lambda_{\mathcal F}P_{\rho(\mathcal F)},
\tag{15}
\]

其中 \(P_{\rho(\mathcal F)}\) 是 spanning forest 的 full-vertex component
averaging projection。因为 spectral norm of an affine map 是凸函数，
\(r_G(Q)\) 不超过这些 full partition projections 中的最大值。

**第二步：隐藏列被准确保留。** 固定一个 full partition projection
\(\widehat P\)。对其每个 component \(D\)，令
\(B=D\cap\Gamma\)，并忽略 \(B=\varnothing\) 的 component。其 port principal
block 为

\[
X=\widehat P[\Gamma,\Gamma]
=\sum_B\rho_B P_B,
\qquad
\rho_B=\frac{|B|}{|D|}\in(0,1].
\tag{16}
\]

由于 \(\widehat P^2=\widehat P\)，有两个精确恒等式

\[
\widetilde C\widehat P\widetilde F= C X F,
\qquad
(\widetilde C\widehat P)(\widetilde C\widehat P)^T=CXC^T.
\tag{17}
\]

所以 full hidden-input error 的平方恰为

\[
\Psi_Q(X)=\lambda_{\max}\!\left(
(\mu^{-1}CXF-Q)(\mu^{-1}CXF-Q)^T
+\mu^{-2}CXC^T
\right).
\tag{18}
\]

没有把第二项错误地写成 \(CX^2C^T\)；hidden columns 的能量正是通过
\(CXC^T\) 保存的。

**第三步：fractional blocks 不会比 0/1 blocks 更坏。** 映射
\(X\mapsto(\mu^{-1}CXF-Q)(\cdot)^T\) 是 matrix convex，因为

\[
\theta AA^T+(1-\theta)BB^T
-(\theta A+(1-\theta)B)(\theta A+(1-\theta)B)^T
=\theta(1-\theta)(A-B)(A-B)^T\succeq0.
\tag{19}
\]

加上线性的 \(CXC^T\) 再取 \(\lambda_{\max}\)，说明 \(\Psi_Q\) 是 convex。
另一方面，(16) 是 partial projections 的凸组合：独立地以概率 \(\rho_B\)
保留每个 block \(B\)，所得随机 partial projection 的期望正是 \(X\)。故

\[
\Psi_Q(X)\le \max_{P\in\mathcal P_q^{\partial}}\Psi_Q(P).
\tag{20}
\]

对 partial projection 有 \(P^2=P\)，所以 (18) 的平方根正好是 (13) 右侧的
operator norm。这证明上界。

**第四步：degree-2 sharp converse。** 对 \(P_{\pi,U}\)：

- 每个 active block \(B\in\pi\) 用一条经过 \(B\) 中全部 ports 的 path，边权取 \(t\to\infty\)。其 normalized resolvent 收敛到 \(P_B\)。
- 每个 inactive port 接一条含 \(h\) 个 hidden vertices 的独立 path，并先令 path conductance 足够大，再令 \(h\to\infty\)。该 component 在 port row 上的可见系数是 \(1/(h+1)\)，整条 hidden-input row 的 \(\ell_2\) norm 是 \(1/\sqrt{h+1}\)，两者都趋于零。

可一次取 \(t_h=\mu h^3\)：长度 \(h+1\) path 的第一非零 Laplacian
eigenvalue 为 \(\Theta(h^{-2})\)，故所有非恒定 mode 的 resolvent 也趋于零。
整个构造是 path forest，最大度 2，且 (10) 的 norm 收敛到指定 partial
projection 场景。下界和定理得证。

### 定理 4.2（异质 grounding 在归一化能量几何下仍然 exact）

固定端口 grounding
\(D_\Gamma=\operatorname{diag}(d_1,\ldots,d_q)\succ0\)。completion 可以加入任意
有限个 hidden vertices、任意正对角 \(D_Z\) 和任意非负 Laplacian edges；为实现
inactive ports，只要求 admissible family 包含总 grounding 可趋于无穷的 hidden
paths（例如可加入任意多个 unit-grounded hidden vertices）。写

\[
D_G=\operatorname{diag}(D_\Gamma,D_Z),\qquad
T_G=D_G^{1/2}(D_G+L_G)^{-1}D_G^{1/2}.
\tag{20a}
\]

在 grounding-normalized coordinates 中令
\(\widetilde C=[C\ 0]\)、
\(\widetilde F=\begin{bmatrix}F\\0\end{bmatrix}\)，并定义

\[
r_G^{(D)}(Q)=
\left\|
\left[
\widetilde C T_G\widetilde F-Q,
-\widetilde C T_G
\right]
\right\|_2.
\tag{20b}
\]

对 active block \(B\subseteq\Gamma\)，令

\[
v_B=D_\Gamma^{1/2}\mathbf1_B,\qquad
P_B^{(D)}=\frac{v_Bv_B^T}{v_B^Tv_B}
=\frac{(D_\Gamma^{1/2}\mathbf1_B)
(D_\Gamma^{1/2}\mathbf1_B)^T}
{\mathbf1_B^TD_\Gamma\mathbf1_B},
\tag{20c}
\]

并令 \(\mathcal P_q^\partial(D_\Gamma)\) 包含对任意 active subset 的每个分割
\(\pi\) 所形成的 \(P_{\pi,U}^{(D)}=\sum_{B\in\pi}P_B^{(D)}\)。则

\[
\boxed{
\sup_G r_G^{(D)}(Q)
=\max_{P\in\mathcal P_q^\partial(D_\Gamma)}
\left\|[CPF-Q,-CP]\right\|_2.}
\tag{20d}
\]

证明与定理 4.1 相同，但需要记录权重。一般 factor/DPP identity 把 \(T_G\) 写成
full-vertex weighted component projections 的凸组合。若 full component \(\mathcal D\) 与
ports 的交为 \(B\)，其 port principal block 为

\[
\frac{(D_\Gamma^{1/2}\mathbf1_B)
(D_\Gamma^{1/2}\mathbf1_B)^T}
{\sum_{j\in\mathcal D}d_j}
=\alpha_B P_B^{(D)},\qquad
\alpha_B=\frac{\sum_{i\in B}d_i}{\sum_{j\in\mathcal D}d_j}\in(0,1].
\tag{20e}
\]

独立地以概率 \(\alpha_B\) 激活每个 port block，又把其 port principal block 写成
weighted partial projections 的凸组合；(17)--(20) 的 matrix-convex argument 不变。
反向构造中，active blocks 用高导通 paths 融合；每个 inactive port 接一条总
grounding 趋于无穷的高导通 hidden path，使 \(\alpha_B\to0\)。故 sharp converse
仍只需 maximum degree 2。

若原始矩阵是 \(K_G=D_G+L_G\)，则 (20b) 对应于
\(C=C_{\rm raw}D_\Gamma^{-1/2}\)、
\(F=D_\Gamma^{-1/2}F_{\rm raw}\) 和 grounding-whitened input block
\(-\widetilde C_{\rm raw}K_G^{-1}D_G^{1/2}\)。uniform
\(D_G=\mu I\) 时全局尺度可吸收到 \(C\)，定理 4.1 正是它的未加权写法。

### 异质 grounding + 未加权 Euclidean RHS：最小 NO-GO 反例

上述输入几何不能静默删除。若第二块仍是原始
\(-\widetilde C_{\rm raw}K_G^{-1}\)，最小反例如下：只有一个 port，port grounding
为 3；加入一个 grounding 为 1 的 hidden vertex，并用 conductance \(t\) 的一条边
相连。取 \(C_{\rm raw}=F=R_i=1,Q=1/4\)，所以
\(C_{\rm raw}=R_iF^T\)。无 hidden vertex 与 inactive endpoint 的平方损失分别为
\(17/144\) 与 \(1/16\)，而

\[
\lim_{t\to\infty}
\left\|
\left[
e_1^TK_t^{-1}e_1-\frac14,
-e_1^TK_t^{-1}
\right]
\right\|_2^2
=\left\|\begin{bmatrix}0&-1/4&-1/4\end{bmatrix}\right\|_2^2
=\frac18>\frac{17}{144},
\qquad
K_t=
\begin{bmatrix}3&0\\0&1\end{bmatrix}
+t\begin{bmatrix}1&-1\\-1&1\end{bmatrix}.
\tag{20f}
\]

只用一个 port 的 raw inactive/no-hidden endpoints 会漏掉这个更坏的
hidden-input 场景。因此 heterogeneous grounding 的 **normalized theorem 是
exact**，但其未归一化 Euclidean-RHS analogue 被 (20f) 严格否定。

### 推论 4.3（state bias + raw noise + score disturbance 的 exact 场景式）

为避免把 score-space theorem 误写成 raw-noise theorem，考虑一个明确的混合输入模型。
保留定理 4.1 的 uniform grounding，并令完整 WLS normal matrix
\(J_G=H_G^TH_G\) 具有 fixed local blocks \(A,E\)；新增 hidden variables 在 \(E\) 中
对应零列。写

\[
\Delta_G(Q)=
\left[
\widetilde C K_G^{-1}\widetilde F-Q,
-\widetilde C K_G^{-1}
\right],\qquad
B_Q=\left[-QA,-QE-C\right].
\tag{20g}
\]

其中 (1g) 保证 \(B_Q\) 是 noiseless local-plus-port state 相对 centralized WLS 的
bias map。给定 state ellipsoid/covariance \(\Pi\succeq0\) 和 score/message disturbance
半径 \(\rho\ge0\)，定义 joint worst-energy central-approximation loss

\[
\mathscr R_G(Q;\Pi,\rho)=
\left\|
\left[
B_Q\Pi^{1/2},\ 
\Delta_G(Q)H_G^T,\ 
\rho\Delta_G(Q)
\right]
\right\|_2.
\tag{20h}
\]

三个 blocks 分别对应 noiseless state variation、whitened raw measurement noise 和在
normal-equation score/message 上的独立 bounded disturbance。由 (1e)，

\[
\mathscr R_G(Q;\Pi,\rho)^2
=\lambda_{\max}\!\left(
B_Q\Pi B_Q^T+QAQ^T
+\widetilde C K_G^{-1}\widetilde C^T
+\rho^2\Delta_G(Q)\Delta_G(Q)^T
\right).
\tag{20i}
\]

对 \(P\in\mathcal P_q^\partial\) 令

\[
\Delta_P(Q)=
\left[\mu^{-1}CPF-Q,-\mu^{-1}CP\right]
\]

以及

\[
\Omega_P(Q)=B_Q\Pi B_Q^T+QAQ^T
+\mu^{-1}CPC^T+\rho^2\Delta_P(Q)\Delta_P(Q)^T.
\]

则有精确等式

\[
\boxed{
\sup_{G\in\mathfrak G_q}\mathscr R_G(Q;\Pi,\rho)^2
=\max_{P\in\mathcal P_q^\partial}
\lambda_{\max}\!\left(\Omega_P(Q)\right).}
\tag{20j}
\]

证明只需在定理 4.1 的 matrix-convex quantity 中加入 fixed PSD term
\(B_Q\Pi B_Q^T+QAQ^T\) 和 linear term
\(\mu^{-1}CXC^T\)；score term乘 \(\rho^2\) 后仍是 matrix convex。degree-2
converse 同时让三个 terms 收敛，所以 (20j) 是 equality。

当 \(\Pi=0,\rho=0\) 时，(20j) 退化为纯 raw-noise comparison，\(Q=0\) 的退化与
(1f) 一致；当 \(\Pi\ne0\) 时 bias--variance tradeoff 已非平凡，而 \(\rho>0\) 又使
connected/partial scenarios 真正影响共同 \(Q\)。这给当前定理一个正确的 state-estimation
应用：它认证 local computation 对 centralized WLS 的 fidelity，并显式包含通信量化、
消息压缩或 score-domain 扰动；它仍不应冒充 local estimator 对 true state 的完整 MSE。

### 凸包的 exactness 与场景不可随意删除

定理 4.1 的证明同时给出如下固定维数的 port-principal convex-hull 等式：

\[
\overline{\operatorname{conv}}
\left\{\mu\bigl[(\mu I+L_G)^{-1}\bigr]_{\Gamma\Gamma}:
G\in\mathfrak G_q\right\}
=\operatorname{conv}\mathcal P_q^{\partial}
\tag{21}
\]

（完整 inverse 的维数会随 hidden-node 数变化，不能把它们直接放进同一个凸包。）
每个 \(P\in\mathcal P_q^{\partial}\) 都是右侧 polytope
的 exposed vertex。事实上，对任意另一个 projection \(P'\)，

\[
\langle 2P-I,P'\rangle
=\operatorname{tr}P-\|P-P'\|_F^2,
\tag{22}
\]

故线性泛函 \(X\mapsto\langle2P-I,X\rangle\) 唯一暴露 \(P\)。所以
\(B_{q+1}\) 个极端场景不是证明过程随意制造的冗余列表。特定的
\((C,F)\) 当然可能只激活其中少数，但不存在对**所有 convex tests** 通用的更小 vertex set。

## 5. 一个有限、精确的 necessary-and-sufficient 条件

先看一般固定 factor family。给定 \(\varepsilon\ge0\)，存在共同
\(Q\in\mathcal Q_{\rm loc}\) 使所有 \(w\ge0\) 的误差不超过 \(\varepsilon\)，
当且仅当

\[
\begin{bmatrix}
\varepsilon I_p&A_S^\Lambda(Q)\\
A_S^\Lambda(Q)^T&\varepsilon I_{k+q}
\end{bmatrix}\succeq0
\quad\text{对 represented matroid 的每个 flat }S.
\tag{22a}
\]

这是定理 2.1 的逐场景 Schur-complement 改写。一般 representable matroid 的 flat
数可能很大；仅从 rank \(q\) 得到的直接枚举界
\(\sum_{j=0}^q\binom mj\) 是 XP 而非 FPT，所以这里不虚报一般 matroid-FPT。
下面 graph-incidence family 才有只依赖 port cut size 的 Bell-number 界。

对 \(P\in\mathcal P_q^{\partial}\) 写

\[
A_P(Q)=\left[\mu^{-1}CPF-Q,\,-\mu^{-1}CP\right].
\tag{23}
\]

给定 \(\varepsilon\ge0\)，存在同一个 admissible local rule 使所有 grounded
completions 的误差不超过 \(\varepsilon\)，当且仅当存在
\(Q\in\mathcal Q_{\rm loc}\) 使

\[
\boxed{
\begin{bmatrix}
\varepsilon I_p&A_P(Q)\\
A_P(Q)^T&\varepsilon I_{k+q}
\end{bmatrix}\succeq0
\quad\text{对每个 }P\in\mathcal P_q^{\partial}.}
\tag{24}
\]

这是 graph-generated completion 类上的 exact iff，不是 sufficient
S-procedure bound。优化 \(\varepsilon\) 即得到 (14)。固定宿主图时只把
\(\mathcal P_q^{\partial}\) 换成 \(\{P_\pi:\pi\in\Pi_c(H)\}\)。
定理 4.2 的归一化异质-grounding 版本也完全相同：只须把 (23) 换成
\([CPF-Q,-CP]\)，并令 \(P\) 遍历
\(\mathcal P_q^\partial(D_\Gamma)\)。因此 weighted hidden theorem 同样给出有限
exact LMI iff，而不只是一个上界。

场景数只依赖 cut size：

- arbitrary hidden topology：\(B_{q+1}=2^{\Theta(q\log q)}\)；
- fixed tree host：\(2^{q-1}\)；
- fixed path host：同为 \(2^{q-1}\)，且可按 cut set 直接生成。

因此标准 SDP 算法可在
\(f(q)\operatorname{poly}(p+k,\log(1/\eta),\text{input bits})\) 时间内求到
精度 \(\eta\)：这是对 separator/cut size 的 exact FPT formulation。
本报告没有声称 Bell-size enumeration 在 \(q\) 上是最优复杂度，也没有证明
separation oracle 的 NP-hardness；该 complexity 下界仍是后续任务。

作为零误差特例，若通常的线性 locality feasible set 满足
\(0\in\mathcal Q_{\rm loc}\)，则因 \(P=0\) 迫使 \(Q=0\)，而 \(P=I\) 的
invisible block 迫使 \(C=0\)，所以

\[
R_{\rm net}=0\iff C=0.
\tag{25}
\]

若 \(R_i=e_s^T\) 选择单个 local target、\(A\) 是由局部图产生的 symmetric
nonsingular M-matrix，且 \(E\le0\) 是无符号抵消的标准负边耦合，则
\((A^{-1})_{sv}>0\) 当且仅当 \(s,v\) 位于同一 irreducible block。因此 \(C=0\)
等价于目标所在的 local component 没有任何 vertex 与 cut ports 耦合，即不存在局部
传递路径。若 \(R_i\) 或 \(E\) 允许 signed cancellation，这个图论等价不再自动成立；
代数 iff 仍是 \(C=0\)。非零容差的完整 iff 则是 (24)。

## 6. 优化共同 \(Q\) 后，中间割仍然必需

只检查“所有 edge weight 为 0”和“所有 edge weight 趋于无穷”很诱人，但即使
宿主图是一条三端口 path，这也不成立。

令宿主 path 的次序为 \(0-2-1\)，\(\mu=1\)，并取

\[
C=\begin{bmatrix}-1&-1&2\end{bmatrix},\qquad
F=\begin{bmatrix}
-1&3&3\\
-1&2&-1\\
2&3&-3
\end{bmatrix}.
\tag{26}
\]

注意 \(F\) 的第一列就是 \(C^T\)，所以取
\(R_i=e_1^T\) 时满足 WLS 必需关系 \(C=R_iF^T\)。该 path 的四个连通分割为

\[
I,\quad P_{\{0,2\}|\{1\}},\quad
P_{\{1,2\}|\{0\}},\quad P_{\{0,1,2\}}.
\]

对 scalar output，令

\[
a_P=CPF,\qquad b_P=CPC^T,\qquad
\|A_P(Q)\|_2^2=\|a_P-Q\|_2^2+b_P.
\tag{27}
\]

四个 \((a_P,b_P)\) 依次为

\[
\begin{array}{c|c|c}
P&a_P&b_P\\\hline
I&(6,1,-8)&6\\
P_{02|1}&(3/2,1,1)&3/2\\
P_{12|0}&(3/2,-1/2,-5)&3/2\\
P_{012}&(0,0,0)&0.
\end{array}
\tag{28}
\]

完整 minimax 的一个最优共同系数是

\[
Q^*=\left(\frac{317}{84},\frac{20}{21},-\frac{157}{42}\right),
\qquad
R_{\rm all\ cuts}^2=\frac{9785}{336}.
\tag{29}
\]

在 \(Q^*\) 处，\(I,P_{02|1},P_{012}\) 三个场景都达到
\(9785/336\)，而 \(P_{12|0}\) 的损失只有 \(3485/336\)。三条 active
quadratics 的 gradients 以正权

\[
\left(\frac{197}{378},\frac{163}{378},\frac1{21}\right)
\tag{30}
\]

加权后为零；权重和为一，所以 convex KKT 条件证明 (29) 是全局最优，而非
数值搜索猜测。

若错误地只检查 \(I\) 和 \(P_{012}\)，最优中心为

\[
Q_{\rm endpoints}=\frac{107}{202}(6,1,-8),
\qquad
R_{\rm endpoints}^2=\frac{11449}{404}
<\frac{9785}{336}.
\tag{31}
\]

所以中间 cut 在**已经优化共同规则以后**仍产生严格 gap。该构造同时是
推论 3.1 中“所有 connected partitions 都可能有用”的 matching
counterexample。

## 7. Petersen/full-block relaxation 的严格 gap

对任意 grounded graph，

\[
0\prec T_G=\mu(\mu I+L_G)^{-1}\preceq I.
\]

若忘掉 Laplacian/forest structure，只保留 \(0\preceq T\preceq I\)，就回到
self-adjoint full-block uncertainty，Petersen lemma 给出的 LMI 对这个**外包集**是
lossless。定理 4.1 表明网络 completion 的 convex extreme points 却只是
coordinate block-averaging projections，而不是任意 subspace projection。

取 \(q=2,\mu=1\)，

\[
C=\begin{bmatrix}1&0\end{bmatrix},\qquad
F=\begin{bmatrix}1&0\\0&\sqrt5\end{bmatrix},\qquad
R_i=\begin{bmatrix}1&0\end{bmatrix},\qquad Q=0.
\tag{32}
\]

这里 \(C=R_iF^T\)，所以该 gap 满足原始 WLS 的代数兼容关系，而不是任取
互不相关的 \(C,F\)。引理 1.1 还给出它的 two-sparse measurement realization。

两端口的五个 partial partition projections 是

\[
0,\quad e_1e_1^T,\quad e_2e_2^T,\quad I,
\quad \frac12\mathbf1\mathbf1^T.
\]

逐个计算给出网络精确值

\[
\sup_{G\in\mathfrak G_2}r_G(0)^2=2,
\qquad
R_{\rm net}=\sqrt2.
\tag{33}
\]

而 full interval 中的 rank-one projection

\[
T_*=\frac14
\begin{bmatrix}3&\sqrt3\\\sqrt3&1\end{bmatrix}
\tag{34}
\]

满足 \(T_*^2=T_*\) 且

\[
\left\|[CT_*F,-CT_*]\right\|_2^2=\frac94,
\qquad
R_{\rm full\ block}=\frac32>\sqrt2.
\tag{35}
\]

一般 rank-one projection \(uu^T\) 令 \(y=u_1^2\)，平方损失为

\[
2y^2+6y(1-y)=6y-4y^2,
\]

在 \(y=3/4\) 达到 \(9/4\)。full interval
\(\{0\preceq T\preceq I\}\) 的 extreme points 是 orthogonal projections；rank 0
与 rank 2 的值分别为 0 和 2，而上式已穷尽所有 rank-one projections，所以 (35)
也是 full interval 的精确值。

这准确说明结构限制破坏了什么：Petersen elimination 没有错；错的是把由
grounded sparse networks 生成的 uncertainty 当成整个 contraction ball。

## 8. 圆盘平面 completion 的拓扑版本

若所有 ports 按固定循环次序位于圆盘边界，completion 必须嵌入圆盘内部且边不交叉，则任何 spanning forest 的不同 components 不可能连接交替出现的 boundary terminals。因此它在 ports 上诱导的 partition 必须 noncrossing。

反过来，每个 noncrossing partition 都能由互不相交的 port paths 实现；inactive
ports 的 dilution paths 可放在各自边界小邻域内。重复定理 4.1 的证明得到

\[
\sup_{G\ {\rm circular\ planar}}r_G(Q)
=\max_{P\in\mathcal P_{q,\rm nc}^{\partial}}
\left\|[\mu^{-1}CPF-Q,-\mu^{-1}CP]\right\|_2,
\tag{36}
\]

其中 \(\mathcal P_{q,\rm nc}^{\partial}\) 只保留 active blocks noncrossing 的
partial partitions。其场景数为

\[
\sum_{u=0}^q{q\choose u}\operatorname{Cat}_u,
\tag{37}
\]

严格小于 \(B_{q+1}\) 从 \(q=4\) 开始。Curtis--Ingerman--Morrow 对 circular
planar response matrices 的 circular-minor characterization 是更完整的经典背景；
(36) 针对的是本文的 resolvent robust performance，而不是重新声称发明该
response-matrix characterization。

这条 planar specialization 很适合论文：它把“拓扑要求更低/更高”直接变成
可数的 worst-completion scenario family，而不是一句模糊的“planarity helps”。

## 9. 与已有理论的边界

### 9.1 DPP mean projection 与 random spanning forests

Dereziński--Khanna--Mahoney 的 Lemma 5 已逐式给出 (2e)：L-ensemble 所选列空间
投影的期望等于 ridge hat matrix。[^8] Kassel--Lévy 回顾并推广了另一条
determinantal mean-projection theorem。[^9] 因此一般 factor family 的 convex-mixture
identity 不是新数学。

在图情形，Chebotarev--Shamis 的 Theorems 5--7 已给出
\((I+L)^{-1}\) 的 rooted-forest entry formula。[^1] 更直接地，
Pilavcı--Amblard--Barthelmé--Tremblay (2021) Proposition 2 定义随机 partition
averaging matrix \(S\)，并明确证明
\(\mathbb E S=(Q+L)^{-1}Q\)。[^10] uniform \(Q=\mu I\) 时，这正是 (6) 的
partition-projection mixture。初版把这一步列为潜在新组合过于乐观，现已纠正。

仍未找到直接先例的是：对全部非负 weights 取 supremum 后的 every-flat/every-connected-
partition sharp converse、保留 hidden RHS 的 partial-partition 极值，以及这些结论与
common Schur minimax 的 exact LMI 组合。

### 9.2 Ratio-cut / connected-subpartition polytopes

矩阵 \(P_\pi=\sum_B|B|^{-1}\mathbf1_B\mathbf1_B^T\) 正是 clustering 文献中的
partition matrix。De Rosa--Khajavirad 研究其 fixed-cluster-number convex hull 的
off-diagonal projection，即 ratio-cut polytope，并由此构造 K-means LP relaxations。[^11]
Moura--Yaman--Leus 研究带宿主图连通约束的 labelled connected
\(k\)-subpartition incidence polytope及其 facets/separation complexity。[^12]

它们不直接给出本报告的 inverse-resolvent supremum 或 hidden-input norm，但术语和
polyhedral geometry 已有成熟归属；不能把“partition polytope”本身称为新对象。
反过来，这些工作也提示 Bell-size scenario enumeration 未必是最终算法。K-means 的
fixed-\(K\) linear optimization 已知 NP-hard，并**不能**自动推出本报告不限制 block 数、
且目标是 operator norm 的 worst-scenario evaluation NP-hard。当前没有完成保持目标形式的
归约，所以复杂性二分仍标为 PENDING。

### 9.3 Kron reduction 与 resistor-network DtN

Dörfler--Bullo 证明 loopy Laplacians 对 Schur/Kron reduction 闭合，并刻画消元后
boundary clique、self-loop 和 connectivity。[^2] 这保证 grounded network 是自然的
Schur completion class，但该文没有给出 (8)、(13) 或共同 local-rule minimax。

Curtis--Ingerman--Morrow 则刻画 circular planar response matrices：对称、行和为零且
所有 circular minors 非负是可实现性的 iff。[^3] 它解释 planarity 为什么比一般
M-matrix 增加真实代数约束；本报告的 noncrossing-partition extremum 与之兼容，
但不是同一个 inverse problem。

### 9.4 Sparse PSD completion

Grone--Johnson--Sá--Wolkowicz 的经典 theorem 说明 chordal specified-entry graph 上的
partial Hermitian positive-definite completion 条件，并刻画 maximum-determinant
completion 的 inverse zero pattern。[^4] 那是“给定部分 entries 是否存在 SPD
completion”的问题；这里是“由 Laplacian network 合成的 inverse family 对 convex
operator loss 的 worst point”。Chordality 不会直接推出分割归约。

### 9.5 Structured singular value \(\mu\)

Braatz--Young--Doyle--Morari 证明 pure real 或 mixed real/complex structured
singular-value recognition 是 NP-hard。[^5] 这警告多个 independent uncertainty
blocks 通常不能期待 Petersen-style exact scalar multiplier。但该 hardness 不应被
移植成本文的 claim：在 uniform grounded-Laplacian model 中，matrix-forest positivity
反而给出了 exact finite partition reduction。本报告没有证明 (14) 对变动 \(q\) 的
NP-hardness。

### 9.6 Robust approximate inverse 与 Petersen lemma

El Ghaoui 已定义 common approximate inverse minimax，并给 structured uncertainty
的 SDP bounds；single full block 的 exact elimination 又属于 Petersen/full-block
S-procedure 谱系。[^6][^7] 所以“写一个 common robust inverse”或“得到一个 full-block
LMI”都不是新意。本报告真正离开这条已知线的地方是 (35)：network-generated set 的
extreme geometry 是 graph partitions，full-block LMI 对它会严格保守。

## 10. 尚未覆盖的边界与 NO-GO 条件

这条结果足以作为候选主定理，但不能越过以下边界宣传。

1. **fixed vertices 与 arbitrary hidden vertices 的输入几何要分开。** 对 fixed factor
   family，定理 2.1 已允许任意 \(\Lambda\succ0\)，包括 heterogeneous unary grounding
   和 fixed SPD affine shift。对 arbitrary hidden nodes，定理 4.2 已在
   grounding-normalized input norm 下给出 weighted partial projections 的 exact theorem，
   前提是 dilution paths 的 hidden total grounding 可无界增长；但未加权 Euclidean RHS
   版本被一端口反例 (20f) 否定。
2. **finite edge upper bound 破坏 sharp converse。** 若 \(0\le w_e\le W_e<\infty\)，
   DPP/Jensen 仍给 flat projections 的上界，但 selected weights 无法趋于无穷。最小严格
   反例取 \(q=1,\mu=C=F=b=1,Q=2,w\in[0,1]\)。令
   \(T=(1+w)^{-1}\in[1/2,1]\)，则
   \[
   \sup_{0\le w\le1}\|[T-2,-T]\|_2^2=\frac52,
   \qquad
   \max_{P\in\{0,1\}}\|[P-2,-P]\|_2^2=4.
   \tag{37a}
   \]
   因而 projection LMI 此时只是严格保守的 sufficient certificate。
3. **signed/general pairwise factors：fixed library 已解决，arbitrary hidden 未解决。**
   定理 2.1 的 DPP proof 不需要 M-matrix，故 fixed two-sparse factor library 可以带任意
   signs。尚未刻画的是允许新增任意 hidden variables 和任意 signed pairwise factors 后，
   port 上究竟出现哪些 subspace/principal-block extremes。
4. **coercivity 不能消失。** fixed \(K_0\succ0\) 可直接作为 \(\Lambda\) whitening；若
   base 仅半正定且 completion 可保持其 kernel，inverse 不存在或 minimax 可发散，不能
   沿用定理 2.1。
5. **定理针对共同线性 Schur correction。** 任意 nonlinear decoder 的跨-completion
   Chebyshev radius 不是由 (24) 自动解决。
6. **FPT 不是 polynomial-in-cut。** \(B_{q+1}=2^{\Theta(q\log q)}\)。要把计算贡献再
   提升一档，应证明 separation 的复杂性，或在 treewidth/planarity 下给更好的 DP。

若论文坚持“每条 edge weight 有固定上界 + heterogeneous hidden grounding 的原始
Euclidean RHS + arbitrary hidden signed factors + 任意 nonlinear decoder”四项同时成立，
本报告的主定理就不能原样使用；在未得到新证明前必须写 **NO-GO**，不能靠术语包装。

## 11. 可量化应用收益与实验判据

### 11.1 把无限 completion 类压缩成 cut-size certificate

定理 4.1 的直接工程含义不是“又多一个 SDP”，而是：只要已知本地子问题与远端通过
\(q\) 个 ports 相连，就可把**任意远端节点数、任意非负边权、任意有环拓扑**的连续
completion 类，精确压缩成 \(B_{q+1}\) 个大小为
\((p+k+q)\times(p+k+q)\) 的 LMI。certificate 的规模与 hidden-node 数和远端边数
无关。固定宿主树时是 \(2^{q-1}\) 个场景；圆盘平面类则降为 (37) 的
noncrossing partial partitions。

这实质上放宽了 2011 模型常见的 acyclic-topology 前提：上界允许所有 finite grounded
graphs，包括任意多 cycles 和 hidden nodes；而 matching lower witnesses 甚至只需
maximum degree 2。一般 fixed factor library 又由定理 2.1 覆盖 arbitrary signed
two-sparse rows 和任意 fixed SPD base，而非只覆盖 incidence/M-matrix rows。

与三个自然基线的取舍是清楚的：

- full-block Petersen 通常只需较少 LMI，但对 network-generated set 可能严格保守；
- 只检查“无耦合/全融合” endpoints 更便宜，却不是可靠上界；
- 有限随机采样真实图可发现坏例子，但不能证明没有更坏的权重或拓扑。

flat/partition formulation 的价值是**非保守 certificate**，不是声称它总比单个
Petersen LMI 更快。若 \(q\) 很小而远端网络很大或未知，它以 cut-parameter
exponential 换掉 network-size dependence；若 \(q\) 很大，则必须进一步做 separation、
column generation 或利用 treewidth/planarity。

### 11.2 两个已经精确量化的 score-space gap

二端口例 (33)--(35) 中，在 (1) 的 normal-equation RHS metric 下，Petersen 外包
相对真实 network completion 的保守性为

\[
\frac{R_{\rm full\ block}^2}{R_{\rm net}^2}-1
=12.5000\%,\qquad
\frac{R_{\rm full\ block}}{R_{\rm net}}-1
=6.0660\%.
\tag{38}
\]

三端口 path 例 (29)--(31) 中，真实最坏值相对 endpoint-only 结果高出

\[
\frac{R_{\rm all\ cuts}^2}{R_{\rm endpoints}^2}-1
=2.7627\%,\qquad
\frac{R_{\rm all\ cuts}}{R_{\rm endpoints}}-1
=1.3719\%.
\tag{39}
\]

它们分别量化“外包过松”和“场景删得过多”的代价。有限权上界反例 (37a) 又说明
相反方向的风险：若物理上有 \(w\le W\)，仍使用无限权 projection certificate 会把
平方值从精确的 \(5/2\) 抬到 4，即在该最小例上高出 60%。

### 11.3 可扩展实验设计

后续 CPU/GPU 实验应同时报告以下量，而不是只画 estimator error：

1. **准确性：**
   \(R_{\rm Petersen}/R_{\rm exact}-1\)、
   \(R_{\rm exact}/R_{\rm endpoints}-1\)，以及 Monte Carlo 最大值到 exact
   certificate 的 sampling gap；
2. **计算：** 场景数、真正 active flats/partitions 数、SDP wall time、peak memory、
   column-generation iterations 与 separation time；
3. **扩展轴：** cut size \(q\)、factor 数 \(m\)、hidden-node 数、host treewidth、
   planarity、grounding condition number 和 weight upper bound \(W\)；
4. **结构收益：** Bell enumeration、tree/path enumeration、noncrossing enumeration、
   factor-flat 去重和 full-block Petersen 的 performance profile；
5. **模型一致性：** 用引理 1.1 生成 raw two-sparse WLS instances，并把
   normalized heterogeneous-input 结果与 (20f) 的 raw Euclidean failure 分开统计。

首轮可在小 \(q\) 上全枚举作为 ground truth，再逐步增加 \(q\) 比较 column generation。
这能同时验证理论 exactness、保守性收益和计算代价，满足“数学主定理必须导出可测应用
收益”的门槛。

## 12. 作为论文主定理的评价与下一步

### 目前可以稳妥写的 claim

- 对 fixed SPD base 与任意 fixed factor library，weight-orthant 上的 worst convex
  Schur loss 精确归约到 representable-matroid flat/nullspace projections；
- uniform-grounded arbitrary hidden graph 的 exact partial-partition theorem，及
  heterogeneous grounding 在 normalized energy geometry 下的 weighted 版本；
- 两个 hidden theorems 的 sharp converse 都只需 maximum-degree-2 path forests；
- common local rule 的 \(\varepsilon\)-可行性有 finite exact LMI iff，graph case 对 cut
  size 是 exact FPT formulation；
- raw two-sparse WLS normal-matrix realization（不声称输入范数等价）、intermediate-cut
  necessity、Petersen strict gap、
  heterogeneous raw-RHS NO-GO 和 bounded-weight NO-GO 都有显式 matching examples。

### 仍不能写的 claim

- “DPP mean-projection、random-forest partition mixture 或 partition polytope 是本文原创”；
- “一般 matroid-flat worst-scenario evaluation 已证明 NP-hard/不可近似”；
- “Bell-size enumeration 在参数 \(q\) 上最优”或“一般 factor case 是 FPT”；
- “bounded weights 仍由 projections exact 表示”；
- “heterogeneous hidden grounding 对未归一化 Euclidean RHS 仍 exact”；
- “normal-equation RHS gain 等同于 raw WLS measurement-noise gain”；
- “arbitrary hidden signed factors 或 nonlinear decoders 已被彻底解决”。

### 建议的 theorem chain 与 GO/NO-GO

论文应把引理 1.1（仅作 realizability）、定理 2.1、定理 4.1--4.2 与
(22a)/(24) 作为一个 theorem chain，
再用 (29)--(31)、(33)--(35)、(20f) 和 (37a) 给 matching boundaries。标题可直接写：

> Exact worst-completion reduction for pairwise-factor graph extensions.

已知 DPP/forest identity 使定理 2.1 的上界成为很短的 Jensen corollary；因此不能只靠
“凸组合”承担原创性。真正可能过主定理门槛的是 sharp every-flat converse、hidden
partial-partition theorem、degree-2 realization、Schur-minimax iff 及 strict application
gaps 的组合。当前判断是：

> **作为 theorem package：GO；作为单独一条顶刊数学新定理：边缘，必须继续逐式查重。**

原创性置信度维持在约 **0.45**。下一步优先做：

1. MathSciNet/zbMATH、cited-by、博士论文和 resistor-network synthesis 文献的逐定理查重；
2. 严格复杂性研究：只在完成保持本报告 operator-norm objective 的 reduction 后声称
   NP-hardness，并与 treewidth/cut-size DP 或 separation algorithm 配对；
3. 刻画 arbitrary hidden signed factors 的 port subspace extremes；
4. 实现 exact SDP + column generation，并按第 11.3 节跑规模化 CPU/GPU benchmark；
5. 若应用最终要求 raw sensor observation norm，而非 normal-equation RHS norm，显式带回
   \(H^TR_z^{-1}\) 的输入映射后重新比较所有 bounds。

## 13. 本地验证

可复现实验：

```powershell
python experiments\theory_search\forest_partition_extrema.py
python experiments\theory_search\hidden_partial_partition_checks.py
```

第一份脚本验证九项：matrix-forest identity、高导通极限、优化后 intermediate-cut
反例、Petersen gap、一般 DPP projection identity、raw pairwise-WLS realization、
raw-measurement metric orthogonalization、heterogeneous-hidden Euclidean obstruction 和
bounded-weight gap。第二份脚本验证五项：uniform/heterogeneous random hidden-graph
upper bounds、mixed state/raw/score corollary，以及两种 degree-2 dilution converses；
所有 hidden-input columns 均被保留。本次运行 **14 项全部通过**。

## Sources

[^1]: Pavel Chebotarev and Elena Shamis, “[Matrix-Forest Theorems](https://arxiv.org/abs/math/0602575),” 2006 manuscript; see Theorems 5--7 for weighted undirected multigraphs and \((I+L)^{-1}\). Earlier journal version: P. Yu. Chebotarev and E. V. Shamis, “[A Matrix-Forest Theorem and Measuring Relations in Small Social Groups](https://www.mathnet.ru/eng/at2672),” *Automation and Remote Control* 58(9) (1997), 1505--1514.

[^2]: Florian Dörfler and Francesco Bullo, “[Kron Reduction of Graphs with Applications to Electrical Networks](https://doi.org/10.1109/TCSI.2012.2215780),” *IEEE Transactions on Circuits and Systems I* 60(1) (2013), 150--163; see Lemma III.3 and Theorem III.4. [Author manuscript](https://motion.me.ucsb.edu/pdf/2011d-db.pdf).

[^3]: Edward B. Curtis, David Ingerman, and James A. Morrow, “[Circular Planar Graphs and Resistor Networks](https://doi.org/10.1016/S0024-3795(98)10087-3),” *Linear Algebra and its Applications* 283 (1998), 115--150; in particular the response-matrix characterization by circular minors.

[^4]: Robert Grone, Charles R. Johnson, Eduardo M. Sá, and Henry Wolkowicz, “[Positive Definite Completions of Partial Hermitian Matrices](https://doi.org/10.1016/0024-3795(84)90207-6),” *Linear Algebra and its Applications* 58 (1984), 109--124.

[^5]: Richard P. Braatz, Peter M. Young, John C. Doyle, and Manfred Morari, “[Computational Complexity of \(\mu\) Calculation](https://doi.org/10.1109/9.284107),” *IEEE Transactions on Automatic Control* 39(5) (1994), 1000--1002. [Author manuscript](https://web.mit.edu/braatzgroup/6_Computational%20complexity%20of%20mu%20calculation.pdf).

[^6]: Laurent El Ghaoui, “[Inversion Error, Condition Number, and Approximate Inverses of Uncertain Matrices](https://doi.org/10.1016/S0024-3795(01)00273-7),” *Linear Algebra and its Applications* 343--344 (2002), 171--193.

[^7]: H. J. van Waarde, M. Kanat Camlibel, J. Eising, and H. L. Trentelman, “[Quadratic Matrix Inequalities with Applications to Data-Based Control](https://doi.org/10.1137/22M1486807),” *SIAM Journal on Control and Optimization* 61(4) (2023), 2251--2281; see Proposition 4.16 for the nonstrict Petersen lemma and Proposition 4.14 for the single-full-block lossless setting.

[^8]: Michał Dereziński, Rajiv Khanna, and Michael W. Mahoney, “[Improved Guarantees and a Multiple-Descent Curve for Column Subset Selection and the Nyström Method](https://proceedings.neurips.cc/paper/2020/hash/342c472b95d00421be10e9512b532866-Abstract.html),” *Advances in Neural Information Processing Systems* 33 (2020); see Lemma 5 and Appendix C. [arXiv:2002.09073](https://arxiv.org/abs/2002.09073).

[^9]: Adrien Kassel and Thierry Lévy, “[On the Mean Projection Theorem for Determinantal Point Processes](https://doi.org/10.30757/ALEA.v20-17),” *ALEA, Latin American Journal of Probability and Mathematical Statistics* 20 (2023), 497--504.

[^10]: Yusuf Yiğit Pilavcı, Pierre-Olivier Amblard, Simon Barthelmé, and Nicolas Tremblay, “[Graph Tikhonov Regularization and Interpolation via Random Spanning Forests](https://doi.org/10.1109/TSIPN.2021.3084879),” *IEEE Transactions on Signal and Information Processing over Networks* 7 (2021), 359--374; see Proposition 2 and equations (19)--(21). [arXiv:2011.10450](https://arxiv.org/abs/2011.10450).

[^11]: Antonio De Rosa and Aida Khajavirad, “[The Ratio-Cut Polytope and K-Means Clustering](https://doi.org/10.1137/20M1348601),” *SIAM Journal on Optimization* 32(1) (2022), 173--203.

[^12]: Phablo F. S. Moura, Hande Yaman, and Roel Leus, “[On the Connected (Sub)partition Polytope](https://doi.org/10.1007/s10107-025-02321-1),” *Mathematical Programming* (accepted/online-first; 2026 bibliographic year). [arXiv:2401.01716](https://arxiv.org/abs/2401.01716).
