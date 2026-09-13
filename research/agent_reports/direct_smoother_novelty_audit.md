# 任意隐藏 completion 下 graph-Tikhonov 局部 smoother：证明、查重与敌对审稿

## 0. 结论先行

候选结论经过独立逐步核验后，结论是 **GO，但必须同时收紧模型表述、承认已知基础、并把设计变量 (Q) 的约束写进主问题**。

最稳妥的数学结论是：固定 (q) 个 ports、固定正定 port base (Lambda\succ0)，令任意有限 completion (G) 的

\[
H_G=\operatorname{diag}(\Lambda,D_Z),\qquad
T_G=H_G^{1/2}(H_G+L_G)^{-1}H_G^{1/2},
\]

其中 (D_Z\succ0) diagonal，(L_G) 是无向非负加权图 Laplacian。若 completion 类允许任意多 hidden vertices、无界总 hidden grounding、无界边导通，并允许自由选择 ports 之间的未知边，则对任意固定 (C\in\mathbb R^{m\times q})、(Q\in\mathbb R^{m\times q})，

\[
\boxed{
\sup_G
\left\|C[T_G]_{\Gamma,:}-[Q,0_Z]\right\|_2
=
\max_{P\in\mathcal P_\Lambda^\partial}
\|CP-Q\|_2 .}
\tag{1}
\]

这里 (\mathcal P_\Lambda^\partial) 是所有 partial partitions 所生成的有限 projector family；若
(\pi=\{B_1,\ldots,B_s\}) 是某个 active subset (U\subseteq\Gamma) 的 partition，令

\[
R_\pi=[\mathbf1_{B_1}\ \cdots\ \mathbf1_{B_s}],\qquad
P_{\pi,U}^{\Lambda}
=\Lambda^{1/2}R_\pi
(R_\pi^T\Lambda R_\pi)^{-1}
R_\pi^T\Lambda^{1/2}.
\tag{2}
\]

它是到
(\operatorname{span}\{\Lambda^{1/2}\mathbf1_B:B\in\pi\})
的 Euclidean orthogonal projection。场景数为

\[
|\mathcal P_\Lambda^\partial|
=\sum_{u=0}^q{q\choose u}B_u=B_{q+1}.
\tag{3}
\]

每个下界场景都可由最大度数 (2) 的 path forests 逼近；若要求每个有限 completion 连通，可用趋于零的弱边把各 component 串成一条路径。因此 (1) 的左侧通常是 **supremum 而非 maximum**。

另外两个结论同样通过了核验：

1. 对完整 hidden-input operator，(1) 的 Schatten-(s) 版本在所有实例上一致成立，当且仅当 (s\ge2)（限于真正的 norms (1\le s\le\infty)）。(1\le s<2) 有一个 (q=1)、一个 hidden vertex 的最小反例。
2. 若 (Q) 在整个 (\mathbb R^{m\times q}) 中自由，则 minimax 设计完全退化：

   \[
   Q^*=C/2,
   \qquad
   \min_Q\max_{P\in\mathcal P_\Lambda^\partial}\|CP-Q\|_2
   =\|C\|_2/2.
   \tag{4}
   \]

   所有 projector 场景在 (C/2) 处甚至严格等距。因而论文应区分：**fixed-(Q) exact certificate 是非平凡定理；自由 (Q) 的 design 没有应用内容；只有真实的 local/shared/unbiased constraint (Q\in\mathcal Q_{\rm loc}) 才形成非平凡优化问题。**

查重方面，random-forest smoother、DPP mean projection、normalized partition matrices、Kron reduction、grove boundary partitions、cactus compactification 都已有强先例。[^1][^2][^3][^4][^5][^6][^7] 截至 2026-09-14，本次对 primary sources 的定向全文核查没有找到同时给出“任意 hidden dimension 的 grounded resolvent port hull + 完整输入 norm equality + degree-2 converse + constrained local minimax”的论文。这个组合的原创性置信度评为 **中等，约 0.60**；它足以支撑应用数学/GSP/network estimation 论文，但绝不能写成“首次发现 forest mixture”或“彻底终结一般 distributed state estimation”。

## 1. 精确量词：论文中不能省略的 completion class

令 ports 为固定有标签集合

\[
\Gamma=\{1,\ldots,q\}.
\]

一个 admissible completion 选择任意有限 hidden set (Z)、任意 positive diagonal (D_Z)，以及
(V=\Gamma\mathbin{\dot\cup}Z) 上任意无向、非负加权、可含环的 graph Laplacian (L_G)。port base
(\Lambda\in\mathbb S_{++}^q) 在所有 completions 中固定。令

\[
E_\Gamma=[I_q\ 0],\qquad
\Delta_G(C,Q)
:=CE_\Gamma T_G-[Q,0_Z].
\tag{5}
\]

式 (1) 是以下量词下的 universal equality：

- (q,m,\Lambda,C,Q) 固定；
- hidden cardinality ( |Z| ) 只要求每个实例有限，但在 supremum 中没有统一上界；
- completion 类至少包含任意长、具有固定正 hidden grounding 的 paths，因此 hidden total grounding 可趋于无穷；
- strong edges 的 conductance 没有统一有限上界；
- 为实现 active blocks，允许直接加入 port--port edges；
- 等式对所有 completions 的 union 取 supremum，不是说每个固定 (G) 的最坏值都由某个 partition 给出。

“graph extension”若是指**保留一个固定可见子图并且只在边界外加 hidden nodes**，则上述 converse 一般不成立。例如固定 port--port edge 不能在 (P=I) 场景中被删除；若禁止新增 port--port edge、同时 hidden grounding 有正下界，则把两个 ports 通过 hidden intermediary 熔合会引入额外 grounding，通常也到不了 (2)。对这类受限 completion，(1) 的右端仍给 universal upper bound，但不能称为 exact uncertainty set。

还应明确：这里的 (D_Z) 或 (H_G) 是 Tikhonov/data-fidelity 的 positive base，不是必须等于 graph degree matrix。若把它绑死为由 (L_G) 决定的 degree，证明与可达性都需重做。

## 2. 第一次 DPP：任意 SPD port base 下的 forest mixture

写

\[
L_G=BW B^T,
\qquad
Y=H_G^{-1/2}BW^{1/2}.
\]

则

\[
T_G=(I+YY^T)^{-1}.
\tag{6}
\]

对 columns of (Y) 定义 (L)-ensemble DPP：

\[
\Pr(S)=
\frac{\det(Y_S^TY_S)}{\det(I+Y^TY)}.
\]

DPP mean-projection identity 给出[^2][^3]

\[
\mathbb E\,\operatorname{Proj}(\operatorname{span}Y_S)
=Y(I+Y^TY)^{-1}Y^T
=I-T_G.
\]

所以

\[
T_G
=\mathbb E\,\Pi_S,
\qquad
\Pi_S:=\operatorname{Proj}(\ker Y_S^T).
\tag{7}
\]

只有 linearly independent edge subsets 才有正概率。incidence columns 的独立集恰是 forests，因此 (7) 确实是 forest projections 的 convex mixture；whitening by (H_G^{-1/2}) 可逆，不改变线性独立性。

若 forest (F) 的 component-indicator matrix 是 (U_F)，则

\[
\ker Y_F^T
=H_G^{1/2}\ker B_F^T
=\operatorname{range}(H_G^{1/2}U_F),
\]

从而

\[
\Pi_F
=H_G^{1/2}U_F
(U_F^TH_GU_F)^{-1}
U_F^TH_G^{1/2}.
\tag{8}
\]

这一层本身不是新的：Pilavcı 等已经对 diagonal grounding 逐式证明
((L+D)^{-1}D) 是 random forest component-wise weighted averaging matrix 的期望；乘以 (D^{1/2}) 与 (D^{-1/2}) 就得到 diagonal 版本的 (7)。[^1]

## 3. 第二次 DPP：forest port block 为什么恰落入 partial-projector hull

丢掉没有 terminal 的 forest components。设剩下的 components 在 ports 上的非空交集为
(B_1,\ldots,B_r)，令

\[
R=[\mathbf1_{B_1}\ \cdots\ \mathbf1_{B_r}],
\qquad
h_j=\sum_{z\in C_j\cap Z}(D_Z)_{zz}.
\]

由 (8)，forest projector 的 terminal compression 是

\[
X_F:=[\Pi_F]_{\Gamma,\Gamma}
=V(V^TV+\operatorname{diag}h)^{-1}V^T,
\qquad V=\Lambda^{1/2}R.
\tag{9}
\]

若所有 (h_j>0)，令 (A=V\operatorname{diag}(h)^{-1/2})。再次应用 mean-projection identity：

\[
X_F
=A(I+A^TA)^{-1}A^T
=\mathbb E_J\operatorname{Proj}(\operatorname{span}V_J).
\tag{10}
\]

正 scaling 不改变 column span。选取 (J\subseteq[r]) 等于从 forest components 中保留若干个，其 terminal intersections 两两不交，因此右侧 projectors 正是 (2)。若某些 (h_j=0)，把它们替为 (h_j+\varepsilon)，再令 (\varepsilon\downarrow0)；因为 (q) 固定且 projector family 有限，其 convex hull 已闭，极限仍在其中。于是

\[
X_F\in\operatorname{conv}\mathcal P_\Lambda^\partial.
\tag{11}
\]

dense (\Lambda) 下有一个容易写错的细节：不同 block vectors
(\Lambda^{1/2}\mathbf1_{B_j}) 一般**不正交**，所以 (2) 不能写成若干 rank-one block projectors 的和。正确对象始终是到它们**联合 span**的 projector。对 diagonal
(\Lambda=D_\Gamma)，支撑不交恢复正交，(2) 简化为

\[
P_{\pi,U}^{D_\Gamma}
=\sum_{B\in\pi}
\frac{D_\Gamma^{1/2}\mathbf1_B\mathbf1_B^TD_\Gamma^{1/2}}
{\mathbf1_B^TD_\Gamma\mathbf1_B}.
\tag{12}
\]

uniform case 就是熟悉的
(\sum_B\mathbf1_B\mathbf1_B^T/|B|)。

由第一次 mixture 与 (11)，立即得到 terminal-level closed hull：

\[
\overline{\operatorname{conv}}
\{[T_G]_{\Gamma,\Gamma}:G\}
=\operatorname{conv}\mathcal P_\Lambda^\partial,
\tag{13}
\]

其中反向 inclusion 由第 5 节的 paths 给出。因而任何只依赖 terminal block 的 continuous convex loss 都可精确归约到有限场景。完整 hidden-input error 则还需要下一节的 Gram 论证，不能只引用 (13)。

## 4. 从 forest mixture 到完整输入 operator norm

令

\[
\bar C=[C,0_Z],\qquad \bar Q=[Q,0_Z].
\]

对 full forest projector (\Pi_F)，定义

\[
A_F=\bar C\Pi_F-\bar Q.
\]

关键是保留 hidden columns 后仍有

\[
\begin{aligned}
A_FA_F^T
&=CX_FC^T-CX_FQ^T-QX_FC^T+QQ^T\\
&=:M_Q(X_F).
\end{aligned}
\tag{14}
\]

这里用了 (\Pi_F^2=\Pi_F)，所以所有 hidden-column contribution 已经包含在
(CX_FC^T) 中，并没有被丢掉。更重要的是，(M_Q(X)) 对 (X) affine。若第二次 DPP 给出

\[
X_F=\sum_j\theta_jP_j,
\qquad P_j\in\mathcal P_\Lambda^\partial,
\]

则

\[
M_Q(X_F)
=\sum_j\theta_j(CP_j-Q)(CP_j-Q)^T.
\tag{15}
\]

谱范数满足

\[
\|A_F\|_2^2
=\lambda_{\max}(M_Q(X_F))
\le \max_j\|CP_j-Q\|_2^2.
\tag{16}
\]

第一次 DPP 又有

\[
\Delta_G(C,Q)=\mathbb E_F A_F,
\]

故由 norm convexity 得到 (1) 的 upper bound。这个两层结构不可互换：第一层把实际 (T_G) 写成 full forest projectors 的均值；第二层只把每个 forest 的 terminal Gram 压缩为 partial scenarios。

## 5. sharp converse、closure 与 degree-2 paths

固定目标 partial partition ((U,\pi))。构造 strong-edge forest：

- 对每个 active block (B\in\pi)，用只经过 (B) 中 ports 的一条 path 把这些 ports 连起来；
- 对每个 inactive port (i\notin U)，从 (i) 接一条含 (h) 个 unit-grounded hidden vertices 的 path；
- 所有 strong edges 使用 conductance (\kappa)。

固定 (h) 并令 (\kappa\to\infty)，由 (6) 的谱分解，(T_G) 趋于该 forest 的 projector (8)。把 active 与 inactive component columns 分块为
(R=[R_A,R_I])，其 component Gram matrix 为

\[
\begin{bmatrix}
R_A^T\Lambda R_A&R_A^T\Lambda R_I\\
R_I^T\Lambda R_A&R_I^T\Lambda R_I+hI
\end{bmatrix}.
\tag{17}
\]

令 (h\to\infty)，block inverse 给出 terminal--terminal block 收敛到

\[
\Lambda^{1/2}R_A(R_A^T\Lambda R_A)^{-1}R_A^T\Lambda^{1/2}
=P_{\pi,U}^{\Lambda}.
\]

terminal--hidden block 的每个 inactive component 是一个 (O(h^{-1})) column 重复 (h) 次，所以其 Frobenius norm 与 operator norm都是 (O(h^{-1/2}))。因此

\[
\left\|
C[T_G]_{\Gamma,:}-[CP_{\pi,U}^{\Lambda},0_Z]
\right\|_{S_s}\longrightarrow0
\quad(s\ge2),
\tag{18}
\]

其中对每个 (h) 再选足够大的有限 (\kappa_h)，形成 diagonal sequence。维数随 (h) 改变不是漏洞：在每个有限维空间里与相同宽度的 padded target 比较，(18) 直接控制 norm values。

若 completion 必须 connected，把上述 path components 按端点用 conductance
(\delta_h\downarrow0) 的弱边串起来。每个有限 graph 都是一条 connected path、maximum degree (2)；先用 inverse continuity 取足够小的 (\delta_h)，再取足够大的 (\kappa_h)。于是每个右端场景均是实际 completions 的 norm-limit，结合 upper bound 得到 equality。

这一 converse 同时说明几个量词：

- projection 通常要求 (\kappa\to\infty)，所以不能把 supremum 写成 attained maximum；
- (P=0) 需要 inactive dilution，若 hidden size 或 hidden total grounding 有界，通常不可达；
- 若 conductance 有统一上界，active fused projectors 通常不可达；
- 若 ports 之间不能加边而 hidden grounding 有固定正下界，active block path 也未必可达。

## 6. Schatten norm 的 sharp threshold：(s\ge2) 当且仅当

对 (2\le s<\infty)，有

\[
\|A_F\|_{S_s}^2
=\|A_FA_F^T\|_{S_{s/2}}.
\]

因为 (s/2\ge1)，Schatten-(s/2) 是 convex norm。由 (15)，

\[
\|A_F\|_{S_s}^2
\le\sum_j\theta_j
\|(CP_j-Q)(CP_j-Q)^T\|_{S_{s/2}}
\le\max_j\|CP_j-Q\|_{S_s}^2.
\tag{19}
\]

第一次 mixture 再用 Schatten norm convexity；第 5 节给 lower bound。因此

\[
\sup_G\|\Delta_G(C,Q)\|_{S_s}
=\max_{P\in\mathcal P_\Lambda^\partial}\|CP-Q\|_{S_s},
\qquad 2\le s\le\infty.
\tag{20}
\]

阈值是 sharp。取一个 port、一个 hidden vertex、(H=I_2)，并令两者高导通融合，则 full projector 极限为

\[
\Pi=\frac12\begin{bmatrix}1&1\\1&1\end{bmatrix}.
\]

取

\[
C=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
Q=\begin{bmatrix}1\\0\end{bmatrix}.
\]

一个 port 的 partial scenarios 只有 (P=0,1)，两端误差分别为 (-Q) 与 (C-Q)，所有 Schatten norms 都等于 (1)。full hidden-limit error 却是

\[
[C,0]\Pi-[Q,0]
=\frac12\begin{bmatrix}-1&1\\1&1\end{bmatrix},
\]

其两个 singular values 都是 (1/\sqrt2)。于是

\[
\left\|[C,0]\Pi-[Q,0]\right\|_{S_s}
=2^{1/s-1/2}>1,
\qquad 1\le s<2.
\tag{21}
\]

所以“对任意输出维数、任意 (C,Q)、任意 admissible completion 的 universal reduction”在 Schatten family 中 iff (s\ge2)。这不是说每个特例在 (s<2) 都失败；例如单输出时矩阵 rank 至多 (1)，所有 Schatten norms 相同。阈值也只针对**完整输入 operator**；若 adversarial input 被限制在 ports，loss 只看 (CT_{\Gamma,\Gamma}-Q)，则 terminal convex-hull theorem 对所有 convex norms 都适用。

该阈值的证明使用的是经典 Schatten norm convexity，纯矩阵不等式本身不新。潜在原创点是它与 hidden-completion/partial-scenario equality 的 sharp 结合。本次定向检索没有找到同型的 random-forest/Tikhonov norm-frontier theorem，但“未找到”不能升级成全球首创证明。

## 7. 自由 (Q) 的退化与真正有价值的 constrained design

因为 (0,I\in\mathcal P_\Lambda^\partial)，对任意 (Q)，

\[
\max_P\|CP-Q\|_2
\ge\max\{\|Q\|_2,\|C-Q\|_2\}
\ge\|C\|_2/2.
\tag{22}
\]

另一方面，取 (Q=C/2)。任意 orthogonal projection (P) 满足
(2P-I) orthogonal，所以

\[
\|CP-C/2\|_2
=\frac12\|C(2P-I)\|_2
=\frac12\|C\|_2.
\tag{23}
\]

这证明 (4)，而且每个 projector 都达到相同距离。相同论证对所有 unitarily invariant norms 也成立。

因此，论文不能把

\[
\min_{Q\in\mathbb R^{m\times q}}
\max_{P\in\mathcal P_\Lambda^\partial}\|CP-Q\|_2
\]

作为主要算法贡献。自然且已有应用依据的约束包括：

- (r)-hop/sparsity mask，使每个 output 只读取实际可达 ports；
- graph-filter/polynomial realization，filter order 直接对应 communication rounds；
- constant/consensus preservation，例如 uniform coordinates 下 (Q\mathbf1=C\mathbf1)；
- shared opcode，如 (Q=\alpha I) 或多个节点共享少量系数；
- quantized coefficients、low rank、有限字长或固定消息维数。

distributed graph-filter 文献已经把 locality、filter order、communication rounds、quantization 与 topology randomness 当作自然设计约束，因此无需创造新术语。[^10][^11][^12][^13]

对给定 convex set (\mathcal Q_{\rm loc})，spectral design 是精确 SDP：

\[
\begin{aligned}
\min_{Q\in\mathcal Q_{\rm loc},\,r}\quad&r\\
\text{s.t.}\quad&
\begin{bmatrix}
rI_m&CP-Q\\
(CP-Q)^T&rI_q
\end{bmatrix}\succeq0,
\qquad P\in\mathcal P_\Lambda^\partial.
\end{aligned}
\tag{24}
\]

它的规模依赖 (q) 与 (B_{q+1})，不依赖 hidden node/edge count。这才是 exact completion theorem 的可测算法收益。

### 7.1 最小 shared-opcode 例：connected fused scenario 改变最优解

取 (q=m=2)、

\[
C=\operatorname{diag}(1,2),
\qquad Q=\alpha I.
\]

两个 ports 的五个 partial scenarios 是

\[
0,quad
\operatorname{diag}(1,0),quad
\operatorname{diag}(0,1),quad
I,quad
P_f=\frac12\mathbf1\mathbf1^T.
\]

若漏掉 fused connected scenario (P_f)，objective 是

\[
\max\{|\alpha|,|1-\alpha|,|2-\alpha|\},
\]

故 (\alpha=1)、(R=1)。加入 (P_f) 后，

\[
CP_f-\alpha I
=\begin{bmatrix}
1/2-\alpha&1/2\\
1&1-\alpha
\end{bmatrix}.
\]

全局 optimum 由其 spectral norm 与 (2-\alpha) 相等得到；行列式方程化为

\[
15\alpha^2-40\alpha+24=0.
\]

取相关的较小根：

\[
\boxed{
\alpha^*=\frac43-\frac{2\sqrt{10}}{15}
=0.91169631198,}
\]

\[
\boxed{
R^*=2-\alpha^*
=\frac{10+2\sqrt{10}}{15}
=1.08830368802.}
\tag{25}
\]

因此忽略 connected fused completion 会把 robust radius 低报 (8.83\%\)，平方风险低报 (18.44\%\)。应用解释很朴素：两个通道具有不同 output scaling，但硬件/协议要求同一个 scalar opcode；未知网络可能把两端强耦合，所以 shared gain 必须为该 completion 留余量。

### 7.2 exact network geometry 相比 generic contraction 可显著降低保守性

取 uniform base、(q=2,m=1)，

\[
C=[0,2],\qquad Q=[a,2-a],
\]

其中 (Q\mathbf1=C\mathbf1=2) 保持常数信号。五个 network scenarios 给

\[
R_{\rm net}(a)
=\max\left\{
\sqrt{a^2+(2-a)^2},
\sqrt2|a|,
\sqrt2|1-a|
\right\}.
\]

故 (a=1)，(R_{\rm net}^*=\sqrt2\)。若不用 graph completion theorem，而把 uncertainty 粗放大到所有
(0\preceq X\preceq I)，其 extreme points 是任意 orthogonal projections；(CX) 的 rank-one 边界在圆心 ([0,1])、半径 (1) 的圆上，于是

\[
R_{\rm full}(a)
=1+\sqrt{a^2+(1-a)^2},
\qquad
R_{\rm full}^*=1+1/\sqrt2.
\]

exact topology geometry 把最优 radius 从 (1.7071) 降到 (1.4142)，相对 generic bound 降低 (17.16\%\)。在 network-optimal (Q=[1,1]) 处，generic set 给 (2)，而 exact value 是 (\sqrt2)，保守高估 (41.42\%\)。这是可直接放入论文的解析收益指标。

## 8. 任意 SPD $\Lambda$ 的数学意义与 locality 代价

两层 DPP 证明说明 diagonal port grounding 不是最大代数范围；任意
(\Lambda\succ0) 都成立。这个 extension 数学上正确，但应用解释必须克制：

- dense (\Lambda) 表示 ports 之间已有 correlated fidelity/base coupling，不再是每个 node 独立的 unary grounding；
- normalized coordinates 要使用 (\Lambda^{1/2})，它通常会混合多个 port inputs；若 ports 不在同一设备或小区域，这一步本身可能非局部；
- 因而 dense-$\Lambda$ 版本适合写成 algebraic extension，主应用最好仍用 diagonal $D_\Gamma$，或者明确把 $q$ 个 ports 视为一个可共同访问的小 interface。

SPD extension 的 proof engine 仍是已知 DPP identity，所以它单独的原创性不强；价值在于给 theorem 找到自然的最大 whitening-invariant范围，并揭示场景是 subspace projectors 而不是逐 block averages。

## 9. 与最接近 primary sources 的逐对象对照

| 方向与原始文献 | 已有对象/公式 | 与候选 theorem 的关键差别 |
|---|---|---|
| Pilavcı et al. 2021[^1] | (K=(L+D)^{-1}D=\mathbb E S_F)，(S_F) 在每个 sampled tree 内按 (d_j) weighted average；研究 Monte Carlo unbiasedness、variance、SURE/LOOCV 与运行时间。 | **覆盖最直接的第一层 forest identity。** 固定一个已知 graph，不取 arbitrary hidden completions；没有 port-only common (Q)、partial-projector hull、worst-completion equality 或 path converse。 |
| Dereziński--Khanna--Mahoney 2020[^2]；Kassel--Lévy 2023[^3] | selected-column/DPP mean projection，例如 (\mathbb EP_S=A(I+A^TA)^{-1}A^T)。 | 两次 DPP 都是已知工具。候选的新意只能是把它们组合到 graph-completion terminal image 与 robust loss，不能声称 mean-projection identity 新。 |
| Dörfler--Bullo 2013[^4] | Kron/Schur reduction 保持 loopy-Laplacian 结构，并研究 boundary clique、effective resistance 与 sensitivity。 | response/Kron object 是消去 hidden variables 后的 Schur complement；候选研究其 grounded inverse 的 normalized terminal block及其跨所有 hidden dimensions 的 convex envelope。没有相同 minimax equality。 |
| Kenyon--Wilson 2006/2011[^5] | planar groves 的 boundary connection probabilities 可由 DtN/response matrix 表示。其“partition projection matrix” (P^{(t)}) 是从 all set partitions 到 planar partitions 的一个 rectangular integer-coefficient linear map。 | **术语极易误撞。** 他们的 (P^{(t)}) 不是 Euclidean orthogonal projector (2)。对象是已知 response \(\to\) grove probabilities；候选是未知 completion \(\to\) grounded resolvent 的 worst convex envelope。planar case只有 noncrossing partitions；候选 arbitrary graph 有全部 partial partitions。 |
| Lam 2014/2018[^6] | 用 cactus networks compactify circular planar electrical networks；把一组 boundary vertices 视为由 infinite conductance identify，并以 noncrossing grove projective coordinates/Grassmannian strata 描述 closure。 | **最危险、必须正文引用的数学近邻。** infinite-conductance boundary identification 与 active-block limit 是同一物理极限。但 Lam 的对象是 planar DtN equivalence/compactification，不是 (H^{1/2}(H+L)^{-1}H^{1/2}) 的 terminal convex hull；没有 hidden dilution/inactive blocks、operator minimax、finite LMI、degree-2 universal witnesses 或 Schatten threshold。 |
| De Rosa--Khajavirad 2022[^7]；Cheng--Zhang--Schuurmans 2013[^8] | full partition 的 normalized equivalence/cluster matrix (Z=\sum_B\mathbf1_B\mathbf1_B^T/|B|) 以及其 convex relaxations/polytope。 | uniform full-partition projector 顶点早已存在。候选的潜在新增是 partial/inactive terminal image、weighted/SPD geometry 与 grounded-network realization，不是“发明 partition projector”。 |
| Moura--Yaman--Leus 2026[^9] | fixed host graph 上 connected (k)-subpartition 的 incidence-vector polytope，研究 facets 与 separation。 | 表示变量是 labeled incidence vectors，host graph 固定，connectivity 也相对该 host；不是 normalized Gram/projector，不是任意 hidden resolvent image。 |
| Foucart--Liao--Veldt 2023[^14] | known graph 上，以 graph smoothness 与 observation error sets 为假设的 local/global worst-case optimal recovery及 regularization parameter selection。 | “worst-case graph signal recovery”主题相近，但 uncertainty 在 signal/noise，不在任意 hidden topology；没有共同 local opcode (Q) 或 partition completion theorem。 |
| Kishida 2026[^15] | uncertain Tikhonov regularization 的 LFT/\(\mu\)-analysis 与 LMI bound；对 fixed vertex/edge set 的 independent bounded edge-weight uncertainty，robust objective 等价于使用 (W_0+R^2) 的 augmented graph。 | **截至 2026 最直接的 robust-Tikhonov 近邻。** 它允许给定 edges 的 bounded parametric uncertainty，不允许任意 hidden nodes、任意新 topology 或 unbounded conductance；其 closed form 也不是 partial-partition scenario reduction。 |
| Segarra et al.、Coutino et al.、Romero et al.、Ben Saad et al.[^10][^11][^12][^13] | known topology 上以 polynomial/node/edge-variant graph filters 实现或逼近 operators；filter order、轮数、量化、随机 link loss的 tradeoff。 | 给出了 (\mathcal Q_{\rm loc}) 的自然工程形式，但未发现对 arbitrary hidden completion 的 exact worst-case set。候选可把这些 realization constraints 接到 (24)，而非替代它们。 |

### 9.1 DtN/groves 与 grounded resolvent 必须逐式分开

Kenyon--Wilson/Lam 的 response matrix 是 Laplacian block的 Schur complement：

\[
\mathcal R_G
=L_{\Gamma\Gamma}
-L_{\Gamma Z}L_{ZZ}^{-1}L_{Z\Gamma},
\tag{26}
\]

它把 boundary voltages 映到 boundary currents，通常 positive semidefinite 且
(\mathcal R_G\mathbf1=0)。候选 terminal block则是 grounded Schur complement 的**逆**：

\[
[T_G]_{\Gamma\Gamma}
=\Lambda^{1/2}
\left(
\Lambda+L_{\Gamma\Gamma}
-L_{\Gamma Z}(D_Z+L_{ZZ})^{-1}L_{Z\Gamma}
\right)^{-1}
\Lambda^{1/2}.
\tag{27}
\]

两者共享 elimination 与 spanning-forest combinatorics，却不是同一个 boundary map。Lam 的 cactus limit说明“infinite conductance identifies boundary vertices”已有明确先例，因此 path converse 的 active half 不宜单独包装成重大新发现；真正不同的是 grounding、inactive dilution、convex mixture和 robust norm consequence。

### 9.2 查重结论与置信度

本次检索覆盖 graph Tikhonov/random spanning forests、DPP mean projections、Kron/DtN、groves/cactus compactification、normalized partition polytopes、connected subpartition polytopes、optimal graph-signal recovery、uncertain graph Tikhonov、distributed/local graph filtering，以及 Schatten/random-projection 关键词组合。对上述 primary papers 的公开全文检查，没有看到以下四项同时出现：

1. hidden cardinality 可变且无上界；
2. normalized grounded inverse 的 terminal partial-partition closed hull；
3. full hidden-input local approximation的 exact spectral/Schatten-(s\ge2) equality；
4. 每个 scenario 的 connected degree-2 realization 与 constrained-(Q) finite SDP。

需要把 **unbounded-completion limit** 与 **finite hidden budget** 明确分开。本报告验证的是前者：partial-partition projectors 是允许 hidden cardinality、总 grounding 与强边 conductance 随序列增长时的 closed limiting scenarios。若统一限制 hidden vertex 数、总 grounding、边权或禁止直接 port--port 边，当前证明只给出一个外层上界，并没有刻画可达集合的精确 convex hull。此次 primary-source 检索同样没有找到这一 finite-budget hull 的现成刻画；但这只能表述为“未发现直接结果”，不能当作完备的不存在证明。有限预算版本因此是一个独立且可能更难的后续问题，而不是主定理已经覆盖的 corollary。

因此给出分层主观判断：

- forest/DPP expectation identity：原创性极低，约 (0.05)；
- partition projector/cluster matrix 对象本身：低，约 (0.10)；
- arbitrary-hidden terminal closed-hull theorem：中等，约 (0.55)；
- fixed-(Q) full-input spectral equality + degree-2 sharp converse：中等偏上，约 (0.65)；
- Schatten “iff (s\ge2)” 在该 completion 问题中的 sharp frontier：中等，约 (0.65)，但基础 norm convexity 是经典的；
- arbitrary SPD (\Lambda) extension 单独看：中低，约 (0.35)；
- 整体 theorem + constrained design + empirical package：中等，约 (0.60)。

这些数字是“在已查范围内没有直接 collision”的信心，不是概率学证明。最安全的论文用语是 “we are not aware of an exact characterization under arbitrary hidden-node completion”，并紧接着列出 Pilavcı、Kenyon--Wilson、Lam、Kishida 等最近先例。

## 10. 可发表的收益指标

这个方向真正能量化的收益不是“比全局通信天然省电”这一句泛论，而是以下四类可复现指标：

1. **certificate exactness**：随机或 adversarial completions 的实际
   (\|\Delta_G\|) 与 (B_{q+1})-scenario radius 的 ratio；degree-2 witness 应从下方逼近 (1)。
2. **reduced conservatism**：与 generic positive-contraction、0/I endpoint-only、有限 random-topology sampling 三种 baseline 比较 robust radius。第 7.2 节已有 exact (17.16\%\) improvement。
3. **constrained-design gain**：在同一 sparsity/hop/shared/quantization budget 下，比较 exact scenario SDP 与 nominal/endpoint/sample-based (Q) 的 realized worst-case error；第 7.1 节给出漏场景时 (8.83\%\) radius、(18.44\%\) squared-risk 低报。
4. **dimension removal**：优化的 LMI count 与 size 只依赖 (q,m)，不依赖 hidden (n)。报告 scenario generation time、solver time、peak RAM、active LMI 数以及相对直接采样/求逆的加速。

Bell growth也必须如实呈现：

| (q) | (B_{q+1}) |
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

所以 exact SDP 是 boundary-size fixed-parameter method，不是对大 (q) 的 polynomial-time miracle。(q\ge9) 应研究 active-scenario generation、symmetry reduction 或 problem-specific separation；当前 theorem 没有自动给出高效 separation oracle。

## 11. 可直接运行的实验设计

### 11.1 Proof-regression 单元测试

- 枚举 (q=1,2,3) 的全部 partial partitions，检查 projector symmetry/idempotence；dense random (\Lambda\) 时用 QR/SVD 构造 (2)，禁止用错误的 rank-one sum。
- 对随机 (C,Q) 验证 q=1 Schatten counterexample，扫描
  (s\in\{1,1.25,1.5,1.75,2,3,4,\infty\})，比较有限 conductance 与极限
  (2^{1/s-1/2})。
- 锁定第 7.1 节的 exact constants，以及第 7.2 节的 (\sqrt2) 与
  (1+1/\sqrt2)。

### 11.2 Random-completion upper-bound stress test

对 (q=2,\ldots,7)，生成以下 hidden topologies：Erdős--Rényi、random geometric、Barabási--Albert、grid、cycle-rich graphs、随机 trees 与长 paths；hidden size 取
(10,10^2,10^3,10^4)，edge weights/groundings 取多个数量级的 log-uniform。每个实例用 sparse Cholesky 或 multi-RHS CG 计算
(CE_\Gamma T_G)，记录

\[
\rho_G=
\frac{\|\Delta_G(C,Q)\|_2}
{\max_P\|CP-Q\|_2}.
\]

正确实现应始终 (\rho_G\le1)（容许 solver tolerance）。任何稳定超界都优先视为 proof/model mismatch，而不是“interesting data”。

### 11.3 Degree-2 sharpness curves

对每个 predicted active scenario，构造第 5 节 path witness，分别扫描

\[
h\in\{1,2,4,8,\ldots,4096\},\qquad
\kappa\in\{10,10^2,\ldots,10^{10}\},
\]

并在 connected 版本扫描 weak link (\delta)。画出 certificate gap
(1-\rho_G) 对 (h,\kappa,\delta) 的 log--log 图，验证 hidden-column norm 的
(O(h^{-1/2})) 趋势。uniform 与 diagonal heterogeneous base 为主，dense SPD 放 supplement。

### 11.4 Constrained (Q) benchmark

依次使用：

- shared scalar (Q=\alpha I)；
- prescribed sparsity mask；
- (Q\mathbf1=C\mathbf1)；
- (r)-hop graph-filter coefficients；
- 8/12/16-bit quantized coefficient grid。

比较 exact scenario SDP、0/I-only、generic (0\preceq X\preceq I)、(10^k) random completions 与 nominal-known-topology design。输出 robust radius、sampled worst error、消息数、rounds、coefficient storage、FLOPs 与 solver wall time。energy proxy 应明确为 transmitted bytes/packets 与 radio model，不要只以“local”二字代替测量。

### 11.5 本地 CPU/GPU 与现有数据

- (q\le7)：CPU/CVXPY 或 MOSEK/Clarabel 做 exact SDP；scenario SVD/LMI generation 可多核并行。
- (q=8,9)：批量 GPU 计算 (\|CP-Q\|_2) 与找 active scenarios，优化仍由 CPU conic solver完成；记录 host--device transfer，避免虚假 GPU speedup。
- sparse completion stress test 可用已归档 PGLib cases 与 SuiteSparse usroads。PGLib 是 topology/physics base，不是现成 SE benchmark；应采用项目已冻结的 measurement/noise generation protocol。usroads 适合大直径、长路径和 boundary dilution压力测试。
- synthetic paths 必须保留为主验证集，因为真实有限网络不太可能自动命中 sharp singular limits。

## 12. 敌对审稿清单与推荐论文定位

提交前至少逐项回答：

1. 是否明确引用 Pilavcı 的 (\mathbb ES=(L+D)^{-1}D)？
2. 是否说明两次 DPP identity 都是已知工具？
3. 是否引用 Kenyon--Wilson/Lam，并区分 DtN/grove coordinates 与 grounded inverse？
4. 是否避免把 Kenyon--Wilson 的“partition projection matrix”误认为 orthogonal projector？
5. 是否把 equality 写成 supremum，并公开 unbounded hidden size/grounding/conductance assumptions？
6. 是否说明 fixed visible subgraph、forbidden port edges、bounded hidden budget 时 converse 会失效？
7. 是否把自由 (Q=C/2) 退化写成 proposition，而不是藏在 supplement？
8. 是否给出真实 (\mathcal Q_{\rm loc})，并把 rounds/messages/energy 纳入实验？
9. 是否限定 Schatten theorem 为完整输入、(1\le s\le\infty) 的 universal statement？
10. 是否把 dense (\Lambda) 的 whitening locality cost 说清？

推荐 headline：

> **Exact robust local approximation of grounded graph smoothers under arbitrary hidden-node completion.**

推荐主贡献顺序：

1. arbitrary-hidden terminal closed hull；
2. full-input spectral/Schatten-(s\ge2) exact robust certificate；
3. connected degree-2 sharp converse；
4. (s=2) sharp norm frontier及 (s<2) 最小反例；
5. constrained local/shared design 的 exact finite SDP 与真实收益。

不推荐的宣传：

- “终结一般 network state estimation”；
- “首次发现 Tikhonov inverse 的 forest representation”；
- “最 relaxed 的所有 local algorithms 充要条件”；
- “任意 convex full-input loss 都可 finite reduction”；
- “free-(Q) minimax 产生新最优 decoder”。

本 theorem 精确解决的是一个范围清楚但相当广的 robust graph-smoother completion 问题。若正文把 constrained design、通信实现和实证收益做扎实，它有明显论文价值；若只剩两次已知 DPP identity 加一个自由 (Q) minimax，则 reviewer 很容易判为漂亮但增量有限。

## 13. 本次新增本地证据

为核查 electrical-network 最近先例，已把以下合法 arXiv OA 版本归档到
`literature/09_network_completion/`：

- `2006_kenyon_wilson_boundary_partitions.pdf`；
- `2018_lam_electroid_varieties.pdf`。

相应 `SOURCE_URLS.md`、`SHA256SUMS`、`research/corpus/09_network_completion/*.txt` 与
`research/corpus/manifest.tsv` 已更新，并已逐项重新核对 SHA-256。

## Sources

[^1]: Yusuf Yiğit Pilavcı, Pierre-Olivier Amblard, Simon Barthelmé, and Nicolas Tremblay, “[Graph Tikhonov Regularization and Interpolation via Random Spanning Forests](https://arxiv.org/abs/2011.10450),” *IEEE Transactions on Signal and Information Processing over Networks* 7 (2021), DOI [10.1109/TSIPN.2021.3084879](https://doi.org/10.1109/TSIPN.2021.3084879), especially Proposition 2 and equations (19)--(21).

[^2]: Michał Dereziński, Rajiv Khanna, and Michael W. Mahoney, “[Improved Guarantees and a Multiple-Descent Curve for Column Subset Selection and the Nyström Method](https://arxiv.org/abs/2002.09073),” NeurIPS 2020, especially Lemma 5.

[^3]: Adrien Kassel and Thierry Lévy, “[On the Mean Projection Theorem for Determinantal Point Processes](https://arxiv.org/abs/2203.04628),” *ALEA* 20 (2023), 497--504, DOI [10.30757/ALEA.v20-17](https://doi.org/10.30757/ALEA.v20-17).

[^4]: Florian Dörfler and Francesco Bullo, “[Kron Reduction of Graphs with Applications to Electrical Networks](https://arxiv.org/abs/1102.2950),” *IEEE Transactions on Circuits and Systems I* 60(1) (2013), 150--163, DOI [10.1109/TCSI.2012.2215780](https://doi.org/10.1109/TCSI.2012.2215780).

[^5]: Richard W. Kenyon and David B. Wilson, “[Boundary Partitions in Trees and Dimers](https://arxiv.org/abs/math/0608422),” *Transactions of the American Mathematical Society* 363(3) (2011), 1325--1364, DOI [10.1090/S0002-9947-2010-04964-5](https://doi.org/10.1090/S0002-9947-2010-04964-5), especially §1.2 and Appendix A.

[^6]: Thomas Lam, “[Electroid Varieties and a Compactification of the Space of Electrical Networks](https://arxiv.org/abs/1402.6261),” arXiv:1402.6261v4 (2018), especially §§1.1--1.4 and 4.1--4.3.

[^7]: Antonio De Rosa and Aida Khajavirad, “[The Ratio-Cut Polytope and K-Means Clustering](https://arxiv.org/abs/2006.15225),” *SIAM Journal on Optimization* 32(1) (2022), 173--203, DOI [10.1137/20M1348601](https://doi.org/10.1137/20M1348601).

[^8]: Hao Cheng, Xinhua Zhang, and Dale Schuurmans, “[Convex Relaxations of Bregman Divergence Clustering](https://arxiv.org/abs/1309.6823),” UAI 2013; see the normalized equivalence matrix (M=Y(Y^TY)^{-1}Y^T).

[^9]: Phablo F. S. Moura, Hande Yaman, and Roel Leus, “[On the Connected (Sub)partition Polytope](https://arxiv.org/abs/2401.01716),” *Mathematical Programming* (2026), DOI [10.1007/s10107-025-02321-1](https://doi.org/10.1007/s10107-025-02321-1).

[^10]: Santiago Segarra, Antonio G. Marques, and Alejandro Ribeiro, “[Distributed Linear Network Operators Using Graph Filters](https://arxiv.org/abs/1510.03947),” arXiv:1510.03947 (2015).

[^11]: Mario Coutino, Elvin Isufi, and Geert Leus, “[Advances in Distributed Graph Filtering](https://arxiv.org/abs/1808.03004),” *IEEE Transactions on Signal Processing* 67(9) (2019), 2320--2333.

[^12]: Daniel Romero, Siavash Mollaebrahim, Baltasar Beferull-Lozano, and César Asensio-Marco, “[Fast Graph Filters for Decentralized Subspace Projection](https://arxiv.org/abs/2011.07579),” arXiv:2011.07579 (2020).

[^13]: Leila Ben Saad, Baltasar Beferull-Lozano, and Elvin Isufi, “[Quantization Analysis and Robust Design for Distributed Graph Filters](https://arxiv.org/abs/2004.06692),” arXiv:2004.06692 (2020).

[^14]: Simon Foucart, Chunyang Liao, and Nate Veldt, “[On the Optimal Recovery of Graph Signals](https://arxiv.org/abs/2304.00474),” 2023 International Conference on Sampling Theory and Applications (SampTA), 2023.

[^15]: Masako Kishida, “[Tikhonov Regularization With Uncertain Graph Laplacians: Application to Graph Signal Processing](https://doi.org/10.1109/LCSYS.2026.3702187),” *IEEE Control Systems Letters* 10 (2026), 985--990, DOI 10.1109/LCSYS.2026.3702187.
