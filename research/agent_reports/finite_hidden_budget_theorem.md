# 有限隐藏节点预算：精确场景定理、完整输入误差与定量层次

## 0. 结论

固定 (q) 个可见端口，每个节点有 unit grounding，completion 至多加入
(h) 个 unit-grounded hidden vertices，并允许任意无向拓扑和非负边权。令

\[
T_G=(I+L_G)^{-1},\qquad X_G=T_G[\Gamma,\Gamma].
\]

本专项得到以下结论。

1. 有限预算的端口闭凸包有一个**精确有限描述**。对端口 full
   partition \(\pi\) 和整数分配 \(k_B\ge0\)、
   \(\sum_{B\in\pi}k_B\le h\)，令

   \[
   X_{\pi,k}
   =\sum_{B\in\pi}\frac{\mathbf1_B\mathbf1_B^T}{|B|+k_B}.
   \]

   则

   \[
   \overline{\operatorname{conv}}\{X_G:|Z|\le h\}
   =\operatorname{conv}\{X_{\pi,k}\}.
   \]

2. “至多 (h)”和“恰好 (h)”给出同一个闭凸包。若要求每个 graph
   connected，闭凸包仍不变；每个场景可由 connected、最大度数 (2) 的
   path completion 逼近。无限大和无限小边权只出现在极限中，所以 closure
   一般不能删除。

3. 场景参数化没有重复，但并非每个场景都是端口 hull 的 vertex。对
   (h\ge1)，vertex 恰好是

   \[
   k=0
   \quad\text{或}\quad
   \sum_Bk_B=h.
   \]

   两类 vertex 全部 exposed。任何
   (0<\sum_Bk_B<h) 的场景都不是 vertex。

4. 对保留全部 hidden-input columns 的 direct smoother，finite-(h)
   精确归约对**所有 Schatten-(p) 范数 (1\le p\le\infty)**成立，特别包括
   nuclear norm。若 (X=X_{\pi,k})，正确场景误差不是
   \(\|CX-Q\|\)，而是

   \[
   \left\|[CX-Q,\ C(X-X^2)^{1/2}]\right\|_{S_p}.
   \]

5. 令有限预算 hull 为 (K_h)，无限隐藏预算的 partial-partition hull 为
   (K_\infty)。operator-norm Hausdorff distance **精确等于**

   \[
   d_H^{\mathrm{op}}(K_h,K_\infty)=\frac q{q+h}.
   \]

6. 有限预算确实能收紧应用 bound。最小非平凡解析例
   (q=2,h=1,C=P_{\rm av}=\mathbf1\mathbf1^T/2) 的自由局部映射最优解是

   \[
   Q_* = \frac23P_{\rm av},\qquad
   R_1^*=\frac{\sqrt2}{3}=0.471404\ldots,
   \]

   而无限预算值为 (R_\infty^*=1/2)。半径下降约 (5.72\%\)，平方风险下降
   (11.11\%\)。(h=0) 时更有 (Q=P_{\rm av}) 和零误差。

数学判决是 **GO**。基础 forest/DPP identity 已知；本次未在检索到的原始文献
中发现“有限 hidden-node budget 的精确 resolvent port hull + path converse +
完整输入 robust approximation”的同一结论。这个阴性检索不是全球原创性证明，
论文仍应把 novelty 写成该组合，而不是 matrix-forest identity 本身。

## 1. 模型与场景

令

\[
\Gamma=\{1,\ldots,q\},\qquad Z\cap\Gamma=\varnothing,qquad |Z|=m\le h.
\]

(G) 是 (V=\Gamma\mathbin{\dot\cup}Z) 上任意有限无向加权图；边权属于
([0,\infty))，允许 cycles 和 disconnected components。令 (L_G) 为 weighted
Laplacian，并定义

\[
T_G=(I_{q+m}+L_G)^{-1},\qquad
E_\Gamma=[I_q\ 0],\qquad
X_G=E_\Gamma T_GE_\Gamma^T.
\tag{1}
\]

对 (Gamma) 的一个 full set partition

\[
\pi=\{B_1,\ldots,B_s\}
\]

以及 weak allocation

\[
k=(k_{B_1},\ldots,k_{B_s})\in\mathbb Z_+^s,
\qquad \sum_{B\in\pi}k_B\le h,
\tag{2}
\]

定义

\[
X_{\pi,k}
:=\sum_{B\in\pi}\frac{\mathbf1_B\mathbf1_B^T}{|B|+k_B}
=\sum_{B\in\pi}\frac{|B|}{|B|+k_B}P_B,
\qquad
P_B=\frac{\mathbf1_B\mathbf1_B^T}{|B|}.
\tag{3}
\]

记所有 (3) 组成的有限集合为 \(\mathcal X_{q,h}\)，并记

\[
K_h:=\operatorname{conv}\mathcal X_{q,h}.
\tag{4}
\]

## 2. 精确有限预算闭凸包

### 定理 2.1

在上述模型下，

\[
\boxed{
\overline{\operatorname{conv}}
\{X_G:G\text{ has }|Z|\le h\}
=K_h.}
\tag{5}
\]

以下三个模型变体给出同一个右端：

- hidden vertex 数量至多 (h)；
- hidden vertex 数量恰好 (h)，但允许 unused hidden-only components；
- 数量恰好 (h)，且要求每个 finite graph connected，等式按 closure 理解。

每个 (X_{\pi,k}\) 都由一列 connected path graphs 逼近，所有 witnesses 的最大
degree 均不超过 (2)。

### 证明

给 (G) 任意定向并写

\[
L_G=AA^T,
\]

其中 (A) 的 edge column 是
(\sqrt{w_e}(e_u-e_v))。edge (L)-ensemble DPP 的
mean-projection identity 给出

\[
(I+AA^T)^{-1}=\mathbb E P_{\ker A_S^T}.
\tag{6}
\]

只有 linearly independent edge sets (S) 有正概率；对 graph incidence matrix，
这些集合恰是 forests。若 forest (S) 的 components 为 (C)，则

\[
P_{\ker A_S^T}
=\sum_C\frac{\mathbf1_C\mathbf1_C^T}{|C|}.
\tag{7}
\]

每个含端口的 component (C) 诱导非空 block
(B=C\cap\Gamma)；这些 blocks 构成 (Gamma) 的 full partition。令

\[
k_B=|C\cap Z|.
\]

则 (7) 的端口主子块正是 (3)，并且

\[
\sum_Bk_B\le |Z|\le h.
\]

hidden-only components 只消耗未使用预算，在端口主子块中没有贡献。对 (6)
取端口主子块，便得 (X_G\in K_h)，从而证明 (5) 的上界。

反向固定 ((\pi,k))。对每个 block (B)，把 (B) 的 ports 与恰好 (k_B)
个 hidden vertices 排成一条 path；剩余 (h-\sum_Bk_B) 个 hidden vertices
另排成 hidden-only paths。先把每个 path component 内部的 edge conductance
(t\to\infty)。由谱分解，

\[
(I+tL_F)^{-1}\longrightarrow P_{\ker L_F},
\tag{8}
\]

其端口块就是 (X_{\pi,k})。若要求 graph connected，把这些 path segments
依次用 conductance \(\varepsilon\downarrow0\) 的弱边连接；整个 graph 仍是一条
path。沿任意对角序列 \(t_n\to\infty,\varepsilon_n\to0\)，continuity 给出同一
极限。这同时证明 at-most/exactly-(h) 和 connected 版本。

一般需要 closure：finite edge weights 下 (T_G\succ0)，而合并两个以上 ports
所得的场景常为 rank deficient，只有 infinite-conductance limit 才达到。

## 3. 场景数量与精确 vertex pruning

### 3.1 没有重复

固定含 (s) 个 blocks 的 partition，满足 (2) 的 weak allocations 数是

\[
{h+s\choose s}.
\]

不同参数产生不同矩阵：(i\ne j) 的 entry 为正当且仅当 (i,j) 属于同一
block，因此矩阵先唯一恢复 partition；随后 block value
(1/(|B|+k_B)) 唯一恢复 (k_B)。所以总场景数**恰好**是

\[
\boxed{
N(q,h)=\sum_{s=1}^q S(q,s){h+s\choose s},}
\tag{9}
\]

其中 (S(q,s)) 是第二类 Stirling number。

### 3.2 哪些场景是 terminal hull 的 vertices

### 定理 3.1

若 (h\ge1)，则 (K_h) 的 vertex 集恰为

\[
\boxed{
\operatorname{vert}(K_h)
=\{X_{\pi,0}\}_\pi
\mathbin{\cup}
\{X_{\pi,k}:\sum_Bk_B=h\}.}
\tag{10}
\]

每个 vertex 都是 exposed。因而 terminal convex loss、spectral squared risk 和
所有可以证明为 (X) 的 convex function 的风险只需检查

\[
\boxed{
N_{\rm vert}(q,h)
=\sum_{s=1}^q S(q,s)
\left[1+{h+s-1\choose s-1}\right]}
\tag{11}
\]

个场景。(h=0) 时 vertex 数为 Bell number (B_q)。

### 证明：未饱和分配不是 vertex

固定 (pi)，假设某个 (k_B>0) 且 (sum_Dk_D<h)。令

\[
K=h-\sum_{D\ne B}k_D.
\]

则 (0<k_B<K)。因为
(1/(|B|+k_B)) 严格位于
(1/|B|) 和 (1/(|B|+K)) 之间，(X_{\pi,k}) 是把该坐标改成
(0) 和 (K) 后两个可行场景的严格 convex combination。因此它不是 vertex。

### 证明：(k=0) exposed

(X_{\pi,0}=P_\pi) 是 orthogonal projection。所有场景均满足
(0\preceq X\preceq I)，而 linear functional

\[
X\longmapsto\langle2P_\pi-I,X\rangle
\]

在整个 positive-contraction interval 上唯一最大于 (P_\pi)。故它在
(K_h) 上 exposed。

### 证明：每个饱和分配 exposed，包括跨 partition 的比较

固定目标 ((\pi,k^*))，其中 (sum_Bk_B^*=h)。构造 symmetric matrix
(H_\pi)：若 (i\in B)，令

\[
(H_\pi)_{ii}=-(|B|-1);
\]

同一目标 block 内的非对角 entry 取 (1)，不同目标 blocks 间的 entry 取
(-1)。对任意非空 (S\subseteq\Gamma)，令
(r_B=|S\cap B|)，则

\[
\mathbf1_S^TH_\pi\mathbf1_S
=\sum_{B\in\pi}r_B(r_B-|B|)
-2\sum_{B<D}r_Br_D\le0.
\tag{12}
\]

等号成立当且仅当 (S) 恰是某一个目标 block。于是

\[
\langle H_\pi,X_{\pi,k}\rangle=0
\quad\text{对所有同 partition 的 }k,
\]

而对任何 (pi'\ne\pi) 及任何 (k')，

\[
\langle H_\pi,X_{\pi',k'}\rangle<0.
\tag{13}
\]

现在在每个目标 block (B) 上取 constant matrix (A_B)，使

\[
\mathbf1_B^TA_B\mathbf1_B
=-\gamma_B,
\qquad
\gamma_B=(|B|+k_B^*)^2,
\]

并令跨 block entries 为零。对同一 partition，目标函数为

\[
f(k)=\langle A,X_{\pi,k}\rangle
=-\sum_B\frac{\gamma_B}{|B|+k_B}.
\tag{14}
\]

给 block (B) 分配从 (r) 到 (r+1) 的一个 hidden vertex，使 (14) 增加

\[
\Delta_B(r)
=\frac{\gamma_B}{(|B|+r)(|B|+r+1)}.
\tag{15}
\]

目标分配中的每个已选 increment (r<k_B^*) 都满足
(Delta_B(r)>1)，每个未选 increment (r\ge k_B^*) 都满足
(Delta_B(r)<1)。所以在 (sum k_B\le h) 下，(k^*) 是 (14) 的唯一
maximizer。候选集有限，故取充分大的 (M)，linear functional

\[
X\longmapsto\langle MH_\pi+A,X\rangle
\]

先由 (13) 排除所有其他 partitions，再由 (14)--(15) 唯一选择 (k^*)。目标
场景 exposed，(10) 得证。

### 一个重要的 pruning 限制

(10) 是 **terminal convex hull** 的最小 vertex 集。对下一节的完整输入误差，
若 (2\le p<\infty)，

\[
X\longmapsto
\|\mathcal E_Q(X)\|_{S_p}^p
=\operatorname{tr}(M_Q(X)^{p/2})
\]

是 convex；(p=\infty) 时其平方是
(\lambda_{\max}(M_Q(X)))，同样 convex。因此 full-input Schatten
(p\ge2) 也只需检查 (11) 的 vertices。若 (1\le p<2)，场景风险不一定是
(X) 的 convex function；不能仅凭“(X) 不是 hull vertex”删除未饱和场景。
对所有 (1\le p\le\infty) 都安全的统一 certificate 应保留 (9) 的全部场景。

## 4. 完整 hidden-input direct smoother

固定 (C\in\mathbb R^{r\times q}) 和一个只能读取 port inputs 的局部映射
(Q\in\mathbb R^{r\times q})。在有 (m) 个 hidden vertices 的 completion 上，
定义完整误差算子

\[
\mathcal A_G(Q)
=CE_\Gamma T_G-QE_\Gamma
\in\mathbb R^{r\times(q+m)}.
\tag{16}
\]

不同 (m) 时可在右侧补零列到 (q+h) 维；Schatten norms 不变。

### 定理 4.1（所有 Schatten norms 的精确有限归约）

对任意固定 (Q) 和任意 (1\le p\le\infty)，

\[
\boxed{
\sup_{G:\,|Z|\le h}\|\mathcal A_G(Q)\|_{S_p}
=\max_{X\in\mathcal X_{q,h}}
\left\|[CX-Q,\ C(X-X^2)^{1/2}]\right\|_{S_p}.}
\tag{17}
\]

更一般地，(17) 对每个 unitarily invariant matrix norm 成立。

### 证明

由 (6)，(T_G=\mathbb E\widehat P)，其中 (widehat P) 是 full forest
orthogonal projection。令

\[
A_{\widehat P}=CE_\Gamma\widehat P-QE_\Gamma.
\]

任意 matrix norm 的 convexity 给出

\[
\|\mathcal A_G(Q)\|
=\|\mathbb EA_{\widehat P}\|
\le\mathbb E\|A_{\widehat P}\|.
\tag{18}
\]

写 (X=E_\Gamma\widehat PE_\Gamma^T)。利用
(widehat P^2=\widehat P)，有 exact Gram identity

\[
\begin{aligned}
A_{\widehat P}A_{\widehat P}^T
&=CXC^T-CXQ^T-QXC^T+QQ^T\\
&=(CX-Q)(CX-Q)^T+C(X-X^2)C^T\\
&=:M_Q(X).
\end{aligned}
\tag{19}
\]

因此 (A_{\widehat P}) 与 canonical matrix

\[
\mathcal E_Q(X)
:=[CX-Q,\ C(X-X^2)^{1/2}]
\tag{20}
\]

有相同 left Gram matrix，故有相同 singular values。forest compression
(X) 必在 (mathcal X_{q,h}) 中，(18) 给出 (17) 的上界。

反向对每个 (X_{\pi,k}) 使用定理 2.1 的 high-conductance path sequence。
完整 projection rows 收敛，其 Gram 收敛到 (19)，故每个右侧场景值都能从
admissible graphs 逼近。这给出下界。证明只使用 norm convexity，而没有使用
(M\mapsto\operatorname{tr}(M^{p/2})) 的 convexity，所以 nuclear norm 和所有
(1\le p<2) 均包含在内。

### 谱范数的 finite LMI

令 \(\rho=\tau^2\)。对一个固定场景 \(X\)，条件
\(\|\mathcal E_Q(X)\|_2^2\le\rho\) 等价于

\[
\boxed{
\begin{bmatrix}
\rho I_r-CXC^T+CXQ^T+QXC^T & Q\\
Q^T&I_q
\end{bmatrix}\succeq0.}
\tag{21}
\]

这是 (19) 的 Schur complement，且对 ((\rho,Q)) affine。因此，若
(Q\in\mathcal Q_{\rm loc}) 是 affine/convex local-communication constraint，
对所有 (11) 中的 terminal-hull vertices 施加 (21) 就是 exact finite SDP；
其余 (N(q,h)-N_{\rm vert}(q,h)) 个 LMI 是这些 vertex LMIs 的 convex
combinations，因而可以删除。这里不要与 nuclear/(p<2) 的 full-input
certificate 混淆：后者在没有额外证明时仍保留 (9) 的全部 candidates。

## 5. 自由 (Q) 不再普遍退化，以及最小收益例

对任何 (X)，取 (Q=C/2) 都有

\[
M_{C/2}(X)=\frac14CC^T.
\tag{22}
\]

所以 \(Q=C/2\) 总给出半径 \(\|C\|_2/2\)，但 finite \(h\) 时这只是一个
universal upper bound，不一定最优。无限预算 family 同时含 (X=0) 和 (X=I)，
两者强迫最优值恰为 \(\|C\|_2/2\)；有限预算没有 \(X=0\)，因此可能严格改善。

取

\[
q=2,\quad h=1,\quad
P=P_{\rm av}=\frac12\mathbf1\mathbf1^T,\quad C=P.
\]

五个场景是

\[
I,\quad
\operatorname{diag}(1/2,1),\quad
\operatorname{diag}(1,1/2),\quad
P,\quad
\frac23P.
\tag{23}
\]

取

\[
Q_*=\frac23P.
\]

它们的 spectral radii 依次是

\[
\frac13,\quad
\frac{\sqrt7}{6},\quad
\frac{\sqrt7}{6},\quad
\frac13,\quad
\frac{\sqrt2}{3}.
\tag{24}
\]

最后一个场景对任意 (Q) 都有不可消除的 hidden-input lower bound

\[
M_Q(2P/3)
\succeq C\left(\frac23P-\frac49P\right)C^T
=\frac29P.
\tag{25}
\]

故 (R\ge\sqrt2/3)，而 (24) 达到该值，证明

\[
Q_* =\frac23P,
\qquad R_1^*=\frac{\sqrt2}{3}.
\tag{26}
\]

与无限预算值 (1/2) 相比，半径相对下降

\[
1-\frac{2\sqrt2}{3}=5.7191\%\ldots,
\]

平方风险从 (1/4) 降到 (2/9)，下降 (1/9=11.11\%\)。若 (h=0)，
(PT_G=P) 对任意 port-only graph 都成立，所以 (Q=P) 给零误差。这形成一个
清楚的 budget--robustness hierarchy，而不是人为制造的约束收益。

## 6. (h\to\infty) 的精确 Hausdorff 速率

令 (K_\infty) 是全部 partial-partition projections 的 convex hull。

### 定理 6.1

对所有 (q\ge1,h\ge0)，

\[
\boxed{
d_H^{\rm op}(K_h,K_\infty)=\frac q{q+h}.}
\tag{27}
\]

### 证明：上界

每个 finite 场景写成

\[
X_{\pi,k}=\sum_B\alpha_BP_B,
\qquad
\alpha_B=\frac{|B|}{|B|+k_B}\in(0,1].
\]

独立地以概率 (alpha_B) retain block (B)，否则 drop 它，就把
(X_{\pi,k}) 写成 partial projections 的 convex combination。因此

\[
K_h\subseteq K_\infty.
\tag{28}
\]

反过来取任意 partial projection (P)。把所有 inactive ports
(I\subseteq\Gamma) 合成一个 block，对 active blocks 分配 (k=0)，对
(I) 分配全部 (h) 个 hidden vertices。所得 finite 场景是

\[
X=P+\frac{|I|}{|I|+h}P_I,
\]

所以

\[
\|X-P\|_2=\frac{|I|}{|I|+h}\le\frac q{q+h}.
\]

对 (K_\infty) 中的 convex combinations 逐 vertex 这样逼近，再用三角不等式，
得到 (27) 的上界。

### 证明：下界与所有 (q) 的 sharpness

(0\in K_\infty)。令 (u=\mathbf1/\sqrt q)。对任一 finite 场景，Cauchy--
Schwarz/Titu Andreescu inequality 给出

\[
\begin{aligned}
u^TX_{\pi,k}u
&=\frac1q\sum_{B\in\pi}\frac{|B|^2}{|B|+k_B}\\
&\ge\frac1q\frac{(\sum_B|B|)^2}{\sum_B(|B|+k_B)}\\
&\ge\frac q{q+h}.
\end{aligned}
\tag{29}
\]

该 inequality 对 convex combinations 保持，所以任意 (X\in K_h) 都满足
\(\|X\|_2\ge u^TXu\ge q/(q+h)\)。另一方面，一整个端口 block 配置 \(k=h\)
给

\[
X=\frac q{q+h}P_\Gamma,
\]

其 operator norm 恰为 (q/(q+h))。于是

\[
\operatorname{dist}_{\rm op}(0,K_h)=\frac q{q+h},
\]

与上界合并即得 (27)。

### fixed-(Q) 谱风险的定量推论

令

\[
g_Q(X)=\lambda_{\max}M_Q(X),
\qquad
L_Q=\|C\|_2^2+2\|C\|_2\|Q\|_2.
\]

因为 (M_Q(X)) 对 (X) affine，

\[
|g_Q(X)-g_Q(Y)|\le L_Q\|X-Y\|_2.
\tag{30}
\]

记 (R_h(Q)) 和 (R_\infty(Q)) 分别为 finite/unbounded completion 的 exact
spectral radius，则

\[
\boxed{
0\le R_\infty(Q)^2-R_h(Q)^2
\le
L_Q\frac q{q+h}.}
\tag{31}
\]

第一项来自 (K_h\subseteq K_\infty) 及 (g_Q) 的 convexity；第二项来自
(27)--(30)。这是一个明确的 (O(q/h)) certification hierarchy。它是通用
upper bound；具体实例可以像第 5 节那样得到更精确值。

## 7. 文献尽调与原创性边界

最接近的已知工具/方向如下。

| 原始文献 | 已解决内容 | 与本定理的差别 |
|---|---|---|
| Chebotarev--Shamis, matrix-forest theorem | ((I+L)^{-1}) 的 rooted-forest 展开 | 是证明起点；没有 finite hidden-node budget 的 port hull、vertex rule 或 minimax |
| Pilavci et al. 2021 | graph Tikhonov smoother 等于 random-forest component averaging 的期望 | 最直接 collision；没有 unknown hidden completion、budget-indexed hull 或 robust local (Q) |
| Dereziński--Khanna--Mahoney 2020；Kassel--Lévy 2023 | 一般 DPP mean-projection identity | 覆盖 (6) 的抽象母题，不覆盖 graph terminal image |
| Dörfler--Bullo 2013 | loopy Laplacian/Kron reduction 的闭合性与网络性质 | 研究 Schur reduction，不给 (5)、(17) 或 (27) |
| De Rosa--Khajavirad 2022 | normalized full-partition matrices/ratio-cut polytope | full (k=0) partition geometry已知；没有 hidden attenuation integers 与 electrical/path realization |
| Lam 2018；Kenyon--Wilson 2011 | planar electrical response、groves、boundary partitions 和 infinite-conductance compactification | planar response/DtN 对象；不是 grounded resolvent port hull，也没有 finite node budget |

本次还以 “hidden/interior vertex budget”、“bounded Steiner vertices”、
“Laplacian resolvent convex hull”、“electrical response with bounded interior nodes”
等组合检索了 primary literature。检索结果主要落在 matrix-forest formulas、planar
electrical response/Grassmannian、Kron reduction 与几何 Steiner network optimization；
没有定位到直接陈述 (5)、(10)、(17) 或 (27) 的论文。应使用如下保守表述：

> 据我们截至 2026 年 9 月的定向检索，尚未发现有限 hidden-node budget 下
> grounded resolvent port compression 的上述精确 hull 和 robust-approximation
> consequence；基础 random-forest projection identity 明确属于既有工作。

不能写成“首次发现 resolvent 是 forest projections 的均值”，也不能把本定理
外推到 directed/signed/block-valued edges、非 unit grounding、有限 edge-weight
区间或任意 nonlinear estimator。

## 8. 适用范围与容易写错之处

- 本报告证明的是 scalar、undirected、nonnegative edge weights、unit grounding。
- 若 edge weights 有统一有限上界，high-conductance converse 不成立，(5) 只剩
  universal upper envelope。
- 若 connected graphs 的 edge weights 还有统一正下界，weak-link decoupling 不成立。
- 对 exactly (h)，unused hidden vertices 必须允许成为 hidden-only components，
  或在 connected closure 中通过趋零弱边连接。
- (X_{\pi,k}) 一般不是 projection；完整输入误差必须保留
  (C(X-X^2)C^T)。
- terminal hull 的非vertex pruning 对 nuclear/(p<2) full-input risk 不自动安全；
  全 Schatten 范数的无条件 finite certificate 使用全部 (N(q,h)) 场景。
- 定理给 fixed boundary size 的 exact finite method，不是关于 (q) 的 polynomial
  algorithm；场景数仍随 Stirling numbers 快速增长。

## 9. 数值复核

独立脚本：

`experiments/theory_search/finite_hidden_budget_checks.py`

它包含六组 deterministic tests：

1. 枚举 edge-DPP subsets，重构 ((I+L)^{-1}) 并逐 forest 核对场景公式；
2. 随机含环 graph 上核对 Schatten (p=1,2,4,\infty) upper bound；
3. connected maximum-degree-2 path 的 high-conductance/weak-link converse；
4. (q=3,h=2) 的 31 个候选和 21 个 predicted vertices 的 LP 核对；
5. Hausdorff construction；
6. (q=2,h=1) 解析收益例及全部五个场景值。

运行结果：

```text
Ran 6 tests in 0.035s
OK
```

数值检查只用于防止 hidden-column Gram、allocation 和极限量词写错；定理依据仍是
第 2--6 节的证明。

## Sources

1. Pavel Chebotarev and Elena Shamis, [“Matrix-Forest Theorems”](https://arxiv.org/abs/math/0602575), arXiv:math/0602575; earlier journal version, *Automation and Remote Control* 58(9), 1997. Local: `literature/09_network_completion/2006_chebotarev_shamis_matrix_forest_theorems.pdf`.
2. Yusuf Yiğit Pilavcı, Pierre-Olivier Amblard, Simon Barthelmé, and Nicolas Tremblay, [“Graph Tikhonov Regularization and Interpolation via Random Spanning Forests”](https://arxiv.org/abs/2011.10450), *IEEE TSIPN* 7 (2021), DOI [10.1109/TSIPN.2021.3084879](https://doi.org/10.1109/TSIPN.2021.3084879). Local: `literature/09_network_completion/2021_pilavci_et_al_graph_tikhonov_random_forests.pdf`.
3. Michał Dereziński, Rajiv Khanna, and Michael W. Mahoney, [“Improved Guarantees and a Multiple-Descent Curve for Column Subset Selection and the Nyström Method”](https://arxiv.org/abs/2002.09073), NeurIPS 2020. Local: `literature/09_network_completion/2020_derezinski_khanna_mahoney_dpp_projection.pdf`.
4. Adrien Kassel and Thierry Lévy, [“On the Mean Projection Theorem for Determinantal Point Processes”](https://doi.org/10.30757/ALEA.v20-17), *ALEA* 20 (2023), 497--504. Local: `literature/09_network_completion/2023_kassel_levy_mean_projection_dpp.pdf`.
5. Florian Dörfler and Francesco Bullo, [“Kron Reduction of Graphs with Applications to Electrical Networks”](https://arxiv.org/abs/1102.2950), *IEEE TCAS-I* 60(1), 2013, DOI [10.1109/TCSI.2012.2215780](https://doi.org/10.1109/TCSI.2012.2215780). Local: `literature/09_network_completion/2013_dorfler_bullo_kron_reduction.pdf`.
6. Antonio De Rosa and Aida Khajavirad, [“The Ratio-Cut Polytope and K-Means Clustering”](https://arxiv.org/abs/2006.15225), *SIAM Journal on Optimization* 32(1), 2022, DOI [10.1137/20M1348601](https://doi.org/10.1137/20M1348601). Local: `literature/09_network_completion/2022_de_rosa_khajavirad_ratio_cut_polytope.pdf`.
7. Thomas Lam, [“Electroid Varieties and a Compactification of the Space of Electrical Networks”](https://arxiv.org/abs/1402.6261), *Advances in Mathematics* 338 (2018), DOI [10.1016/j.aim.2018.09.014](https://doi.org/10.1016/j.aim.2018.09.014). Local: `literature/09_network_completion/2018_lam_electroid_varieties.pdf`.
8. Richard W. Kenyon and David B. Wilson, [“Boundary Partitions in Trees and Dimers”](https://arxiv.org/abs/math/0608422), *Transactions of the AMS* 363(3), 2011, DOI [10.1090/S0002-9947-2010-04964-5](https://doi.org/10.1090/S0002-9947-2010-04964-5). Local: `literature/09_network_completion/2006_kenyon_wilson_boundary_partitions.pdf`.
