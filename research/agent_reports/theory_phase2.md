# 第二阶段数学专项：多数节点、固定半径、仅局部模型的 Gaussian/WLS 估计

> 检索与核验截止：2026-09-13。本文不用新术语包装已有概念。
>
> 状态标记：**[已知]** 表示可在所引原始来源中直接找到；**[可直接证明]** 表示下面给出完整等式或足以复核的证明，不代表已经通过全球文献穷尽证明其原创性；**[猜想/缺口]** 表示尚缺关键证明或反例。

## 0. 结论先行

1. **固定全局模型、允许离线知道中心系数时，最小半径误差已经有精确答案。** 若中心线性估计的第 \(i\) 个块行为 \(T_i\)，本地数据抽取为 \(S_{i,r}\)，输入统一约束为 \(\|z\|_2\le 1\)，则任意（甚至非线性）本地解码器的最坏误差恰为

   \[
   \|T_i(I-S_{i,r}^{*}S_{i,r})\|_{2\to2}.
   \]

   这是一条精确充要条件，但本质上是信息集/行支撑结论，不是深结构定理。

2. **“算法只知道局部模型”不能仍用单个固定模型的行尾来定义。** 正确的线性 minimax 量是

   \[
   a_r(v)=\inf_{\ell}\sup_{(\theta,i):\mathcal V_r(i,\theta)=v}
   \|T_i(\theta)-\ell S_{i,r}^{\theta}\|_{2\to2}. \tag{0.1}
   \]

   它确实是共同局部线性系数存在的充要条件；但对任意非线性算法，正确对象是每个局部数据纤维上可能输出集合的 Chebyshev 半径，不能把 (0.1) 冒充任意算法的充要条件。

3. **在线性方程 \(Jx=b\) 的原始正规方程右端模型中，Schur 补给出了真正有内容的双向桥梁。** 对球 \(B\ni i\)，局部 principal/Dirichlet 解与全局最优行截断误差 \(g_i(B)\) 满足

   \[
   g_i(B)\le d_i(B)\le
   \sqrt{1+\|J_{OB}J_{BB}^{-1}\|^2}\;g_i(B). \tag{0.2}
   \]

   右侧因子只涉及球内 principal block 与割边。若全族 \(mI\preceq J\preceq MI\)，它由仅依赖 \(\kappa=M/m\) 的常数控制。因此，“离线知道全局系数的多数节点局部性”和“同一个仅看局部模型的 Dirichlet 规则”在这种模型类上等价到常数因子；无环完全不需要。

4. **原始观测 \(z\) 的 WLS 估计 \(Kz=J^{-1}H^*R^{-1}z\) 不能只研究 \(J^{-1}\)。** \(J^{-1}\) 的长程项可以与 \(H^*R^{-1}\) 精确抵消。只有当 \(G=H^*R^{-1}\) 具有有限传播且有有限传播、范数受控的右逆时，\(K\) 的测量行尾与 \(J^{-1}\) 的状态行尾才有带半径平移的双向界。

5. **对标量 node-plus-relative-measurement，即 grounded Laplacian/M-matrix，存在最干净的结构解释。** 写

   \[
   J=S(I-P),\qquad J^{-1}=(I-P)^{-1}S^{-1},\qquad (I-P)^{-1}=\sum_{t\ge0}P^t,
   \]

   并把输入归一成 \(c=S^{-1}b\)。在 \(\|c\|_\infty\le1\) 下：固定模型的最优球外误差恰是 killed random walk 在球外的期望总占用时间；局部 Dirichlet 误差恰是首次出球后、被 kill 前的全部剩余生存时间；\(r\) 阶 Neumann 算法误差恰是第 \(r+1\) 步后的剩余生存时间。这些表示本身是经典 Green/potential theory，不可声称原创；把三者、局部模型 minimax 与多数节点分位数放在同一 theorem 中，尚未检索到直接先例。

6. **最重要的 no-go：** 若允许同一局部视图的任意规模 SPD/M-matrix completion，球外既无统一 resolvent/lifetime 上界，输入又只约束 \(\|b\|\)，那么任何有非零边界传递的视图通常都有 \(a_r(v)=\infty\)。在 grounded Laplacian 中，给球外接一个任意长、无 ground 的“悬挂支路”，度数和权重都可保持有界，而根节点的 \(\ell_2\) 行尾按 \(\sqrt N\) 增长、\(\ell_\infty\) 对偶行尾按 \(N\) 增长。允许少量节点失败并不能自动消除此障碍：单锚点路径只有 \(o(n)\) 个特殊节点，却能让几乎所有节点依赖远端锚点。

7. **截至 2026 年的撞车判断：** Gaussian screening、walk-sums、稀疏 SPD 逆衰减、local linear solvers、GFF/potential theory、measured coarse geometry 都覆盖了重要组成部分；但未找到同时处理“共同局部系数 + 任意全局 completion + 节点测度 \(1-\delta\) + WLS/Gaussian”的现成严格 theorem。最可发表的方向不是再写定义式 iff，而是证明第 2 节的 Schur-oracle sandwich 及第 3 节的 killed-walk completion dichotomy，并把适用范数和失败量词写死。

## 1. 先固定量词：否则 \(a_r\) 会被错误地写成无穷或 tautology

### 1.1 模型与数据必须分开

令 \(\Theta\) 是有限网络模型类。一个模型 \(\theta\in\Theta\) 包含：通信图、节点/边系数、坐标维数、数据归属规则，以及中心线性算子

\[
T_i(\theta):\mathcal Z_\theta\longrightarrow \mathcal Y_i.
\]

这里 **\(\theta\) 不包含本次实现的数值数据 \(z\)**。数据另取 \(\|z\|_{\mathcal Z_\theta}\le1\)。若把全局数据也塞进 \(\omega\)，再对不受统一范数约束的远端数据取 supremum，只要存在一个非零远端系数，误差当然就是 \(+\infty\)；那不是结构结论，只是尺度未归一化。

取正交直和

\[
\mathcal Z_\theta=\bigoplus_{u\in V(\theta)}\mathcal Z_u,
\]

令 \(S_{i,r}^{\theta}:\mathcal Z_\theta\to\mathcal Z_v\) 抽取 \(B_r(i)\) 内数据。半径 \(r\) 的带标记 rooted view 记为 \(v=\mathcal V_r(i,\theta)\)。标记须包含块维数、端口、可见系数和目标坐标类型；同构视图通过固定的等距映射识别其本地输入、输出空间。不同网络规模的球外空间不必互相识别，因为每个 operator norm 先在自己的全局空间中计算，再取 supremum。

以下 LOCAL 存在性结果默认每条消息大小和本地计算不受限：\(r\) 轮后根可收集完整 \(r\)-view。它不自动给出低 bit-hop 或低计算复杂度。

### 1.2 固定模型的最小不可避免误差

令 \(P=S^*S\) 是球内坐标投影，\(Q=I-P\)。

**定理 1 [可直接证明，固定模型 oracle]。** 对固定 \((\theta,i)\)，

\[
\inf_{f}\sup_{\|z\|_2\le1}
\|T_i z-f(Sz)\|_2
=\|T_iQ\|_{2\to2}, \tag{1.1}
\]

其中 infimum 允许所有确定性函数 \(f\)。线性规则 \(f(u)=T_iS^*u\) 已达到上界。

证明只有两行。上界是直接截断。下界取任意球外单位向量 \(q\)，则 \(Sq=S(-q)=0\)，而两个中心答案为 \(\pm T_iq\)；任意同一个 \(f(0)\) 对二者至少一个误差不小于 \(\|T_iq\|\)。对 \(q\) 取 supremum 即得。

因此固定模型的精确 \(r\)-轮条件是 \(T_iQ=0\)，误差 \(\varepsilon\) 条件是 \(\|T_iQ\|\le\varepsilon\)。这一定理允许图有任意环。

### 1.3 仅局部模型：共同系数的精确 minimax

对每个 view \(v\)，定义父任务建议的量

\[
a_r^{\rm lin}(v):=
\inf_{\ell:\mathcal Z_v\to\mathcal Y_v}
\sup_{(\theta,i):\mathcal V_r(i,\theta)=v}
\|T_i(\theta)-\ell S_{i,r}^{\theta}\|_{2\to2}. \tag{1.2}
\]

**定理 2 [可直接证明，共同局部线性系数]。** 在确定性线性 LOCAL 模型中，view \(v\) 上存在一个共同 \(r\)-轮规则、对所有 completion 和所有 \(\|z\|_2\le1\) 误差至多 \(\varepsilon\)，当且仅当 (1.2) 的 infimum 可由某个 \(\ell\) 以值不超过 \(\varepsilon\) 达到。若本地空间有限维，且该 view 的 completion 满足 \(\sup\|T_i(\theta)\|<\infty\)，则 infimum 一定达到，故条件简化为 \(a_r^{\rm lin}(v)\le\varepsilon\)。

理由是：一个 uniform 线性 LOCAL 算法在同构 view 上必须使用同一个 \(\ell\)，反之任意 \(\ell_v\) 都定义一个 \(r\)-LOCAL 规则。达到性来自目标函数连续且

\[
\sup_\theta\|T_i(\theta)-\ell S\|\ge \|\ell\|-\sup_\theta\|T_i(\theta)\|,
\]

所以它在有限维空间上 coercive。

若 completion 为单点，(1.2) 精确退化为 (1.1)。若把 \(T_i=[A_\theta,B_\theta]\) 按本地/球外输入分块，则

\[
a_r^{\rm lin}(v)=\inf_\ell\sup_\theta
\|[A_\theta-\ell,\ B_\theta]\|, \tag{1.3}
\]

并有不可忽略的 completion ambiguity 下界

\[
a_r^{\rm lin}(v)\ge
\max\left\{
\sup_\theta\|B_\theta\|,
\frac12\sup_{\theta,\theta'}\|A_\theta-A_{\theta'}\|
\right\}. \tag{1.4}
\]

只有当本地块 \(A_\theta\) 已由 view 唯一确定时，才有

\[
a_r^{\rm lin}(v)=\sup_\theta\|B_\theta\|. \tag{1.5}
\]

因此，“所有 completion 的远端尾都小”一般仍不够；全局 completion 还可能改变中心算子的本地系数。

**任意非线性规则。** 对本地数据 \(u\) 定义输出不确定集

\[
\mathcal C_v(u)=\left\{
T_i(\theta)z:\mathcal V_r(i,\theta)=v,\ Sz=u,\ \|z\|_2\le1
\right\}.
\]

忽略可测选择的技术问题时，任意解码器的精确 minimax 是

\[
a_r^{\rm any}(v)=
\sup_{\|u\|\le1}\operatorname{rad}\mathcal C_v(u), \tag{1.6}
\]

其中 \(\operatorname{rad}C=\inf_y\sup_{x\in C}\|x-y\|\)。跨 completion 时最优中心可随 \(u\) 非线性变化，所以通常只能说 \(a_r^{\rm any}\le a_r^{\rm lin}\)，不能把 (1.2) 宣称成任意算法的 iff。

### 1.4 “允许 \(\delta\) 节点失败”有三个不同版本

令 \(\mu_\theta\) 是节点概率测度，并定义分位数

\[
Q_{1-\delta}^{\mu}(e):=
\inf\{t:\mu\{i:e_i\le t\}\ge1-\delta\}. \tag{1.7}
\]

对 view-indexed 线性规则 \(\ell=(\ell_v)_v\)，令

\[
e_i(\theta;\ell)=
\|T_i(\theta)-\ell_{\mathcal V_r(i,\theta)}S_{i,r}^\theta\|.
\]

则真正的“每个好节点对所有单位输入都好”的 class-level 数值是

\[
A^{\rm loc}_{r,\delta}:=
\inf_{(\ell_v)}\sup_{\theta\in\Theta}
Q_{1-\delta}^{\mu_\theta}
\big(e_i(\theta;\ell)\big). \tag{1.8}
\]

**[可直接证明]** 若 infimum 达到，存在一个仅看本地模型的线性规则，使每个模型至少 \(1-\delta\) 测度节点对所有 \(\|z\|\le1\) 误差不超过 \(\varepsilon\)，iff \(A^{\rm loc}_{r,\delta}\le\varepsilon\)。一般情形严格地说是：任意 \(\eta>0\) 可达 \(A^{\rm loc}_{r,\delta}+\eta\)。

不能用简单的

\[
\sup_\theta\mu_\theta\{i:a_r^{\rm lin}(\mathcal V_r(i,\theta))>\varepsilon\}\le\delta \tag{1.9}
\]

代替 (1.8)。式 (1.9) 是充分条件，但一般不必要：不同 view 的最坏 completion 可能互不相容，无法在同一个全局模型里同时出现。只有模型类具有足够强的 disjoint-union/pasting 闭包时，才可能把逐 view supremum 与每实例分位数交换。

还必须区分下列较弱版本：

* **每个输入只要求多数节点同时正确：**

  \[
  \inf_\ell\sup_{\theta,\|z\|\le1}
  \mu_\theta\{i:\|(T_i-\ell_iS_i)z\|>\varepsilon\}.
  \]

  它不等于行 operator norm 的分位数。例：误差算子为 \(I_n\)，每一行最坏误差都是 1；但任意 \(\|z\|_2\le1\) 中超过 \(\varepsilon\) 的坐标至多 \(1/\varepsilon^2\) 个，所以固定 \(\delta,\varepsilon\) 时该版本随 \(n\) 增大可以成功。

* **随机数据/Bayes 失败率：** 好坏还对数据分布平均，见第 2.5 节。

* **实现后按实际残差挑好节点：** 这是事后保证，弱于事前的 nodewise robust 保证。

后文的“多数节点”默认采用 (1.8)，因为它最接近“允许一组节点失败，其余节点有一致误差保证”。

对固定模型的 oracle 行尾

\[
g_i(r)=\|T_i(I-S_{i,r}^*S_{i,r})\|,
\]

条件

\[
\forall\varepsilon,\delta>0\ \exists r<\infty:\quad
\sup_\theta\mu_\theta\{i:g_i(r)>\varepsilon\}\le\delta \tag{1.10}
\]

就是随机根 \(I\sim\mu_\theta\) 的行尾对模型类一致地依概率趋零。若 \(g_i(r)^2\) 一致可积，且输出块维数一致有界，则 (1.10) 等价于加权 normalized Frobenius 尾趋零。它不是全网 \(\ell_2\to\ell_2\) operator norm 局部化。

## 2. 线性 WLS：行尾、Schur 边界影响和局部有限子问题

### 2.1 先研究 primitive normal RHS：\(Jx=b\)

设 \(J=J^*\succ0\)，目标是根块 \(x_i=P_ix=P_iJ^{-1}b\)。取含根的球 \(B\)，外部为 \(O\)，分块写成

\[
J=\begin{bmatrix}A&E\\E^*&D\end{bmatrix},\qquad
T:=D-E^*A^{-1}E\succ0. \tag{2.1}
\]

固定模型、只使用 \(b_B\) 的最佳 oracle 误差为

\[
g_i(B):=\|P_iJ^{-1}P_O\|. \tag{2.2}
\]

局部 principal/零 Dirichlet 规则为

\[
L_i^{D}b:=P_iA^{-1}b_B,
\qquad
d_i(B):=\|P_iJ^{-1}-P_iA^{-1}P_B\|. \tag{2.3}
\]

**定理 3 [可直接证明，精确 Schur 因子分解]。** 令

\[
Y:=P_iJ^{-1}P_O=-P_iA^{-1}ET^{-1},\qquad
F:=E^*A^{-1}.
\]

则

\[
P_iJ^{-1}-P_iA^{-1}P_B=Y[-F,\ I_O], \tag{2.4}
\]

从而

\[
g_i(B)\le d_i(B)\le
\gamma_i(B)g_i(B),\qquad
\gamma_i(B):=\sqrt{1+\|F\|^2}. \tag{2.5}
\]

这里 \(\|[-F,I]\|=\sqrt{1+\|F\|^2}\) 是等式。下界只需把输入限制为 \([0,b_O]\)。另一个等价恒等式是

\[
P_BJ^{-1}-A^{-1}P_B=-A^{-1}EP_OJ^{-1}. \tag{2.6}
\]

这精确说明：Dirichlet 算法的额外误差来自球外扰动通过边界进入，以及已经离开球的信息重新影响球内系数。

若 radius-view 包含 \(A\) 以及 cut coupling \(E\)（标准端口模型中边界节点知道其 incident edge；否则多取一圈 halo），则 \(P_iA^{-1}\) 和 \(\gamma_i(B)\) 都只由局部模型决定。

### 2.2 completion minimax 与多数节点的双向 theorem

考虑所有共享同一 closed view、因而共享 \(A,E\) 的 completion。由 block inverse，中心根行可写成

\[
P_iJ^{-1}=[P_iA^{-1}-YF,\ Y]. \tag{2.7}
\]

因此该 view 的共同线性系数误差满足

\[
\sup_{\theta\succ v}\|Y_\theta\|
\le a_r^{\rm lin}(v)
\le \gamma(v)\sup_{\theta\succ v}\|Y_\theta\|. \tag{2.8}
\]

下界是每个 completion 的固定模型 oracle 尾；上界由同一个 Dirichlet 系数 \(P_iA^{-1}\) 达到。式 (2.8) 比“Chebyshev 半径定义”更有结构：在 primitive RHS 的 SPD 问题中，本地系数 ambiguity 与外部行尾不会完全脱钩，二者由 Schur 补强制绑定。

定义离线 oracle 的多数节点值

\[
G^{\rm off}_{r,\delta}:=
\sup_{\theta}Q_{1-\delta}^{\mu_\theta}\big(g_i(B_r(i))\big). \tag{2.9}
\]

**定理 4 [可直接证明，多数节点 Dirichlet–oracle sandwich]。** 若所有节点的 \(\gamma_i(B_r(i))\le\Gamma\)，则

\[
G^{\rm off}_{r,\delta}
\le A^{\rm loc}_{r,\delta}
\le \Gamma G^{\rm off}_{r,\delta}. \tag{2.10}
\]

左侧因为任何 local-model rule 在固定模型上也不可能优于 (2.2)；右侧使用同一个 view-to-Dirichlet 映射和 (2.5)。因此在 \(\Gamma\) 一致有界的模型类上，oracle 行尾按节点测度消失，**当且仅当** 存在仅看局部模型的 Dirichlet 规则按节点测度消失；二者半径只差可选的一圈 coefficient halo，误差只差常数。

不必要求所有节点的 \(\gamma_i\) 都好。阈值形式更稳妥：若某模型中至少 \(1-\delta_1\) 节点有 \(g_i\le a\)，至少 \(1-\delta_2\) 节点有 \(\gamma_i\le b\)，则至少 \(1-\delta_1-\delta_2\) 节点的 Dirichlet 误差不超过 \(ab\)。这允许少量局部病态边界区域，但 **\(g_i\) 本身仍是全局 exterior influence，不能仅靠球内谱隙替代。**

若存在统一谱窗

\[
mI\preceq J_\theta\preceq MI,\qquad \kappa=M/m, \tag{2.11}
\]

则 principal block \(A\) 和外部 Schur 补 \(T\) 也落在 \([m,M]\)。又因

\[
E=P_B\left(J-\frac{M+m}{2}I\right)P_O,
\]

有更尖的 \(\|E\|\le(M-m)/2\)，故

\[
\gamma_i(B)\le
\sqrt{1+\left(\frac{\kappa-1}{2}\right)^2}
\le\sqrt{1+\kappa^2}. \tag{2.12}
\]

这给出不依赖网络规模、环数或树宽的 (2.10)。

还可定义完全局部的 boundary-transfer 数值

\[
\tau_i(B):=\|P_iA^{-1}E\|. \tag{2.13}
\]

由 \(Y=-P_iA^{-1}ET^{-1}\) 得

\[
\frac{\tau_i(B)}{M}\le g_i(B)\le\frac{\tau_i(B)}m,
\qquad
d_i(B)\le\gamma_i(B)\frac{\tau_i(B)}m. \tag{2.14}
\]

注意上界仍用了 exterior Schur 补的统一 lower bound \(m\)。若没有任何 exterior stability promise，\(\tau_i\) 小并不排除 \(T^{-1}\) 极大。

**与 Demko/Chebyshev 的关系 [已知 + 直接推论]。** 若 \(J\) 的传播半径为 \(s\)，(2.11) 成立，取 \(1/x\) 在 \([m,M]\) 上的 \(k\) 次最佳多项式 \(p_k\)，则 \(p_k(J)\) 的传播半径至多 \(ks\)，且 Demko–Moss–Smith Proposition 2.1（PDF p. 492）给出

\[
\|J^{-1}-p_k(J)\|\le C_0q^{k+1},\quad
q=\frac{\sqrt\kappa-1}{\sqrt\kappa+1},\quad
C_0=\frac{(1+\sqrt\kappa)^2}{2m\kappa}. \tag{2.15}
\]

故每个 \(g_i(B_{ks}(i))\le C_0q^{k+1}\)，再由 (2.10) 得 only-local Dirichlet 保证。指数衰减本身是 1984 年以来的经典充分条件，不是新贡献。

### 2.3 从“principal 子问题”到“只保留球内测量”的差别

实际本地 WLS 常把跨边界因子直接丢掉，得到 \(\widetilde A\)，而不是使用全局正规矩阵的 principal block \(A=J_{BB}\)。二者不能混写。恒等式

\[
A^{-1}-\widetilde A^{-1}
=A^{-1}(\widetilde A-A)\widetilde A^{-1} \tag{2.16}
\]

只给出 triangle upper bound；除非有 Loewner 单调性、符号或 range 对齐假设，没有反向下界，因为 boundary error 与 model-dropping error 可以抵消。要实现定理 3 的规则，必须让 boundary nodes 报告构成 \(J_{BB}\) 的全部对角/内部项；如果一个跨边界测量对 \(J_{BB}\) 有贡献，就需要把该因子信息纳入 halo。

Cheng–Jiang–Sun 的 localized finite-system 结果与此最接近：在 polynomial-growth/doubling 图和 \(J_\alpha\) 衰减类中，其 Theorems 6.1–6.2（本地 PDF pp. 15–16）给出全局 stability 与所有足够大 quasi-main restrictions 统一 stability 的条件化双向桥梁；Proposition 7.1（pp. 17–18）给出局部有限系统在半球内的误差率。它并不覆盖任意 SPD completion 或 \(\delta\)-比例异常根。

### 2.4 原始观测 WLS：必须保留 \(G=H^*R^{-1}\)

写

\[
J=H^*R^{-1}H+Q\succ0,
\qquad G=H^*R^{-1},
\qquad K=J^{-1}G. \tag{2.17}
\]

若本地测量集合为 \(M\)，局部 principal 规则取

\[
L_{i,B,M}=P_iA^{-1}P_BG P_M.
\]

则有精确分解

\[
P_iK-L_{i,B,M}
=\underbrace{(P_iJ^{-1}-P_iA^{-1}P_B)G}_{\text{边界/Schur 误差}}
+\underbrace{P_iA^{-1}P_BG(I-P_M)}_{\text{未收集 RHS 贡献}}. \tag{2.18}
\]

两个项可能抵消，因此一般只有三角上界。若 \(P_BG(I-P_M)=0\)，第二项消失；这要求所有能贡献到 \(b_B=Gz\) 的测量都已收集。

设 \(G\) 从 measurement-owner metric 到 state metric 的传播半径至多 \(h\)。定义

\[
j_i(t)=\|P_iJ^{-1}Q^{x}_{\overline B_t(i)}\|,\qquad
k_i(r)=\|P_iKQ^{z}_{\overline B_r(i)}\|.
\]

则

\[
k_i(r)\le\|G\|j_i(r-h). \tag{2.19}
\]

若还存在右逆 \(L_G\) 满足 \(GL_G=I\)、\(\|L_G\|\le\Lambda\)，且传播半径至多 \(h'\)，则

\[
j_i(r+h')\le\Lambda k_i(r). \tag{2.20}
\]

所以只有“\(G\) 局部 + 有局部受控右逆”时，两种行尾才双向等价。\(H\) 局部并不自动意味着 \(G\) 局部：稠密噪声精度 \(R^{-1}\) 会直接产生长程耦合。

### 2.5 Gaussian/Bayes 版本的精确风险等式

令 \((X_i,Z)\) 联合 Gaussian，中心后验均值为 \(Y=T_iZ=\mathbb E[X_i\mid Z]\)。把 \(Z=(Z_B,Z_O)\)，协方差写成 \(\Sigma_Z\)。本地 Bayes 最优估计是

\[
\mathbb E[X_i\mid Z_B]=\mathbb E[Y\mid Z_B].
\]

**定理 5 [标准 Gaussian 投影 + 可直接证明]。** 局部相对中心的 excess Bayes risk 恰为

\[
\begin{aligned}
\Delta_i(B)
&=\mathbb E\|Y-\mathbb E[Y\mid Z_B]\|^2\\
&=\operatorname{tr}\!\left(
T_{i,O}\Sigma_{O\mid B}T_{i,O}^*
\right),\\
\Sigma_{O\mid B}
&=\Sigma_{OO}-\Sigma_{OB}\Sigma_{BB}^{-1}\Sigma_{BO}.
\end{aligned} \tag{2.21}
\]

而且由条件期望的正交投影恒等式，\(\Delta_i(B)\) 正是“仅本地数据估计 \(X_i\)”与“全部数据估计 \(X_i\)”的 MSE 之差。

若 \(cI\preceq\Sigma_Z\preceq CI\)，则其 Schur 补也满足同一谱窗，故

\[
c\|T_{i,O}\|_F^2
\le\Delta_i(B)\le
C\|T_{i,O}\|_F^2. \tag{2.22}
\]

标量输出时 Frobenius 行尾就是 \(\ell_2\) 行尾；固定输出块维数时二者只差维数常数。没有协方差谱窗，operator coefficient tail 与 Bayes 价值不等价：远端数据可几乎由本地数据预测，或其条件方差可随规模爆炸。

Stein 的 screening effect 正是 (2.21) 的近亲，但其判据通常是

\[
\operatorname{MSE}(B\cup O)/\operatorname{MSE}(B)\to1,
\]

即相对 Bayes risk，而不是 adversarial operator tail；研究的是连续空间、单预测点和特定渐近采样几何。2011 年论文只证明其一般猜想的两个特殊情形，2026 年 Meysami–Lotfi 仍是在更广 spectral 条件下给充分条件和速率，并未给本报告的 graph/majority/local-completion iff。

## 3. Grounded Laplacian/M-matrix：精确 killed-walk 表示

### 3.1 正规化

考虑标量相对测量与 node anchors：

\[
J=L_W+\operatorname{diag}(\kappa),\qquad
s_i=\sum_jw_{ij}+\kappa_i,
\qquad S=\operatorname{diag}(s_i), \tag{3.1}
\]

\[
P_{ij}=w_{ij}/s_i,\qquad
J=S(I-P). \tag{3.2}
\]

\(P\) 是 substochastic transition matrix，节点 \(i\) 每步以 \(\kappa_i/s_i\) 的概率进入 cemetery。若每个 connected component 都能到达正 killing，\(\rho(P)<1\)，且

\[
J^{-1}=(I-P)^{-1}S^{-1},\qquad
N:=(I-P)^{-1}=\sum_{t\ge0}P^t. \tag{3.3}
\]

必须说明输入归一化。令 \(c=S^{-1}b\)，则 \(x=Nc\)。下面先约束 \(\|c\|_\infty\le1\)。若约束原始 \(\|b\|_\infty\le1\)，所有 occupation 都带 reward \(1/s_j\)；仅当 \(s_j\) 有统一上下界时两种范数才等价到常数。

令 \(X_t\) 是从 \(i\) 出发的 killed chain，\(\zeta=\inf\{t:X_t=\dagger\}\)。于是

\[
N_{ij}=\mathbb E_iL_j,\qquad
L_j:=\sum_{t=0}^{\zeta-1}\mathbf1\{X_t=j\}. \tag{3.4}
\]

这是 absorbing Markov chain/离散 potential theory 的标准 fundamental matrix 解释；Malioutov–Johnson–Willsky 的 walk-sums、Bendito–Carmona–Encinas 的 finite-network Green/Poisson kernels，以及 Sheffield 的 GFF covariance 都属于同一经典结构。

### 3.2 三种 \(r\)-局部误差的精确概率公式

取 \(B=B_r(i)\)，令

\[
\tau_B=\inf\{t\ge0:X_t\in V\setminus B\},
\]

若被 kill 前未出球则 \(\tau_B=\infty\)。

**定理 6 [经典表示的直接组合]。** 对标量根输出、\(\|c\|_\infty\le1\)：

1. 固定全局模型、允许 oracle 选择最优球内系数时，精确误差是

   \[
   o_i(r):=\sum_{j\notin B}N_{ij}
   =\mathbb E_i\sum_{t<\zeta}\mathbf1\{X_t\notin B\}. \tag{3.5}
   \]

   即被 kill 前在球外的期望总占用时间。

2. 局部 Dirichlet 解 \(x_i^D=e_i^T(I-P_{BB})^{-1}c_B\) 的精确误差是

   \[
   d_i(r)=
   \mathbb E_i\!\left[(\zeta-\tau_B)\mathbf1\{\tau_B<\zeta\}\right]. \tag{3.6}
   \]

   它计算首次出球以后直到被 kill 的所有 live time，包括离开后又返回球内的路径。因此它通常严格大于 oracle 球外行尾。

3. \(r\) 阶 Neumann 规则 \(\sum_{t=0}^{r}P^tc\) 的精确误差是

   \[
   n_i(r)=e_i^T\sum_{t=r+1}^{\infty}P^t\mathbf1
   =\mathbb E_i[(\zeta-r-1)_+]. \tag{3.7}
   \]

逐路径比较给出

\[
o_i(r)\le d_i(r)\le n_i(r), \tag{3.8}
\]

因为走出 \(r\)-球至少需要 \(r+1\) 步。证明也可直接把 \(N\) 分成从未出球的 paths 与已经出球的 paths；前者恰是 \((I-P_{BB})^{-1}\)。

这组公式的用途不是宣称 Green 函数新表示，而是把三件常被混淆的事分开：全局 oracle 行截断、仅局部模型的 finite Dirichlet solve、以及有限时间 Neumann/message passing。

### 3.3 逃逸概率何时成为局部证书

定义在 kill 前出球的概率

\[
h_i(r):=\Pr_i(\tau_B<\zeta)
=e_i^T(I-P_{BB})^{-1}P_{BO}\mathbf1. \tag{3.9}
\]

若 boundary outgoing probabilities 是 view 的一部分，\(h_i(r)\) 可完全由局部模型计算。每次成功出球至少贡献一次球外访问，故

\[
h_i(r)\le o_i(r)\le d_i(r). \tag{3.10}
\]

若从可能的 exit states 出发的剩余期望寿命统一有界

\[
\sup_{y\in\partial^+B}\mathbb E_y\zeta\le L, \tag{3.11}
\]

强 Markov 性给出

\[
h_i(r)\le o_i(r)\le d_i(r)\le Lh_i(r). \tag{3.12}
\]

所以在具有统一 lifetime envelope 的模型类上，“多数根的 escape probability 趋零”与 oracle/Dirichlet 误差按节点测度趋零是等价到常数的充要条件。它不要求一个控制全网所有方向的 Euclidean condition number；若 (3.11) 只对绝大多数根的 exit distribution 成立，仍可用好集交集得到 \(\delta_1+\delta_2\) 的失败率。

但 (3.11) 绝不是可以省去的小技术条件。仅知道 \(h_i(r)\) 很小，不知道逃出去后会活多久，\(d_i(r)\) 可任意大。全局每步 killing 至少 \(p_0\) 是一个强但简单的充分条件，此时 \(L\le1/p_0\)。

类似地，令 survival probability

\[
s_i(r)=e_i^TP^{r+1}\mathbf1=\Pr_i(\zeta>r+1).
\]

若所有可能的 time-\(r+1\) states 余寿命至多 \(L\)，则

\[
s_i(r)\le n_i(r)\le Ls_i(r). \tag{3.13}
\]

### 3.4 completion-rich 模型类的 no-go

**命题 7 [可直接证明，抽象 SPD/M-matrix completion]。** 固定一个局部分块 \(A,E\)。若根到 cut 的 transfer \(C=P_iA^{-1}E\ne0\)，且模型类允许选择 exterior block

\[
D_\eta=E^*A^{-1}E+\eta I,\qquad\eta\downarrow0,
\]

则每个 completion 都 SPD，但 Schur 补 \(T_\eta=\eta I\)，所以

\[
\|P_iJ_\eta^{-1}P_O\|=\|C\|/\eta\to\infty.
\]

因此同一 view 的 \(a_r^{\rm lin}(v)=\infty\)。若要同时留在 symmetric M-matrix 类，可取一维 exterior、令 \(A\) 为 nonsingular M-matrix 且 \(E\) 为非正单列；此时 \(D_\eta\) 是正标量，整个矩阵仍有非正 off-diagonal。多维 exterior 的上式一般会在 \(D_\eta\) 中产生正 off-diagonal，不能不加说明地称为 M-matrix。

**命题 8 [可直接证明，grounded-Laplacian 有界度版本]。** 即使限制为 \(J=L_W+\operatorname{diag}\kappa\)、度数和非零权重一致有界，任意有 open boundary 的固定球后接一条任意长、球外无 killing 的悬挂支路，也可使根的球外行尾无界。若球内至少有一个 ground，整个有限网络仍 SPD。对支路上任意源点 \(j\)，由 reciprocity：在根 \(i\) 注入单位电流时，无 sink 的悬挂支路不承载净电流，整条支路电位等于 attachment 电位；故

\[
(J^{-1})_{ij}=c_i>0
\]

与支路长度和 \(j\) 无关。于是含 \(N\) 个球外节点时

\[
\|e_i^TJ^{-1}P_O\|_2\ge c_i\sqrt N,
\qquad
\|e_i^TJ^{-1}P_O\|_1\ge c_iN. \tag{3.14}
\]

这说明：

* “only-local-model + 任意规模 completion + 仅约束输入范数”没有非平凡、纯局部、统一误差证书；
* bounded degree、bounded coefficients、甚至每个固定半径 principal ball 都良态，仍不够；
* 必须加入某种 class-level resolvent/lifetime/tightness、规模上界、completion 分布，或改用相对误差；
* 允许 \(\delta\) 节点失败只有在病态 exterior 对大多数根不可达或影响尾确实小的时候才有用，不能只数“病态节点”占比。

### 3.5 \(\ell_2\) 输入：两条独立 killed walks 的交叉局部时间

对归一化 \(c\) 的 \(\ell_2\) 单位球，固定模型 oracle 行尾为

\[
o_{i,2}(r)^2
=\sum_{j\notin B}N_{ij}^2. \tag{3.15}
\]

取两条从 \(i\) 出发、相互独立的 killed walks，局部时间分别为 \(L_j^{(1)},L_j^{(2)}\)，则

\[
o_{i,2}(r)^2
=\mathbb E_i\sum_{j\notin B}L_j^{(1)}L_j^{(2)}. \tag{3.16}
\]

这是“球外 intersection local time”的精确表示。它可把 Green-row \(\ell_2\) locality 转化成两条 walk 的相交估计，但仅是经典 occupation identity 的直接推论。逃逸概率不能单独双向控制 (3.15)：相同总 occupation 可集中在一个节点或分散到巨大体积，\(\ell_1\) 与 \(\ell_2\) 会给完全不同答案。对原始 \(b\) 再乘 reward \(s_j^{-2}\)。

### 3.6 为什么一般 block 系统不能直接概率化

标量 M-matrix 的关键是 \(N_{ij}\ge0\)，所以 operator \(\ell_\infty\) 行范数就是行和，路径不会抵消。一般 block SPD 的 off-diagonal 是矩阵，路径乘积有旋转和符号抵消：

* 若存在共同不变正锥、系统是 cone-preserving block M-matrix，可在扩展状态空间保留正概率解释；
* 更一般地，若 \(J=D(I-R)\)，可定义标量 majorant \(Q_{ij}=\|R_{ij}\|\)。在 \(\rho(Q)<1\) 或 block row diagonal dominance 下，\(\|(I-R)^{-1}_{ij}\|\le(I-Q)^{-1}_{ij}\)，给出充分衰减界；
* 该 majorant 不是必要条件，因为 matrix path cancellations 可能使真实 Green kernel 很小。这与 Gaussian walk-summability \(\rho(|R|)<1\) 在一般 loopy 图上只充分、不必要的限制相同。

因此 killed-walk occupation 的精确必要充分结论应先限定在标量 attractive/M-matrix 类；对 generic blocks，可靠主工具仍是第 2 节的 Schur 因子分解。

## 4. 必须通过的反例审计

### 4.1 稀有锚点与近零全局模态

取路径状态 \(x_1,\ldots,x_n\)，测量

\[
z_0=x_1,\qquad z_k=x_{k+1}-x_k,\quad k=1,\ldots,n-1,
\]

单位噪声。测量矩阵 \(H\) 是可逆三角型，故全局 WLS 唯一；

\[
J=H^*H=L_{P_n}+e_1e_1^*.
\]

\(\lambda_{\min}(J)=\Theta(n^{-2})\)：取缓慢线性增长的试验向量即得 \(O(n^{-2})\) 上界，离散 Poincaré 给匹配下界。模型上只有一个 anchor，是 \(o(n)\) 个特殊节点，但它控制全局 gauge。

令所有 relative data 为 0，而 anchor data 分别为 \(+1\) 和 \(-1\)。两实例的中心答案分别是 \(x_i\equiv+1\) 和 \(x_i\equiv-1\)。距 anchor 超过 \(r\) 的节点拥有完全相同的局部模型和局部数据。对每个这种节点，同一个算法输出不可能同时距 \(\pm1\) 都小于 1；把两实例的失败节点数相加，至少是远端节点数，所以至少一个实例让约一半远端节点失败。该例同时击穿：

* “坏模型节点只有 \(o(n)\)，所以其余节点必然可局部估计”；
* “每个固定半径球都可逆/看起来良态，所以全局 completion 不重要”；
* 允许少量节点失败可自动绕过单个锚点。

Pirani–Sundaram 对 grounded Laplacian 的 Theorem 1 及后续图类结果给出系统的小特征值界，并证明单 ground 的随机 \(d\)-regular 图有 \(\lambda_{\min}=\Theta(d/n)\)；路径甚至更差。其对象是谱，不是本报告的 LOCAL iff，但支持“稀有 ground 可产生全局脆弱模态”。

### 4.2 同一局部视图、不同远端 completion

最小矩阵例取

\[
J_d=\begin{bmatrix}1&-c\\-c&d\end{bmatrix},\qquad d>c^2.
\]

根只知道 \(1\) 和 cut coefficient \(c\)，不知道外部 \(d\)。而

\[
e_1^TJ_d^{-1}=\frac1{d-c^2}[d,c]. \tag{4.1}
\]

当 \(d\downarrow c^2\) 时本地系数和球外系数同时爆炸。这是 (2.7)–(2.8) 与命题 7 的最小见证，也说明仅以一个 local principal singular value 作证书不够。

### 4.3 \(J^{-1}\) 与 observation map 的精确抵消

取任意稀疏 SPD \(J\)，令 \(H=I,\ R^{-1}=J,\ Q=0\)。则正规矩阵仍是 \(J\)，但

\[
K=(H^*R^{-1}H)^{-1}H^*R^{-1}=J^{-1}J=I. \tag{4.2}
\]

即使 \(J^{-1}\) 极长程，WLS 估计器也严格 0-hop。任何声称“\(K\) 局部 iff \(J^{-1}\) 局部”的 theorem 都被此例否定；必须加入第 2.4 节的 \(G\) 条件。

### 4.4 稠密噪声精度与范数偷换

在任意 bounded-degree 通信图上取 \(H=I,\ Q=\lambda I\)，但令

\[
R^{-1}=W=I+\alpha P_n,\qquad
P_n=\frac1n\mathbf1\mathbf1^*.
\]

则

\[
K=(W+\lambda I)^{-1}W=aI+bP_n,
\]

\[
a=\frac1{1+\lambda},\qquad
b=\frac{\alpha\lambda}{(1+\lambda)(1+\alpha+\lambda)}. \tag{4.3}
\]

固定半径外的每行 \(\ell_2\) 尾约为 \(b/\sqrt n\to0\)，所以 nodewise \(\ell_2\)-input 误差对所有节点都消失；但 \(\ell_\infty\)-input 的对偶 \(\ell_1\) 行尾趋于 \(b\)，全网 spectral norm 的 rank-one 部分也为 \(b\)。同时 \(R^{-1}\) 本身稠密，违反局部模型可获得性。该例要求论文必须同时写清：噪声 precision 的传播、输入范数、逐节点还是全向量误差。

### 4.5 block/path cancellation

对 generic block matrices，两个大路径族或式 (2.18) 的两个误差项可精确相消。不存在只从 \(\|A\|\)、\(\|B\|\) 推出 \(\|A+B\|\) 的非零下界；也不存在从 block-norm walk majorant 反推真实 inverse tail 的一般必要条件。若要双向 theorem，必须增加正锥、range/singular-value 对齐，或局部右逆条件。

## 5. 与已有领域的准确碰撞边界

### 5.1 已经直接解决的片段

1. **稀疏 SPD 逆衰减。** Demko–Moss–Smith Proposition 2.1、2.2 和 Theorem 2.4（PDF pp. 492–494）已给带状/稀疏正定算子逆的指数衰减，速率由谱窗控制。Benzi–Razouk 等将多项式函数逼近推广到一般 sparse matrix functions。这里是强充分条件，不是最宽必要条件。

2. **Gaussian walk-sums/GaBP。** Malioutov–Johnson–Willsky Proposition 1（JMLR p. 2040）给 walk-summability 的等价刻画 \(\rho(|R|)<1\iff I-|R|\succ0\) 等；Proposition 21（p. 2049）给 GaBP 收敛充分条件。论文 p. 2051 明确指出一般图上 walk-summability 不是 GaBP 收敛必要条件。它研究指定迭代，不是任意 \(r\)-LOCAL estimator。

3. **本地线性系统算法。** Chung–Simpson 用 Dirichlet Green/heat-kernel PageRank 解给定边界子问题；Lee/Ozdaglar/Shah（后以 Ozdaglar–Shah–Yu 名义发表）用 Neumann series 和 weighted walks 异步近似一个坐标，在 \(\rho(|G|)<1\) 下收敛；Andoni–Krauthgamer–Pogrow Theorem 1.1 对稀疏 SDD、小 condition number 给单坐标 sublinear 算法，而 Theorem 1.2 对一般 bounded-sparsity、condition number \(\le3\) 的 PSD 系统仍给 \(n^{\Omega(1)}\) probes 下界。它们的 access model、相对 \(\|x\|_\infty\) 误差和随机查询复杂度与本报告的 fixed-radius completion minimax 不同。

4. **局部 finite-section stability。** Cheng–Jiang–Sun Theorems 5.3、6.1、6.2 和 Proposition 7.1 已覆盖 doubling/polynomial-growth + Jaffard decay 类中的 inverse-closed、local/global stability 和 finite solve 误差。不能声称一般 local/global stability 是新 idea；可能的新处只在按节点测度允许异常根、且不要求处处统一 Jaffard/spectral 条件。

5. **Gaussian screening。** Stein 2011 Theorems 1–2 只在一、二维非可微 stationary processes 和特定 near/far geometry 下证明 screening；论文自身把一般版列为 Conjecture 1。Stein 2015 给不成立的谱反例；Meysami–Lotfi 2026 扩到 O-regular variation、临界指数、非平稳/各向异性和 Delone sampling 的充分条件/速率。它们与 (2.21) 的 Bayes 语义直接相邻，但没有 local-model completion 或多数图节点 theorem。

6. **potential theory/GFF。** \(N_{ij}=\mathbb E_iL_j\)、Dirichlet Green kernel、Poisson kernel、GFF covariance 与 domain Markov property 都是经典。Bendito–Carmona–Encinas 2005 系统处理 finite Schrödinger networks；Sheffield 2007 是 GFF 标准综述；Carmona–Encinas–Jiménez–Martín 2024 明确研究 symmetric M-matrix 关联的 random walks。式 (3.5)–(3.7) 的每个概率对象都不应单独报原创。

### 5.2 measured coarse geometry 与 \(\delta\)-node criterion 不是同一个 theorem

Li–Špakula–Zhang 的对象是 sparse coarse union 上的 block-rank-one projection

\[
P_n=\xi_n\xi_n^*,\qquad m_n(x)=|\xi_n(x)|^2.
\]

其 Proposition 4.8（本地 PDF pp. 15–16）把 \(P\) quasi-local 精确等价为 measured asymptotic expansion；Theorem 6.1（pp. 26–27）在 bounded geometry 下把该 rank-one \(P\) 的 quasi-locality 等价到 Roe membership；Corollary 4.21（pp. 22–23）的 \(1-\alpha_k\) 高测度 core 是 measured expander core。

它与本报告的

\[
\mu_\theta\{i:\|P_i(T-L_r)\|>\varepsilon\}\le\delta \tag{5.1}
\]

有四个不可交换之处：

* \(m_n(x)=|\xi_n(x)|^2\) 是 rank-one 振幅质量，未必是均匀节点比例；
* 高测度 core 是 expansion 的结构证人，不是“允许估计错误的节点集”；
* 结论控制整个 operator norm，而 (5.1) 是块行范数的分位数；
* 证明高度依赖 rank one，一般 WLS/Kalman 增益秩随规模增长。

因此 measured coarse geometry 对 \(\delta\)-node 问题目前是可借用的测度/随机根语言，不是直接答案。逐行分位数

\[
\inf_{S:\mu(S)\ge1-\delta}\sup_{i\in S}\|P_i(T-L_r)\|
\]

甚至不是通常的 \(C^*\)-operator norm；一般不具备把它直接变成 Roe ideal 所需的代数性质。

### 5.3 Ozawa 反例准确限制什么

Ozawa 的 Theorems A–B 与 Corollary C（预印本 pp. 1–2）证明：含 expander sequence 的 uniformly locally finite 空间上可有 quasi-local operator 不属于 uniform Roe algebra，因而不能在全网 \(\ell_2\to\ell_2\) operator norm 中由有限传播算子逼近。Špakula–Zhang Theorem 3.3 则说明 Property A 下 quasi-locality可推出 Roe membership。

这否定了“任意 bounded-degree 图族上，远隔集合作用趋零 iff finite-round 全网 operator approximation”的无条件版本。但 Ozawa 的算子未被证明是某个 uniformly sparse SPD Gaussian/WLS 的 \(J^{-1}G\)，也未研究删除 \(\delta\) 比例输出行后的 row quantile。因此它 **不能** 直接否定第 2 节的特殊 SPD theorem 或第 3 节的 M-matrix theorem；同样也不能被引用来声称多数节点版本已经解决。

## 6. 最可能可发表的两个 theorem 版本

### 版本 A：SPD/WLS 的多数节点 Dirichlet–oracle 等价

**建议正式命题。** 对一族 finite block-SPD normal systems \(J_\theta x=b\)：

* \(J\) 的非零块受通信图控制，\(b_j\) 是节点本地初始数据；
* radius-view 包含 principal block 与 cut coefficients；
* 目标是根状态块；节点测度任意；
* 性能采用 (1.8) 的 nodewise robust \(\ell_2\)-input 误差。

证明定理 3–4，并给出两种 corollary：

1. 统一谱窗时，所有 \(r,\delta\) 上有显式 \(\kappa\)-常数 sandwich (2.10)–(2.12)；
2. 无统一 condition number 时，用大多数节点的局部 \(\gamma_i\) 分位数与 oracle \(g_i\) 分位数做 \(\delta_1+\delta_2\) 组合，从而允许少量局部病态区域。

**证明路线已闭合：** block inverse \(\to\) 精确 factorization (2.4) \(\to\) fixed-model optimality (1.1) \(\to\) quantile monotonicity。无需树、junction tree 或每个 \(A_i\) 满列秩。

**真正需要新增的内容：**

* 给 \(\Gamma\) 常数的 sharpness 例子，或改进为角度/有效边界维数常数；
* 处理 raw WLS \(z\) 时，把 finite-propagation \(G\) 和 local right inverse 的半径平移写进 theorem；
* 把收集 view、形成 principal system、求解的 rounds/bits/flops 分开；
* 和 Cheng–Jiang–Sun finite-section theorem、Chung–Simpson local Dirichlet solver 做正面对比。

**撞车风险：中。** Schur identities 都是标准线性代数，单独不足以发论文；尚未找到 (2.10) 这种“固定模型 oracle 最优误差 vs 同一 only-local-model Dirichlet 规则、再取节点分位数”的现成 theorem。文章贡献应放在准确 minimax、sharp constants、completion lower bounds 和 resource implications 的组合上。

### 版本 B：grounded Laplacian 的 occupation criterion 与 completion dichotomy

**建议正式命题。** 对 \(J=L_W+\operatorname{diag}\kappa\)、输入 \(c=S^{-1}b\)、\(\ell_\infty\) 性能：

1. 证明 (3.5)–(3.8) 的 oracle/Dirichlet/Neumann 精确三分法；
2. 对任意节点测度，把三种误差的 \(1-\delta\) 分位数写成 killed-walk occupation/survival 的分位数；
3. 在 rootwise exterior lifetime envelope 下证明 escape probability 的双向证书 (3.12)，故按节点测度消失 iff；
4. 对允许 arbitrary-size ungrounded tentacles 的 completion-rich 类证明命题 8：任何非零 open boundary view 的 robust completion minimax 为无穷；
5. 给 \(\ell_2\) 的 two-walk intersection identity (3.16)，并明确它不能由单 walk escape probability 替代。

**证明路线已基本闭合：** Neumann path expansion + strong Markov property + electrical reciprocity。它给出一个相当干净的结论：最 relaxed 的固定模型条件不是无环、全局条件数或局部满秩，而是相应 Green occupation tail；从 fixed model 升级到 arbitrary completion 时，必须加入余寿命/tightness，且命题 8 说明这类假设在量级上不可省。

**关键未解引理/边界：**

* **[猜想/缺口]** 对 generic block-SPD，能否用 exterior Schur complement 的某种 capacity/transfer singular values 得到与 occupation theorem 同等尖锐、允许 path cancellation 的 completion dichotomy；简单 block-norm majorant 只充分。
* **[猜想/缺口]** 对不是“任意 completion”而是 Benjamini–Schramm/随机生成模型类，局部可观测的 escape/capacity 加何种 uniform integrability 恰好推出 majority Green-row locality；局部弱收敛本身在谱靠近 0 时不够。
* **[猜想/缺口]** 是否可将 rootwise lifetime envelope 替换成一个更局部、可抽样验证的 effective-resistance/capacity 条件，同时仍对所有允许 completion 成立。命题 8 表明若 completion 类仍可自由添加长 tentacle，答案必为否。

**撞车风险：中高（若只写随机游走表示），低到中（若包含 completion dichotomy + 节点分位数 minimax）。** Green/occupation/GFF 都是经典；论文必须把新意落在“哪一种局部算法误差对应哪一种 stopping-time functional”、only-local completion no-go 和多数节点 sharp theorem 上。

### 版本 A、B 都未自动解决计算与能量

上述 theorem 刻画的是信息半径。把 \(B_r(i)\) 完整汇集后做 dense Dirichlet solve，若球内状态维数为 \(m_i(r)\)，朴素存储/计算分别是 \(O(m_i(r)^2)\)、\(O(m_i(r)^3)\)；bounded degree 也可能有 \(m_i(r)\asymp(\Delta-1)^r\)。稀疏 Cholesky 的代价进一步依赖局部 treewidth。Neumann/message-passing 则可流式做 \(r\) 次 sparse matvec，通信轮数正好 \(r\)，但误差是较大的 \(n_i(r)\) 而非 \(o_i(r)\)。因此后续资源论文应以 (3.8) 作为“更强局部计算换更小误差”的第一条 Pareto 曲线，而不能从存在一个 \(r\)-local 映射直接推出低能耗。

## 7. 最终判断：哪些已结束，哪些还没结束

### 已经可以严谨宣布“结束”的部分

* 固定模型、固定输入范数、任意本地解码器：式 (1.1) 是最宽充要条件。
* only-local-model 的线性规则：式 (1.2) 是精确充要条件，但属于 minimax 定义展开；式 (2.8) 才是 SPD primitive-RHS 的结构化双向化。
* 多数节点的 nodewise robust 版本：式 (1.8) 是精确量词；不能与 per-input majority 或 Bayes failure 混用。
* grounded Laplacian、\(\ell_\infty\) 归一输入、固定模型：式 (3.5) 是最优误差的精确 potential-theoretic 条件。
* arbitrary-size completion 且无 stability/tightness：命题 7–8 给出非平凡局部证书通常不可能的 no-go。

### 还不能宣布“彻底终结”的部分

* 对所有 generic block Gaussian/WLS，在不加正锥、局部右逆或谱/tightness 条件时，把 completion minimax 化成一个纯图论、纯局部、易验证的 iff；反例表明这种 theorem 很可能根本不存在。
* 让坏节点只占 \(\delta\) 后，仅由“坏区节点数少”推出好节点的 Green tails 小；单锚点路径明确为假，正确条件必须涉及 capacity、可达性或 occupation tail。
* 同时最优 rounds、bits、energy、flops 的唯一 theorem；这些目标依赖通信与硬件模型，需另设 Pareto/下界问题。

## 8. 主来源、定理位置与本地路径

1. Stephen Demko, William F. Moss, Philip W. Smith, *Decay Rates for Inverses of Band Matrices*, Math. Comp. 43 (1984), 491–499. Proposition 2.1, Proposition 2.2, Theorem 2.4, pp. 492–494. DOI: [10.1090/S0025-5718-1984-0758197-9](https://doi.org/10.1090/S0025-5718-1984-0758197-9). 本地：literature/03_sparse_inverse_graph_filters/1984_demko_moss_smith_decay_inverse_band_matrices.pdf。

2. Dmitry M. Malioutov, Jason K. Johnson, Alan S. Willsky, *Walk-Sums and Belief Propagation in Gaussian Graphical Models*, JMLR 7 (2006), 2031–2064. Proposition 1 p. 2040; Proposition 21 p. 2049; necessity discussion pp. 2051–2052. [JMLR](https://jmlr.org/papers/v7/malioutov06a.html). 本地：literature/02_graphical_models_local_inference/2006_malioutov_johnson_willsky_walk_sums_gaussian_bp.pdf。

3. Cheng Cheng, Yingchun Jiang, Qiyu Sun, *Spatially Distributed Sampling and Reconstruction*, Appl. Comput. Harmon. Anal. 47 (2019), 109–148. Theorem 5.3 manuscript p. 14; Theorems 6.1–6.2 pp. 15–16; Proposition 7.1 pp. 17–18. DOI: [10.1016/j.acha.2017.07.007](https://doi.org/10.1016/j.acha.2017.07.007). 本地：literature/03_sparse_inverse_graph_filters/2015_cheng_jiang_sun_spatially_distributed_sampling_reconstruction.pdf。

4. Fan Chung, Olivia Simpson, *Solving Local Linear Systems with Boundary Conditions Using Heat Kernel PageRank*, Internet Mathematics 11 (2015), 449–471. DOI: [10.1080/15427951.2015.1009522](https://doi.org/10.1080/15427951.2015.1009522), [arXiv:1503.03157](https://arxiv.org/abs/1503.03157). 当前未见规范本地副本。

5. Christina Lee Yu, Asuman Ozdaglar, Devavrat Shah, *Asynchronous Approximation of a Single Component of the Solution to a Linear System*, IEEE Trans. Network Sci. Eng. 7(3) (2020), 975–986; preprint [arXiv:1411.2647](https://arxiv.org/abs/1411.2647). 其核心条件为 \(\rho(|G|)<1\)，并用 Neumann/walk expansion；当前未见规范本地副本。

6. Alexandr Andoni, Robert Krauthgamer, Yosef Pogrow, *On Solving Linear Systems in Sublinear Time*, ITCS 2019, Theorems 1.1–1.2；PSD lower-bound proof见 Proposition 3.3. DOI: [10.4230/LIPIcs.ITCS.2019.3](https://doi.org/10.4230/LIPIcs.ITCS.2019.3), [arXiv:1809.02995](https://arxiv.org/abs/1809.02995). 当前未见规范本地副本。

7. Enrique Bendito, Ángeles Carmona, Andrés M. Encinas, *Potential Theory for Schrödinger Operators on Finite Networks*, Rev. Mat. Iberoam. 21 (2005), 771–818. DOI: [10.4171/RMI/435](https://doi.org/10.4171/RMI/435). 当前未见规范本地副本。

8. Ángeles Carmona, Andrés M. Encinas, María José Jiménez, Àngela Martín, *Random Walks Associated with Symmetric M-matrices*, Linear Algebra Appl. 693 (2024), 324–338. DOI: [10.1016/j.laa.2023.10.009](https://doi.org/10.1016/j.laa.2023.10.009). 当前未见规范本地副本。

9. Scott Sheffield, *Gaussian Free Fields for Mathematicians*, Probab. Theory Relat. Fields 139 (2007), 521–541. DOI: [10.1007/s00440-006-0050-1](https://doi.org/10.1007/s00440-006-0050-1), [arXiv:math/0312099](https://arxiv.org/abs/math/0312099). 当前未见规范本地副本。

10. Michael L. Stein, *2010 Rietz Lecture: When Does the Screening Effect Hold?*, Ann. Statist. 39(6) (2011), 2795–2819. Conjecture 1 and Theorems 1–2, article pp. 2807–2809 附近（arXiv HTML §2–3）. DOI: [10.1214/11-AOS909](https://doi.org/10.1214/11-AOS909), [arXiv:1203.1801](https://arxiv.org/abs/1203.1801). 当前未见规范本地副本。

11. Michael L. Stein, *When Does the Screening Effect Not Hold?*, Spatial Statistics 11 (2015), 65–80. DOI: [10.1016/j.spasta.2014.12.003](https://doi.org/10.1016/j.spasta.2014.12.003). 当前未见规范本地副本。

12. Mohammad Meysami, Ali Lotfi, *Screening Effects in Gaussian Random Fields under Generalized Spectral Conditions*, Results in Applied Mathematics 29 (2026), 100681. DOI: [10.1016/j.rinam.2025.100681](https://doi.org/10.1016/j.rinam.2025.100681). 当前未见规范本地副本。

13. Mohammad Pirani, Shreyas Sundaram, *On the Smallest Eigenvalue of Grounded Laplacian Matrices*, IEEE Trans. Automat. Control 61(2) (2016), 509–514. Theorem 1 及 random graph corollaries. DOI: [10.1109/TAC.2015.2444191](https://doi.org/10.1109/TAC.2015.2444191), [arXiv:1406.2271](https://arxiv.org/abs/1406.2271). 当前未见规范本地副本。

14. Ján Špakula, Jiawen Zhang, *Quasi-locality and Property A*, J. Funct. Anal. 278(1) (2020), 108299. Theorem 3.3, PDF pp. 9–13. DOI: [10.1016/j.jfa.2019.108299](https://doi.org/10.1016/j.jfa.2019.108299). 本地：literature/03_sparse_inverse_graph_filters/2019_spakula_zhang_quasi_locality_property_a.pdf。

15. Kang Li, Ján Špakula, Jiawen Zhang, *Measured Asymptotic Expanders and Rigidity for Roe Algebras*, IMRN 2023(17), 15102–15154. Proposition 4.8, Corollary 4.21, Theorem 6.1. DOI: [10.1093/imrn/rnac242](https://doi.org/10.1093/imrn/rnac242), [arXiv:2010.10749](https://arxiv.org/abs/2010.10749). 本地：literature/03_sparse_inverse_graph_filters/2020_li_spakula_zhang_measured_asymptotic_expanders_roe.pdf。

16. Narutaka Ozawa, *Embeddings of Matrix Algebras into Uniform Roe Algebras and Quasi-local Algebras*, Theorems A–B, Corollary C, pp. 1–2. 截至核验日为预印本：[arXiv:2310.03677](https://arxiv.org/abs/2310.03677). 本地：literature/03_sparse_inverse_graph_filters/2023_ozawa_matrix_embeddings_quasi_local_not_uniform_roe.pdf。

## 9. 给主项目的直接建议

不要以“找一个比无环+局部满秩更弱的纯拓扑条件”为主问题；一般情形不存在这样的单一条件。建议把论文对象固定为下列二选一：

1. **block-SPD primitive sufficient-statistic 模型：** 主打定理 3–4、completion ambiguity 与多数节点 quantile；再用 raw WLS 的 local-\(G\)/local-right-inverse corollary 接回传感器观测。
2. **标量 grounded Laplacian/M-matrix 模型：** 主打定理 6、escape/lifetime 双向证书和 arbitrary-completion tentacle no-go。这条线更有直观图论/概率论解释，也最接近“允许少量节点失败但不要求全局 condition number”。

无论选哪条，都应先在正文第一页声明：输入范数、好节点是否对所有输入统一、模型是否固定还是对 completion robust、节点是否知道全局谱常数。否则所谓“最 relaxed iff”会因量词变化得到互相矛盾但都正确的答案。
