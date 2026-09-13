# 局部状态估计定理包 v1

> 数学核验日期：2026-09-13  
> 目的：把“有限轮局部通信何时能逼近中心 Gaussian/WLS 估计”写成可逐行审计的定义、定理、证明与反例。本文不以新术语重命名经典结果。

## 0. 结论状态与使用边界

本文使用四种状态标记：

- **[标准恒等式]**：经典线性代数、条件期望或 Markov 链恒等式；不可主张原创。
- **[直接推论]**：由定义或标准恒等式数行推出；即使未找到相同表述，也不能仅凭此主张原创。
- **[可能新的组合]**：各组成部分经典，但“共同局部规则、completion、节点分位数”这一组合在第二阶段检索中未发现完全相同的主定理；仍需正式 novelty search。
- **[仍有缺口]**：当前证明不能在所声称的一般性下闭合。

本文闭合的核心结论是：

1. 固定模型时，任意确定性局部解码器的最小最坏误差就是中心算子的球外块行范数。
2. 只知道局部模型时，线性算法的精确对象是“同一局部视图必须共用一个系数”的 minimax，而不是单个全局模型的行尾。
3. 对原始正规方程右端 $Jx=b$，Schur 补把最优球外行尾与一个真正只用局部系数的 principal/Dirichlet 解连接起来；对多数节点可得到双向常数夹逼。
4. 对 grounded Laplacian/M-matrix，最优截断、Dirichlet 解和有限 Neumann 迭代分别等于三种 killed-walk 路径泛函。
5. 若允许同一局部视图外接任意长、近奇异的 completion，则不存在不附带稳定性或余寿命控制的非平凡统一局部保证。

这里的“充要条件”始终相对于明确写出的模型类、输入范数、信息视图和失败语义。不存在脱离这些量词的单一“最宽拓扑条件”。

---

## 1. 统一模型与局部信息

### 1.1 网络、模型与数据

一个有限网络模型记为

\[
\theta=(G_\theta,\{\mathsf X_u\},\{\mathsf Z_a\},o_\theta,\lambda_\theta,T_\theta).
\]

其中：

- $G_\theta=(V_\theta,E_\theta)$ 是通信图；
- 节点 $u$ 的状态空间是有限维实或复 Hilbert 空间 $\mathsf X_u$；
- 数据块索引为 $a$，数据空间为 $\mathsf Z_a$，其所有者 $o_\theta(a)\in V_\theta$；
- $\lambda_\theta$ 是初始可见的模型标签，例如块维数、端口、边权、测量矩阵块和噪声精度块；
- 中心线性估计算子为

  \[
  T_\theta:\mathsf Z_\theta:=\bigoplus_a\mathsf Z_a
  \longrightarrow
  \mathsf Y_\theta:=\bigoplus_{u\in V_\theta}\mathsf Y_u.
  \]

根 $i$ 的中心目标行写成 $T_{\theta,i}=P_iT_\theta:\mathsf Z_\theta\to\mathsf Y_i$。模型 $\theta$ **不包含本次实现的数据值** $z$；数据另行归一化为 $\|z\|_2\le1$，或在第 6 节明确改用 $\ell_\infty$ 范数。

把数据值混进模型，再对无范数约束的远端数据取 supremum，会把任何非零远端系数机械地变成无穷误差；本文排除这种无尺度的表述。

### 1.2 $r$ 轮视图与局部规则

令 $B_r(i)$ 是通信图中的闭球。半径 $r$ 的数据限制算子

\[
S_{\theta,i,r}:\mathsf Z_\theta\to
\mathsf U_{\theta,i,r}:=
\bigoplus_{a:o_\theta(a)\in B_r(i)}\mathsf Z_a
\tag{1.1}
\]

是坐标抽取，故 $S S^*=I_{\mathsf U}$，而 $P:=S^*S$ 是全局数据空间上的正交投影。记 $Q:=I-P$。

根标记**模型**视图 $\mathcal V_r(i,\theta)$ 包括：根、半径 $r$ 球的图结构、端口、块空间、节点初始模型标签和数据所有权标签，但不包括本次实现的数据值。算法在运行时看到的是二元组
\((\mathcal V_r(i,\theta),S_{\theta,i,r}z)\)。两个模型视图同构时，用预先固定的等距映射识别其本地输入和目标空间。不同 completion 的球外空间不必彼此识别；每个算子范数先在各自空间计算，再取 supremum。

本文采用确定性 LOCAL 信息模型，并允许消息大小与本地计算不受限。因此 $r$ 轮后根可以重建完整的标记 $r$-view；反过来，任何只依赖该 view 与 $Sz$ 的函数都可视为一个 $r$ 轮信息规则。这个约定只刻画轮数/信息半径，不自动给出低 bit complexity 或低 flops。

### 1.3 闭合系数视图与 halo

Schur 定理需要一个比“诱导子图”稍强、但仍局部的信息约定。对 $B=B_r(i)$，称 view 对正规矩阵 $J$ **系数闭合**，若它确定

\[
A=J_{BB},\qquad E=J_{BO},
\tag{1.2}
\]

其中 $E$ 只需按边界端口识别其非零列。一个标准实现是：每个节点初始知道自己 incident 的全部 $J$-块，边界节点可在 $r$ 轮内把这些块和端口报告给根。若测量因子的所有权使边界节点不知道跨边界因子，必须显式多收集一圈 halo；后文所有“$r$ 轮”结论都允许这一固定半径平移。

只保留球内测量因子得到的矩阵 $\widetilde A$ 通常不等于 principal block $A=J_{BB}$。本文的 Dirichlet 定理使用 $A$，不把二者偷换。

### 1.4 三种互不等价的“允许 $\delta$ 失败”

每个模型给定节点概率测度 $\mu_\theta$；均匀节点比例是特例。对有限节点上的非负函数 $e_i$，定义

\[
Q_{1-\delta}^{\mu_\theta}(e)
:=\inf\{t:\mu_\theta\{i:e_i\le t\}\ge1-\delta\},
\qquad 0\le\delta<1.
\tag{1.3}
\]

本小节先用线性规则写最紧凑的公式。设线性局部规则 $L$ 在根 $i$ 的误差算子为

\[
E_{\theta,i}(L)=T_{\theta,i}-L_{\mathcal V_r(i,\theta)}S_{\theta,i,r}.
\]

下列三种语义必须分开。

若规则 \(f\) 非线性，就把下文的行算子范数替换为
\[
e_{\theta,i}(f)
:=\sup_{\|z\|_2\le1}
\|T_{\theta,i}z-f_{\mathcal V_r(i,\theta)}(S_{\theta,i,r}z)\|_2.
\]
三种失败语义的量词顺序不变。

**(R) 节点鲁棒多数。** 存在一个与数据无关的好节点集，测度至少 $1-\delta$，且每个好节点对所有 $\|z\|_2\le1$ 都满足误差界。其数值是

\[
\mathcal A^{\rm R}_{r,\delta}
=\inf_L\sup_{\theta}
Q_{1-\delta}^{\mu_\theta}
\bigl(\|E_{\theta,i}(L)\|_{2\to2}\bigr).
\tag{1.4}
\]

**(I) 每个输入的多数。** 好节点集可以依赖本次输入：

\[
\mathcal A^{\rm I}_{r,\delta}
=\inf_L\sup_{\theta}\sup_{\|z\|_2\le1}
Q_{1-\delta}^{\mu_\theta}
\bigl(\|E_{\theta,i}(L)z\|_2\bigr).
\tag{1.5}
\]

总有 $\mathcal A^{\rm I}_{r,\delta}\le\mathcal A^{\rm R}_{r,\delta}$，反向一般为假。例如误差算子是 $I_n$ 时，每行算子范数均为 1；但对任意单位 $\ell_2$ 输入，超过 $\varepsilon$ 的坐标至多 $1/\varepsilon^2$ 个。

**(B) 随机数据/Bayes 多数。** 给定数据律 $\mathbb P_\theta$，一种明确语义是联合失败概率

\[
\mathcal F^{\rm B}_{r}(L,\varepsilon)
=\sup_\theta
(\mu_\theta\otimes\mathbb P_\theta)
\{(i,Z):\|E_{\theta,i}(L)Z\|_2>\varepsilon\}.
\tag{1.6}
\]

也可另定义“节点 MSE 的 $1-\delta$ 分位数”，但必须单独写出；它不等于 (1.6)。第 5.3 节给出 Gaussian MSE 的精确公式。除非特别声明，本文的多数节点定理一律采用最强的 (R)。

---

## 2. 固定模型：任意解码器的精确充要条件

### 背景引理 2.1（固定模型的最小不可避免误差）

**状态：[直接推论]。** 这是坐标信息结构的精确恒等式，不是图结构新定理。

固定 $(\theta,i,r)$，令 $T=T_{\theta,i}$、$S=S_{\theta,i,r}$、$P=S^*S$、$Q=I-P$。允许 $f:\mathsf U\to\mathsf Y_i$ 为任意确定性函数，则

\[
\boxed{
\inf_f\sup_{\|z\|_2\le1}\|Tz-f(Sz)\|_2
=\|TQ\|_{2\to2}.}
\tag{2.1}
\]

线性规则 $f_0(u)=TS^*u$ 达到该值。

#### 证明

对任意 $z$，

\[
Tz-f_0(Sz)=Tz-TS^*Sz=TQz,
\]

所以上界为 $\|TQ\|$。

反过来，取任意 $q\in\ker S$、$\|q\|_2\le1$。输入 $q$ 和 $-q$ 的本地数据都等于 0。令 $a=f(0)$、$v=Tq$。三角不等式给出

\[
2\|v\|=\|v-(-v)\|
\le\|v-a\|+\|-v-a\|,
\]

故两输入中至少一个的误差不小于 $\|Tq\|$。对 $q$ 取 supremum 得下界 $\|TQ\|$。有限维时范数达到；不依赖达到性时取逼近极值序列即可。证毕。

#### 充要条件

- 精确 $r$ 轮恢复当且仅当 $TQ=0$，即目标行完全因子化为 $T=(TS^*)S$。
- 最坏误差不超过 $\varepsilon$ 当且仅当 $\|TQ\|\le\varepsilon$。
- 图可以有任意环；无环既非必要，也没有出现在条件中。

同一证明适用于任何使球内/球外坐标投影为收缩的 $\ell_p$ 直和范数。标量输出、输入为 $\ell_\infty$ 时，右侧就是球外系数的 $\ell_1$ 行和。

---

## 3. 只知道局部模型：共同线性系数 minimax

固定一个同构视图 $v$。其所有 rooted completion 组成集合

\[
\mathfrak C(v)=\{\alpha=(\theta,i):\mathcal V_r(i,\theta)\cong v\}.
\]

经同构识别后，它们有共同本地输入空间 $\mathsf U_v$ 和输出空间 $\mathsf Y_v$，但全局数据空间 $\mathsf Z_\alpha$ 可以不同。记 $S_\alpha:\mathsf Z_\alpha\to\mathsf U_v$、$T_\alpha:\mathsf Z_\alpha\to\mathsf Y_v$。定义

\[
a_r^{\rm lin}(v)
:=\inf_{\ell\in\mathcal L(\mathsf U_v,\mathsf Y_v)}
\sup_{\alpha\in\mathfrak C(v)}
\|T_\alpha-\ell S_\alpha\|_{2\to2}.
\tag{3.1}
\]

### 定理 3.1（共同局部线性规则的 iff）

**状态：[直接推论]。** 这是 LOCAL 一致性量词的展开；式 (3.1) 本身不能作为结构创新。

给定 $\varepsilon\ge0$，view $v$ 上存在同一个线性 $r$ 轮规则，对每个 completion 和每个单位输入误差均不超过 $\varepsilon$，当且仅当 (3.1) 的 infimum 由某个 $\ell$ 达到且值不超过 $\varepsilon$。不要求达到时，严格结论是：对每个 $\eta>0$，存在共同规则达到 $a_r^{\rm lin}(v)+\eta$。

若 $\mathsf U_v,\mathsf Y_v$ 有限维、$\mathfrak C(v)\ne\varnothing$，且

\[
\sup_{\alpha\in\mathfrak C(v)}\|T_\alpha\|<\infty,
\tag{3.2}
\]

则 infimum 必定达到。

#### 证明

在同构 view 上，确定性线性 LOCAL 规则必须采用同一个 $\ell$；反之任何 $\ell$ 都定义一个该 view 上的规则。因此误差陈述与 (3.1) 逐字等价。

只需证明达到性。令

\[
F(\ell)=\sup_\alpha\|T_\alpha-\ell S_\alpha\|.
\]

因为每个 $S_\alpha$ 是 coisometry，$\|\ell S_\alpha\|=\|\ell\|$。固定任意 $\alpha_0$，有

\[
F(\ell)\ge\|T_{\alpha_0}-\ell S_{\alpha_0}\|
\ge\|\ell\|-\|T_{\alpha_0}\|,
\]

所以 $F$ coercive。又

\[
|F(\ell)-F(m)|\le\|\ell-m\|,
\]

且由 (3.2) 知 $F(0)<\infty$。有限维连续 coercive 函数达到最小值。证毕。

### 引理 3.2（远端尾与本地系数歧义）

**状态：[直接推论]。**

按 $\mathsf Z_\alpha=\mathsf U_v\oplus\mathsf W_\alpha$ 写

\[
T_\alpha=[A_\alpha,B_\alpha],\qquad S_\alpha=[I,0].
\]

则

\[
a_r^{\rm lin}(v)=
\inf_\ell\sup_\alpha
\|[A_\alpha-\ell,B_\alpha]\|,
\tag{3.3}
\]

并且

\[
a_r^{\rm lin}(v)\ge
\max\left\{
\sup_\alpha\|B_\alpha\|,
\frac12\sup_{\alpha,\beta}\|A_\alpha-A_\beta\|
\right\}.
\tag{3.4}
\]

若 $A_\alpha=A_v$ 对所有 completion 相同，则

\[
a_r^{\rm lin}(v)=\sup_\alpha\|B_\alpha\|.
\tag{3.5}
\]

#### 证明

式 (3.3) 是分块。把输入限制在 $\mathsf W_\alpha$ 得第一项下界。对任意 $\ell$，

\[
\|A_\alpha-A_\beta\|
\le\|A_\alpha-\ell\|+\|A_\beta-\ell\|,
\]

得到第二项。若所有本地块相同，取 $\ell=A_v$ 给出与第一项下界相同的上界。证毕。

### 备注 3.3（任意非线性 completion 规则）

对跨 completion 的任意解码器，(3.1) 不再是精确 iff。固定本地数据 $u$，令

\[
\mathcal C_v(u)=
\{T_\alpha z:\alpha\in\mathfrak C(v),\ S_\alpha z=u,\ \|z\|\le1\}.
\]

忽略可测选择问题时，精确 minimax 是

\[
\sup_{\|u\|\le1}\operatorname{rad}\mathcal C_v(u),
\qquad
\operatorname{rad}C:=\inf_y\sup_{x\in C}\|x-y\|.
\tag{3.6}
\]

其最优中心可以随 $u$ 非线性变化。因此后文 completion sandwich 明确针对线性规则；固定模型的引理 2.1 则已经允许任意非线性规则。

---

## 4. block-SPD 原始正规方程右端：Schur 精确分解

本节先研究

\[
Jx=b,\qquad J=J^*\succ0,
\tag{4.1}
\]

并把 $b=\bigoplus_{u\in V}b_u$ 当作节点本地的 primitive sufficient statistic。原始传感器观测 $z$ 与 $b=Gz$ 的差别留到第 5 节。

### 4.1 分块记号

固定根 $i$ 和包含它的球 $B$，令 $O=V\setminus B$。按
$\mathsf X=\mathsf X_B\oplus\mathsf X_O$ 写

\[
J=
\begin{bmatrix}
A&E\\ E^*&D
\end{bmatrix},
\qquad
\Sigma:=D-E^*A^{-1}E.
\tag{4.2}
\]

这里 $A\succ0$ 且 Schur 补 $\Sigma\succ0$。令

\[
R_i:\mathsf X_B\to\mathsf X_i
\]

为根块抽取，并定义

\[
C:=R_iA^{-1}E:\mathsf X_O\to\mathsf X_i,
\qquad
F:=E^*A^{-1}:\mathsf X_B\to\mathsf X_O.
\tag{4.3}
\]

固定模型、只观察 $b_B$ 时的 oracle 尾和 principal/Dirichlet 误差分别是

\[
g_i(B):=\|[R_i,0]J^{-1}[0,I_O]^*\|
\tag{4.4}
\]

和

\[
d_i(B):=
\|[R_i,0]J^{-1}-[R_iA^{-1},0]\|.
\tag{4.5}
\]

所有范数均为 Hilbert 空间的 $2\to2$ 算子范数。

### 引理 4.1（Schur 逆公式）

**状态：[标准恒等式]。**

\[
J^{-1}=
\begin{bmatrix}
A^{-1}+A^{-1}E\Sigma^{-1}E^*A^{-1}
&
-A^{-1}E\Sigma^{-1}\\[2mm]
-\Sigma^{-1}E^*A^{-1}
&
\Sigma^{-1}
\end{bmatrix}.
\tag{4.6}
\]

#### 证明

作块消元

\[
J=
\begin{bmatrix}I&0\\E^*A^{-1}&I\end{bmatrix}
\begin{bmatrix}A&0\\0&\Sigma\end{bmatrix}
\begin{bmatrix}I&A^{-1}E\\0&I\end{bmatrix},
\]

逐项求逆并相乘即得 (4.6)。直接乘回 $J$ 也可验证。证毕。

### 背景引理 4.2（根行误差的精确 Schur 因子分解）

**状态：[标准恒等式 + 直接推论]。** 因子分解来自经典块逆；oracle/Dirichlet 解释是直接推论。

令

\[
Y:=[R_i,0]J^{-1}[0,I_O]^*
=-C\Sigma^{-1}.
\tag{4.7}
\]

则

\[
\boxed{
[R_i,0]J^{-1}-[R_iA^{-1},0]
=Y[-F,I_O].}
\tag{4.8}
\]

因此

\[
\boxed{
g_i(B)\le d_i(B)\le
\gamma_i(B)g_i(B),\qquad
\gamma_i(B):=\sqrt{1+\|F\|^2}.}
\tag{4.9}
\]

#### 证明

由 (4.6)，根行的球内块与球外块是

\[
\bigl[
R_iA^{-1}+C\Sigma^{-1}F,\;
-C\Sigma^{-1}
\bigr].
\]

减去 $[R_iA^{-1},0]$，并用 $Y=-C\Sigma^{-1}$，得到

\[
[C\Sigma^{-1}F,-C\Sigma^{-1}]
=[-YF,Y]=Y[-F,I_O],
\]

即 (4.8)。

把 (4.8) 的输入限制为 $[0,b_O]$，得到 $d_i(B)\ge\|Y\|=g_i(B)$。另一方面，

\[
d_i(B)\le\|Y\|\,\|[-F,I_O]\|.
\]

而

\[
[-F,I_O][-F,I_O]^*=FF^*+I_O,
\]

故 $\|[-F,I_O]\|=\sqrt{1+\|F\|^2}$。证毕。

### 命题 4.3（常数 $\gamma$ 的精确地位）

**状态：[直接推论]。** 这是对上界常数的 sharpness 审计，不声称文献原创。

对一个固定 completion，若 $g_i(B)>0$，真正的比值是

\[
\gamma_{\rm eff}(J;i,B)
:=\frac{\|C\Sigma^{-1}[-F,I_O]\|}
{\|C\Sigma^{-1}\|}
\le\sqrt{1+\|F\|^2}.
\tag{4.10}
\]

它可能严格小于 $\gamma_i(B)$，但依赖未知的 exterior Schur 补 $\Sigma$。

若固定局部 $(A,E,R_i)$、$C\ne0$，并允许所有正定 Schur 补
$\Sigma\succ0$ 作为 completion，则

\[
\sup_{\Sigma\succ0}\gamma_{\rm eff}
=\sqrt{1+\|F\|^2}.
\tag{4.11}
\]

所以在不限制 exterior 的模型类中，(4.9) 的 $\gamma$ 是
**Dirichlet 误差相对 fixed-model oracle 尾**的最小局部统一乘法常数；一般不能无条件改小。这里没有声称它也是所有共同局部算法相对
\(\bar g\) 的 minimax sharp factor；命题 4.5 明确展示共同系数可优于 Dirichlet。

#### 证明

上界已由引理 4.2 证明。只证反向。取单位向量 $v_k\in\mathsf X_O$，使

\[
\|F^*v_k\|\to\|F\|,
\qquad Cv_k\ne0.
\]

这样的序列总能选到：若某个最大奇异向量落在 $\ker C$，用不在
$\ker C$ 的向量作任意小扰动。对固定 $v=v_k$，令

\[
M_t=I+t\,vv^*\succ0,\qquad \Sigma_t=M_t^{-1}.
\]

则

\[
CM_t=C+t(Cv)v^*.
\]

除以 $t\|Cv\|$ 并令 $t\to\infty$，有

\[
\frac{\|CM_t[-F,I]\|}{\|CM_t\|}
\longrightarrow
\|v^*[-F,I]\|
=\sqrt{1+\|F^*v\|^2}.
\]

再令 $k\to\infty$ 得 (4.11)。每个 $\Sigma_t$ 都可通过
$D=E^*A^{-1}E+\Sigma_t$ 实现为一个 SPD completion。若 $C=0$，则所有 completion 都有 $g_i=d_i=0$，不存在 sharpness 问题。证毕。

若 completion 另受 $mI\preceq\Sigma\preceq MI$ 等限制，精确的最小常数是

\[
\sup_{\Sigma\ {\rm allowed}}
\frac{\|C\Sigma^{-1}[-F,I]\|}{\|C\Sigma^{-1}\|},
\tag{4.12}
\]

它可小于 $\gamma$。在一般 block 情形，尚未得到一个仅用少数局部标量就闭式计算 (4.12) 的公式；把 (4.10) 冒充局部证书会循环使用未知 exterior。

### 4.2 completion minimax sandwich

考虑共享同一个系数闭合 view $v$ 的 SPD completions。技术上，若外部维数变化，view 给出一个固定边界端口空间 $\mathsf C_v$ 和
$E_v:\mathsf C_v\to\mathsf X_B$；每个 completion 把 $\mathsf C_v$ 等距嵌入其
$\mathsf X_O$，并令 $E$ 在其余外部坐标上为零。因此 $A$、非零 cut columns 和 $\|F\|$ 都由 view 唯一确定。

令 $a_{r,{\rm SPD}}^{\rm lin}(v)$ 是 (3.1) 对映射
$b\mapsto R_ix$ 的值，并令

\[
\bar g(v):=\sup_{\alpha\in\mathfrak C(v)}g_\alpha(B).
\]

### 定理 4.4（共同 Dirichlet 规则与 completion 尾）

**状态：[可能新的组合]。** 两边证明只用标准 Schur 公式和固定模型信息下界；“only-local completion minimax”的组合表述尚未发现直接先例，但不能把组成恒等式报为原创。

\[
\boxed{
\bar g(v)
\le a_{r,{\rm SPD}}^{\rm any}(v)
\le a_{r,{\rm SPD}}^{\rm lin}(v)
\le\gamma(v)\bar g(v),}
\tag{4.13}
\]

其中 $a^{\rm any}$ 允许任意确定性共同解码器，
$\gamma(v)=\sqrt{1+\|E^*A^{-1}\|^2}$。最后一个上界由同一个 view-to-rule 映射

\[
\ell_v^D(b_B)=R_iA^{-1}b_B
\tag{4.14}
\]

达到。

#### 证明

固定任一 completion。即使解码器离线知道完整 completion，引理 2.1 也说明其误差至少是球外行尾 $g_\alpha(B)$。只知道共同 view 不会更容易，所以任意共同解码器的最坏误差至少为
$\sup_\alpha g_\alpha(B)=\bar g(v)$。

任意规则的 infimum 不大于线性规则的 infimum。最后，$A$ 由 view 决定，故 (4.14) 对所有 completion 是同一个合法规则；引理 4.2 对每个 completion 给出误差
$d_\alpha(B)\le\gamma(v)g_\alpha(B)$。取 supremum 得最右上界。证毕。

这条定理还说明：在 primitive RHS 的 SPD 类中，completion 改变本地最优系数的歧义不会与外部行尾完全脱钩；共同 Dirichlet 规则把二者控制在局部因子 $\gamma$ 内。

### 命题 4.5（单输出、单割边、已知 Schur 谱窗的精确 robust 系数）

**状态：[精确受约束标量特例；一般 robust approximate-inverse 框架已知；该闭式仍待专项查重]。**

El Ghaoui (2002) 已把 uncertain-matrix 的共同 approximate inverse 定义为

\[
\inf_X\sup_{\Delta\in\Delta_\rho}
\|A(\Delta)^{-1}-X\|_2
\tag{4.14-ELG}
\]

（原文 Eq. (1.4) 另乘常数 \(1/\rho\)），并在 Theorem 6.2/Eq. (6.3)
对 structured LFR uncertainty 给 SDP 上界；unstructured
\(\Delta=\mathbb R^{p\times q}\) 时该界精确，Eq. (6.5) 给出解析
approximate inverse。因此，“为不确定逆矩阵族选一个共同中心”以及
“用 SDP 处理 structured uncertainty”已经直接撞车，不能作为本文创新。

本命题的目标比 (4.14-ELG) 多一个通信诱导的仿射可行集。只允许根修改
**可见列**，不可见列必须为 0：

\[
\inf_q\sup_s
\left\|
\underbrace{[R_iA^{-1}+csF,\,-cs]}_{\text{完整根行}}
-
\underbrace{[R_iA^{-1}+q,\,0]}_{\text{局部规则}}
\right\|_2.
\tag{4.14-LOC}
\]

El Ghaoui 的 \(X\) 是逼近整个 inverse 的自由满矩阵；Eq. (6.5) 也只针对
unstructured perturbation 与自由 \(X\)。它没有直接给出
(4.14-LOC) 中“选定根行 + 不可见列强制为 0”的分段公式。若解除这个
局部列约束，当前标量不确定集只是一条线段，其自由 Chebyshev center
就是两端点中点，已完全落入经典 approximate-inverse 思路。

所以，下面的计算至多保留为一个**局部可行集下的 exact scalar reduction**；
它不是 robust approximate inverse 新框架。当前核验未发现 Eq. (6.5)
直接包含该分段式，但尚不能据此确认公式原创。

设根输出为标量、外部接口为一维。于是

\[
C=c\in\mathbb R,\qquad F=f\,u^*,\qquad \|u\|=1,\quad f=\|F\|\ge0.
\]

假设允许的 exterior Schur 补是所有标量

\[
\Sigma\in[m,M],\qquad 0<m<M,
\]

并记

\[
\alpha:=1/M,\qquad \beta:=1/m.
\]

完整根行的球内系数为 \(R_iA^{-1}+csF\)，球外系数为
\(-cs\)，其中 \(s=\Sigma^{-1}\in[\alpha,\beta]\)。共同局部规则写成
\(\ell=R_iA^{-1}+q\)。则精确 robust minimax 为

\[
a_{\rm sc}
=\inf_q\sup_{s\in[\alpha,\beta]}
\sqrt{\|csF-q\|^2+c^2s^2}.
\tag{4.14a}
\]

若 \(c=0\)，则 \(a_{\rm sc}=0\)。若 \(c\ne0\)，最优 \(q\) 平行于
\(F\)，并有：

\[
\boxed{
a_{\rm sc}=
\begin{cases}
|c|\beta,
&
f^2\le\dfrac{\alpha+\beta}{\beta-\alpha},\\[3mm]
|c|\sqrt{
\beta^2+
\dfrac{\bigl(f^2(\beta-\alpha)-(\alpha+\beta)\bigr)^2}{4f^2}},
&
f^2\ge\dfrac{\alpha+\beta}{\beta-\alpha}.
\end{cases}}
\tag{4.14b}
\]

第一种情形可取

\[
q^*=c\beta F,
\tag{4.14c}
\]

第二种情形可取

\[
q^*=
c\,\frac{(f^2+1)(\alpha+\beta)}{2f^2}\,F.
\tag{4.14d}
\]

特别地，completion 中最坏 oracle 尾是
\(\bar g=|c|\beta\)。当
\(f^2\le(\alpha+\beta)/(\beta-\alpha)\) 时，尽管 completion 会改变本地系数，共同规则仍精确达到 \(\bar g\)；而 Dirichlet 规则 \(q=0\) 的最坏误差是
\(\sqrt{1+f^2}\,\bar g\)。

#### 证明

把 \(q\) 分解成 \(\operatorname{span}\{u^*\}\) 内外两部分。正交部分在每个
\(s\) 的平方误差中只增加同一个非负量，所以最优解必有
\(q=\rho u^*\)。除以 \(|c|\) 并吸收符号，问题化成

\[
\inf_{\nu\in\mathbb R}
\max_{s\in[\alpha,\beta]}
\bigl[(fs-\nu)^2+s^2\bigr].
\tag{4.14e}
\]

对固定 \(\nu\)，括号是 \(s\) 的严格凸二次函数，所以区间最大值在
\(\alpha,\beta\) 两端之一。

若取 \(\nu=f\beta\)，端点 \(\beta\) 的值为 \(\beta^2\)，端点
\(\alpha\) 的值为
\[
f^2(\beta-\alpha)^2+\alpha^2.
\]
它不超过 \(\beta^2\) 当且仅当
\[
f^2\le\frac{\alpha+\beta}{\beta-\alpha}.
\]
由于任何规则在 \(s=\beta\) 都不能消除球外项，minimax 至少为
\(\beta^2\)，故第一种情形最优。

在相反情形，最优点必须让两个端点损失相等。解

\[
(f\alpha-\nu)^2+\alpha^2
=(f\beta-\nu)^2+\beta^2
\]

得到

\[
\nu^*=\frac{(f^2+1)(\alpha+\beta)}{2f}.
\]

此时
\[
f\alpha<\nu^*\le f\beta.
\]
两个 active convex quadratics 对 \(\nu\) 的导数分别非负、非正，故 0
属于其 active subgradient 的凸包；因此该交点确为 convex maximum
的全局最小点，而不只是一个 equalizer。

代回 \(\beta\) 端点即得 (4.14b) 第二行；由
\(q=\rho u^*\)、\(F=fu^*\) 得 (4.14d)。证毕。

若 \(m=M\)，则 completion 实际唯一，取
\(q=c(1/m)F\) 后误差恰为 \(|c|/m\)，即 oracle 尾。

对一般 block cut，robust common coefficient 可精确写成凸问题

\[
\inf_Q\sup_{mI\preceq\Sigma\preceq MI}
\|[C\Sigma^{-1}F-Q,\,-C\Sigma^{-1}]\|.
\tag{4.14f}
\]

对有限个候选 \(\Sigma\)，它是标准 operator-norm epigraph SDP；连续谱窗给出 semi-infinite robust problem。El Ghaoui 的 LFR/S-procedure/SDP 已经提供处理一般 uncertain inverse 的已知框架；把 affine locality constraints 加到 \(X\) 上仍保持凸性，并可直接得到同类 structured-SDP **上界**。因此，“把 (4.14f) 写成 robust minimax/SDP”本身不再是创新候选。

尚未由该文直接解决的窄问题是：对 Schur-complement uncertainty 和
\([Q,0]\) 这种局部列约束，何时 structured SDP 上界精确；非交换谱窗的最坏
\(\Sigma\) 是什么；能否得到类似 (4.14b) 的 exact closed form。El Ghaoui
明确只在 unstructured uncertainty 下声称其条件必要且 Eq. (6.5) 精确，
不能把 Theorem 6.2 的 structured SDP 上界误报成 (4.14f) 的精确解。

### 4.3 多数节点 sandwich

对模型 $\theta$ 与根 $i$，取 $B=B_r(i)$，定义 $g_{\theta,i}(r)$、
$d_{\theta,i}(r)$、$\gamma_{\theta,i}(r)$ 如上。令

\[
\mathcal G_{r,\delta}
:=\sup_\theta Q_{1-\delta}^{\mu_\theta}
\bigl(g_{\theta,i}(r)\bigr),
\tag{4.15}
\]

并令 $\mathcal A^{\rm any,R}_{r,\delta}$、$\mathcal A^{\rm lin,R}_{r,\delta}$
分别是所有共同 view-based 确定性规则与线性规则在节点鲁棒语义 (R) 下的最优值。

### 定理 4.6（多数节点的 oracle–local sandwich）

**状态：[可能新的组合]。** 逐节点不等式是直接推论；把它与共同 view 规则和节点测度量词组合，是本项目最值得保留的 theorem 形式之一，但原创性仍需同行检索确认。

若所有相关 closed views 可在 $r$ 轮加固定 halo 内获得，且

\[
\sup_{\theta,i}\gamma_{\theta,i}(r)\le\Gamma_r<\infty,
\tag{4.16}
\]

则

\[
\boxed{
\mathcal G_{r,\delta}
\le\mathcal A^{\rm any,R}_{r,\delta}
\le\mathcal A^{\rm lin,R}_{r,\delta}
\le\Gamma_r\mathcal G_{r,\delta}.}
\tag{4.17}
\]

#### 证明

对任意模型、任意根、任意局部解码器，固定模型引理 2.1 给出误差至少
$g_{\theta,i}(r)$。因此逐节点误差支配 $g$；分位数、对模型 supremum 和对算法 infimum 都保持此下界。

共同 Dirichlet 规则 (4.14) 是线性的，并逐节点满足
$d_{\theta,i}\le\Gamma_rg_{\theta,i}$。有限节点分位数的单调性给出

\[
Q_{1-\delta}(d)\le\Gamma_rQ_{1-\delta}(g).
\]

再对模型取 supremum 得最右上界。中间不等式来自规则类包含关系。证毕。

不要求处处有统一 $\Gamma$ 时，以下阈值形式更准确。

### 推论 4.7（允许局部病态节点）

**状态：[直接推论]。**

若对每个模型

\[
\mu_\theta\{i:g_{\theta,i}>a\}\le\delta_1,
\qquad
\mu_\theta\{i:\gamma_{\theta,i}>b\}\le\delta_2,
\tag{4.18}
\]

则同一个 Dirichlet 局部规则在至少
$1-\delta_1-\delta_2$ 测度节点上，对所有单位输入误差不超过 $ab$。

#### 证明

在两好集的交上，$d_i\le\gamma_ig_i\le ab$；坏集包含于两个坏集之并，使用 union bound。证毕。

注意 (4.18) 仍含全局量 $g_i$；$\gamma_i$ 是局部可算的，但它本身不能替代 exterior influence。

### 推论 4.8（统一谱窗）

**状态：[Kantorovich/antieigenvalue 型标准谱事实 + 直接推论]。**

若整个模型类满足

\[
mI\preceq J_\theta\preceq MI,\qquad 0<m\le M<\infty,
\tag{4.19}
\]

则

\[
\boxed{
\|E^*A^{-1}\|
\le\frac{M-m}{2\sqrt{mM}}
=\frac{\kappa-1}{2\sqrt\kappa},}
\tag{4.20}
\]

从而

\[
\boxed{
\gamma_{\theta,i}(r)
\le
\frac{M+m}{2\sqrt{mM}}
=\frac{\kappa+1}{2\sqrt\kappa},
\qquad \kappa=M/m.}
\tag{4.21}
\]

#### 证明

谱演算给出

\[
(J-mI)(MI-J)\succeq0.
\]

把它压缩到 $B$ 块，并使用

\[
(J^2)_{BB}=A^2+EE^*,
\]

得到

\[
EE^*\preceq(M+m)A-A^2-mMI.
\tag{4.22}
\]

左右同乘 $A^{-1}$：

\[
\begin{aligned}
A^{-1}EE^*A^{-1}
&\preceq
(M+m)A^{-1}-I-mMA^{-2}.
\end{aligned}
\tag{4.23}
\]

由于 $mI\preceq A\preceq MI$，右侧最大特征值不超过标量函数

\[
\max_{\lambda\in[m,M]}
\left(\frac{M+m}{\lambda}-1-\frac{mM}{\lambda^2}\right).
\]

最大点是

\[
\lambda_*=\frac{2mM}{m+M},
\]

最大值为 $(M-m)^2/(4mM)$。左侧是
\((E^*A^{-1})^*(E^*A^{-1})\)，开平方得到 (4.20)；代入
\(\gamma=\sqrt{1+\|F\|^2}\) 得 (4.21)。证毕。

常数是 universal sharp。取二维矩阵

\[
J=U\operatorname{diag}(m,M)U^*,
\]

并选择旋转使 $A=2mM/(m+M)$，等价地令相应权重
\(\cos^2\vartheta=M/(m+M)\)、\(\sin^2\vartheta=m/(m+M)\)。此时

\[
\frac{|E|}{A}=\frac{M-m}{2\sqrt{mM}},
\]

所以 (4.20)–(4.21) 取等。

较直接的估计
\[
\|E\|\le(M-m)/2,\qquad \|A^{-1}\|\le1/m
\]
只给
\(\|E^*A^{-1}\|\le(M-m)/(2m)\)。它是正确但不 sharp 的粗界；第二阶段报告中的对应常数应按此理解，而非错误结论。

这不要求无环、有限树宽或逐节点测量矩阵满列秩。它是强的全局稳定性充分条件，不是必要条件。

### 4.4 一个局部边界数值及其局限

令

\[
\tau_i(B):=\|R_iA^{-1}E\|=\|C\|.
\tag{4.24}
\]

若额外知道 exterior Schur 补满足
$mI\preceq\Sigma\preceq MI$，则

\[
\frac{\tau_i(B)}{M}\le g_i(B)\le\frac{\tau_i(B)}m.
\tag{4.25}
\]

**状态：[直接推论]。** 因为 $C=(C\Sigma^{-1})\Sigma$ 给出左界，而
$\|C\Sigma^{-1}\|\le\|C\|/m$ 给出右界。没有 exterior lower bound 时，
$\tau_i$ 再小也可能被 $\Sigma^{-1}$ 任意放大；第 7 节给出精确反例。

---

## 5. 原始 WLS/Gaussian 数据：必须审计 \(G\) 映射

### 5.1 从传感器观测到正规方程右端

设

\[
J=H^*R^{-1}H+Q\succ0,\qquad
R\succ0,\quad Q\succeq0,
\tag{5.1}
\]

并令

\[
G:=H^*R^{-1},\qquad K:=J^{-1}G.
\tag{5.2}
\]

中心 WLS/MAP 估计是 \(x=Kz\)，而第 4 节研究的是给定
\(b=Gz\) 后的 \(x=J^{-1}b\)。二者不能在没有条件时互换。

令 \(P_B^x:\mathsf X\to\mathsf X_B\) 是状态限制，
\(P_M^z:\mathsf Z\to\mathsf Z_M\) 是已收集测量限制；在公式中把
\((P_M^z)^*P_M^z\) 简写为同名正交投影 \(P_M^z\)。局部 principal 规则为

\[
L_{i,B,M}
=R_iA^{-1}P_B^xG P_M^z.
\tag{5.3}
\]

### 命题 5.1（raw WLS 的精确误差分解）

**状态：[标准代数恒等式]。**

\[
\boxed{
[R_i,0]K-L_{i,B,M}
=
\bigl([R_i,0]J^{-1}-R_iA^{-1}P_B^x\bigr)G
+R_iA^{-1}P_B^xG(I-P_M^z).}
\tag{5.4}
\]

#### 证明

在右侧加减 \(R_iA^{-1}P_B^xG\)，再用 \(K=J^{-1}G\)。证毕。

第一项是 Schur/边界误差经过 \(G\) 后的作用；第二项是形成球内
\(b_B= P_B^xGz\) 时漏掉的测量贡献。两项可能相消，所以一般只有

\[
\|[R_i,0]K-L_{i,B,M}\|
\le d_i(B)\|G\|
+\|R_iA^{-1}P_B^xG(I-P_M^z)\|.
\tag{5.5}
\]

若 \(P_B^xG(I-P_M^z)=0\)，第二项才严格消失。换言之，局部算法必须收齐所有会贡献到 \(b_B\) 的测量。

### 命题 5.2（\(K\) 与 \(J^{-1}\) 行尾的半径平移）

**状态：[直接推论]。** 这是有限传播支持关系；双向结论必须显式假设局部右逆。

对根 \(i\)，令 \(Q_i^x(t)\) 投影到距离 \(i\) 大于 \(t\) 的状态块，
\(Q_i^z(r)\) 投影到所有者距离 \(i\) 大于 \(r\) 的测量块，并定义

\[
j_i(t):=\|[R_i,0]J^{-1}Q_i^x(t)\|,
\qquad
k_i(r):=\|[R_i,0]KQ_i^z(r)\|.
\tag{5.6}
\]

若

\[
G_{u,a}=0\quad\text{只要}\quad
d_G(u,o(a))>h,
\tag{5.7}
\]

则对 \(r\ge h\)

\[
\boxed{k_i(r)\le\|G\|\,j_i(r-h).}
\tag{5.8}
\]

若还存在 \(L_G:\mathsf X\to\mathsf Z\) 满足

\[
GL_G=I_{\mathsf X},\qquad
\|L_G\|\le\Lambda,\qquad
(L_G)_{a,u}=0\ \text{当}\ d_G(o(a),u)>h',
\tag{5.9}
\]

则

\[
\boxed{j_i(r+h')\le\Lambda\,k_i(r).}
\tag{5.10}
\]

#### 证明

有限传播 (5.7) 给出支持恒等式

\[
GQ_i^z(r)=Q_i^x(r-h)GQ_i^z(r).
\]

左乘根行并取范数得到 (5.8)。同理，(5.9) 给出

\[
L_GQ_i^x(r+h')=Q_i^z(r)L_GQ_i^x(r+h').
\]

于是

\[
\begin{aligned}
[R_i,0]J^{-1}Q_i^x(r+h')
&=[R_i,0]J^{-1}GL_GQ_i^x(r+h')\\
&=[R_i,0]KQ_i^z(r)L_GQ_i^x(r+h'),
\end{aligned}
\]

取范数即得 (5.10)。证毕。

因此，“\(J^{-1}\) 的状态行尾小”推出“\(K\) 的测量行尾小”需要
\(G\) 局部；反向还需要一个传播有限、范数受控的右逆。仅有 \(H\) 局部不够，因为稠密噪声精度 \(R^{-1}\) 可使 \(G=H^*R^{-1}\) 立即变成全局算子。

### 5.2 principal block 与只保留球内因子

若实际算法丢弃跨边界测量，得到 \(\widetilde A\)，则

\[
A^{-1}-\widetilde A^{-1}
=A^{-1}(\widetilde A-A)\widetilde A^{-1}.
\tag{5.11}
\]

**状态：[标准 resolvent 恒等式]。** 它给出 triangle upper bound；没有符号、Loewner 单调性或 range 对齐时不存在一般反向界，因为该误差可以与 (5.4) 的 Schur 误差抵消。故第 4 节主定理不能直接套给“只解球内因子”的实现，除非证明它确实形成 \(J_{BB}\) 或另行控制 (5.11)。

### 5.3 Gaussian/Bayes 风险

设 \((X_i,Z)\) 是中心化联合 Gaussian，且中心后验均值

\[
Y:=\mathbb E[X_i\mid Z]=T_iZ.
\]

把 \(Z=(Z_B,Z_O)\)，并设 \(\Sigma_{BB}\succ0\)。

### 命题 5.3（局部数据相对中心估计的 excess Bayes risk）

**状态：[标准 Gaussian 投影恒等式]。**

局部 Bayes 最优估计为 \(\mathbb E[X_i\mid Z_B]\)，且

\[
\boxed{
\Delta_i(B):=
\mathbb E\|Y-\mathbb E[Y\mid Z_B]\|^2
=
\operatorname{tr}\!\left(
T_{i,O}\Sigma_{O\mid B}T_{i,O}^*
\right),}
\tag{5.12}
\]

其中

\[
\Sigma_{O\mid B}
=\Sigma_{OO}-\Sigma_{OB}\Sigma_{BB}^{-1}\Sigma_{BO}.
\tag{5.13}
\]

此外，\(\Delta_i(B)\) 正好等于只用 \(Z_B\) 与使用全部 \(Z\) 估计
\(X_i\) 的 MSE 之差。

#### 证明

塔式性质给出

\[
\mathbb E[X_i\mid Z_B]
=\mathbb E[\mathbb E[X_i\mid Z]\mid Z_B]
=\mathbb E[Y\mid Z_B].
\]

联合 Gaussian 的线性条件期望给出

\[
Y-\mathbb E[Y\mid Z_B]
=T_{i,O}\bigl(Z_O-\mathbb E[Z_O\mid Z_B]\bigr),
\]

括号内条件协方差为 (5.13)，取平方期望得到 (5.12)。最后，条件期望是
\(L_2\) 正交投影；嵌套信息空间的 Pythagoras 恒等式给出 MSE 差。证毕。

若

\[
cI\preceq\Sigma_Z\preceq CI,
\tag{5.14}
\]

则 Schur 补也满足 \(cI\preceq\Sigma_{O\mid B}\preceq CI\)，从而

\[
c\|T_{i,O}\|_F^2
\le\Delta_i(B)
\le C\|T_{i,O}\|_F^2.
\tag{5.15}
\]

这是 adversarial \(\ell_2\) 行尾与 Bayes 价值之间的条件化桥梁。没有 (5.14) 时，远端数据可几乎被本地数据预测，或条件方差可异常放大，二者不等价。

---

## 6. grounded Laplacian/M-matrix：killed-walk 精确公式

### 6.1 归一化与 Markov 链

设

\[
J=L_W+\operatorname{diag}(\kappa),\qquad
w_{uv}=w_{vu}\ge0,\quad \kappa_u\ge0,
\tag{6.1}
\]

其中每个连通分量都能到达某个 \(\kappa_u>0\) 的节点。令

\[
s_u=\kappa_u+\sum_vw_{uv},\qquad
S=\operatorname{diag}(s_u),\qquad
P_{uv}=\frac{w_{uv}}{s_u}.
\tag{6.2}
\]

于是

\[
J=S(I-P),\qquad
N:=(I-P)^{-1}=\sum_{t=0}^{\infty}P^t,
\qquad
J^{-1}=NS^{-1}.
\tag{6.3}
\]

这里 \(P\) 是 substochastic 矩阵；从 \(u\) 出发，每步以
\(\kappa_u/s_u\) 的概率进入 cemetery 状态 \(\dagger\)。有限性和“每个分量可到达 killing”保证 \(\rho(P)<1\)。

把输入归一成

\[
c=S^{-1}b,\qquad x=Nc,
\tag{6.4}
\]

并在本节先约束 \(\|c\|_\infty\le1\)。若约束的是
\(\|b\|_\infty\le1\)，下面每次访问节点 \(j\) 的 reward 要改成 \(1/s_j\)；只有当 \(s_j\) 有统一上下界时，两种输入归一化才等价到常数。

令 \((X_t)\) 为上述 killed chain，

\[
\zeta:=\inf\{t\ge0:X_t=\dagger\},\qquad
L_j:=\sum_{t=0}^{\zeta-1}\mathbf 1\{X_t=j\}.
\]

### 引理 6.1（fundamental matrix 的占用解释）

**状态：[标准恒等式]。**

\[
\boxed{N_{ij}=\mathbb E_iL_j.}
\tag{6.5}
\]

#### 证明

\((P^t)_{ij}=\Pr_i(X_t=j,\ t<\zeta)\)。对 \(t\ge0\) 求和，并用非负项的 Tonelli 定理：

\[
N_{ij}
=\sum_{t\ge0}(P^t)_{ij}
=\mathbb E_i\sum_{t<\zeta}\mathbf1\{X_t=j\}.
\]

证毕。

这是 absorbing Markov chain/离散 potential theory 的经典 fundamental matrix 表示，不是本文的新结果。

### 6.2 oracle、Dirichlet 与 Neumann 三种误差

固定根 \(i\)，取 \(B=B_r(i)\)。令

\[
\tau_B:=\inf\{t\ge0:X_t\notin B\},
\]

并约定被 kill 前未出球时 \(\tau_B=\infty\)。记
\(P_{BB}\) 为 \(P\) 的 principal block。

### 背景命题 6.2（三种局部误差的精确路径公式）

**状态：[标准表示的直接组合]。** 每个 Green/Dirichlet/Neumann 路径展开都是经典的；把三者并列用于区分三种信息算法是直接组合，不应把路径表示本身报为原创。

对标量根输出和单位 \(\ell_\infty\) 输入：

1. 固定完整模型、允许 oracle 使用全局系数但只读取 \(c_B\) 时，任意解码器的最小最坏误差是

   \[
   \boxed{
   o_i(r)
   :=\sum_{j\notin B}N_{ij}
   =\mathbb E_i\sum_{t=0}^{\zeta-1}
   \mathbf1\{X_t\notin B\}.}
   \tag{6.6}
   \]

2. 只用局部转移块的 Dirichlet 规则

   \[
   x_i^D=e_i^*(I-P_{BB})^{-1}c_B
   \tag{6.7}
   \]

   的精确最坏误差是

   \[
   \boxed{
   d_i(r)
   =\mathbb E_i\!\left[
   (\zeta-\tau_B)\mathbf1\{\tau_B<\zeta\}
   \right].}
   \tag{6.8}
   \]

3. 有限时间 Neumann/message-passing 规则

   \[
   x_i^{N,r}=e_i^*\sum_{t=0}^{r}P^tc
   \tag{6.9}
   \]

   的精确最坏误差是

   \[
   \boxed{
   n_i(r)
   =e_i^*\sum_{t=r+1}^{\infty}P^t\mathbf1
   =\mathbb E_i[(\zeta-r-1)_+].}
   \tag{6.10}
   \]

并且逐路径有

\[
\boxed{o_i(r)\le d_i(r)\le n_i(r).}
\tag{6.11}
\]

#### 证明

由引理 2.1 的 \(\ell_\infty\) 版本，oracle 误差是根行在球外坐标上的
\(\ell_1\) 范数。因为 \(N_{ij}\ge0\)，它等于球外行和。使用引理 6.1 并交换有限/非负求和，得到 (6.6)。

\((I-P_{BB})^{-1}=\sum_{t\ge0}P_{BB}^t\) 只计算在 kill 或首次出球以前始终留在 \(B\) 的路径。因此全局 Green 行减去 (6.7) 的行，恰由“已经出球”的路径组成，所有系数非负。其 \(\ell_\infty\to\mathbb R\) 范数等于系数行和；逐条样本路径上，被删去的 live time 是
\(\zeta-\tau_B\)（若先被 kill 则为 0），得到 (6.8)。

Neumann 余项的系数也非负，其行和为

\[
\sum_{t=r+1}^{\infty}\Pr_i(t<\zeta)
=\mathbb E_i\#\{t:r+1\le t<\zeta\}
=\mathbb E_i[(\zeta-r-1)_+],
\]

得到 (6.10)。最后，出球至少需要 \(r+1\) 步；每条路径的球外访问数不超过出球后的 live time，而出球后的 live time不超过第 \(r+1\) 步后的 live time。因此 (6.11) 逐路径成立。证毕。

这里 oracle 与 Dirichlet 的区别很重要：oracle 对本地 \(c_B\) 使用的是完整
\(N\) 的本地列，其中包含“先出球再返回”的影响；Dirichlet 规则不知道这些 completion-dependent 系数，所以把首次出球后的整段路径全部舍去。

### 6.3 escape probability 与条件剩余寿命：精确式和可检验桥梁

定义 kill 前出球概率

\[
h_i(r):=\Pr_i(\tau_B<\zeta)
=e_i^*(I-P_{BB})^{-1}P_{BO}\mathbf1.
\tag{6.12}
\]

若闭合 view 告知边界节点的总 outgoing transition probability，
\(h_i(r)\) 完全由球内系数计算，不需要知道球外如何 completion。

在 \(h_i(r)>0\) 时定义两个条件均值

\[
\begin{aligned}
m_i^{\rm out}(r)
&:=
\mathbb E_i\!\left[
\sum_{t<\zeta}\mathbf1\{X_t\notin B\}
\ \middle|\ \tau_B<\zeta
\right],\\
\ell_i^{\rm rem}(r)
&:=
\mathbb E_i[
\zeta-\tau_B\mid\tau_B<\zeta].
\end{aligned}
\tag{6.13}
\]

若 \(h_i=0\)，约定两者为 0。

### 命题 6.3（精确乘积分解）

**状态：[标准条件期望恒等式]。** 这是精确 iff 的正确量词，但本身只是对误差的重写，不是新的结构判据。

\[
\boxed{
o_i(r)=h_i(r)m_i^{\rm out}(r),\qquad
d_i(r)=h_i(r)\ell_i^{\rm rem}(r).}
\tag{6.14}
\]

当 \(h_i>0\) 时

\[
1\le m_i^{\rm out}(r)\le\ell_i^{\rm rem}(r).
\tag{6.15}
\]

#### 证明

球外占用数和出球后 live time 在未出球事件上均为 0。对出球事件条件化即得 (6.14)。成功出球至少在球外访问一次，而球外访问次数不超过出球后 live time，得到 (6.15)。证毕。

因此，对任意模型族和节点测度，Dirichlet 误差按节点测度趋零的**完全精确**条件是

\[
h_i(r)\ell_i^{\rm rem}(r)\longrightarrow0
\quad\text{一致地依 }\mu_\theta\text{-概率},
\tag{6.16}
\]

即对每个 \(\varepsilon>0\)，

\[
\sup_\theta
\mu_\theta\{i:h_i(r)\ell_i^{\rm rem}(r)>\varepsilon\}
\longrightarrow0.
\tag{6.17}
\]

但 (6.16) 与 \(d_i\to0\) 由 (6.14) 完全同义，不能把它包装成一个新的可检验定理。真正的结构问题是：何时能去掉未知 exterior 的
\(\ell_i^{\rm rem}\)，只检查局部可算的 \(h_i\)。

### 定理 6.4（escape-only 桥梁的准确假设）

**状态：[直接推论；可能有用的组合]。** 条件期望与 union bound 是标准的；“多数根 + 条件余寿命 tightness”的组合可作为论文假设框架，但 novelty 未确认。

总有

\[
h_i(r)\le o_i(r)\le d_i(r).
\tag{6.18}
\]

进一步：

1. 若对所有模型、根和半径

   \[
   \ell_i^{\rm rem}(r)\le L,
   \tag{6.19}
   \]

   则

   \[
   h_i(r)\le o_i(r)\le d_i(r)\le Lh_i(r).
   \tag{6.20}
   \]

   因而 \(h_i(r)\to0\) 与 \(d_i(r)\to0\) 在 all-node 或节点测度语义下等价到常数。

2. 更一般地，若对每个模型

   \[
   \mu_\theta\{i:h_i(r)>a\}\le\delta_1,\qquad
   \mu_\theta\{i:\ell_i^{\rm rem}(r)>L\}\le\delta_2,
   \tag{6.21}
   \]

   则

   \[
   \mu_\theta\{i:d_i(r)>aL\}\le\delta_1+\delta_2.
   \tag{6.22}
   \]

3. 若 \(\ell_i^{\rm rem}(r)\) 对 \((\theta,r)\) 在节点测度下**一致 tight**，即

   \[
   \lim_{L\to\infty}
   \sup_{\theta,r}
   \mu_\theta\{i:\ell_i^{\rm rem}(r)>L\}=0,
   \tag{6.23}
   \]

   则

   \[
   h_i(r)\to0\text{ 一致依节点测度}
   \quad\Longleftrightarrow\quad
   d_i(r)\to0\text{ 一致依节点测度}.
   \tag{6.24}
   \]

#### 证明

(6.18) 来自 (6.15)。在 (6.19) 下用 (6.14) 得 (6.20)。
在 (6.21) 的两个好集交上，\(d=h\ell\le aL\)，所以 (6.22) 是 union bound。

对 (6.24)，\(d\to0\Rightarrow h\to0\) 来自 \(h\le d\)。反向固定
\(\varepsilon,\eta>0\)。由 tightness 选 \(L\) 使第一项小于 \(\eta/2\)，再由 \(h\to0\) 取足够大的 \(r\) 使第二项小于 \(\eta/2\)：

\[
\mu_\theta\{h\ell>\varepsilon\}
\le
\mu_\theta\{\ell>L\}
+\mu_\theta\{h>\varepsilon/L\}.
\]

对 \(\theta\) 取 supremum 即得。证毕。

一个更强但易陈述的 sufficient promise 是：对所有可能的首次出球点
\(y\)，从 \(y\) 开始的期望剩余寿命不超过 \(L\)。强 Markov 性立即推出
(6.19)。例如每个 live 状态每步 killing 概率至少 \(p_0>0\) 时，
\(L\le1/p_0\)。

必须强调：\(\ell_i^{\rm rem}\) 本身依赖 exterior，不是局部证书。
(6.14) 是精确但同义的 iff；(6.19)/(6.23) 才把局部可算的 \(h_i\) 变成可用桥梁，而它们是模型类的全局 promise。没有这类 promise，escape probability 很小但逃出后寿命极长，乘积仍可很大。

### 6.4 Neumann 的对应乘积

令

\[
s_i(r):=\Pr_i(\zeta>r+1)=e_i^*P^{r+1}\mathbf1,
\]

并在 \(s_i(r)>0\) 时定义

\[
\ell_i^{N}(r):=
\mathbb E_i[\zeta-r-1\mid\zeta>r+1].
\]

### 命题 6.5（survival–remainder 分解）

**状态：[标准条件期望恒等式]。**

\[
n_i(r)=s_i(r)\ell_i^N(r).
\tag{6.25}
\]

因此 survival probability 单独控制 Neumann 误差也需要条件余寿命有界或在节点测度下一致 tight；证明与定理 6.4 完全相同。

### 6.5 \(\ell_2\) 输入的两条独立 walk 表示

对归一化 \(c\) 的单位 \(\ell_2\) 球，固定模型的 oracle 尾为

\[
o_{i,2}(r)^2=\sum_{j\notin B}N_{ij}^2.
\tag{6.26}
\]

### 命题 6.6（球外交叉局部时间）

**状态：[标准占用恒等式的直接推论]。**

取两条从 \(i\) 出发、相互独立的 killed walks，其局部时间为
\(L_j^{(1)},L_j^{(2)}\)。则

\[
\boxed{
o_{i,2}(r)^2
=\mathbb E_i
\sum_{j\notin B}L_j^{(1)}L_j^{(2)}.}
\tag{6.27}
\]

#### 证明

独立性和引理 6.1 给出

\[
\mathbb E[L_j^{(1)}L_j^{(2)}]
=(\mathbb E L_j)^2=N_{ij}^2.
\]

求和即可。证毕。

这说明 \(\ell_2\) locality 对应两条 walk 的 intersection local time，而不是单条 walk 的总 occupation。单独的 escape probability 一般不能双向控制 (6.26)：相同 \(\ell_1\) occupation 可集中在少数节点，也可分散到巨大体积。

---

## 7. completion-rich 模型类的精确 no-go

### 定理 7.1（任意 SPD Schur completion 的零—无穷二分）

**状态：[直接推论；可能新的 completion-minimax 表述]。**
证明只是 Schur 公式，但它给出了一个真正的结构结论：若对 exterior 不施加任何稳定性下界，则“有限统一误差”与“根到割边精确解耦”等价，不存在中间情形。

固定局部 \(A\succ0\)、\(E\) 和根抽取 \(R_i\)，并令

\[
C=R_iA^{-1}E.
\]

考虑完整 completion 类

\[
\mathfrak J(A,E)
=\left\{
\begin{bmatrix}
A&E\\E^*&E^*A^{-1}E+\Sigma
\end{bmatrix}
:\ \Sigma\succ0
\right\}.
\tag{7.1}
\]

对 primitive RHS \(b\)，记

\[
\mathcal E(J,f):=
\sup_{\|b\|_2\le1}
\|[R_i,0]J^{-1}b-f(b_B)\|.
\]

则

\[
\boxed{
\sup_{J\in\mathfrak J(A,E)}
\inf_{f_J}\mathcal E(J,f_J)
=
\inf_{\text{共同 }f}
\sup_{J\in\mathfrak J(A,E)}
\mathcal E(J,f)
=
\begin{cases}
0,&C=0,\\
+\infty,&C\ne0.
\end{cases}}
\tag{7.2}
\]

左边允许 \(f_J\) 知道各自完整 completion；中间要求所有 completion
共用同一个局部 \(f\)。两种信息条件仍给出相同的零—无穷二分。

#### 证明

若 \(C=0\)，由块逆公式，根行的球外块
\(-C\Sigma^{-1}=0\)，球内修正
\(C\Sigma^{-1}E^*A^{-1}=0\)。所以对所有 completion，

\[
[R_i,0]J^{-1}=[R_iA^{-1},0],
\]

共同 Dirichlet 规则精确。

若 \(C\ne0\)，取 \(\Sigma_\eta=\eta I\)、\(\eta\downarrow0\)。每个
\(J_\eta\) 都 SPD，而固定模型 oracle 尾为

\[
g_i(B)=\|C\Sigma_\eta^{-1}\|=\frac{\|C\|}{\eta}\to\infty.
\]

引理 2.1 对任意解码器给出同一下界，故 supremum 为无穷。证毕。

式 (7.2) 对“任意 SPD completion”是完整 iff，但它也说明这个模型类太宽：除了精确解耦外，没有非平凡局部鲁棒性。要得到渐近衰减，必须缩小 completion 类，例如要求
\(\Sigma\succeq mI\)、限制余寿命/规模，或给定随机生成律。

### 命题 7.2（M-matrix 内的二维最小见证）

**状态：[直接构造]。**

取

\[
J_d=
\begin{bmatrix}
1&-c\\-c&d
\end{bmatrix},
\qquad c>0,\quad d>c^2.
\tag{7.3}
\]

它是 SPD Z-matrix，因而是 nonsingular symmetric M-matrix，且

\[
e_1^*J_d^{-1}
=\frac1{d-c^2}[d,c].
\tag{7.4}
\]

根的局部块 \(A=1\) 与 cut coefficient \(E=-c\) 固定，而
\(d\downarrow c^2\) 时本地系数和球外系数同时爆炸。它证明 no-go 并非 generic block 的符号抵消伪影；一维 exterior 已足够。

更一般地，若 \(A\) 是 nonsingular M-matrix、\(E\le0\) 是单列且
\(R_iA^{-1}E\ne0\)，取

\[
D_\eta=E^*A^{-1}E+\eta
\]

仍得到 SPD Z-matrix。多维 exterior 时
\(E^*A^{-1}E\) 可能有正 off-diagonal，不能不加条件地把定理 7.1 的整个 completion 类称为 M-matrix。

### 定理 7.3（有界度 grounded-Laplacian 的 dangling-tentacle no-go）

**状态：[直接构造；可能新的 LOCAL completion 用法]。**
电网络 reciprocity 与逆 M-matrix 正性都是标准事实；把它们用于固定视图 minimax no-go 是直接组合。

固定一个有限局部 core，含根 \(i\)、边界节点 \(a\)，并在 closed view 中声明一条从 \(a\) 到球外端口 \(u_0\) 的固定正权边。假设 completion 类允许：

1. 在 core 内或 \(u_0\) 处放置一个固定正 ground，使全图 SPD；
2. 从 \(u_0\) 再接一条无 ground 的路径
   \(u_0-u_1-\cdots-u_N\)，路径边权固定为 \(w>0\)；
3. \(N\) 任意增大，且 \(u_0\) 以外的信息不出现在根的半径 \(r\) view。

令

\[
J_N=L_{W_N}+\operatorname{diag}(\kappa_N).
\]

则存在与 \(N\) 无关的常数 \(c_i>0\)，使

\[
(J_N^{-1})_{i,u_k}=c_i,\qquad k=1,\ldots,N.
\tag{7.5}
\]

因而

\[
\boxed{
\|e_i^*J_N^{-1}P_O\|_2\ge c_i\sqrt N,\qquad
\|e_i^*J_N^{-1}P_O\|_1\ge c_iN.}
\tag{7.6}
\]

图的最大度、非零边权上下界和 ground 权重都可与 \(N\) 无关。

#### 证明

在根 \(i\) 注入单位电流并由 ground 吸收，电位向量是
\(v=J_N^{-1}e_i\)。无 ground 的 dangling path 没有任何电流 sink。由叶端 Kirchhoff 方程开始反推，每条路径边电流都为 0，因此

\[
v_{u_1}=\cdots=v_{u_N}=v_{u_0}.
\]

路径对 \(u_0\) 的净电流也为 0，所以 \(v_{u_0}\) 与 \(N\) 无关，等于把 dangling path 删除后的端口电位。整个 grounded 网络连通时，nonsingular irreducible M-matrix 的逆严格为正，故
\(c_i:=v_{u_0}>0\)。

由 \(J_N^{-1}\) 对称，

\[
(J_N^{-1})_{i,u_k}=(J_N^{-1})_{u_k,i}=v_{u_k}=c_i.
\]

对这 \(N\) 个球外坐标求 \(\ell_2\) 或 \(\ell_1\) 行范数，得到 (7.6)。路径内部度最多 2，所有新增权重固定，所以有界度/有界系数声明成立。证毕。

若使用第 6 节的归一输入 \(c=S^{-1}b\)，中心算子为
\(N=J_N^{-1}S\)。路径节点的 \(s_j\) 在固定正区间内，所以相同的
\(\sqrt N\) 与 \(N\) 发散仍成立，只改变常数。

### 推论 7.4（local-only uniform certificate 不可能）

**状态：[直接推论]。**

在定理 7.3 的 completion-closed 类中，任意固定半径 \(r\)：

- 对单位 \(\ell_2\) primitive input，固定模型 oracle 误差随 \(N\) 至少按
  \(\sqrt N\) 增长；
- 对单位 \(\ell_\infty\) input，误差至少按 \(N\) 增长；
- 因此即使算法知道全局规模 \(N\)，任何只读取球内数据的规则也没有统一有限误差；若 \(N\) 不属于节点初始标签，同一个局部 view 的 completion minimax 直接为无穷。

证明只需把 (7.6) 代入引理 2.1。允许少量节点失败不会自动消除此障碍：必须证明这种 exterior 对绝大多数根不可达或其 Green 行尾确实小，而不能只数“病态节点”本身占比。

---

## 8. 最小候选主定理及其发表性判断

### 定理 8.1（共同局部规则、节点分位数与 completion 二分）

**状态：[可能新的组合；数学证明闭合，原创性与发表性未闭合]。**

设 \(\Theta\) 是有限 block-SPD primitive systems

\[
J_\theta x=b,\qquad J_\theta=J_\theta^*\succ0
\]

的模型类，并满足：

1. 状态块和 \(b\)-块由通信节点持有；
2. 根在 \(r\) 轮加一个固定 halo 后取得系数闭合 view，故该 view 确定
   \(A=J_{BB}\)、cut ports 与 \(E=J_{BO}\)；
3. 同构 view 使用同一个确定性局部规则；
4. 性能采用单位 \(\ell_2\) primitive input 和节点鲁棒多数语义 (R)；
5. 每个模型带任意节点概率测度 \(\mu_\theta\)。

定义

\[
\begin{aligned}
\mathcal G_{r,\delta}
&=\sup_\theta Q_{1-\delta}^{\mu_\theta}(g_{\theta,i}),\\
\mathcal D_{r,\delta}
&=\sup_\theta Q_{1-\delta}^{\mu_\theta}(d_{\theta,i}),\\
\mathcal A^{\rm any}_{r,\delta}
&=\inf_{\text{共同 view-based }f}
\sup_\theta Q_{1-\delta}^{\mu_\theta}
\left(
\sup_{\|b\|_2\le1}
\|[R_i,0]J_\theta^{-1}b-f_{\mathcal V_r}(b_B)\|
\right),
\end{aligned}
\tag{8.1}
\]

并类似定义只允许线性规则的
\(\mathcal A^{\rm lin}_{r,\delta}\)。则：

**(a) 无附加谱条件的基本链**

\[
\boxed{
\mathcal G_{r,\delta}
\le\mathcal A^{\rm any}_{r,\delta}
\le\mathcal A^{\rm lin}_{r,\delta}
\le\mathcal D_{r,\delta}.}
\tag{8.2}
\]

**(b) 局部边界因子控制**

若 \(\gamma_{\theta,i}(r)\le\Gamma_r\) 对所有模型和根成立，则

\[
\boxed{
\mathcal G_{r,\delta}
\le\mathcal A^{\rm any}_{r,\delta}
\le\mathcal A^{\rm lin}_{r,\delta}
\le\mathcal D_{r,\delta}
\le\Gamma_r\mathcal G_{r,\delta}.}
\tag{8.3}
\]

若 \(\sup_r\Gamma_r<\infty\)，则对每个固定 \(\delta\)，下列三件事等价：

\[
\mathcal G_{r,\delta}\to0,\qquad
\mathcal A^{\rm any}_{r,\delta}\to0,\qquad
\mathcal A^{\rm lin}_{r,\delta}\to0.
\tag{8.4}
\]

实现 sufficiency 的是同一个、完全由 view 决定的 principal/Dirichlet 规则，而不是为每个 completion 单独挑选的 oracle 系数。

**(c) 大多数边界良态而非处处良态**

若 \(g_i\le a_r\) 至少在 \(1-\delta_1\) 测度节点上成立，而
\(\gamma_i\le b_r\) 至少在 \(1-\delta_2\) 测度节点上成立，则同一个
Dirichlet 规则在至少 \(1-\delta_1-\delta_2\) 测度节点上误差不超过
\(a_rb_r\)。

**(d) 无 exterior stability 时的二分**

对任一固定 view，若 completion 类包含 (7.1) 的所有 SPD Schur 补，则其共同局部 minimax 有限当且仅当

\[
R_iA^{-1}E=0;
\tag{8.5}
\]

此时误差实际为 0。若不等于 0，minimax 为 \(+\infty\)。在单输出、单 cut 且只允许
\(m\le\Sigma\le M\) 时，精确的非平凡中间值和最优共同系数由
(4.14b)–(4.14d) 给出。

#### 证明

(a) 的第一项是引理 2.1 对每个固定模型和根的下界；任意规则类包含线性规则类的关系给出中间项；共同 Dirichlet 规则的实际误差就是
\(d_{\theta,i}\)，给出最后一项。

(b) 由逐节点 Schur 分解
\(d_{\theta,i}\le\gamma_{\theta,i}g_{\theta,i}\) 和分位数单调性得到。
(8.4) 由 sandwich theorem。  
(c) 是两个好集取交并用 union bound。  
(d) 是定理 7.1 与命题 4.5。证毕。

### 8.2 这一定理比逐项拼接多了什么

严格来说，(8.2)–(8.4) 的证明就是四个标准步骤的组合；它没有产生一个新的深谱定理。它提供的实质主要在量词而不在代数：

1. **算法选择在 completion supremum 之前。** 同一局部 view 必须选同一个规则，不能把每个固定模型的 oracle 截断系数偷偷传给节点。
2. **necessary 与 constructive sufficient 使用同一个量 \(g_i\)。**
   \(g_i\) 是任何算法的固定模型下界，而 principal solve 是只看本地模型的统一构造；两者只差一个可审计的局部边界放大因子。
3. **异常节点可以局部化。** 不要求全网统一 condition number，而允许
   \(g_i\) 与 \(\gamma_i\) 的坏集分别占小测度，再精确累计失败比例。
4. **假设边界被 no-go 钉死。** 若 exterior Schur 补没有任何 lower stability，结论不是“常数变差”，而是除精确解耦外 minimax 直接无穷。
5. **标量 cut 的 locality-constrained 最优系数不是 Dirichlet。**
   命题 4.5 给出 spectrum-window 内、不可见列强制为 0 时的共同 minimax，
   并量化简单 Dirichlet 上界可能浪费的幅度。

其中第 1–4 点作为一个 theorem package 具有澄清价值，但每一步都很短。第 5 点是目前唯一直接求解了 locality-constrained common rule、而非仅把逐点界取分位数的新增计算；但 robust approximate inverse 的总框架已经由 El Ghaoui (2002) 建立，因此这里只剩“局部列约束下的特定闭式是否已有”待核，不能再称为新的 robust-minimax 思路。

### 8.3 发表性如实判断

**当前版本尚不足以稳妥宣称一篇独立理论论文的主定理已经完成。**

原因是：

- fixed-model radius-of-information 是经典信息半径；
- Schur 因子分解接近离散 Dirichlet-to-Neumann/透明边界与 geometric resolvent identity；
- Green occupation 三公式是经典 potential theory 的直接推论；
- 从逐节点 inequality 取 \(1-\delta\) 分位数本身只是单调性和 union bound；
- 零—无穷 completion 二分虽尖锐，但证明是一行 Schur scaling；
- El Ghaoui (2002), Eq. (1.4), Theorem 6.2/Eq. (6.3) 已直接覆盖
  uncertain inverse 的共同 minimax 中心和 structured-SDP 上界；
- 标量 cut 的命题 4.5 只比该已知框架多“选定根行、不可见列强制为 0”
  的特定仿射约束。Eq. (6.5) 不直接给出我们的分段式，但该闭式尚未与
  constrained optimal recovery/robust approximation 文献完成专项查重。

一个更可信的最小发表包应至少再完成以下之一：

1. 不再把一般 robust SDP formulation 当贡献，而是证明
   locality-constrained Schur uncertainty 下某个 **exact** finite-SDP
   reduction/duality theorem，明确超出 El Ghaoui 对 structured uncertainty
   只给 SDP 上界的范围，并刻画最坏 Schur 补；
2. 对一个有自然统计生成机制、但无全局 condition-number 界的图族，给出可由局部抽样验证的多数节点 exterior-stability 条件；
3. 把信息半径 theorem 与 rounds/bits/flops 的 matching lower/upper bounds 结合，而不只证明“收集球后存在一个映射”。

在这些增强完成前，定理 8.1 最适合被称为“严谨理论骨架/主 lemma package”，而不是已经终结该领域的原创主定理。

### 8.4 一个不能进行的量词交换

对每个 view 的 worst-completion 数值 \(a_r(v)\) 很方便，但一般不能用

\[
\sup_\theta
\mu_\theta\{i:a_r(\mathcal V_r(i,\theta))>\varepsilon\}
\tag{8.6}
\]

替代 (8.1) 的 class-level minimax。式 (8.6) 小是一个充分条件；一般不必要，因为不同 view 的最坏 completion 可能互不相容，不能同时出现在同一个全局模型中。只有证明模型类对 disjoint union/pasting 有足够闭包性，才可能交换“逐 view supremum”和“单模型节点分位数”。本文没有假设或证明这种闭包。

---

## 9. 反例套件

### 反例 9.1（单个 anchor 可使几乎所有节点依赖远端）

**状态：[直接构造]。**

取路径状态 \(x_1,\ldots,x_n\)，单位权测量

\[
z_0=x_1,\qquad
z_k=x_{k+1}-x_k,\quad k=1,\ldots,n-1.
\tag{9.1}
\]

令 \(H\) 为该测量矩阵。它是可逆下三角差分矩阵，所以全局 WLS 唯一，

\[
J=H^*H=L_{P_n}+e_1e_1^*\succ0.
\tag{9.2}
\]

考虑两组单位数据

\[
z^+=(1,0,\ldots,0),\qquad z^-=(-1,0,\ldots,0).
\]

中心解分别为

\[
x^+=\mathbf1,\qquad x^-=-\mathbf1.
\tag{9.3}
\]

对距节点 1 超过 \(r\) 的根，它们在两组数据下看到完全相同的本地数据
（全为 0），而中心答案相差 2。因此任意同一本地输出对两答案至少一个误差不小于 1。

更强地，从 \(z_0\) 到任一 \(x_i\) 的中心系数等于 1；当 anchor 在球外时，引理 2.1 直接给出 \(g_i(r)\ge1\)。所以在节点鲁棒语义 (R) 下，除 \(O(r)\) 个靠近 anchor 的节点外全部是坏节点。anchor 自身只占 \(1/n\)，但这不意味着其影响只污染 \(1/n\) 的节点。

这个例子同时说明：

- 有向/无向拓扑是否无环不是固定轮局部估计的充分判据；该图本身就是树；
- 全局 \(H\) 满秩不意味着每个远端根可在固定半径估计；
- “允许少量病态节点失败”不能只按 anchor 数量计数，必须按 Green/估计行尾计数；
- 局部弱收敛或每个固定球看起来正常，不能控制靠近 0 的全局 gauge 模态。

对每个远端节点，两输入中至少一个失败；把失败数在 \(z^+,z^-\) 上相加，还可推出至少一个实际输入使一半以上远端节点失败。这是语义 (I) 下的配套 lower bound，但弱于前述 (R) 结论。

### 反例 9.2（同一 \((A,E)\)，远端 Schur 补近奇异）

**状态：[直接构造]。**

二维族 (7.3) 满足

\[
e_1^*J_d^{-1}
=\frac1{d-c^2}[d,c].
\]

固定 \(A=1,E=-c\)，让 \(d\downarrow c^2\)。根的 principal ball 完全不变且良态，但中心算子的本地块和远端块一起爆炸。它击穿任何只检查
\(\lambda_{\min}(A)\)、局部度数或 cut weight，而不限制 exterior Schur complement 的 uniform certificate。

### 反例 9.3（\(J^{-1}\) 长程，但原始 WLS 增益严格 0-hop）

**状态：[直接构造]。**

取任意稀疏 SPD 矩阵 \(J\)，令

\[
H=I,\qquad R^{-1}=J,\qquad Q=0.
\]

则正规矩阵仍为 \(J\)，但

\[
K=(H^*R^{-1}H)^{-1}H^*R^{-1}
=J^{-1}J=I.
\tag{9.4}
\]

即使 \(J^{-1}\) 没有任何有用的固定半径衰减，原始观测估计器也严格本地。这击穿无条件的

\[
K\text{ 局部}\quad\Longleftrightarrow\quad J^{-1}\text{ 局部}.
\]

第 5 节反向结论中的局部受控右逆不是装饰性条件：本例的自然右逆
\(G^{-1}=J^{-1}\) 正是非局部的。

### 反例 9.4（输入范数与全网算子范数不可偷换）

**状态：[直接构造]。**

在任意 bounded-degree 通信图上，令

\[
H=I,\qquad Q=\lambda I,\qquad
R^{-1}=W=I+\alpha P_n,\qquad
P_n=\frac1n\mathbf1\mathbf1^*,
\tag{9.5}
\]

其中 \(\alpha,\lambda>0\)。谱分解给出

\[
K=(W+\lambda I)^{-1}W=aI+bP_n,
\tag{9.6}
\]

\[
a=\frac1{1+\lambda},\qquad
b=\frac{\alpha\lambda}
{(1+\lambda)(1+\alpha+\lambda)}.
\tag{9.7}
\]

删除任意固定半径球后，rank-one 部分每行的
\(\ell_2\) 尾为 \(b/\sqrt n+o(n^{-1/2})\)，故对单位 \(\ell_2\) 输入每个根都渐近局部；但其
\(\ell_\infty\) 对偶 \(\ell_1\) 行尾趋于 \(b\)，全网
\(2\to2\) operator norm 也等于 \(b\)。此外 \(R^{-1}\) 是稠密的，根不能从 bounded-degree 通信 view 获得它。

所以任何 theorem 必须同时写清：

1. 输入是 \(\ell_2\) 还是 \(\ell_\infty\)；
2. 误差是逐根块行、节点分位数还是全网 operator norm；
3. \(R^{-1}\) 或 \(G\) 是否真由局部模型持有。

### 反例 9.5（per-input majority 严格弱于 robust-node majority）

**状态：[直接构造]。**

令全局误差算子为 \(E=I_n\)，节点均匀测度。每个块行范数都是 1，所以当
\(\varepsilon<1\) 时，语义 (R) 的坏节点比例是 1。另一方面，对任意
\(\|z\|_2\le1\)，

\[
\#\{i:|z_i|>\varepsilon\}
\le\frac1{\varepsilon^2},
\]

故语义 (I) 的坏节点比例至多 \(1/(n\varepsilon^2)\to0\)。这证明
\(\sup_z Q_i(\|E_iz\|)\) 与 \(Q_i(\sup_z\|E_iz\|)\) 不能交换。

### 反例 9.6（block/path cancellation 阻止一般必要 walk 条件）

**状态：[标准线性代数警告；未提出新 theorem]。**

在一般 block-SPD 或 raw WLS 中，矩阵路径乘积有符号和旋转，(5.4) 的两个误差项也可相消。仅从
\(\|A\|\)、\(\|B\|\) 不可能推出 \(\|A+B\|\) 的正下界：取
\(B=-A\) 即可。因而用 block norm 构造的正 walk majorant 至多给充分条件；除非增加正锥、M-matrix 符号、range/singular-vector 对齐或局部右逆，不能反推真实 inverse/WLS 行尾的必要条件。

---

## 10. 证明审计清单与尚未闭合之处

### 10.1 已闭合

- **范数：** 第 2–5、7–9 节默认 Hilbert \(2\to2\)；第 6 节明确切换为标量
  \(\ell_\infty\to|\cdot|\)，其对偶为 \(\ell_1\) 行和；第 6.5 节另列
  \(\ell_2\) 版本。
- **投影：** 本地数据抽取 \(S\) 是 coisometry，\(S^*S\) 是坐标投影；
  WLS 中状态投影和测量投影分开记为 \(P^x,P^z\)。
- **Schur 方向：** \(E:\mathsf X_O\to\mathsf X_B\)，
  \(F=E^*A^{-1}:\mathsf X_B\to\mathsf X_O\)，
  \(\Sigma=D-E^*A^{-1}E\)，根球外行
  \(Y=-R_iA^{-1}E\Sigma^{-1}\)；(4.8) 的左右维数一致。
- **达到性：** fixed-model 线性截断直接达到；共同线性 minimax 只有在本地空间有限维且 completion 算子统一有界时保证达到，否则陈述保留任意
  \(\eta>0\)。
- **量词：** 数据不属于模型；算法选择在模型/completion supremum 之前；节点鲁棒 supremum 在节点分位数内部。逐 view worst completion 与单实例分位数没有被交换。
- **halo：** principal \(A\) 与 cut \(E\) 必须由闭合系数 view 决定；因子所有权不满足时明确增加固定 halo。
- **raw WLS：** \(K=J^{-1}G\) 始终保留 \(G=H^*R^{-1}\)；双向 locality 只在
  \(G\) 局部且有局部受控右逆时成立。
- **grounded walk：** 输入是 \(c=S^{-1}b\)；若改回 \(b\)，reward 必须乘
  \(1/s_j\)。
- **\(\gamma\) sharpness：** 对单 completion 使用 global
  \(\gamma_{\rm eff}\)；对任意 SPD Schur completion，局部常数
  \(\sqrt{1+\|F\|^2}\) 对 Dirichlet/oracle 比值渐近可达。没有声称它是所有共同算法 minimax 的 sharp factor。
- **escape：** 精确对象是 \(h_i\ell_i^{\rm rem}\)；escape-only iff 明确附带条件余寿命有界或 tightness。

### 10.2 仍有缺口

1. **locality-constrained structured robust inverse 的精确解。**
   一般 uncertain inverse 的 common minimax 与 structured-SDP 框架已由
   El Ghaoui (2002) 解决；(4.14f) 作为 robust formulation/SDP 方向不再构成缺口或创新。真正未闭合的是：给 approximant 加
   \([Q,0]\) 的局部列约束、并把 uncertainty 限为 SPD Schur 谱窗后，
   能否得到 exact finite SDP、dual certificate、最坏 \(\Sigma\) 和一般
   block 闭式。其 Theorem 6.2 对 structured uncertainty 给的是上界，
   Eq. (6.5) 的 exactness 只属于 unstructured case。
2. **主定理原创性。** common-rule + Schur + quantile + no-go 的完全相同组合尚未找到，但 common approximate inverse 这一核心优化量词已有明确先例，现有其余证明又大多是标准恒等式的短组合，不能仅以“未找到同句 theorem”推断足够原创。
3. **可局部验证的 exterior stability。** \(\gamma_i\) 和 escape \(h_i\) 可由 closed view 算，但
   \(g_i\)、Schur lower bound、\(\ell_i^{\rm rem}\) 都含 exterior。尚未得到既纯局部可抽样、又对广泛 completion 类必要充分的替代；定理 7.3 表明对任意 tentacle-closed 类这种替代不可能。
4. **block 概率化。** 没有共同正锥时，matrix path cancellation 使 scalar killed-walk iff 不能推广成一般 block-SPD 的必要条件。
5. **计算与能量。** 本文只刻画信息半径。收集球后 dense solve 的存储/计算可达
   \(O(m_r^2)\)/\(O(m_r^3)\)，且球体积可能指数增长；尚无 matching rounds–bits–energy–flops Pareto theorem。
6. **Bayes 与 adversarial 的统一。** (5.12) 精确，但只有协方差谱窗下才与 coefficient tail 双向等价。随机 completion/Benjamini–Schramm 极限下的多数节点定理未证明。

### 10.3 经典来源定位

以下只用于标明哪些成分已有经典来源，不构成 novelty 声明：

1. Demko–Moss–Smith, *Decay Rates for Inverses of Band Matrices*, Math. Comp. 43 (1984), 491–499，Propositions 2.1–2.2、Theorem 2.4；DOI
   [10.1090/S0025-5718-1984-0758197-9](https://doi.org/10.1090/S0025-5718-1984-0758197-9)。本地：
   literature/03_sparse_inverse_graph_filters/1984_demko_moss_smith_decay_inverse_band_matrices.pdf。
2. Malioutov–Johnson–Willsky, *Walk-Sums and Belief Propagation in Gaussian Graphical Models*, JMLR 7 (2006), 2031–2064，尤其 Proposition 1 与 Proposition 21。本地：
   literature/02_graphical_models_local_inference/2006_malioutov_johnson_willsky_walk_sums_gaussian_bp.pdf。
3. Bendito–Carmona–Encinas, *Potential Theory for Schrödinger Operators on Finite Networks*, Rev. Mat. Iberoam. 21 (2005), 771–818；DOI
   [10.4171/RMI/435](https://doi.org/10.4171/RMI/435)。
4. Carmona–Encinas–Jiménez–Martín, *Random Walks Associated with Symmetric M-matrices*, Linear Algebra Appl. 693 (2024), 324–338；DOI
   [10.1016/j.laa.2023.10.009](https://doi.org/10.1016/j.laa.2023.10.009)。
5. Cheng–Jiang–Sun, *Spatially Distributed Sampling and Reconstruction*, Appl. Comput. Harmon. Anal. 47 (2019), 109–148，Theorems 5.3, 6.1–6.2 与 Proposition 7.1；DOI
   [10.1016/j.acha.2017.07.007](https://doi.org/10.1016/j.acha.2017.07.007)。本地：
   literature/03_sparse_inverse_graph_filters/2015_cheng_jiang_sun_spatially_distributed_sampling_reconstruction.pdf。
6. Laurent El Ghaoui, *Inversion Error, Condition Number, and Approximate Inverses of Uncertain Matrices*, Linear Algebra Appl. 343–344 (2002), 171–193。Eq. (1.4) 定义 common approximate inverse；Theorem 6.2/Eq. (6.3) 给 structured uncertainty 的 SDP 上界；Eq. (6.5) 给 unstructured case 的解析最优 inverse；§7.4（article p. 186）解释 additive case 与 total least squares 的关系。DOI
   [10.1016/S0024-3795(01)00273-7](https://doi.org/10.1016/S0024-3795(01)00273-7)，[作者页 PDF](https://people.eecs.berkeley.edu/~elghaoui/Pubs/InvErr_LAA02.pdf)。当前未在规范 literature 目录归档。

Schur 边界公式应视为离散 Dirichlet-to-Neumann/透明边界和 geometric resolvent identity 的近亲；fixed-model (2.1) 是经典 radius-of-information 类型结论；第 6 节的 Green/occupation 公式属于标准 potential theory；common uncertain inverse minimax 则已由 El Ghaoui 明确定义并系统处理。这四条碰撞边界均已确认。
