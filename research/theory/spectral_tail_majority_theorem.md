# 多数节点上的谱尾充要条件：共同多项式、Richardson 与 WLS 输入几何

> 数学与文献核验日期：2026-09-13  
> 结论状态：第 3 节的多数节点谱尾 iff 严格成立，但核心是谱定理、低谱截断和 Richardson 残差的短推论，**不建议单独宣称为原创主定理**。第 10 节进一步严格闭合了 \(\mathbb Z^d\) 上的任意局部算法下界、Christoffel 最佳多项式上界与临界维数；这是明显更强的主定理骨架，但因 Green 渐近、Christoffel 渐近和 GFF/membrane 临界维数分别都是经典成分，仍须继续查重与推广后才能作原创性声明。

## 1. 先给结论

对一个没有统一谱下界的 SPD 网络矩阵族，决定“一个共同有限轮多项式能否在至少 \(1-\delta\) 的节点上逼近逆”的，不是无环性，也不是全局最小特征值，而是小特征值在**各节点谱测度中携带多少逆矩能量**。

更准确地，令估计目标为

\[
x=J^{-1}G\xi,
\]

并假设输入几何满足统一 Loewner 比较

\[
c_-J^\nu\preceq GG^T\preceq c_+J^\nu,
\qquad 0<c_-\le c_+<\infty,
\qquad 0\le \nu\le2.
\tag{1.1}
\]

记

\[
s:=2-\nu\in[0,2].
\tag{1.2}
\]

那么，对整个模型族共同的有限次多项式存在，当且仅当对任意
\(a>0,0<\delta<1\)，都能统一选一个低谱阈值 \(\eta>0\)，使至少 \(1-\delta\) 的节点满足

\[
\int_{(0,\eta)}\lambda^{-s}\,d\mu_{J,i}(\lambda)\le a.
\tag{1.3}
\]

而且不需要另找某个神秘多项式：最简单的 Richardson/Neumann 多项式

\[
p_t(\lambda)=\frac{1-(1-\lambda/M)^t}{\lambda}
=\frac1M\sum_{k=0}^{t-1}(1-\lambda/M)^k
\tag{1.4}
\]

已经在定性上达到这个充要条件。

两个最重要的特例不能混为一谈：

| 输入与目标 | \(G G^T\) | \(s\) | 正确的低谱尾 |
|---|---:|---:|---:|
| 原始正规方程右端 \(Jx=b,\ \|b\|_2\le1\) | \(I=J^0\) | \(2\) | \(\int\lambda^{-2}d\mu_{J,i}\) |
| whitened WLS 数据 \(x=J^{-1}G\xi,\ \|\xi\|_2\le1\) | \(J\) | \(1\) | \(\int\lambda^{-1}d\mu_{J,i}\) |

因此，把 RHS 模型的 \(\lambda^{-2}\) 条件直接搬到 WLS measurement-input 模型，会多要求一整个逆幂，是实质性的范数/输入几何错误。

在 exact whitened WLS 中，必要性还可越过 polynomial 假设：只要估计器写成
\(K G\xi\) 且 \(K\) 在模型族上一致稳定，即使 \(K\) 非对称、与 \(J\) 不交换并依赖实例，也必须满足定理 3.5 的低谱投影不等式。另一方面，在
\(\mathbb Z^d\) 上可以得到真正的 matching round theorem：任意半径
\(r\) 的局部解码与最佳 degree-\(r\) 多项式的误差都为
\(\Theta(r^{-(d-2q)/2})\)，而 Richardson 只有
\(\Theta(r^{-(d-2q)/4})\)。这给出明确的 latency/标准同步通信成本收益，但该格点组合结论的原创性目前仍标为 PENDING。

---

## 2. 完整模型与量词

### 2.1 模型族

令 \(\mathfrak F\) 是一族有限模型。每个 \(\theta\in\mathfrak F\) 包含：

- 有限状态节点集 \(V_\theta\)，标量状态总维数 \(n_\theta=|V_\theta|\)；
- 实对称正定矩阵 \(J_\theta\in\mathbb R^{n_\theta\times n_\theta}\)；
- 输入算子 \(G_\theta\in\mathbb R^{n_\theta\times m_\theta}\)；
- 统一谱上界
  \[
  0\prec J_\theta\preceq MI
  \qquad(\theta\in\mathfrak F),
  \tag{2.1}
  \]
  其中 \(M<\infty\) 与 \(\theta\) 无关；
- 统一输入几何 (1.1)。

不假设统一正下界 \(J_\theta\succeq mI\)。每个有限矩阵自身可逆，但

\[
\inf_{\theta\in\mathfrak F}\lambda_{\min}(J_\theta)
\]

可以为零。

### 2.2 节点谱测度

对 \(i\in V_\theta\)，定义概率测度

\[
\mu_{\theta,i}(B)
:=e_i^T\mathbf 1_B(J_\theta)e_i,
\qquad B\subset(0,M]\text{ Borel}.
\tag{2.2}
\]

等价地，若 \(J_\theta u_k=\lambda_k u_k\) 且 \(\{u_k\}\) 正交归一，则

\[
\mu_{\theta,i}
=\sum_{k=1}^{n_\theta}|u_k(i)|^2\delta_{\lambda_k}.
\tag{2.3}
\]

低谱尾记为

\[
\tau_{\theta,i}^{(s)}(\eta)
:=\int_{(0,\eta)}\lambda^{-s}\,d\mu_{\theta,i}(\lambda),
\qquad 0<\eta<M.
\tag{2.4}
\]

### 2.3 多数节点与“对整个族统一”

对任意节点量 \(Z_{\theta,i}\)，定义最坏模型坏节点比例

\[
\mathsf{Bad}_{\mathfrak F}(Z>a)
:=\sup_{\theta\in\mathfrak F}
\frac1{|V_\theta|}
\big|\{i\in V_\theta:Z_{\theta,i}>a\}\big|.
\tag{2.5}
\]

本文中“至少 \(1-\delta\) 节点”始终是

\[
\mathsf{Bad}_{\mathfrak F}(Z>a)\le\delta.
\tag{2.6}
\]

这比只对每个单独矩阵各选参数更强：多项式 \(p\)、阈值 \(\eta\) 和迭代次数 \(t\) 都不能依赖于具体 \(\theta\)。它们可以依赖于 \(\varepsilon,\delta,M,c_\pm,\nu\) 以及整个模型类的先验。

若研究的是一个渐近序列并希望忽略任意有限前缀，可把 (2.5) 换成

\[
\limsup_{N\to\infty}
\sup_{\substack{\theta\in\mathfrak F\\|V_\theta|\ge N}}
\frac{|\{i:Z_{\theta,i}>a\}|}{|V_\theta|}.
\tag{2.7}
\]

下面全部证明逐字适用于这个渐近版本。正文采用更强的 class-uniform 定义 (2.5)。

### 2.4 误差量

对实多项式 \(p\)，节点 \(i\) 的单位输入最坏误差为

\[
\mathcal E_{\theta,i}(p)
:=\left\|
e_i^T\big(J_\theta^{-1}-p(J_\theta)\big)G_\theta
\right\|_2.
\tag{2.8}
\]

它等于

\[
\sup_{\|\xi\|_2\le1}
\left|
e_i^T\big(J_\theta^{-1}-p(J_\theta)\big)G_\theta\xi
\right|.
\tag{2.9}
\]

如果 \(\xi\sim N(0,I)\)，它也正好是该节点均方误差的平方根。这个随机解释依赖 whitened covariance；不能把任意 raw measurement 的欧氏单位球与它混同。

---

## 3. 一个输入感知的多数节点充要定理

### 引理 3.1（精确谱误差与上下夹逼）

令

\[
r_p(\lambda):=1-\lambda p(\lambda).
\tag{3.1}
\]

则对每个 \(\theta,i,p\)，

\[
c_-\int_{(0,M]}\lambda^{-s}|r_p(\lambda)|^2
\,d\mu_{\theta,i}(\lambda)
\le
\mathcal E_{\theta,i}(p)^2
\le
c_+\int_{(0,M]}\lambda^{-s}|r_p(\lambda)|^2
\,d\mu_{\theta,i}(\lambda).
\tag{3.2}
\]

若 \(GG^T=J^\nu\)，两边取等号且 \(c_-=c_+=1\)。

#### 证明

记 \(f_p(J)=J^{-1}-p(J)=J^{-1}r_p(J)\)。由 (1.1) 作合同变换，

\[
c_- f_p(J)J^\nu f_p(J)
\preceq
f_p(J)GG^Tf_p(J)
\preceq
c_+ f_p(J)J^\nu f_p(J).
\tag{3.3}
\]

因为 \(f_p(J)\) 与 \(J^\nu\) 都是 \(J\) 的函数，

\[
f_p(J)J^\nu f_p(J)
=J^{\nu-2}r_p(J)^2=J^{-s}r_p(J)^2.
\tag{3.4}
\]

以 \(e_i\) 左右夹取二次型，再用 (2.2) 的谱积分表示即得。证毕。

### 定理 3.2（共同多项式、低谱尾与 Richardson 的等价性）

在第 2 节假设下，下列三件事等价。

**(A) 共同有限多项式可在多数节点逼近。** 对任意
\(\varepsilon>0\) 和 \(0<\delta<1\)，存在一个实多项式
\(p=p_{\varepsilon,\delta}\)，使

\[
\mathsf{Bad}_{\mathfrak F}
\big(\mathcal E(p)>\varepsilon\big)\le\delta.
\tag{3.5}
\]

**(B) 小特征值的节点逆矩尾统一消失。** 对任意 \(a>0\) 和
\(0<\delta<1\)，存在 \(0<\eta<M\)，使

\[
\mathsf{Bad}_{\mathfrak F}
\big(\tau^{(s)}(\eta)>a\big)\le\delta.
\tag{3.6}
\]

等价地，对每个固定 \(a>0\)，

\[
\lim_{\eta\downarrow0}
\mathsf{Bad}_{\mathfrak F}
\big(\tau^{(s)}(\eta)>a\big)=0.
\tag{3.7}
\]

**(C) 固定的 Richardson 族已经足够。** 对任意
\(\varepsilon>0\) 和 \(0<\delta<1\)，存在整数 \(t\ge1\)，使 (1.4) 满足

\[
\mathsf{Bad}_{\mathfrak F}
\big(\mathcal E(p_t)>\varepsilon\big)\le\delta.
\tag{3.8}
\]

因此，在这个共同多项式模型内，(B) 是最宽的定性充要条件；无环、逐节点局部测量满列秩、统一条件数等都不是必要条件。

#### 证明：\((C)\Rightarrow(A)\)

\(p_t\) 是次数至多 \(t-1\) 的实多项式，所以直接成立。

#### 证明：\((A)\Rightarrow(B)\)

固定 \(a>0\) 和 \(0<\delta<1\)。在 (A) 中取

\[
\varepsilon_0:=\frac12\sqrt{a c_-}.
\tag{3.9}
\]

得到一个共同多项式 \(p\)，使至少 \(1-\delta\) 的节点满足
\(\mathcal E_{\theta,i}(p)\le\varepsilon_0\)。令

\[
B_p:=\max_{0\le\lambda\le M}|p(\lambda)|<\infty.
\tag{3.10}
\]

若 \(B_p=0\)，任取 \(\eta<M\)；否则取

\[
0<\eta\le\min\left\{\frac M2,\frac1{2B_p}\right\}.
\tag{3.11}
\]

于是对 \(0<\lambda<\eta\)，

\[
|r_p(\lambda)|=|1-\lambda p(\lambda)|\ge\frac12.
\tag{3.12}
\]

由引理 3.1，任一上述好节点满足

\[
\varepsilon_0^2
\ge \mathcal E_{\theta,i}(p)^2
\ge\frac{c_-}{4}
\tau_{\theta,i}^{(s)}(\eta),
\tag{3.13}
\]

故 \(\tau_{\theta,i}^{(s)}(\eta)\le a\)。同一个 \(\eta\) 对整个族有效，得到 (B)。

#### 证明：\((B)\Rightarrow(C)\)

固定 \(\varepsilon>0,0<\delta<1\)，再固定任意 \(0<\vartheta<1\)。在 (B) 中取

\[
a:=\frac{\vartheta\varepsilon^2}{c_+},
\tag{3.14}
\]

得到共同阈值 \(\eta\in(0,M)\)。Richardson 残差恰为

\[
r_{p_t}(\lambda)=(1-\lambda/M)^t.
\tag{3.15}
\]

对每个满足 \(\tau_{\theta,i}^{(s)}(\eta)\le a\) 的节点，把谱积分在
\((0,\eta)\) 与 \([\eta,M]\) 分开。因 \(s\ge0\)，有

\[
\begin{aligned}
\mathcal E_{\theta,i}(p_t)^2
&\le c_+
\int_{(0,M]}\lambda^{-s}(1-\lambda/M)^{2t}
\,d\mu_{\theta,i}(\lambda)\\
&\le c_+\left[
a+\eta^{-s}(1-\eta/M)^{2t}
\right].
\end{aligned}
\tag{3.16}
\]

选 \(t\) 使第二项不超过
\((1-\vartheta)\varepsilon^2/c_+\)，便有
\(\mathcal E_{\theta,i}(p_t)\le\varepsilon\)。原来至少
\(1-\delta\) 的节点满足尾条件，所以 (C) 成立。证毕。

### 推论 3.3（显式迭代次数）

若对至少 \(1-\delta\) 的节点已经有

\[
\tau_{\theta,i}^{(s)}(\eta)
\le\frac{\vartheta\varepsilon^2}{c_+},
\tag{3.17}
\]

则以下整数足够：

\[
t\ge
\max\left\{
1,
\left\lceil
\frac{
\log_+\!\left(
\dfrac{c_+\eta^{-s}}
{(1-\vartheta)\varepsilon^2}
\right)}
{-2\log(1-\eta/M)}
\right\rceil
\right\},
\tag{3.18}
\]

其中 \(\log_+(x)=\max\{0,\log x\}\)。较简单但稍松的充分界为

\[
t\ge
\max\left\{
1,
\left\lceil
\frac{M}{2\eta}
\log_+\!\left(
\frac{c_+}{(1-\vartheta)\varepsilon^2\eta^s}
\right)
\right\rceil
\right\}.
\tag{3.19}
\]

特别取 \(\vartheta=1/2\)：

- raw RHS（\(s=2,c_+=1\)）需要
  \[
  t\ge\left\lceil\frac{M}{2\eta}
  \log_+\frac{2}{\varepsilon^2\eta^2}\right\rceil;
  \tag{3.20}
  \]
- exact whitened WLS（\(s=1,c_+=1\)）需要
  \[
  t\ge\left\lceil\frac{M}{2\eta}
  \log_+\frac{2}{\varepsilon^2\eta}\right\rceil.
  \tag{3.21}
  \]

两式都还要与 \(1\) 取最大值。

### 推论 3.4（双向的坏节点定量比较）

对任意 \(0<\eta<M\) 和整数 \(t\ge1\)，若

\[
b_t(\eta,\varepsilon)
:=\frac{\varepsilon^2}{c_+}
-\eta^{-s}(1-\eta/M)^{2t}>0,
\tag{3.22}
\]

则

\[
\mathsf{Bad}_{\mathfrak F}
(\mathcal E(p_t)>\varepsilon)
\le
\mathsf{Bad}_{\mathfrak F}
(\tau^{(s)}(\eta)>b_t(\eta,\varepsilon)).
\tag{3.23}
\]

反过来，若一个多项式 \(p\) 在 \((0,\eta)\) 上满足

\[
\sup_{0<\lambda<\eta}|\lambda p(\lambda)|\le\rho<1,
\tag{3.24}
\]

则

\[
\mathsf{Bad}_{\mathfrak F}
\left(
\tau^{(s)}(\eta)>
\frac{\varepsilon^2}{c_-(1-\rho)^2}
\right)
\le
\mathsf{Bad}_{\mathfrak F}
(\mathcal E(p)>\varepsilon).
\tag{3.25}
\]

这两式是定理 3.2 的定量核心。它们也说明：只给“(B) 成立”而不给尾函数随 \(\eta\) 的模量，就不可能推出一个只依赖
\((M,\varepsilon,\delta)\) 的统一轮数率；任意慢的尾消失都会产生任意慢的共同收敛。

### 定理 3.5（WLS 中任意有界、非交换局部算子的必要性）

多项式假设可以在 exact whitened WLS 中大幅放宽。本定理只假设

\[
G_\theta G_\theta^T=J_\theta.
\tag{3.26}
\]

对每个模型，允许一个任意矩阵
\(K_\theta\in\mathbb R^{n_\theta\times n_\theta}\)；它可以非对称，可以不与
\(J_\theta\) 交换，也可以随实例变化。估计器只要先形成 WLS sufficient statistic \(G_\theta\xi\)，再输出

\[
\widehat x=K_\theta G_\theta\xi.
\tag{3.27}
\]

其节点误差记为

\[
\mathcal E_{\theta,i}[K]
:=
\left\|e_i^T(J_\theta^{-1}-K_\theta)G_\theta\right\|_2.
\tag{3.28}
\]

则对任意 \(0<\eta<M\)，只要
\(\|K_\theta\|_{2\to2}\le B\)，就有逐节点的必要不等式

\[
\boxed{
\sqrt{\tau_{\theta,i}^{(1)}(\eta)}
\le
\mathcal E_{\theta,i}[K]+B\sqrt\eta.}
\tag{3.29}
\]

因此，假设 \(J_\theta,G_\theta\) 都有统一有限传播范围。以下两件事等价：

**(A\(_{\rm nc}\))** 对每个 \(\varepsilon>0,0<\delta<1\)，存在有限常数
\(B_{\varepsilon,\delta}\)、有限传播半径
\(R_{\varepsilon,\delta}\)，以及一族可依赖实例的矩阵
\(\{K_\theta\}\)，使

\[
\sup_\theta\|K_\theta\|\le B_{\varepsilon,\delta},
\qquad
\operatorname{prop}(K_\theta)\le R_{\varepsilon,\delta},
\tag{3.30}
\]

且

\[
\mathsf{Bad}_{\mathfrak F}
\big(\mathcal E[K]>\varepsilon\big)
\le\delta.
\tag{3.31}
\]

**(B\(_1\))** 定理 3.2(B) 在 \(s=1\) 时成立，即对任意
\(a>0,0<\delta<1\)，存在共同 \(\eta>0\)，使

\[
\mathsf{Bad}_{\mathfrak F}
\left(
\int_{(0,\eta)}\lambda^{-1}d\mu_{\theta,i}(\lambda)>a
\right)
\le\delta.
\tag{3.32}
\]

换句话说，在“先形成 \(G\xi\)，再施加一个稳定的有限范围状态算子”这个大类中，\(\lambda^{-1}\) 多数节点谱尾仍是充要条件；这个必要性不要求
\(K_\theta=p(J_\theta)\)。
事实上，必要方向 (3.29) 连有限范围都没有使用：它对任意全局但统一稳定的
\(K_\theta\) 也成立。有限范围只是在充分方向把存在性解释成有限通信轮数所需。

#### 证明

固定 \(\theta\)，省略下标。令

\[
U:=J^{-1/2}G,
\qquad UU^T=I,
\tag{3.33}
\]

并记

\[
P_\eta:=\mathbf1_{(0,\eta)}(J),
\qquad
Q_\eta:=U^TP_\eta U.
\tag{3.34}
\]

由 \(UU^T=I\) 可验证 \(Q_\eta^2=Q_\eta=Q_\eta^T\)，所以它是输入空间上的正交投影。对真实 WLS 算子 \(A:=J^{-1}G=J^{-1/2}U\)，

\[
AQ_\eta=J^{-1/2}P_\eta U.
\tag{3.35}
\]

因为右乘 \(U\) 保持状态行范数，

\[
\|e_i^TAQ_\eta\|_2^2
=e_i^TP_\eta J^{-1}P_\eta e_i
=\tau_i^{(1)}(\eta).
\tag{3.36}
\]

另一方面，由 \(G=J^{1/2}U\)，

\[
KGQ_\eta=KJ^{1/2}P_\eta U,
\tag{3.37}
\]

所以

\[
\|e_i^TKGQ_\eta\|_2
\le \|K\|\,\|J^{1/2}P_\eta\|
\le B\sqrt\eta.
\tag{3.38}
\]

再用 \(Q_\eta\) 是投影，

\[
\|e_i^T(A-KG)Q_\eta\|_2
\le\mathcal E_i[K].
\tag{3.39}
\]

对 (3.36) 用三角不等式并代入 (3.38)–(3.39)，得 (3.29)。

现证等价。若 (A\(_{\rm nc}\)) 成立，给定 \(a,\delta\)，在其中取
\(\varepsilon=\sqrt a/4\)，得到某个共同有限 \(B\)。再取

\[
0<\eta<\min\left\{M,\frac{a}{16B^2}\right\}
\tag{3.40}
\]

（\(B=0\) 时忽略第二项）。每个误差至多
\(\sqrt a/4\) 的节点都由 (3.29) 满足
\(\tau_i^{(1)}(\eta)\le a/4<a\)，因而得到
(B\(_1\))。

反过来，若 (B\(_1\)) 成立，定理 3.2 的 Richardson 方向直接给出

\[
K_\theta=p_t(J_\theta),
\qquad
\|K_\theta\|
\le \max_{0\le\lambda\le M}p_t(\lambda)
=\frac tM,
\tag{3.41}
\]

且传播半径至多 \((t-1)r_J\)。它满足
(A\(_{\rm nc}\))。证毕。

#### 边界：哪些“任意 local linear estimator”尚未包含

定理 3.5 允许 \(K\) 非交换、node-varying 且 instance-dependent，但它仍要求估计器因子分解为
\(KG\)。若允许任意 measurement-to-state 算子
\(C_\theta:\mathbb R^{m_\theta}\to\mathbb R^{n_\theta}\)，而不要求它只通过 sufficient statistic \(G\xi\) 作用，上述证明不再给出
\(B\sqrt\eta\)：对低谱输入投影只能粗略估计
\(\|C Q_\eta\|\le\|C\|\)。局部性是否能在充分一般的 holder 模型中恢复一个足够强的必要不等式，目前是 **PENDING**，不应写成已解决。

另外，全族统一有界性不能删掉。它表示噪声增益/数值稳定性受控；如果允许
\(\|K_\theta\|\to\infty\)，例如对角 \(J_\theta\) 上的零通信逐点除法就会绕开 (3.29) 的统一控制。

### 推论 3.6（误差–稳定性预算的必要关系）

定理 3.5 还给出直接可量化的坏节点包含：对任意
\(B,\eta,\varepsilon>0\)，

\[
\left\{i:
\tau_{\theta,i}^{(1)}(\eta)>
(\varepsilon+B\sqrt\eta)^2
\right\}
\subseteq
\left\{i:\mathcal E_{\theta,i}[K]>\varepsilon\right\}.
\tag{3.42}
\]

特别地，若某组节点在一段低谱上有下界

\[
\tau_{\theta,i}^{(1)}(\eta)
\ge c\eta^\kappa,
\qquad 0<\eta\le\eta_0,
\tag{3.43}
\]

则任意 \(\|K\|\le B\) 的 sufficient-statistic 估计器都必须满足

\[
\mathcal E_{\theta,i}[K]
\ge
\sup_{0<\eta\le\eta_0}
\left[
\sqrt c\,\eta^{\kappa/2}-B\sqrt\eta
\right]_+.
\tag{3.44}
\]

当 \(0<\kappa<1\) 时还能显式优化这个 supremum。若

\[
B\ge \kappa\sqrt c\,
\eta_0^{-(1-\kappa)/2},
\tag{3.45}
\]

则驻点落在允许区间内，并给出

\[
\mathcal E_{\theta,i}[K]
\ge
(1-\kappa)\kappa^{\kappa/(1-\kappa)}
c^{1/[2(1-\kappa)]}
B^{-\kappa/(1-\kappa)}.
\tag{3.46}
\]

因此在这个区间里，要达到 \(\mathcal E_{\theta,i}[K]\le\varepsilon\)，稳定性预算至少按
\(B=\Omega(\varepsilon^{-(1-\kappa)/\kappa})\) 增长。这是可量化的
noise-gain/系数幅值代价；例如三维格点 WLS 的
\(\kappa=(d-2)/2=1/2\)，故有 \(B=\Omega(\varepsilon^{-1})\)。它仍不等于轮数或通信能耗。

这是一个噪声增益/稳定性预算下界，但它本身不是通信轮数下界：一个半径为 \(r\) 的算子可以有非常大的 \(B\)，而仅有低谱质量也不告诉我们 Green 行在半径 \(r\) 外的空间能量。要得到 matching round lower bound，仍需要第 7 节或第 10 节那样的 spatial indistinguishability/Green-tail 论证。

---

## 4. WLS 特例：为什么指数从 2 变成 1

### 4.1 whitened measurement input

考虑线性测量

\[
z=Hx+w,\qquad \operatorname{Cov}(w)=R\succ0,
\tag{4.1}
\]

且全局信息矩阵

\[
J=H^TR^{-1}H\succ0.
\tag{4.2}
\]

令 whitened 数据

\[
\xi:=R^{-1/2}z,
\qquad
G:=H^TR^{-1/2}.
\tag{4.3}
\]

这里更精确地说，\(R^{-1/2}w\) 的协方差是 \(I\)；完整的
\(\xi=R^{-1/2}z\) 还含有确定性（或另有先验分布的）信号均值，不能仅凭
“whitened” 就宣称 \(\xi\sim N(0,I)\)。下面的单位球结论是确定性的；白噪声 MSE 解释只施加在噪声/扰动分量上。

则 WLS 解为

\[
\widehat x=J^{-1}G\xi,
\tag{4.4}
\]

并且

\[
GG^T=H^TR^{-1}H=J.
\tag{4.5}
\]

因此，对任意多项式 \(p\)，

\[
\begin{aligned}
\left\|e_i^T(J^{-1}-p(J))G\right\|_2^2
&=e_i^T(J^{-1}-p(J))J(J^{-1}-p(J))e_i\\
&=\int_{(0,M]}
\lambda^{-1}|1-\lambda p(\lambda)|^2
\,d\mu_{J,i}(\lambda).
\end{aligned}
\tag{4.6}
\]

这就是 \(\lambda^{-1}\) 而不是 \(\lambda^{-2}\) 的来源。

### 4.2 三种范数语义必须分开

1. **Whitened 单位球。** \(\|\xi\|_2\le1\) 等价于 raw 数据满足
   \(z^TR^{-1}z\le1\)。此时 (4.6) 是精确 deterministic operator error。
2. **白噪声均方误差。** 若算法误差所作用的扰动分量
   \(R^{-1/2}w\sim N(0,I)\)，(4.6) 也是该节点由白噪声引起的误差方差；它不是对含信号均值的完整 \(\xi\) 的无条件协方差陈述。
3. **Raw Euclidean 单位球。** 若直接规定 \(\|z\|_2\le1\)，输入算子是
   \(G_{\rm raw}=H^TR^{-1}\)，而
   \(G_{\rm raw}G_{\rm raw}^T=H^TR^{-2}H\)，通常不等于 \(J\)。若有统一
   \(r_-I\preceq R\preceq r_+I\)，则
   \[
   r_+^{-1}J
   \preceq G_{\rm raw}G_{\rm raw}^T
   \preceq r_-^{-1}J,
   \tag{4.7}
   \]
   所以仍可用 \(\nu=1\)，但常数变为
   \(c_-=r_+^{-1},c_+=r_-^{-1}\)。没有统一噪声谱界时，不能自动使用同一结论。

### 4.3 先验也可以并入测量

若 WLS/MAP 目标还含有正定二次先验，只要把先验写成额外的虚拟测量行并一并 whiten，就仍有某个扩展 \(G\) 满足
\(GG^T=J\)。若只把先验加到 \(J\) 而不加入输入算子，则 (4.5) 不再成立，应回到 (1.1) 检查 Loewner 比较，而不能靠记号猜测指数。

---

## 5. 多项式次数何时真等于通信半径

定理 3.2 本身只用谱定理，不需要图。图的作用是把多项式次数解释成局部通信轮数。

设状态图距离为 \(d_V\)，且

\[
(J_\theta)_{ij}=0
\quad\text{whenever}\quad
d_V(i,j)>r_J
\tag{5.1}
\]

对所有模型有共同 \(r_J<\infty\)。再设每个输入坐标 \(a\) 有一个 holder
\(h(a)\)，并且

\[
(G_\theta)_{ja}=0
\quad\text{whenever}\quad
d(h(a),j)>r_G.
\tag{5.2}
\]

那么次数 \(k\) 的 \(p(J)G\) 只有在

\[
d(h(a),i)\le r_G+k r_J
\tag{5.3}
\]

时才可能有非零 \((i,a)\) 元素。换言之，它可由局部初始聚合加 \(k\) 次
\(J\)-matvec 实现。

Richardson 的 \(p_t\) 次数为 \(t-1\)。从 \(x_0=0\) 写成

\[
x_{\ell+1}=x_\ell+M^{-1}(G\xi-Jx_\ell)
\tag{5.4}
\]

时，字面上进行了 \(t\) 次 update；第一次 \(Jx_0=0\) 不携带信息，可消去。最终 \(x_t=p_t(J)G\xi\) 的数据依赖半径至多

\[
r_G+(t-1)r_J.
\tag{5.5}
\]

必须同时写明以下限制：

- holder 图与状态稀疏图可能不是同一个图；
- \(H\) 局部不自动推出任意 chosen \(J\)-matvec 在同一网络上一轮完成，需检查
  \(J=H^TR^{-1}H\) 的实际 support；
- 消息位数、量化、异步和丢包不包含在次数结论内；
- 一个一般的本地算法可读取局部系数并使用 node-varying 更新，它远比“全族共享同一个标量多项式 \(p\)”更宽。

最后一点尤其重要。比如对任意正对角 \(J\)，节点知道自己的 \(J_{ii}\) 时可零通信精确输出 \(b_i/J_{ii}\)，即使不存在一个共同低次多项式在全族所有极小对角值上逼近 \(1/\lambda\)。所以定理 3.2 **不是所有 LOCAL opcode 的充要条件**；它是共同 polynomial/Richardson filter 类的精确充要条件。

---

## 6. localized 与 delocalized 小特征值：相同 \(\lambda_{\min}\)，相反结论

下面两个族都有 \(\lambda_{\min}=1/n\)，却表现完全不同。这说明全局最小特征值或 condition number 不能刻画多数节点可局部估计性。

### 6.1 一个局部化坏模态：多数节点可以完全忽略

令

\[
J_n^{\rm loc}=\operatorname{diag}(n^{-1},1,\ldots,1)
\in\mathbb R^{n\times n}.
\tag{6.1}
\]

则

\[
\mu_{n,1}=\delta_{1/n},
\qquad
\mu_{n,i}=\delta_1\quad(i\ge2).
\tag{6.2}
\]

对 raw RHS，节点 1 的低谱逆二阶矩可达 \(n^2\)，但它只占 \(1/n\)；其余节点对任意 \(\eta<1\) 的尾都为零。严格的 class-uniform 量词也成立：给定 \(\delta\)，只需对有限多个
\(n<1/\delta\) 选 \(\eta\) 小于它们的 \(1/n\)；更大的矩阵允许把唯一坏节点计入失败集合。

因此 (B) 成立，继而 (A)、(C) 成立，尽管

\[
\inf_n\lambda_{\min}(J_n^{\rm loc})=0.
\tag{6.3}
\]

这个例子还表明 (B) 严格弱于平均谱测度的 uniform integrability：若 \(s=2\)，

\[
\frac1n\operatorname{Tr}
\left[(J_n^{\rm loc})^{-2}
\mathbf1_{(0,\eta)}(J_n^{\rm loc})\right]
=n
\tag{6.4}
\]

沿适当 \(n\) 发散，但坏能量集中在一个可丢弃节点。

WLS 版本可取本地自测量因子

\[
G_n^{\rm loc}=(J_n^{\rm loc})^{1/2},
\tag{6.5}
\]

它是对角且满足 \(GG^T=J\)。相同的多数节点结论成立，只是坏节点的权重从 \(n^2\) 变为 \(n\)。

### 6.2 一个离域坏模态：每个节点都失败

令 \(C_n\) 是 \(n\)-cycle，\(L(C_n)\) 是组合 Laplacian，并取

\[
J_n^{\rm del}=n^{-1}I+L(C_n).
\tag{6.6}
\]

这是 one-hop SPD 矩阵，且 \(J_n^{\rm del}\preceq5I\)。常数向量
\(u=n^{-1/2}\mathbf1\) 是特征值 \(1/n\) 的单位特征向量，所以对**每个**节点

\[
\mu_{n,i}(\{1/n\})=|u(i)|^2=\frac1n.
\tag{6.7}
\]

于是，只要 \(1/n<\eta\)，

\[
\tau_{n,i}^{(2)}(\eta)\ge n
\tag{6.8}
\]

对 raw RHS 的每个节点成立；而对 whitened WLS，

\[
\tau_{n,i}^{(1)}(\eta)\ge1
\tag{6.9}
\]

也对每个节点成立。因此 raw RHS 的 (B) 失败；WLS 中取任意 \(a<1\)，(B) 同样失败。不是一个小比例的节点坏，而是同一个低模态以 \(1/n\) 权重铺到了全部节点。

WLS 可用完全局部的测量实现 (6.6)：令 \(B_n\) 是 cycle 的有向 incidence matrix，

\[
H_n=
\begin{bmatrix}
n^{-1/2}I\\ B_n
\end{bmatrix},
\qquad
G_n=H_n^T,
\tag{6.10}
\]

则 \(G_nG_n^T=H_n^TH_n=J_n^{\rm del}\)。第一组是节点自测量，第二组是边测量，holder 半径至多一。

### 6.3 结论

两个族都有相同量级的 \(\lambda_{\min}\)，也都有有限范围实现。差别只在小特征向量的空间分布。节点谱测度正是把“特征值多小”和“该模态在节点上有多大”合在一起的最小数学对象。

---

## 7. 一个 algorithm-independent 的通信轮数下界

定理 3.2 的定性等价不能说明 Richardson 是否轮数最优。事实上它通常不是。cycle 给出一个可精确计算的 separation。

### 7.1 cycle Green kernel

固定 \(0<\alpha\le1\)，令

\[
J_{n,\alpha}=\alpha I+L(C_n).
\tag{7.1}
\]

定义

\[
q_\alpha
:=\frac{\alpha+2-\sqrt{\alpha(\alpha+4)}}2
=e^{-\vartheta_\alpha},
\qquad
\vartheta_\alpha
=\operatorname{arcosh}(1+\alpha/2),
\tag{7.2}
\]

以及

\[
c_\alpha:=\frac1{\sqrt{\alpha(\alpha+4)}}.
\tag{7.3}
\]

若 \(d=0,1,\ldots,n-1\) 是有向 cyclic displacement，则直接解二阶差分方程可得

\[
(J_{n,\alpha}^{-1})_{0d}
=c_\alpha
\frac{q_\alpha^d+q_\alpha^{n-d}}
{1-q_\alpha^n}.
\tag{7.4}
\]

### 定理 7.1（任意 \(k\)-轮局部算法的下界）

在 raw RHS 模型 \(J_{n,\alpha}x=b,\ \|b\|_2\le1\) 中，假设每个节点初始只持有自己的 \(b_i\)，通信沿 cycle 边进行，消息大小和本地计算均不受限。若

\[
k+1<\frac n2,
\tag{7.5}
\]

则任意确定性 \(k\)-轮算法在每个节点的最坏绝对误差至少为

\[
\boxed{
c_\alpha q_\alpha^{k+1}.}
\tag{7.6}
\]

所以若误差不超过 \(\varepsilon<c_\alpha\)，必有

\[
k\ge
\frac{\log(c_\alpha/\varepsilon)}
{-\log q_\alpha}-1.
\tag{7.7}
\]

#### 证明

\(k\) 轮后，节点 0 只能看到 \(B_k(0)\) 内的输入。对固定已知模型，单位球上线性目标行的最优任意解码误差，至少是不可见坐标上的行范数；更直接地，取两个输入在球内相同、只在某个球外坐标上符号相反，即得同一结论。

由 (7.5)，距离 \(k+1\) 的坐标不可见。故最坏误差至少是单个逆矩阵元素

\[
|(J_{n,\alpha}^{-1})_{0,k+1}|.
\]

公式 (7.4) 给出

\[
|(J_{n,\alpha}^{-1})_{0,k+1}|
\ge c_\alpha q_\alpha^{k+1}.
\]

cycle 的平移对称性把结论推广到每个节点。证毕。

这个下界不是只针对 polynomial filter；它针对相同信息半径下的任意确定性本地解码器。

### 7.2 Chebyshev 上界与量级匹配

\(J_{n,\alpha}\) 的谱位于 \([\alpha,\alpha+4]\)。令

\[
\kappa_\alpha:=\frac{\alpha+4}{\alpha},
\qquad
\rho_\alpha:=
\frac{\sqrt{\kappa_\alpha}-1}
{\sqrt{\kappa_\alpha}+1}.
\tag{7.8}
\]

缩放的 Chebyshev residual 给出一个次数 \(k\) 的多项式
\(p_k^{\rm Ch}\)，满足

\[
\max_{\lambda\in[\alpha,\alpha+4]}
|1-\lambda p_k^{\rm Ch}(\lambda)|
\le2\rho_\alpha^{k+1}.
\tag{7.9}
\]

所以

\[
\left\|J_{n,\alpha}^{-1}
-p_k^{\rm Ch}(J_{n,\alpha})\right\|_2
\le\frac2\alpha\rho_\alpha^{k+1}.
\tag{7.10}
\]

次数 \(k\) 只需 \(k\) hop。由此，以下次数足够：

\[
k\ge
\frac{\log(2/(\alpha\varepsilon))}
{\log(1/\rho_\alpha)}-1.
\tag{7.11}
\]

当 \(\alpha\downarrow0\) 时，

\[
-\log q_\alpha\sim\sqrt\alpha,
\qquad
\log(1/\rho_\alpha)\sim\sqrt\alpha,
\qquad
c_\alpha\sim\frac1{2\sqrt\alpha}.
\tag{7.12}
\]

因此在 \(n\) 足够大使 (7.5) 不截断传播时，上下界共同给出

\[
k_*(\alpha,\varepsilon)
=\Theta\!\left(
\alpha^{-1/2}
\log\frac1{\varepsilon\alpha}
\right)
\tag{7.13}
\]

（对 \(0<\varepsilon\le1\)，上下对数只差至多常数因子）。这是 local information radius 的 matching order，而不是某个特定迭代法的分析。

### 7.3 Richardson 可以慢一个条件数平方根

Richardson 在常数特征向量上的贡献给出每个节点

\[
\mathcal E_i(p_t)^2
\ge
\frac1n\frac{(1-\alpha/M)^{2t}}{\alpha^2}.
\tag{7.14}
\]

故只要 \(1/(\varepsilon\alpha\sqrt n)>1\)，必要条件为

\[
t\ge
\frac{
\log(1/(\varepsilon\alpha\sqrt n))
}{-\log(1-\alpha/M)}.
\tag{7.15}
\]

取 \(\alpha=1/n\) 和 \(M=5\)，Richardson 需要

\[
t=\Omega\!\left(n\log\frac{\sqrt n}{\varepsilon}\right),
\tag{7.16}
\]

而 Chebyshev 只需

\[
k=O\!\left(\sqrt n\log\frac n\varepsilon\right).
\tag{7.17}
\]

所以 (C) 的意义是“Richardson 在存在性层面 universal”，绝不是“Richardson 在轮数层面最优”。一个可信的强主定理必须保留这一区分。

---

## 8. block nodes

### 8.1 正确对象是矩阵值谱测度

令节点 \(i\) 的状态块维数为 \(d_i\)，\(P_i:\mathbb R^N\to\mathbb R^{d_i}\) 是坐标投影。定义正半定矩阵值测度

\[
\Omega_{J,i}(B)
:=P_i\mathbf1_B(J)P_i^T
\in\mathbb R^{d_i\times d_i}.
\tag{8.1}
\]

其总质量为

\[
\Omega_{J,i}((0,M])=I_{d_i}.
\tag{8.2}
\]

定义 block 尾

\[
\mathcal T_{J,i}^{(s)}(\eta)
:=\left\|
\int_{(0,\eta)}\lambda^{-s}
\,d\Omega_{J,i}(\lambda)
\right\|_{\rm op}.
\tag{8.3}
\]

block operator error 为

\[
\mathcal E_{J,i}^{\rm blk}(p)
:=\|P_i(J^{-1}-p(J))G\|_{2\to2}.
\tag{8.4}
\]

### 定理 8.1（block 版）

把定理 3.2 中的 \(\tau_{J,i}^{(s)}\) 换成
\(\mathcal T_{J,i}^{(s)}\)，把 \(\mathcal E_{J,i}\) 换成
\(\mathcal E_{J,i}^{\rm blk}\)，三条件仍然等价，显式 \(t\) 界完全不变。

#### 证明

引理 3.1 的标量二次型升级为 Loewner 夹逼

\[
\begin{aligned}
c_-\int\lambda^{-s}|r_p(\lambda)|^2d\Omega_{J,i}(\lambda)
&\preceq
P_i f_p(J)GG^Tf_p(J)P_i^T\\
&\preceq
c_+\int\lambda^{-s}|r_p(\lambda)|^2d\Omega_{J,i}(\lambda).
\end{aligned}
\tag{8.5}
\]

中间矩阵的 operator norm 正是
\((\mathcal E_{J,i}^{\rm blk}(p))^2\)。低谱上
\(|r_p|^2\ge1/4\) 给出必要性；高谱上
\(\lambda^{-s}|r_{p_t}|^2\le\eta^{-s}(1-\eta/M)^{2t}\)，再用 (8.2) 给出充分性。其余多数节点量词不变。证毕。

### 8.2 trace 版不能无条件替代 operator 版

若关心 block Frobenius error，可使用标量 trace measure

\[
\mu_{J,i}^{\rm tr}(B):=\operatorname{Tr}\Omega_{J,i}(B).
\tag{8.6}
\]

但若关心任意 block 方向的最坏误差，(8.3) 才是准确条件。对 PSD
\(A\in\mathbb R^{d_i\times d_i}\)，

\[
\|A\|_{\rm op}\le\operatorname{Tr}A
\le d_i\|A\|_{\rm op}.
\tag{8.7}
\]

因此，当 \(d_i\le d_{\max}<\infty\) 一致有界时，trace 与 operator 尾在定性上等价，只差固定常数；若 block 维数无界，normalized trace 可以把一个极坏方向稀释掉，不能推出 operator guarantee。

### 推论 8.2（定理 3.5 的 block 版）

在 exact whitened WLS 条件 \(GG^T=J\) 下，令

\[
\mathcal E_{J,i}^{\rm blk}[K]
:=\|P_i(J^{-1}-K)G\|_{2\to2}.
\tag{8.8}
\]

则任意（可以非对称且不与 \(J\) 交换的）\(K\) 只要满足
\(\|K\|_{2\to2}\le B\)，便有

\[
\boxed{
\left\|
\int_{(0,\eta)}\lambda^{-1}
\,d\Omega_{J,i}(\lambda)
\right\|_{\rm op}^{1/2}
\le
\mathcal E_{J,i}^{\rm blk}[K]+B\sqrt\eta.}
\tag{8.9}
\]

证明与 (3.33)–(3.39) 相同，只需把 \(e_i^T\) 换成
\(P_i\)，并在三角不等式后取 operator norm。因此，若以
block operator error 计量，定理 3.5 的 class-uniform 多数节点充要条件也逐字成立；不要求 block 维数有界。维数一致有界只是在希望把 operator 条件替换成 trace 条件时才需要。

---

## 9. 从小谱质量上界得到显式轮数率

抽象尾条件 (B) 可以很弱，因而不自带统一 rate。若模型类还能给出小谱质量的幂律上界，就能直接算出 Richardson 的多数节点轮数。

### 定理 9.1（节点小谱质量 \(\Rightarrow\) Richardson rate）

固定 \(\beta>s\)、\(C<\infty\) 和 \(0<\delta<1\)。假设对每个
\(\theta\in\mathfrak F\)，至少 \(1-\delta\) 的节点同时满足

\[
F_{\theta,i}(x)
:=\mu_{\theta,i}((0,x])
\le Cx^\beta
\qquad(0<x\le M).
\tag{9.1}
\]

则对这些节点，

\[
\boxed{
\mathcal E_{\theta,i}(p_t)^2
\le
c_+ C\beta\Gamma(\beta-s)
\left(\frac{M}{2t}\right)^{\beta-s}.}
\tag{9.2}
\]

特别地，

\[
t\ge
\frac M2
\left(
\frac{c_+C\beta\Gamma(\beta-s)}
{\varepsilon^2}
\right)^{1/(\beta-s)}
\tag{9.3}
\]

足以使至少 \(1-\delta\) 的节点误差不超过 \(\varepsilon\)。

#### 证明

由

\[
(1-\lambda/M)^{2t}\le e^{-2t\lambda/M}
\tag{9.4}
\]

和引理 3.1，只需估计

\[
\int_0^M \lambda^{-s}e^{-a\lambda}\,dF_{\theta,i}(\lambda),
\qquad a:=2t/M.
\tag{9.5}
\]

把 \(F(x)\) 在 \(x>M\) 延拓为 \(1\)。因 (9.1) 在 \(M\) 处蕴含
\(1\le CM^\beta\)，故 \(F(x)\le Cx^\beta\) 对所有 \(x>0\) 仍成立。又因
\(\beta>s\)，\(x^{-s}F(x)\to0\) 当 \(x\downarrow0\)。Stieltjes 分部积分给出

\[
\begin{aligned}
\int_0^\infty \lambda^{-s}e^{-a\lambda}\,dF(\lambda)
&=\int_0^\infty F(\lambda)
\left(
s\lambda^{-s-1}
+a\lambda^{-s}
\right)e^{-a\lambda}\,d\lambda\\
&\le C\left[
s\Gamma(\beta-s)a^{-(\beta-s)}
+\Gamma(\beta-s+1)a^{-(\beta-s)}
\right]\\
&=C\beta\Gamma(\beta-s)a^{-(\beta-s)}.
\end{aligned}
\tag{9.6}
\]

代回 \(a=2t/M\) 即得。证毕。

若把常见记号“谱维数”写成 \(d_{\rm sp}=2\beta\)，这个充分条件要求

\[
d_{\rm sp}>2s.
\tag{9.7}
\]

于是 exact whitened WLS（\(s=1\)）的阈值是 \(d_{\rm sp}>2\)，raw RHS
（\(s=2\)）的阈值是 \(d_{\rm sp}>4\)。在欧氏格点/平移不变 Laplacian
上，这分别对应 Green function 与 squared Green function 的经典临界维数。该解释很有用，但临界数值本身不是新的。

block 版本也成立：只需把 (9.1) 换成 Loewner 上界

\[
\Omega_{J,i}((0,x])\preceq Cx^\beta I_{d_i},
\tag{9.8}
\]

同一个 Stieltjes 计算给出矩阵 Loewner 上界。

---

## 10. 更强的 sharp theorem：格点上的任意局部算法与最佳多项式同阶

第 3 节的 iff 只回答“是否存在有限次多项式”，而且证明很短。在规则格点上，可以把它升级为一个有匹配上下界的通信轮数定理。这一节的结论不再只针对 Richardson，也不再只针对 polynomial filter；下界允许算法在已见数据上使用任意非线性解码。

### 10.1 模型和两个输入几何

在 \(\ell_2(\mathbb Z^d)\) 上定义组合 Laplacian

\[
(Lu)(x)=\sum_{j=1}^d
\big(2u(x)-u(x+e_j)-u(x-e_j)\big).
\tag{10.1}
\]

令前向差分

\[
(Du)(x,j)=u(x+e_j)-u(x),
\qquad D^*D=L.
\tag{10.2}
\]

下面用 \(q\in\{1,2\}\) 表示低谱逆矩指数：

\[
\begin{array}{c|c|c|c}
q&\text{输入空间 }\mathcal H_q&G_q&根节点目标\ T_q\xi\\ \hline
1&\ell_2(\mathbb Z^d\times[d])&D^*&
\big(L^{-1}D^*\xi\big)(0)\\
2&\ell_2(\mathbb Z^d)&I&
\big(L^{-1}\xi\big)(0).
\end{array}
\tag{10.3}
\]

\(q=1\) 就是 whitened WLS/边噪声输入；因为
\(G_1G_1^*=L\)，目标协方差是 \(L^{-1}\)，也就是 massless
GFF。\(q=2\) 是 raw right-hand side；若输入是 i.i.d. Gaussian，目标协方差是 \(L^{-2}\)，对应 membrane/bilaplacian field。

对 \(q=2\)，数据坐标 \(x\) 由节点 \(x\) 持有。对 \(q=1\)，定向边数据
\((x,j)\) 指边 \(x\to x+e_j\)，约定由尾节点 \(x\) 持有。记这些 holder 中距根至多 \(r\) 的集合为 \(S_{q,r}\)。允许一个 \(r\)-local 算法为任意确定性函数

\[
\Phi_r:\mathbb R^{S_{q,r}}\longrightarrow\mathbb R.
\tag{10.4}
\]

它预先知道整个格点模型和所有系数，本地计算与消息长度不受限。因此下面的下界是一个 oracle 下界，不是由某个具体迭代法造成的。定义

\[
\mathsf E_q^{\rm loc}(r)
:=
\inf_{\Phi_r}
\sup_{\substack{\xi\in\mathcal H_q\\
\|\xi\|_2\le1}}
\left|T_q\xi-
\Phi_r(\xi|_{S_{q,r}})\right|.
\tag{10.5}
\]

另外定义 degree-限制的共同多项式误差

\[
\mathsf E_q^{\rm poly}(n)
:=
\inf_{\deg p\le n-1}
\left\|
e_0^*\big(L^{-1}-p(L)\big)G_q
\right\|_2,
\qquad n\ge1.
\tag{10.6}
\]

这里的 \(L^{-1}\) 先在 Fourier 积分有定义的自然稠密 core 上用 multiplier
\(1/\ell(\theta)\) 解释，也可等价地从 \((L+mI)^{-1}\) 的
\(m\downarrow0\) 极限出发。定理会精确说明何时它的根节点行能唯一延拓为有界 \(\ell_2\) 线性泛函；在低维 recurrent 情形，不能预先假装该乘子对每个有限支撑输入都有有限值。

### 定理 10.1（\(\mathbb Z^d\) 上的临界维数和匹配轮数率）

对 \(q\in\{1,2\}\)，有以下结论。

1. \(T_q\) 能延拓为 \(\mathcal H_q\) 上的有界根节点泛函，当且仅当
   
   \[
   d>2q.
   \tag{10.7}
   \]
   
   因此 WLS 的临界维数是 \(2\)，raw RHS 的临界维数是 \(4\)。当 \(d\le2q\) 时，该根节点泛函没有有界 \(\ell_2\) 延拓，所以不存在有限最坏误差的局部方法。在 Green 行逐点有定义的维数，可直接用其有限空间截断得到发散的单位输入；更低的 recurrent 维数则连某些自然输入的 massless 极限都发散。

2. 若 \(d>2q\)，令
   
   \[
   \gamma_{d,q}:=\frac{d-2q}{2}>0.
   \tag{10.8}
   \]
   
   则存在只依赖 \(d,q\) 的常数
   \(0<c_{d,q}\le C_{d,q}<\infty\)，使对所有 \(r\ge1,n\ge1\)，
   
   \[
   c_{d,q}(r+1)^{-\gamma_{d,q}}
   \le \mathsf E_q^{\rm loc}(r)
   \le C_{d,q}(r+1)^{-\gamma_{d,q}},
   \tag{10.9}
   \]
   
   且
   
   \[
   c_{d,q}(n+1)^{-\gamma_{d,q}}
   \le \mathsf E_q^{\rm poly}(n)
   \le C_{d,q}(n+1)^{-\gamma_{d,q}}.
   \tag{10.10}
   \]
   
   因为 (10.6) 中的滤波在 raw RHS 中至多传播
   \(n-1\) hop，在边数据 WLS 中至多传播 \(n\) hop，还有更强的同阶夹逼
   
   \[
   c_{d,q}(n+1)^{-\gamma_{d,q}}
   \le \mathsf E_q^{\rm loc}(n)
   \le \mathsf E_q^{\rm poly}(n)
   \le C_{d,q}(n+1)^{-\gamma_{d,q}}.
   \tag{10.11}
   \]
   
   换句话说，最佳多项式在轮数指数上已经追平了知道全部模型且可用任意非线性本地计算的 oracle。

   由于格点平移对称，同一个下界在每个节点成立。因此在本文采用的
   “逐节点先取单位输入 supremum、再计算坏节点比例”语义下，允许任意固定
   \(\delta<1\) 的失败节点并不能改善 (10.9) 的指数；这不是只在一个特殊根上的偶然下界。

3. 多项式误差有精确的 Christoffel 表示。若 \(\mu_d\) 是 \(L\) 在根的谱测度，并定义
   
   \[
   d\omega_{d,q}(\lambda):=\lambda^{-q}d\mu_d(\lambda),
   \tag{10.12}
   \]
   
   则
   
   \[
   \boxed{
   \big(\mathsf E_q^{\rm poly}(n)\big)^2
   =
   \inf_{\substack{R\in\mathcal P_n\\R(0)=1}}
   \int R(\lambda)^2d\omega_{d,q}(\lambda)
   =\Lambda_n(\omega_{d,q},0),}
   \tag{10.13}
   \]
   
   其中 \(\Lambda_n\) 是端点 Christoffel function。最佳残差为
   
   \[
   R_n^*(\lambda)
   =\frac{K_n(\lambda,0)}{K_n(0,0)},
   \qquad
   p_n^*(\lambda)=\frac{1-R_n^*(\lambda)}{\lambda},
   \tag{10.14}
   \]
   
   \(K_n\) 是 \(L_2(\omega_{d,q})\) 中 degree \(\le n\) 的 reproducing kernel。

4. Richardson 在定性上能收敛，却严格慢一个平方级的轮数指数。对任意固定 \(M\ge4d\)，
   
   \[
   \left\|
   e_0^*\big(L^{-1}-p_t(L)\big)G_q
   \right\|_2
   \asymp
   t^{-(d-2q)/4},
   \tag{10.15}
   \]
   
   而 (10.10) 的最佳 degree-\(t\) 率是
   \(t^{-(d-2q)/2}\)。因此 Richardson 可以见证“存在有限轮数”，但不能作为这个无谱隙问题的最优通信算法。

#### 证明

**第一步：Fourier 谱测度和临界维数。** Fourier 变换把 \(L\) 变成乘子

\[
\ell(\theta)=2d-2\sum_{j=1}^d\cos\theta_j,
\qquad \theta\in[-\pi,\pi]^d.
\tag{10.16}
\]

因此

\[
\int f(\lambda)d\mu_d(\lambda)
=\frac1{(2\pi)^d}
\int_{[-\pi,\pi]^d}f(\ell(\theta))d\theta.
\tag{10.17}
\]

在 \(\theta=0\) 附近，
\(\ell(\theta)=|\theta|^2+O(|\theta|^4)\)，所以

\[
\mu_d([0,\lambda])
\sim \frac{v_d}{(2\pi)^d}\lambda^{d/2},
\qquad
\rho_d(\lambda):=\frac{d\mu_d}{d\lambda}
\sim
\frac{\lambda^{d/2-1}}
{2^d\pi^{d/2}\Gamma(d/2)}.
\tag{10.18}
\]

根行平方范数为

\[
\|T_q\|^2
=\int_0^{4d}\lambda^{-q}d\mu_d(\lambda).
\tag{10.19}
\]

用 (10.18) 检查零点可积性，(10.19) 有限当且仅当
\(d/2-q>0\)，即 (10.7)。如果该积分发散，谱截断立即说明该 multiplier 泛函没有有界 \(\ell_2\) 延拓；在相应 Green/gradient-Green 系数逐点有定义时，再对系数行作越来越大的有限空间截断，就得到一列有限支撑的单位输入，其根目标值趋于无穷。对剩余 recurrent 情形，\((L+mI)^{-1}\) 的根响应随 \(m\downarrow0\) 已发散，结论只会更强。

**第二步：任意局部解码的精确 oracle 误差。** 设 \(a_q\) 是泛函 \(T_q\) 的 \(\ell_2\) 系数行，分解为已见坐标与不可见坐标

\[
a_q=a_{q,S_{q,r}}+a_{q,S_{q,r}^c}.
\tag{10.20}
\]

截断线性解码
\(\Phi_r(\xi)=\langle a_{q,S_{q,r}},\xi\rangle\)
的最坏误差是
\(\|a_{q,S_{q,r}^c}\|_2\)。反过来，任意甚至非线性的 \(\Phi_r\) 在两个输入

\[
\xi_\pm
=\pm\frac{a_{q,S_{q,r}^c}}
{\|a_{q,S_{q,r}^c}\|_2}
\tag{10.21}
\]

上看到的本地数据完全一样（都为零），但两个目标值相反。因此它在二者中至少一个上的误差不小于该尾范数，从而

\[
\boxed{
\mathsf E_q^{\rm loc}(r)
=\|a_{q,S_{q,r}^c}\|_2.}
\tag{10.22}
\]

这是一个精确的 minimax 等式，不仅是下界。

允许私有随机性也不能降低这个量级：令
\(a:=\|a_{q,S_{q,r}^c}\|_2\)。在 \(\xi_+\) 与
\(\xi_-\) 上，算法看到同一局部输入，所以输出具有同一分布 \(Y\)；而
\(\max\{\mathbb E|a-Y|,\mathbb E|-a-Y|\}\ge a\)，平方损失下也有
\(\max\{\mathbb E(a-Y)^2,\mathbb E(-a-Y)^2\}\ge a^2\)。因此 (10.22)
同时是 worst-case expected absolute error 和 root-MSE 的 randomized lower bound；线性截断上界仍为确定性的。

令

\[
g(x):=(L^{-1})_{0x}.
\tag{10.23}
\]

格点 Green function 的经典渐近式是

\[
g(x)=a_d|x|^{2-d}+O(|x|^{-d}),
\qquad a_d>0.
\tag{10.24}
\]

故 raw RHS 的精确尾为

\[
\big(\mathsf E_2^{\rm loc}(r)\big)^2
=\sum_{|x|_1>r}|g(x)|^2
\asymp r^{4-d}.
\tag{10.25}
\]

对 WLS，数据 \((x,j)\) 的行系数是
\(g(x+e_j)-g(x)\)（正负号取决于 incidence 方向，不影响范数）。对 (10.24) 作一次差分得

\[
\sum_{j=1}^d|g(x+e_j)-g(x)|^2
\asymp |x|^{2-2d},
\tag{10.26}
\]

因而

\[
\big(\mathsf E_1^{\rm loc}(r)\big)^2
=\sum_{|x|_1>r}\sum_{j=1}^d
|g(x+e_j)-g(x)|^2
\asymp r^{2-d}.
\tag{10.27}
\]

式 (10.25)–(10.27) 证明 (10.9)。注意这里不是用谱维数猜指数；下界真正使用了 Green 行在空间中的尾分布。

**第三步：多项式误差就是 Christoffel function。** 令

\[
R(\lambda)=1-\lambda p(\lambda).
\tag{10.28}
\]

由谱定理和 \(G_qG_q^*=L^{2-q}\)（分别是
\(L\) 与 \(I\)），

\[
\left\|e_0^*(L^{-1}-p(L))G_q\right\|_2^2
=\int_0^{4d}\lambda^{-q}R(\lambda)^2d\mu_d(\lambda).
\tag{10.29}
\]

如果 \(\deg p\le n-1\)，则 \(R\in\mathcal P_n,R(0)=1\)。反之，任意满足这两个条件的 \(R\) 都使
\((1-R(\lambda))/\lambda\) 成为次数至多 \(n-1\) 的多项式。这证明 (10.13)。Christoffel 的再生核变分原理给出 (10.14)。

现在给出一个不需要猜测最佳多项式的自包含阶数计算。对 \(\alpha>-1\)，记

\[
\Lambda_n^{(\alpha,M)}
:=
\inf_{\substack{R\in\mathcal P_n\\R(0)=1}}
\int_0^M R(\lambda)^2\lambda^\alpha d\lambda.
\tag{10.30}
\]

用 shifted Jacobi 多项式
\(P_k^{(0,\alpha)}(2\lambda/M-1)\) 并把端点值与平方范数代入 Christoffel–Darboux kernel，可得到精确式

\[
\left(\Lambda_n^{(\alpha,M)}\right)^{-1}
=\frac1{M^{\alpha+1}}
\sum_{k=0}^n
(2k+\alpha+1)
\left[
\frac{\Gamma(k+\alpha+1)}
{\Gamma(\alpha+1)\Gamma(k+1)}
\right]^2.
\tag{10.31}
\]

Gamma 比值渐近式随即给出

\[
\Lambda_n^{(\alpha,M)}
\sim
(\alpha+1)M^{\alpha+1}
\Gamma(\alpha+1)^2
n^{-2\alpha-2}.
\tag{10.32}
\]

对当前问题取

\[
\alpha=\frac d2-q-1>-1.
\tag{10.33}
\]

由 (10.18)，存在 \(\eta,c,C>0\)，使

\[
c\lambda^\alpha
\le \lambda^{-q}\rho_d(\lambda)
\quad(0<\lambda<\eta),
\qquad
\lambda^{-q}\rho_d(\lambda)
\le C\lambda^\alpha
\quad(0<\lambda<4d).
\tag{10.34}
\]

第二个全区间上界可以如下核验：每个一维色散
\(2-2\cos\Theta\) 的特征函数是一个 Bessel 函数乘相位，大频率下为
\(O(|t|^{-1/2})\)；\(d\ge3\) 个独立和的特征函数可积，因而 \(\rho_d\) 有界连续。再把零点附近用 (10.18)、远离零点处用有界性即得 (10.34)。

把 (10.34) 代入变分定义，得

\[
c\Lambda_n^{(\alpha,\eta)}
\le \Lambda_n(\omega_{d,q},0)
\le C\Lambda_n^{(\alpha,4d)}.
\tag{10.35}
\]

由 (10.32)–(10.33)，两端都是
\(\Theta(n^{-(d-2q)})\)。开平方得 (10.10)。有限传播性结合 (10.9) 得 (10.11)。

**第四步：Richardson 的严格慢率。** Richardson 残差为
\((1-\lambda/M)^t\)，所以由 (10.18) 和端点 Laplace 渐近，若
\(\beta=d/2-q>0\)，则

\[
\begin{aligned}
&\left\|e_0^*(L^{-1}-p_t(L))G_q\right\|_2^2\\
&\qquad=
\int_0^{4d}\lambda^{-q}(1-\lambda/M)^{2t}d\mu_d(\lambda)\\
&\qquad\sim
\frac{\Gamma(\beta)}
{2^d\pi^{d/2}\Gamma(d/2)}
\left(\frac{M}{2t}\right)^\beta.
\end{aligned}
\tag{10.36}
\]

开平方即得 (10.15)。定理证毕。

### 推论 10.2（轮数、延迟与一个标准通信能耗模型）

设目标是使每个节点的单位输入最坏误差不超过
\(\varepsilon\)。在定理 10.1 的格点模型上，任意局部算法的最优通信半径/轮数为

\[
r_q^*(\varepsilon)
=\Theta\!\left(
\varepsilon^{-2/(d-2q)}
\right),
\qquad d>2q.
\tag{10.37}
\]

最佳 Christoffel 多项式达到这个界。相比之下，Richardson 需要

\[
r_q^{\rm Rich}(\varepsilon)
=\Theta\!\left(
\varepsilon^{-4/(d-2q)}
\right).
\tag{10.38}
\]

例如：

| 模型 | 维数 | 最优轮数 | Richardson 轮数 |
|---|---:|---:|---:|
| whitened WLS | \(d=3\) | \(\Theta(\varepsilon^{-2})\) | \(\Theta(\varepsilon^{-4})\) |
| whitened WLS | \(d=4\) | \(\Theta(\varepsilon^{-1})\) | \(\Theta(\varepsilon^{-2})\) |
| raw RHS | \(d=5\) | \(\Theta(\varepsilon^{-2})\) | \(\Theta(\varepsilon^{-4})\) |

若使用最简单的全边同步 graph-filter 实现：每条有向边每轮传送一个固定精度标量花费 \(E_{\rm msg}\)，并且全网 \(N\) 个节点同时估计，则 \(d\)-维环面每轮有 \(\Theta(dN)\) 次传输。最佳 polynomial 因而达到

\[
E_{\rm comm}^{\rm poly}(\varepsilon)
=\Theta\!\left(
dN E_{\rm msg}
\varepsilon^{-2/(d-2q)}
\right),
\tag{10.39}
\]

而同一种实现下 Richardson 的通信能耗指数加倍。同样的幂律是同步网络的 **minimax latency** 幂律。最佳 polynomial 可用三项递推或 Horner 形式在每节点
\(O(r_q^*)\) 标量算术和 \(O(1)\) 个工作状态向量内完成；正交多项式递推系数可在格点模型上离线预计算。

式 (10.39) 是一个可实现的能耗上界及与 Richardson 的同模型比较，**不是**任意消息传递协议上的 minimax 能耗下界：轮数下界本身不能推出每条边每轮都必须发消息。若要把星号“最优”也放到能耗上，必须固定消息位数/聚合规则并另证 information-flow 或 cut lower bound。它也不等于“局部通信在所有无线物理层上必然比集中式更省电”。若考虑启动开销、路径损耗、量化位数、冲突与睡眠调度，需要再加一层通信模型才能宣称实际 Joule 节省。

### 10.3 有限 SPD 网络并不逃离这个下界

定理 10.1 用无限格点上的 massless \(L\) 来写，是为了不让边界和有限尺寸参数掩盖主要指数。它可以严格嵌回本文的有限 SPD 模型。

令 \(L_N\) 是 \(d\)-维离散环面 \(\mathbb T_N^d\) 上的 Laplacian，并取

\[
J_{N,m}=L_N+mI\succ0.
\tag{10.40}
\]

对 raw RHS 取 \(G_{2,N,m}=I\)；对 WLS 取完全局部的因子

\[
G_{1,N,m}=\begin{bmatrix}D_N^*&\sqrt m I\end{bmatrix},
\qquad
G_{1,N,m}G_{1,N,m}^*=J_{N,m}.
\tag{10.41}
\]

对每个固定的局部半径 \(r\) 和多项式次数 \(n\)，先令
\(N\to\infty\)，再令 \(m\downarrow0\)，相应根行、局部 oracle 误差与多项式误差分别收敛到定理 10.1 中的量。理由是：固定
\(m>0\) 时 Fourier Riemann sums 收敛；再用
\(d>2q\) 时的 \(\lambda^{-q}\) 可积性作 dominated convergence。对 WLS 中新增的自测量行，

\[
m\int(\lambda+m)^{-2}d\mu_d(\lambda)\longrightarrow0,
\tag{10.42}
\]

因为
\(m/(\lambda+m)^2\le1/(4\lambda)\)。

因此可用对角化选择 \(N_r\uparrow\infty,m_r\downarrow0\)，使一列真正的有限、one-hop、有环 SPD 网络在半径 \(r\) 上保留 (10.11) 的上下界常数。这说明临界维数与局部轮数下界不是无限维奇性的产物，也完全不需要无环假设。

### 10.4 能否只用“谱维数”把定理推到一般图？

不能无条件地这样推。如果根谱密度只满足

\[
d\mu_o(\lambda)\asymp
\lambda^{d_{\rm sp}/2-1}d\lambda
\quad(\lambda\downarrow0),
\tag{10.43}
\]

那么 Christoffel 部分确实预测最佳 degree-\(n\) 多项式误差

\[
\mathsf E_q^{\rm poly}(n)
\asymp n^{-(d_{\rm sp}-2q)/2},
\qquad d_{\rm sp}>2q,
\tag{10.44}
\]

但谱测度没有说明 Green 行的能量在图距离上如何分布。任意 \(r\)-local 算法的最佳误差由 (10.22) 的**空间尾**决定，不由一个标量谱维数单独决定。

在异常扩散图上，这个区别不只是技术细节。若体积维数为
\(d_f\)、walk dimension 为 \(d_w\)，则常见
\(d_{\rm sp}=2d_f/d_w\)。raw Green 核的空间幂律会同时依赖
\(d_f,d_w\)，而不是只依赖 \(d_{\rm sp}\)。因此，要把 (10.11) 推广到一大类图，至少还需要双边 heat-kernel/Green annulus 估计；WLS 还需要梯度 Green 估计。“只假设 spectral dimension 就得到所有 local algorithms 的 matching lower bound”是不正确的。

### 10.5 这个强定理的创新性风险

定理 10.1 比定理 3.2 强得多：它给出临界维数、精确的 oracle 变分式、最佳轮数指数，并证明一个可分布实现的多项式类追平了任意局部解码。但它的每个单独部件都与经典结果紧密相邻：

- Green 核 \(|x|^{2-d}\) 及其差分渐近是经典 random-walk potential theory；
- residual polynomial 与谱测度上的正交多项式是 Krylov/CG 的经典结构；
- Christoffel 端点幂律是成熟的近似论；
- GFF 的 \(d>2\) 和 membrane field 的 \(d>4\) 临界维数早已知道。

截至本次检索，没找到一篇 primary source 把这些组件合成
“WLS/raw 输入几何 \(+\) Christoffel 最佳分布式滤波 \(+\) 任意 LOCAL 解码的 Green-tail minimax 下界”并得到 (10.11)。这使它成为一个**有希望的主定理骨架**，但仍不能在未做更全面查重前宣称原创。如果论文只做 \(\mathbb Z^d\) 标量情形，reviewer 仍可能认为这是三个经典事实的优雅组合。要明显降低这个风险，建议最终论文至少再完成一项：

1. 推到具有统一 ellipticity 的 periodic/random geometric 网络，并保留多数节点失败量词；
2. 在有体积增长、heat-kernel 与梯度估计的非齐次图上给出可验证的一般 matching theorem；
3. 设计不需要预知全局谱密度的自适应三项递推，并证明它在节点失败概率 \(\delta\) 下达到同一 oracle rate。

---

## 11. 与概率、图极限和算子代数的准确关系

### 11.1 它不是通常意义下的平均 uniform integrability

对每个模型先均匀随机选节点 \(I\)，再条件于 \(I=i\) 从
\(\mu_{J,i}\) 抽取 \(\Lambda\)。则

\[
\tau_{J,i}^{(s)}(\eta)
=\mathbb E\!\left[
\Lambda^{-s}\mathbf1_{\{\Lambda<\eta\}}
\mid I=i
\right].
\tag{11.1}
\]

(B) 说的是这些**条件尾期望**在随机节点上依概率、且对模型族一致地趋于零：

\[
\sup_{J\in\mathfrak F}
\mathbb P_I\left(
\mathbb E[\Lambda^{-s}\mathbf1_{\{\Lambda<\eta\}}\mid I]>a
\right)\longrightarrow0.
\tag{11.2}
\]

通常的 averaged uniform integrability 要求更强的

\[
\sup_{J\in\mathfrak F}
\frac1{|V_J|}\operatorname{Tr}
\left[J^{-s}\mathbf1_{(0,\eta)}(J)\right]
\longrightarrow0.
\tag{11.3}
\]

(11.3) 由 Markov 不等式推出 (B)，反向不成立，(6.1)–(6.4) 就是反例。因此若论文使用“uniform integrability”一词，应明确是 (11.2) 还是 (11.3)，不要把二者混写。

### 11.2 Benjamini–Schramm/local weak convergence 只管 bounded tests

有界度标记图的 local weak convergence 能自然产生随机根处的谱测度，并控制多项式、resolvent（离开实谱）或其他有界连续 test function。函数
\(\lambda^{-s}\) 在零点无界；仅有弱谱收敛不能传递逆矩。尾条件 (B) 正是额外需要的控制之一。

但 (B) 本身不要求图列收敛，也不由 Benjamini–Schramm convergence 自动推出。若想把它变成拓扑定理，需要另外证明：某个可由有限邻域检测的几何/锚定条件当且仅当保证 (B)。这一步目前尚未完成。

### 11.3 与 measure topology 的相似和差别

有限 von Neumann 代数中，affiliated unbounded operators 与 convergence in measure 为“没有统一谱隙时逼近逆”提供了成熟语言。共同多项式逼近 \(J^{-1}\) 的思想在该语言下并不突兀。

不过，本文误差是

\[
i\mapsto\|e_i^T(J^{-1}-p(J))G\|_2,
\]

即先取 row \(L_2\) 能量，再对**坐标对角代数**做节点分位数；标准 noncommutative measure topology 通常按算子的谱投影/广义奇异值控制。两者有关但不相同。没有证明二者等价时，不应声称定理 3.2 只是某个标准 operator-algebra theorem 的逐字实例；也不应因为术语不同就反向宣称原创。

### 11.4 与 inverse-closed / uniform Roe 结果的边界

另一个容易混淆的已知方向是“有限传播算子的逆是否仍可由有限传播算子逼近”。空间衰减 Banach algebra 的 inverse-closedness（例如 Motee–Sun）以及 Property A 空间上的 quasi-local/uniform Roe 结果（例如 Špakula–Zhang）都给出很强的 **operator-norm** 局部化结论。若 \(J\) 作为有界算子还有统一正谱隙，那么 \(J^{-1}\) 本身有界，这些理论与 Demko 型结论处在同一个经典制度内。

本文刻意处理的是 \(\inf\lambda_{\min}(J_\theta)=0\)，此时极限逆可能是 unbounded affiliated operator，因而通常不是 uniform Roe algebra 中的元素；同时我们只要求丢弃 \(\delta\) 比例节点后的 row error，而不是全局 operator norm。Roe 文献中的 ghost projection 和 measured asymptotic expander 现象说明“逐元素很小/多数位置很小”与 operator-norm 可逼近确实可能分离，但 Li–Špakula–Zhang 的 rigidity 定理并不直接等价于 (3.6)。因此这些工作是重要的边界和警告，不是定理 3.2 的现成证明。

---

## 12. 截至 2026 年的 primary-source 查重

### 12.1 已明确撞上的组成部分

| 本文成分 | 已有主来源 | 判断 |
|---|---|---|
| 节点/根谱测度，把谱权重写成 \(|u_k(i)|^2\) | Godsil–Mohar (1988)；Bordenave–Lelarge (2010) | 经典谱图论对象 |
| 用多项式 \(p(J)\) 得到有限传播/分布式图滤波 | Hammond–Vandergheynst–Gribonval (2011)；Shuman et al. (2018) | 经典 graph-filter machinery |
| 用多项式逼近矩阵逆、由谱区间给指数衰减 | Demko–Moss–Smith (1984)；Benzi–Razouk (2007) | 经典近似论/矩阵函数结论 |
| 分布式 inverse graph filtering | Jiang–Tay (2019)；Emirov et al. (2022) | 直接相邻的算法文献 |
| Richardson/Landweber 残差 \(r_t(\lambda)=(1-\omega\lambda)^t\) 的谱积分 | Landweber/迭代正则化与 source-condition 文献；Andreev (2015) 给出 Landweber rate characterizations | 核心谱演算经典 |
| uniformly stable、非交换的 \(K\) 所满足的 WLS 低谱投影不等式 (3.29) | Hohage–Weidling (2017) 等 converse/source-condition 文献是最近的抽象邻居 | 未找到逐节点、\(KG\) sufficient-statistic 形式的同一定理；但证明只有一次极分解式等距与投影三角不等式，应按 lemma-level 贡献定位 |
| 在顶点可传图上，\(L_2\) 可积谱 multiplier 由有限半径多项式逼近 | Backhausz–Virág (2017), Theorem 2 与 Lemma 20 | 与固定齐次图上的定性谱尾条件直接碰撞 |
| 瞬态顶点可传图上 GFF 是 linear factor of i.i.d. | Backhausz–Virág (2017), Proposition 26 | WLS \(q=1\) 的随机场/单根定性版本实质已覆盖 |
| 端点 Christoffel function 在 power-type weight 下的 sharp 幂律 | Danka–Totik (2018)；Danka (2017) | 定理 10.1 的最佳 polynomial rate 所需近似论是成熟的 |
| 以谱测度加权、在 \(R(0)=1\) 下选择最小二次误差 residual polynomial | Daniel (1967) 的 Hilbert-space CG 极值表述；Knockaert (1987) 的谱正交多项式解释 | (10.13) 的 Krylov/正交多项式骨架是经典的；把特定格点权重写成 endpoint Christoffel function 不是单独的原创点 |
| \(\mathbb Z^d\) Green function 的 Newtonian-potential 渐近 | Uchiyama (1998) 及经典 random-walk potential theory | 定理 10.1 的空间尾指数有成熟基础 |
| bilaplacian/membrane field 的临界维数 \(4\) | Cipriani (2013) 及 membrane-model 文献 | raw RHS \(q=2\) 的临界数值不新 |
| local weak limit 推出谱测度收敛，但零附近无界函数需要额外可积性 | Bordenave–Lelarge (2010)；Abért–Thom–Virág 的 pointwise spectral-measure 工作 | 图极限中的标准障碍 |
| 无统一谱隙时把逆视为可测/affiliated operator | Nelson (1974) 及 Murray–von Neumann algebra 文献 | 成熟抽象框架 |
| 有谱隙/有界逆的空间衰减 algebra 与 finite-propagation closure | Motee–Sun (2017)；Špakula–Zhang (2020) | operator-norm 制度已有强理论；不覆盖 unbounded inverse 的多数节点 row topology |
| ghost projections、measured asymptotic expanders 与 Roe algebra | Li–Špakula–Zhang (2023) | 说明“局部系数小”与全局 norm locality 可分离；其 rigidity 问题不是本文 iff |
| 2025–2026 的 inverse polynomial / entrywise inversion | Embree et al. (2025)；Ghadiri–Nguyen–Yang (SODA 2026) | 说明该方向仍活跃，但误差模型不同 |

### 12.2 Backhausz–Virág (2017) 的逐式碰撞审计

这篇文章不能只在 related work 中轻描淡写。它对定性版的碰撞很直接。以下命题号按该文的公开 v2/期刊排版核对：

1. **Theorem 2.** 在无限顶点可传图上，一个有限测度是 linear factor of i.i.d. 的谱测度，当且仅当它对图的根谱测度 \(\nu\) 绝对连续。
2. **Lemma 20(a).** 它把闭空间
   \(\overline{\{p(A)\delta_o:p\text{ polynomial}\}}\subset\ell_2(V)\)
   与 \(L_2(\nu)\) 等距同构；关键等式正是
   
   \[
   \|p(A)\delta_o\|_2^2
   =\int |p(x)|^2d\nu(x).
   \tag{12.1}
   \]
   
   证明明确使用了多项式在 \(L_2(\nu)\) 中稠密。
3. **Lemma 20(b,c).** 若 \(h\in L_2(\nu)\)，则有一个 spherical linear factor，其谱测度密度是 \(h^2\)；用多项式
   \(p_n\to h\) 的 \(L_2(\nu)\) 收敛就是有限半径 linear block factors 的均方逼近。
4. **Definition 25 与 Proposition 26.** 对 \(d_0\)-regular 图，GFF 协方差的 adjacency 谱密度是
   
   \[
   \sum_{k\ge0}(x/d_0)^k
   =\frac{d_0}{d_0-x}.
   \tag{12.2}
   \]
   
   瞬态性保证 (12.2) 对 \(\nu\) 可积，Proposition 26 因此得出 GFF 是 linear factor of i.i.d.

把

\[
J=I-A/d_0
\tag{12.3}
\]

代入，(12.2) 就是 \(J^{-1}\) 的谱密度。其平方根 multiplier
\(h(\lambda)=\lambda^{-1/2}\) 属于 \(L_2(\mu_o)\) 当且仅当

\[
\int\lambda^{-1}d\mu_o(\lambda)<\infty.
\tag{12.4}
\]

所以，若把 WLS 看成“从 i.i.d. Gaussian 产生一个协方差为
\(J^{-1}\) 的随机场，并在一个固定顶点可传图的根上作均方逼近”，那么 \(q=1\) 的定性 iff 实质已经被覆盖。同理，在 Lemma 20 中取
\(h(\lambda)=\lambda^{-1}\)，也直接得到 \(q=2\) 的固定图定性条件
\(\int\lambda^{-2}d\mu_o<\infty\)。

本文仍有四个**真实但需克制表述**的差别：

- Backhausz–Virág 处理一个固定无限顶点可传图；定理 3.2 处理有限模型族的共同多项式与失败节点分位数。
- 他们的 GFF 结论只要生成相同 Gaussian law；我们的 WLS 误差固定了原始 measurement holder 与输入算子 \(D^*\)。两者协方差误差一致，但不是逐数据坐标的同一算子表示。
- 他们证明存在某个 polynomial approximating sequence，没有证明本文的固定 Richardson 序列在 class-uniform 多数节点意义下 universal。
- 他们没有给出定理 10.1 的 Christoffel 最佳轮数、任意局部解码下界或 Richardson 的严格慢率。

前三个差别更像 uniform/finite-network/WLS packaging，未必足以支撑强原创性；第四个才是更值得追求的数学提升。

### 12.3 没有检索到的精确组合

截至 2026-09-13 的定向检索中，没有找到一篇 primary source 同时陈述：

1. 一族有限 SPD finite-range matrices，无统一谱下界；
2. 一个对整个族共同的 polynomial filter；
3. 允许每个网络有至多 \(\delta\) 比例失败节点；
4. 用节点条件逆矩尾 (3.6) 给出 iff；
5. 同时区分 raw RHS 的 \(s=2\) 与 whitened WLS 的 \(s=1\)；
6. 证明固定 Richardson 家族在**存在性**意义下 universal。

这只能支持“未找到相同表述”，不能支持“已证明原创”。原因是定理 3.2 的证明只有三个标准动作：谱定理、低/高谱切分、以及多项式在零附近有界。即使完整句子没人写过，reviewer 仍很可能把它视为 dominated/monotone-convergence 风格的直接推论。

### 12.4 对主定理强度的实话判断

**单独作为论文唯一主定理：不够强。**

优点是 scope 大、量词干净、纠正了 WLS 指数，并明确允许少数节点失败；它很适合作为全文的“总判据”。但原创数学的证明难度仍偏低，而且 common polynomial filter 不是所有 local opcode。

**作为更强主定理的第一层：很有价值。** 第 10 节已在 \(\mathbb Z^d\) 上闭合了“任意局部算法下界 \(=\) Green 尾 \(\asymp\) 最佳多项式上界”，因而比只有第 3 节时更接近 headline result。但要将经典成分的组合提升为说服力足够的原创主定理，最好再闭合以下一项：

1. **可局部验证的图/测量 iff。** 从 bounded degree、锚定密度、局部观测几何、heat-kernel/spectral dimension 等有限邻域可见量，精确刻画何时 (B) 成立，并给 converse completion。
2. **超越齐次格点的 instance-optimal lower bound。** 将定理 10.1 从 \(\mathbb Z^d\) 推到非齐次、可随机失败的节点，并保留同一个可局部执行的上界。
3. **自适应且无需全局谱尾先验的算法。** 节点只凭本地 certificates 选择停止时间，在允许 \(\delta\) 失败下达到 oracle tail-rate，并给出匹配通信下界。
4. **topology–spectrum 的 sharp threshold。** 例如在一大类非齐次网络上证明 WLS 的临界低频维数为 2、raw RHS 为 4，包含 sharp 临界/log 修正及必要性，而不仅是规则格点的经典现象。
5. **共同 completion 与多数节点的耦合。** 证明同一个全局稀疏 completion 能同时实现所需坏谱尾，而不是逐节点分别挑选最坏矩阵。

第 7 节的有谱隙 cycle 下界与第 10 节的无谱隙格点幂律共同说明：充要条件的存在性可由 Richardson 见证，但匹配轮数需要真正适配低谱几何的多项式。真正新结果更可能藏在“尾模量 \(\leftrightarrow\) 拓扑 \(\leftrightarrow\) 最优通信轮数”的闭环，而不是 (A)–(C) 的定性等价本身。

### 12.5 当前决策：GO / PENDING / NO-GO

| 命题或论文角色 | 决策 | 理由 |
|---|---|---|
| 定理 3.5：WLS sufficient statistic 后的任意 uniformly bounded finite-range \(K\)，允许非交换 | **GO（作为支撑定理）** | (3.29) 严格、量词完整，明显宽于 common polynomial；但它仍是短的投影不等式，单独不足以承担原创主线 |
| 定理 3.2：class-uniform 多数节点谱尾 iff 与 Richardson 存在性 | **GO（作为统一判据）** | 正确且对 WLS 指数纠错有价值；Backhausz–Virág 与谱逼近文献使其不适合作唯一主定理 |
| 定理 10.1：\(\mathbb Z^d\) 上 all-local minimax \(\asymp\) 最佳 polynomial，并给出临界维数和轮数指数 | **PENDING（强备选主定理）** | 数学上已闭合，且 (10.37)–(10.39) 给出 minimax latency 与标准同步实现的能耗收益；但 Green/Christoffel/GFF 组成均经典，需更广非齐次推广或更完整查重才适合宣称原创 |
| 仅用 spectral dimension 就在一般图上推出 all-local matching round lower bound | **NO-GO** | 谱密度不决定 Green 行的空间尾；还需要 volume/walk dimension、Green annulus，WLS 还需要梯度估计 |
| 将定理 3.5 直接宣称为任意 measurement-to-state local operator \(C\) 的 iff | **NO-GO** | \(C=KG\) 因子约束是 (3.29) 产生 \(B\sqrt\eta\) 的关键；删掉后现有证明失效 |
| 宣称 Richardson 在轮数/能耗上最优 | **NO-GO** | cycle 上它可慢一个 condition-number 平方根；格点无谱隙情形中误差幂指数也只有最佳多项式的一半 |

综合判断：这份谱尾工作适合收口为一条独立备选，不必为了主线强行拔高。它已经提供了一个稳固的 qualitative iff，以及一个有量化应用收益的齐次格点强版。前者不够当原创主定理；后者有潜力，但当前应标为 PENDING 而不是已确认的 headline claim。

---

## 13. 可直接用于论文的安全表述

在未进一步闭合 novelty gap 前，建议使用：

> We establish an input-aware spectral-tail criterion for uniform polynomial locality on a fraction of nodes. The result unifies primitive right-hand-side error and whitened WLS measurement error, and shows that Richardson filters are qualitatively universal under the criterion. The equivalence follows from spectral calculus and is used as an organizing theorem; our originality claim is reserved for the subsequent topology/rate characterization.

对定理 3.5，可以安全地表述为：

> For exact whitened WLS, the necessity direction extends from shared polynomial filters to arbitrary instance-dependent, noncommuting state-space operators that are uniformly stable and act on the normal-equation sufficient statistic. A low-spectral input projection yields the pointwise inequality (3.29).

对定理 10.1，在更全面查重前可以写成“我们将证明/候选贡献”，不要写成已确认的世界首创：

> On lattice benchmarks, the local minimax error is exactly the invisible tail norm of the Green row, while the best degree-constrained distributed filter is an endpoint Christoffel problem. Their asymptotics match, yielding the critical dimensions (2) and (4) and the optimal round/error exponent (r^{-(d-2q)/2}).

不建议使用：

> We completely solve local state estimation for arbitrary local algorithms.

因为定理 3.2 只处理 shared scalar polynomial functional calculus；定理 3.5 虽放宽到 noncommuting \(K\)，仍要求估计器通过 sufficient statistic \(G\xi\) 因子分解。定理 10.1 确实针对任意局部解码，但当前只在齐次格点/其有限 SPD 逼近上闭合，不是对所有网络的完整刻画。

---

## Sources

1. C. D. Godsil and B. Mohar, “[Walk generating functions and spectral measures of infinite graphs](https://doi.org/10.1016/0024-3795(88)90245-5),” *Linear Algebra and its Applications* 107 (1988), 191–206. DOI: 10.1016/0024-3795(88)90245-5.
2. C. Bordenave and M. Lelarge, “[Resolvent of Large Random Graphs](https://doi.org/10.1002/rsa.20313),” *Random Structures & Algorithms* 37(3) (2010), 332–352; [author preprint](https://www.di.ens.fr/~lelarge/papiers/2009/res.pdf). DOI: 10.1002/rsa.20313.
3. D. K. Hammond, P. Vandergheynst, and R. Gribonval, “[Wavelets on Graphs via Spectral Graph Theory](https://doi.org/10.1016/j.acha.2010.04.005),” *Applied and Computational Harmonic Analysis* 30 (2011), 129–150; [arXiv](https://arxiv.org/abs/0912.3848).
4. D. I. Shuman, P. Vandergheynst, D. Kressner, and P. Frossard, “[Distributed Signal Processing via Chebyshev Polynomial Approximation](https://doi.org/10.1109/TSIPN.2018.2824239),” *IEEE Transactions on Signal and Information Processing over Networks* 4(4) (2018), 736–751; [arXiv](https://arxiv.org/abs/1111.5239).
5. S. Demko, W. F. Moss, and P. W. Smith, “[Decay Rates for Inverses of Band Matrices](https://doi.org/10.1090/S0025-5718-1984-0758197-9),” *Mathematics of Computation* 43(168) (1984), 491–499; [AMS PDF](https://www.ams.org/mcom/1984-43-168/S0025-5718-1984-0758197-9/S0025-5718-1984-0758197-9.pdf).
6. M. Benzi and N. Razouk, “[Decay Bounds and O(n) Algorithms for Approximating Functions of Sparse Matrices](https://etna.ricam.oeaw.ac.at/vol.28.2007/pp16-39.dir/pp16-39.pdf),” *Electronic Transactions on Numerical Analysis* 28 (2007), 16–39.
7. J. Jiang and D. B. Tay, “[Decentralised signal processing on graphs via matrix inverse approximation](https://doi.org/10.1016/j.sigpro.2019.07.010),” *Signal Processing* 165 (2019), 292–302. DOI: 10.1016/j.sigpro.2019.07.010.
8. N. Emirov, C. Cheng, J. Jiang, and Q. Sun, “[Polynomial graph filters of multiple shifts and distributed implementation of inverse filtering](https://doi.org/10.1007/s43670-021-00019-x),” *Sampling Theory, Signal Processing, and Data Analysis* 20 (2022), article 2; [arXiv](https://arxiv.org/abs/2003.11152).
9. R. Andreev, “[Tikhonov and Landweber convergence rates: characterization by interpolation spaces](https://doi.org/10.1088/0266-5611/31/10/105007),” *Inverse Problems* 31 (2015), 105007; [arXiv](https://arxiv.org/abs/1503.05742).
10. T. Hohage and F. Weidling, “[Characterizations of Variational Source Conditions, Converse Results, and Maxisets of Spectral Regularization Methods](https://doi.org/10.1137/16M1067445),” *SIAM Journal on Numerical Analysis* 55(4) (2017), 1711–1738.
11. E. Nelson, “[Notes on non-commutative integration](https://doi.org/10.1016/0022-1236(74)90014-7),” *Journal of Functional Analysis* 15 (1974), 103–116. DOI: 10.1016/0022-1236(74)90014-7.
12. F. Chung and O. Simpson, “[Solving Local Linear Systems with Boundary Conditions Using Heat Kernel Pagerank](https://arxiv.org/abs/1503.03157),” 2015 preprint/extended work on local Dirichlet Laplacian systems.
13. Q. Yang, Z. Zhang, and M. Fu, “[Distributed Weighted Least-squares Estimation for Networked Systems with Edge Measurements](https://doi.org/10.1016/j.automatica.2020.109091),” *Automatica* 120 (2020), 109091; [arXiv](https://arxiv.org/abs/2002.11221).
14. M. Embree, J. A. Henningsen, J. Jackson, and R. B. Morgan, “[Polynomial Approximation to the Inverse of a Large Matrix](https://doi.org/10.1137/24M1677599),” *SIAM Journal on Scientific Computing*, published online 2025.
15. M. Ghadiri, H.-A. Nguyen, and J. Yang, “[Entrywise Approximation for Matrix Inversion and Linear Systems](https://doi.org/10.1137/1.9781611978971.188),” *Proceedings of SODA 2026*, 5199–5210.
16. Á. Backhausz and B. Virág, “[Spectral measures of factor of i.i.d. processes on vertex-transitive graphs](https://doi.org/10.1214/16-AIHP790),” *Annales de l'Institut Henri Poincaré, Probabilités et Statistiques* 53(4) (2017), 2260–2278; [author/repository PDF](https://real.mtak.hu/74276/1/1505.07412v2.pdf). DOI: 10.1214/16-AIHP790.
17. T. Danka and V. Totik, “[Christoffel Functions with Power Type Weights](https://doi.org/10.4171/JEMS/776),” *Journal of the European Mathematical Society* 20(3) (2018), 747–796; [arXiv](https://arxiv.org/abs/1504.03968). DOI: 10.4171/JEMS/776.
18. T. Danka, “[Universality limits for generalized Jacobi measures](https://doi.org/10.1016/j.aim.2017.06.026),” *Advances in Mathematics* 316 (2017), 613–666; [arXiv](https://arxiv.org/abs/1605.04275). DOI: 10.1016/j.aim.2017.06.026.
19. K. Uchiyama, “[Green's Functions for Random Walks on \(\mathbb Z^N\)](https://doi.org/10.1112/S0024611598000458),” *Proceedings of the London Mathematical Society* 77(1) (1998), 215–240. DOI: 10.1112/S0024611598000458.
20. A. Cipriani, “[High points for the membrane model in the critical dimension](https://doi.org/10.1214/EJP.v18-2750),” *Electronic Journal of Probability* 18 (2013), paper 86, 1–17; [arXiv](https://arxiv.org/abs/1303.6792). DOI: 10.1214/EJP.v18-2750.
21. L. Knockaert, “[A note on the relationship between the conjugate gradient method and polynomials orthogonal over the spectrum of a linear operator](https://doi.org/10.1109/TAP.1987.1144227),” *IEEE Transactions on Antennas and Propagation* 35(9) (1987), 1089–1091. DOI: 10.1109/TAP.1987.1144227.
22. M. Abért, A. Thom, and B. Virág, “[Benjamini–Schramm convergence and pointwise convergence of the spectral measure](https://users.renyi.hu/~abert/luckapprox.pdf),” preprint (2011/2013 versions circulated). This is a primary preprint rather than a journal DOI source.
23. J. W. Daniel, “[The Conjugate Gradient Method for Linear and Nonlinear Operator Equations](https://doi.org/10.1137/0704002),” *SIAM Journal on Numerical Analysis* 4(1) (1967), 10–26. DOI: 10.1137/0704002. Proposition 1.2.1 already formulates the optimal residual polynomial through a spectral-measure integral and orthogonality.
24. N. Motee and Q. Sun, “[Sparsity and Spatial Localization Measures for Spatially Distributed Systems](https://doi.org/10.1137/15M1049294),” *SIAM Journal on Control and Optimization* 55(1) (2017), 200–235; [arXiv](https://arxiv.org/abs/1402.4148). DOI: 10.1137/15M1049294.
25. J. Špakula and J. Zhang, “[Quasi-Locality and Property A](https://doi.org/10.1016/j.jfa.2019.108299),” *Journal of Functional Analysis* 278(1) (2020), 108299; [arXiv](https://arxiv.org/abs/1809.00532). DOI: 10.1016/j.jfa.2019.108299.
26. K. Li, J. Špakula, and J. Zhang, “[Measured Asymptotic Expanders and Rigidity for Roe Algebras](https://doi.org/10.1093/imrn/rnac242),” *International Mathematics Research Notices* 2023(17), 15102–15154; [arXiv](https://arxiv.org/abs/2010.10749). DOI: 10.1093/imrn/rnac242.
