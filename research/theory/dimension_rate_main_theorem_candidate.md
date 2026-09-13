# 局部状态估计的临界维数与最优轮数率：候选主定理

> 状态：研究草案，2026-09-13。本文记录一条比定性谱尾 iff 更强的候选路线；尚未完成全部有限体积一致常数证明与专项原创性审计，不能作为既成结果引用。

## 1. 为什么考虑这条路线

仅证明“存在某个有限度多项式”还没有回答通信轮数是否最优。一个够强的主定理应在同一模型内同时给出：

1. 什么维度/低频结构下局部估计可能；
2. 半径为 \(r\) 时任何局部算法都不能超过的误差下界；
3. 一个每条边每轮只传常数个标量的算法达到同阶上界；
4. 失败时给出有限网络序列上的发散反例。

欧氏型网络上，空间 Green 核尾给出信息论下界，根谱测度的 Christoffel 极值问题给出多项式消息传递上界。两端预言相同的临界维数与误差指数。

## 2. 两种输入模型

令 \(L\) 是 \(\mathbb Z^d\) 上的最近邻图 Laplacian，根为 \(o\)。先在有限支撑输入上定义目标，再问根坐标是否能连续延拓到单位 \(\ell_2\) 输入球。

### 2.1 普通正规方程右端（\(q=2\)）

给定 \(b\in\ell_2(V)\)，目标根坐标为

\[
x_o=\langle L^{-1}\delta_o,b\rangle .
\]

其代表核是 Green 行 \(g_o=L^{-1}\delta_o\)，故完整根坐标在 \(\ell_2\) 输入上有界，当且仅当

\[
\|g_o\|_2^2
=\int_{(0,4d]}\lambda^{-2}\,d\nu_o(\lambda)<\infty .
\tag{2.1}
\]

这里 \(\nu_o\) 是 \(L\) 在根的谱测度。

### 2.2 白化 WLS 测量（\(q=1\)）

取局部测量算子 \(H\) 满足 \(H^*H=L\)，例如带定向的边—点 incidence 算子。对白化输入 \(\xi\in\ell_2(E)\)，目标为

\[
x_o=\langle HL^{-1}\delta_o,\xi\rangle .
\]

代表核是 Green 函数的边梯度，且

\[
\|HL^{-1}\delta_o\|_2^2
=\int_{(0,4d]}\lambda^{-1}\,d\nu_o(\lambda).
\tag{2.2}
\]

所以两种情形可用 \(q\in\{1,2\}\) 统一：根误差的谱权重为 \(\lambda^{-q}\)。

## 3. 任意局部算法的精确信息下界

令 \(\mathcal I_r(o)\) 是根在 \(r\) 轮后可见的输入坐标。设目标线性泛函的 Riesz 核为 \(k_o\)。即使允许任意大消息、任意本地计算和非线性输出，只要算法只能观察 \(\mathcal I_r(o)\)，就有

\[
\inf_{\phi}
\sup_{\|u\|_2\le1}
\bigl|\langle k_o,u\rangle-\phi(u_{\mathcal I_r(o)})\bigr|
=\|k_{o,\mathcal I_r(o)^c}\|_2.
\tag{3.1}
\]

下界取两个不可区分输入
\(u=\pm k_{o,\mathcal I_r(o)^c}/\|k_{o,\mathcal I_r(o)^c}\|_2\)；上界由截断线性泛函达到。若输入是独立标准 Gaussian，则所有可测局部估计器的最小 MSE 同样是

\[
\inf_\phi \mathbb E
\bigl|\langle k_o,u\rangle-\phi(u_{\mathcal I_r(o)})\bigr|^2
=\|k_{o,\mathcal I_r(o)^c}\|_2^2,
\tag{3.2}
\]

因为条件期望就是可见坐标上的截断线性泛函。

式 (3.1)--(3.2) 本身是经典 information radius/条件期望事实；候选创新不能建立在这一步上。它的作用是让后面的 rate converse 覆盖**所有** \(r\)-local 算法，而不只覆盖多项式滤波器。

## 4. 多项式算法的精确谱极值

度数至多 \(r-1\) 的共同滤波器输出 \(p(L)b\)（普通右端）或 \(p(L)H^*\xi\)（白化 WLS）。它在根的平方误差为

\[
\mathcal R_{r,q}(p)^2
=\int |1-\lambda p(\lambda)|^2\lambda^{-q}
\,d\nu_o(\lambda).
\tag{4.1}
\]

置 \(s(\lambda)=1-\lambda p(\lambda)\)。则 \(\deg s\le r\) 且 \(s(0)=1\)。定义有限测度

\[
d\rho_q(\lambda)=\lambda^{-q}d\nu_o(\lambda).
\]

当 \(\rho_q\) 有限时，最优 \(r\)-轮多项式误差满足精确恒等式

\[
\boxed{
\inf_{\deg p\le r-1}\mathcal R_{r,q}(p)^2
=\inf_{\substack{\deg s\le r\\s(0)=1}}
\int s(\lambda)^2d\rho_q(\lambda)
=\frac1{K_r^{\rho_q}(0,0)} .}
\tag{4.2}
\]

\(K_r^{\rho_q}\) 是次数不超过 \(r\) 的 Christoffel--Darboux 核。最优残差和滤波器可显式写成

\[
s_r^*(\lambda)
=\frac{K_r^{\rho_q}(\lambda,0)}{K_r^{\rho_q}(0,0)},
\qquad
p_{r-1}^*(\lambda)=\frac{1-s_r^*(\lambda)}{\lambda}.
\tag{4.3}
\]

式 (4.2) 是经典正交多项式极值公式；候选新意只能来自它与局部信息下界的匹配，以及网络状态估计中的临界维数/资源结论。

## 5. 候选主定理：欧氏型周期网络上的 matching rate

### 5.1 预期正式陈述

对 \(q\in\{1,2\}\)，考虑 \(d\) 维最近邻格点，或更一般的固定胞元、有限作用距离、均匀椭圆的周期网络。假设最低 Bloch 谱带在唯一的零频点具有非退化二次极小值。则：

1. 根目标泛函在单位 \(\ell_2\) 输入球上有界，当且仅当
   \[
   d>2q.
   \tag{5.1}
   \]
2. 当 \(d>2q\) 时，任意 \(r\)-轮局部算法的最坏情形误差和 Gaussian RMS 误差均满足
   \[
   \mathcal E^{\rm all}_{r,q}
   \ge c\,r^{-(d-2q)/2}.
   \tag{5.2}
   \]
3. 存在同一个度数 \(r-1\) 的线性多项式消息传递算法，每条边每轮只传常数个标量，并满足
   \[
   \mathcal E^{\rm poly}_{r,q}
   \le C\,r^{-(d-2q)/2}.
   \tag{5.3}
   \]
4. 因为多项式算法属于全部局部算法，(5.2)--(5.3) 给出
   \[
   \boxed{
   \mathcal E^{\rm all}_{r,q}
   \asymp
   \mathcal E^{\rm poly}_{r,q}
   \asymp r^{-(d-2q)/2}.}
   \tag{5.4}
   \]
5. 当 \(d\le2q\) 时，可取越来越大的周期网络和趋零的统一 grounding，使任何固定 \(r\) 的局部误差无统一有限上界；临界维数处为对数发散，低维为幂次发散。

特别地：

- 白化 WLS（\(q=1\)）：临界维数是 \(2\)，最优误差为
  \(\Theta(r^{-(d-2)/2})\)；
- 普通 RHS（\(q=2\)）：临界维数是 \(4\)，最优误差为
  \(\Theta(r^{-(d-4)/2})\)。

### 5.2 下界为何给出这个指数

格点 Green 核及其梯度满足

\[
G(o,x)\asymp |x|^{2-d},
\qquad
|\nabla G(o,x)|\asymp |x|^{1-d}
\tag{5.5}
\]

（方向性下界在正式证明中需按锥区或平均平方陈述，不能逐边无条件写成双边点态界）。因此普通 RHS 的不可见平方尾为

\[
\sum_{|x|>r}|G(o,x)|^2
\asymp r^{4-d},
\tag{5.6}
\]

而白化 WLS 的不可见边核平方尾为

\[
\sum_{e:\,d(o,e)>r}|\nabla_eG(o,\cdot)|^2
\asymp r^{2-d}.
\tag{5.7}
\]

开平方后正是 (5.2)。在 \(d=4\) 或 \(d=2\) 时相应和式按 \(\log(N/r)\) 增长。

### 5.3 上界为何给出同一个指数

周期椭圆网络的根谱密度在零附近满足

\[
\frac{d\nu_o}{d\lambda}
\sim c\lambda^{d/2-1}.
\tag{5.8}
\]

故 \(\rho_q\) 的端点权重指数是

\[
a=d/2-1-q>-1.
\]

端点 Christoffel 函数的标准渐近给出

\[
\frac1{K_r^{\rho_q}(0,0)}
\asymp r^{-2a-2}
=r^{-(d-2q)}.
\tag{5.9}
\]

结合 (4.2)，开平方得到 (5.3)。对于纯权重 \(d\rho(\lambda)=\lambda^a d\lambda\)（区间缩放略去），还可直接核对

\[
K_r^{\rho}(0,0)
=\frac1{a+1}
\left[
\frac{\Gamma(r+a+2)}{\Gamma(a+1)\Gamma(r+1)}
\right]^2,
\tag{5.10}
\]

从而 \(K_r(0,0)^{-1}\sim\Gamma(a+1)\Gamma(a+2)r^{-2a-2}\)。

## 6. 与“允许少量节点失败”的连接

周期网络每个同类型根具有相同风险，所以允许固定比例失败不会改变 (5.1)--(5.4)。它提供一个不能靠牺牲少数节点绕开的 matching lower bound。

对异质网络族，已有的节点谱尾 iff 可作为外层定性结论：若除至多 \(\delta\) 比例节点外，局部谱密度在零点具有统一指数 \(\alpha>q\)，则共同多项式在好节点上的预期 rate 为

\[
r^{-(\alpha-q)}
\tag{6.1}
\]

（这里 \(\alpha=d_s/2\)，式 (6.1) 是误差而非平方误差）。要把 (6.1) 升级为异质网络上对全部局部算法的 sharp rate，还需要相同好节点上的空间 Green 尾双边界；仅有谱测度不够推出任意局部算法的空间下界。

## 7. 资源含义

截断 oracle 达到信息下界，但若每个根独立泛洪整个 \(r\)-球，需搬运随球体积增长的数据。式 (4.3) 的多项式可用三项递推或 Horner 形式实现：每轮每条边传常数个标量，节点只保留常数个工作向量。于是 (5.4) 的意义不只是误差率相同：一个受限得多的标量消息算法达到了允许无限消息的 LOCAL 下界。

若每条边每个标量的能耗记为常数，则达到误差 \(\varepsilon\) 所需轮数满足

\[
r(\varepsilon)
=\Theta\!\left(
\varepsilon^{-2/(d-2q)}
\right),
\tag{7.1}
\]

总 bit-hop/能耗还需乘以活跃边数与量化位宽。式 (7.1) 是候选 theorem 可支撑的通信—精度 Pareto law；实际无线能耗不能只由图轮数推出，必须另设物理层模型。

## 8. 当前证明缺口

在把第 5 节称为定理前必须关闭：

1. 用有限 \(d\)-torus 加 \(\tau_N\asymp N^{-2}\) 或严格 mean-zero gauge 写出无歧义的有限问题，并证明 \(N\to\infty\)、\(r\to\infty\) 的量词顺序；
2. 对白化 WLS 明确输入坐标的通信半径，并处理 \(H=[B;\sqrt\tau I]\) 的 grounding 列；
3. 证明有限体积 Green/gradient tail 的统一上下界，临界维数给出对数项；
4. 证明同一个 Christoffel/Jacobi 多项式对有限周期网络的离散谱和极限谱具有所需一致上界；
5. 给出稳定三项递推及有限精度/位宽版本，否则只能声称实数 LOCAL；
6. 将周期 block-WLS 的 Bloch 符号条件写成可核验假设，避免把“欧氏型”当作未定义词。

## 9. 原创性碰撞清单

以下部分已知，不能单独报新：

- information radius 与 Gaussian 条件方差给出 (3.1)--(3.2)；
- Christoffel 函数给出 (4.2)；
- Green 核在格点上的衰减及 GFF 的临界维数是经典概率论；
- Backhausz--Virág 已证明 transient vertex-transitive graph 上 Gaussian free field 是 linear factor of iid，并用谱密度和多项式逼近；这会直接碰撞 \(q=1\) 的定性存在性；
- polynomial graph filtering、Richardson/Landweber 与分布式 inverse filtering 已有成熟文献。

真正需要专项确认的是：是否已有工作在**同一个局部状态估计/线性求解信息模型**中，把任意 \(r\)-local 算法的空间下界、常数消息多项式上界、\(q=1/2\) 两个临界维数及 exact round exponent 同时证明。只有该 matching theorem 未被现成结果直接包含，才有资格作为主定理。

## 10. 当前判决

这条路线在数学形态上比单纯的多数节点谱尾 iff 更强，因为它给出 all-local lower bound 与可执行算法的 matching rate。但每个证明模块都与成熟理论相邻，原创性风险仍高。

当前标签：**强主定理候选，正确性与原创性均 PENDING。** 若有限体积一致界或文献审计不能关闭，则退回为格点示例/背景，不进入论文主张。

