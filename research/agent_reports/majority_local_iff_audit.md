# 允许少量节点失败时的局部估计充要条件：审计报告

**日期：** 2026-09-13  
**范围：** 共同 \(r\)-轮局部规则、模型 completion 不可区分性、任意/线性解码器、\((1-\delta)\)-节点语义，以及一个针对有限传播 SPD 逆问题的谱条件。  
**原创性口径：** “没有检索到同句定理”不等于原创；下文逐条标为标准事实、直接推论、可能的新组合或未闭合。

## 0. 结论先行

这次审计得到两个层次完全不同的结论。

1. **对任意模型族、任意局部规则追求“最宽”充要条件，结论是 NO-GO。** 精确答案可以写成约束族删点后的共同可行性，有限情形又等价于带每个模型删点预算的冲突超图横截。欧氏输出或线性系数下，Helly 定理把最小冲突大小压到有限维数加一。但这基本是 tolerant Helly、maximum feasible subsystem、hitting set 和 optimal recovery 的直接组合，不足以当论文主定理。更严重的是，一个无数据、有限字母的特例就精确等于 **Closest String**；因此只看逐节点信息半径、局部 view 的边际 completion 集或 view 频率，原则上不可能给出充要条件。

2. **在“共同多项式图滤波器逼近 SPD 逆”这一自然但更窄的类上，出现了一个真正的结构 iff 候选。** 对统一有界、有限传播的 SPD 矩阵族 \(J_\theta\)，节点 \(i\) 的最坏输入行误差满足

   \[
   \|e_i^\top(J_\theta^{-1}-p(J_\theta))\|_2^2
   =\int |\lambda^{-1}-p(\lambda)|^2\,d\mu_{\theta,i}(\lambda).
   \]

   存在同一个有限次数多项式 \(p\)，使任意给定精度下除任意小比例节点外都准确，**当且仅当** 零附近的节点谱质量

   \[
   X_{\theta,i}(\eta)
   :=\int_{(0,\eta)}\lambda^{-2}\,d\mu_{\theta,i}(\lambda)
   \]

   对均匀随机节点在整个模型族上一致依概率趋于零。更强的是，不必搜索最优多项式：固定步长 Richardson/Neumann 多项式

   \[
   p_t(\lambda)=\frac{1-(1-\lambda/\Lambda)^t}{\lambda}
   \]

   已经对所有满足该条件的模型族定性普适。必要性则对每个固定多项式在足够小的特征值上给出匹配下界。它允许条件数发散和少量局部病态区，不需要全局谱隙。对白化 WLS 输入，相同 iff 把 \(\lambda^{-2}\) 换成 \(\lambda^{-1}\)。

   这个定理比 Chebyshev radius 或 quantile 恒等式强：它给出谱结构条件、实际有限轮算法和匹配反例。但它只刻画**共同多项式滤波器**，不是所有可能的 view-dependent 局部算法；截至本次定向检索，没有找到完全相同的“节点谱测度 \(\lambda^{-q}\) 尾 iff 多数节点局部逆逼近”先例。其各个证明部件均经典，所以目前应标为 **主定理形态达到、原创性待专项查重**，不能直接宣称原创。

## 1. 普通语言下的统一模型

### 1.1 网络、模型与 \(r\)-轮 view

令 \(\Theta\) 是一族有限网络模型。每个 \(\theta\in\Theta\) 包含：

- 有限节点集 \(V_\theta\) 和节点概率测度 \(m_\theta\)；均匀节点语义取 \(m_\theta(i)=1/|V_\theta|\)；
- 带端口、方向、边/节点参数及本地维数的通信图；
- 允许的数据集合 \(\mathcal X_\theta\)；
- 每个节点 \(i\) 的目标映射 \(F_{\theta,i}:\mathcal X_\theta\to Y_v\)。

一个发生项记为 \(\alpha=(\theta,i)\)。固定轮数 \(r\) 后，\(v_\alpha\) 表示节点在 \(r\) 轮内能够得到的**模型 view**：即完整确定性 transcript 的等价类，但不包含尚未给定的运行时数据。无界消息的标准 LOCAL 模型里，它就是带根、带端口和本地参数的半径 \(r\) 标记球。相同 view 必须同时给出输入空间 \(U_v\) 与输出空间 \(Y_v\) 的规范识别。

节点实际收到的本地数据写成

\[
q_\alpha:\mathcal X_\theta\longrightarrow U_{v_\alpha}.
\]

“所有 completion 使用同一个局部算法”的精确含义是：对每个 view 类 \(v\) 只选一个解码器

\[
f_v:U_v\to Y_v,
\]

发生项 \(\alpha\) 必须输出 \(f_{v_\alpha}(q_\alpha(x))\)。若模型允许唯一 ID，ID 必须是 view 的一部分；若研究目标是拓扑/系数局部性，则应采用匿名、端口标号或 order-invariant 版本，否则 ID 会人为拆散本来相同的 view。相同 \(r\)-view 导致相同确定性输出，是 LOCAL 模型的基本不可区分性，而不是 WLS 的特殊性质。[^1]

### 1.2 三种不能混用的 \(\delta\) 语义

取损失 \(d_v\) 和阈值 \(\varepsilon\)。定义发生项的逐节点 robust 误差

\[
E_\alpha(f)
:=\sup_{x\in\mathcal X_\theta}
d_v\bigl(f_{v_\alpha}(q_\alpha x),F_{\theta,i}(x)\bigr).
\tag{1.1}
\]

以下三种陈述不等价。

**R：固定失败节点、对数据稳健。**

\[
\exists f\ \forall\theta\ \exists B_\theta\subseteq V_\theta:
m_\theta(B_\theta)\le\delta,
\quad
\forall i\notin B_\theta\ \forall x\in\mathcal X_\theta,\ 
d(f(q_i x),F_i x)\le\varepsilon.
\tag{1.2R}
\]

**I：失败节点可随数据实例变化。**

\[
\exists f\ \forall\theta\ \forall x\ \exists B_{\theta,x}:
m_\theta(B_{\theta,x})\le\delta,
\quad
\forall i\notin B_{\theta,x},\
d(f(q_i x),F_i x)\le\varepsilon.
\tag{1.2I}
\]

**B：Bayes/联合平均。** 给定数据律 \(P_\theta\)，要求

\[
(m_\theta\otimes P_\theta)
\{(i,x):d(f(q_i x),F_i x)>\varepsilon\}\le\delta.
\tag{1.2B}
\]

有 \((1.2R)\Rightarrow(1.2I)\)，后者一般再强于平均语义；反向均不成立。本报告第 2--7 节默认最强的 R 语义。I 语义可以把 \((\theta,x)\) 当作场景组重新写同样的有限约束，但其坏节点预算按数据实例变化；B 语义则是风险积分问题，不能直接套用 R 的删点结论。

## 2. 任意解码器：精确条件及其信息含量

### 2.1 精确共同可行性公式

令 \(\mathcal F=\prod_v Y_v^{U_v}\) 是所有集合意义下的共同解码器。对发生项 \(\alpha\) 定义

\[
C_\alpha(\varepsilon)
:=\{f\in\mathcal F:E_\alpha(f)\le\varepsilon\}.
\tag{2.1}
\]

**命题 2.1（精确但定义性的 iff；标准约束改写）。** R 语义成立，当且仅当

\[
\bigcap_{\theta\in\Theta}
\ \bigcup_{\substack{B\subseteq V_\theta\\m_\theta(B)\le\delta}}
\ \bigcap_{i\in V_\theta\setminus B}
C_{(\theta,i)}(\varepsilon)
\neq\varnothing.
\tag{2.2}
\]

**证明。** 若 \(f\) 满足 R，取每个模型的坏节点集 \(B_\theta\)，则 \(f\) 属于右侧对应的每个交集。反之，右侧任一元素同时给出共同规则及每个模型可删的坏节点集。∎

式 (2.2) 是量词无歧义的“最宽充要条件”，但没有提供比原问题更多的结构，不能作为主定理。

### 2.2 有限模型目录的冲突超图

以下假设发生项集合

\[
\Omega=\{(\theta,i):\theta\in\Theta,\ i\in V_\theta\}
\]

有限。把 \(A\subseteq\Omega\) 称为冲突，如果 \(\bigcap_{\alpha\in A}C_\alpha=\varnothing\)；只保留按包含关系极小的冲突作为超边，得超图 \(\mathcal H_\varepsilon\)。均匀节点时令 \(k_\theta=\lfloor\delta|V_\theta|\rfloor\)。

**定理 2.2（删点 iff；直接的 hitting-set 表述）。** 存在共同规则，使每个模型至多 \(k_\theta\) 个节点失败，当且仅当存在 \(B\subseteq\Omega\)，满足

\[
|B\cap(\{\theta\}\times V_\theta)|\le k_\theta
\quad\text{对每个 }\theta,
\tag{2.3}
\]

并且 \(B\) 与 \(\mathcal H_\varepsilon\) 的每条超边相交。

**证明。** 若共同规则 \(f\) 存在，令 \(B\) 收集所有不满足约束的发生项。任何完全落在 \(\Omega\setminus B\) 的冲突边都会与 \(f\) 同时满足其所有约束矛盾，所以 \(B\) 必击中每条边。反之，若 \(B\) 击中所有极小冲突，则幸存约束族不含任何不可行子族；有限性保证其全交非空，取交中的 \(f\) 即得所需规则。∎

加权节点只需把 (2.3) 换成每组权重预算。无限目录下，(2.2) 仍精确；“所有有限子族可行推出全族可行”还需要紧性、系数有界性或其他有限交性质，不能无条件援引。

这一定理已经保留 completion 的**联合共现结构**，因而不是逐节点 quantile。但它也是 tolerant feasibility 的标准横截形式；Montejano--Oliveros 的 tolerance Helly 定理已经明确使用超图横截刻画删去 \(t\) 个集合后的公共交。[^2]

## 3. Helly 给出的有限冲突证书

### 3.1 任意非线性规则，有限维输出

假定 \(Y_v=\mathbb R^{p_v}\)，损失来自任意范数。对发生项 \(\alpha\) 和本地输入值 \(u\in U_v\)，定义

\[
D_\alpha(u)
:=\bigcap_{\substack{x\in\mathcal X_\theta\\q_\alpha(x)=u}}
\overline B\bigl(F_\alpha(x),\varepsilon\bigr),
\tag{3.1}
\]

空 fiber 时约定 \(D_\alpha(u)=Y_v\)。它是闭凸集。

**命题 3.1（冲突秩上界；Helly 的直接推论）。** 对同一 view 的有限发生项集 \(A\)，有

\[
\bigcap_{\alpha\in A}C_\alpha\neq\varnothing
\quad\Longleftrightarrow\quad
\forall u\in U_v,\
\bigcap_{\alpha\in A}D_\alpha(u)\neq\varnothing.
\tag{3.2}
\]

因此每个极小冲突都只含同一 view 的发生项，且大小至多 \(p_v+1\)。

**证明。** 一个规则满足 \(A\) 的全部约束，当且仅当对每个 \(u\)，它的值 \(f_v(u)\) 同时落入所有 \(D_\alpha(u)\)。若每个交集非空，可逐 \(u\) 选择一点得到集合意义下的规则；若不可行，则存在某个 \(u\) 的凸集交为空。有限维 Helly 定理给出至多 \(p_v+1\) 个集合已经交空。不同 view 使用 \(\mathcal F\) 的不同直积分量，不能共同形成极小冲突。∎

这里的“构造”在无限 \(U_v\) 上可能依赖选择公理，未自动保证可测、连续或可计算。有限 \(U_v\) 且有凸可行性 oracle 时才直接成为算法。经典 fractional/tolerant Helly 文献已经覆盖“局部小冲突证书 + 容忍删除”这一几何骨架。[^2][^3]

### 3.2 线性规则与 WLS

令 \(U_v,Y_v\) 为有限维实 Hilbert 空间，规则限制为 \(L_v\in\operatorname{Hom}(U_v,Y_v)\)。对发生项 \(\alpha=(\theta,i)\)，令全局数据空间为 \(Z_\theta\)，本地可见算子为 \(S_\alpha:Z_\theta\to U_v\)，中心目标块行为 \(T_\alpha:Z_\theta\to Y_v\)。在输入统一归一化为 \(\|z\|_2\le1\) 时，约束为

\[
K_\alpha(\varepsilon)
=\{L_v:\|T_\alpha-L_vS_\alpha\|_{2\to2}\le\varepsilon\}.
\tag{3.3}
\]

**定理 3.2（线性共同规则的有限冲突 iff；直接凸几何推论）。** 令

\[
d_v=\dim_{\mathbb R}\operatorname{Hom}(U_v,Y_v)
=\dim U_v\cdot\dim Y_v.
\]

只把同一 view 中大小不超过 \(d_v+1\)、且 (3.3) 共同不可行的发生项集列为超边。则定理 2.2 的分组预算横截条件与线性共同规则存在性等价。

**证明。** 每个 \(K_\alpha\) 是系数空间中的闭凸集。若一个同-view 约束族不可行，Helly 给出至多 \(d_v+1\) 个已经不可行；之后逐字应用定理 2.2。复数系数必须按实维数计，即通常加倍。∎

谱范数约束可写成 LMI

\[
\begin{bmatrix}
\varepsilon I_{Y_v}&T_\alpha-L_vS_\alpha\\
(T_\alpha-L_vS_\alpha)^*&\varepsilon I_{Z_\theta}
\end{bmatrix}\succeq0.
\tag{3.4}
\]

因此有限目录中，小冲突可由 SDP 检查；随后解带分组容量的 hitting set，再按 view 解幸存凸约束，就构造出运行时相同的局部系数。总删除数为 \(K\)、最大边大小为 \(h\) 时，对一条未击中的边逐顶点分支给出深度 \(K\)、分支数至多 \(h\) 的参数化算法。不过边枚举、hitting set 和全目录综合都是**离线全局计算**，不是从单个 \(r\)-view 可验证的证书。

在线性 WLS 中，必须区分原始观测与法方程右端：

\[
J_\theta=H_\theta^*R_\theta^{-1}H_\theta+J_{0,\theta},
\qquad
T_\theta=J_\theta^{-1}H_\theta^*R_\theta^{-1}.
\tag{3.5}
\]

若把 \(b=H^*R^{-1}z\) 当作 primitive input，则目标是 \(J^{-1}b\)；若输入是原始 \(z\)，(3.3) 中必须使用 (3.5) 的整行，并让 \(S_\alpha\) 选择一致的原始测量坐标。抽象 iff 不需要“图无环”或“每个 \(A_i\) 满列秩”；这些只是某篇算法保证可行的强充分条件。

跨不确定模型选择一个共同近似逆/共同线性估计器本身已有成熟 robust minimax 与 SDP 框架，尤其 El Ghaoui 2002 的 Eq. (1.4)、Theorem 6.2/Eq. (6.3) 和 Eq. (6.5)。[^4] 本节增加的是 view 绑定和节点删点，但数学机器并非新机器。

## 4. 为什么多数节点版本通常不是逐节点半径的分位数

定义发生项自己的最优半径

\[
\rho_\alpha:=\inf_g E_\alpha(g),
\tag{4.1}
\]

以及同一 view 跨所有 completion 的半径

\[
a(v):=\inf_g\sup_{\alpha:v_\alpha=v}E_\alpha(g).
\tag{4.2}
\]

单个 information fiber 的 arbitrary-decoder 最优误差等于目标像的 Chebyshev radius，是 optimal recovery 的经典结论；本节关心的是共同规则和每模型删点预算造成的额外量词耦合。[^18]

对节点函数 \(h_i\)，取上分位数

\[
Q_{1-\delta}^{m_\theta}(h)
:=\inf\{t:m_\theta\{i:h_i>t\}\le\delta\}.
\tag{4.3}
\]

共同规则的最优多数节点误差是

\[
A_{r,\delta}:=\inf_f\sup_\theta
Q_{1-\delta}^{m_\theta}\bigl(E_{(\theta,i)}(f)\bigr).
\tag{4.4}
\]

**命题 4.1（总是成立的 sandwich；直接推论）。** 在允许任意共同规则且各 infimum 可用任意小松弛逼近时，

\[
\sup_\theta Q_{1-\delta}(\rho_{\theta,i})
\ \le\ A_{r,\delta}\ \le
\sup_\theta Q_{1-\delta}(a(v_{\theta,i})).
\tag{4.5}
\]

**证明。** 任意 \(f\) 对每个发生项都有 \(E_\alpha(f)\ge\rho_\alpha\)，给出左界。对每个 view 选择 \(a(v)+\gamma\) 最优的 \(g_v\)，则所有该 view 发生项误差不超过 \(a(v)+\gamma\)，给出右界；令 \(\gamma\downarrow0\)。∎

左端允许每个 completion/node 单独选规则，右端要求一个 view 在所有 completion 上同时最坏；真实答案还要考虑每个**全局模型内**哪些 view 的困难 completion 能否同时出现，所以一般严格夹在中间。

### 4.1 两个 completion 类的同边际反例

取两个 view 类型 \(a,b\)，每个模型各占一半节点；没有运行时数据；标量目标在 \(\{-1,+1\}\)，要求 \(\varepsilon=0\)。

- 相关 completion 类
  \[
  \Theta_{\rm corr}=\{(+,+),(-,-)\}.
  \]
  共同规则 \(f(a)=+1,f(b)=-1\) 在两个模型里都恰好成功一半节点，所以 \(\delta=1/2\) 可行。
- 矩形 completion 类
  \[
  \Theta_{\rm rect}=\{(+,+),(+,-),(-,+),(-,-)\}.
  \]
  对任意共同输出对，总有逐坐标相反的 completion，使两个节点都失败；所以 \(\delta=1/2\) 不可行。

两个类具有完全相同的 view 频率，并且每个 view 的边际目标集合都是 \(\{-1,+1\}\)，故 \(\rho_\alpha\)、\(a(v)\) 和任何只看单个 view 边际 completion 集的统计都相同。差别仅在**不同 view 的 completion 是否联合共现**。这证明：

> 不保留 completion 的联合耦合，就不可能得到一般类上的多数节点 iff。

### 4.2 何时分位数公式才精确：强拼接假设

对规则 \(f\) 写

\[
w_f(v):=\sup_{\alpha:v_\alpha=v}E_\alpha(f_v).
\]

设模型族具有如下“同 profile 同时实现最坏 completion”性质：对每个可实现的 view 分布 \(\nu\)、每个共同规则 \(f\) 和每个 \(\gamma>0\)，存在一个仍有 profile \(\nu\) 的模型，使匹配到类型 \(v\) 的每个节点误差都至少 \(w_f(v)-\gamma\)。

**命题 4.2（强拼接下的精确公式；可直接证明的充分结构，不是一般必要条件）。** 在上述性质下，

\[
A_{r,\delta}
=\sup_{\nu\ {\rm achievable}}
Q_{1-\delta}^{\nu}(a(v)).
\tag{4.6}
\]

**证明。** 上界是命题 4.1。对任意 \(f\) 和可实现 profile，拼接性质给出一个模型，其误差分位数至少为 \(Q_{1-\delta}(w_f(v))-\gamma\)，又因 \(w_f(v)\ge a(v)\)，得到反向不等式；再对 \(f\) 取 inf、令 \(\gamma\downarrow0\)。∎

这一假设很强：不同 view 的各自最坏远端 completion 必须能在同一个合法网络里同时实现。闭合于不交并只能在允许 root-weighted 模型或均匀密度拼接开销可忽略时帮助；连通、固定几何或物理 WLS 模型一般不满足。故 (4.6) 只能解释“何时 quantile 化成立”，不能作为最宽结构定理。

## 5. 计算下界：一般问题已包含成熟难题

### 5.1 Closest String 是精确特例

取 \(q\) 个等权 view、无数据、输出字母表 \(\Sigma\)，一个 completion 是字 \(c\in\Sigma^q\)，共同规则是字 \(a\in\Sigma^q\)，损失为是否字符相等。则

\[
\text{每个 completion 至少 }(1-\delta)q\text{ 个节点成功}
\quad\Longleftrightarrow\quad
\min_{a\in\Sigma^q}\max_{c\in\mathcal C}d_H(a,c)\le\delta q.
\tag{5.1}
\]

右侧正是 Closest String/Hamming 1-center。该问题即使二元字母表也已知 NP-hard，并有成熟 FPT 文献。[^5] 所以一般 completion-majority 问题不可能靠一个只依赖拓扑的简单标量条件被“彻底终结”；联合 completion 代码本身就能编码组合难度。

### 5.2 线性规则包含 smallest \(k\)-enclosing ball

取一个 view，\(U=\mathbb R^d,Y=\mathbb R\)。一个模型有 \(n\) 个节点，第 \(i\) 个目标为 \(t_i^\top z_i\)，共同线性规则为 \(\ell^\top z_i\)。则节点 operator error 是

\[
\|t_i-\ell\|_2.
\]

允许 \(k\) 个节点失败等价于寻找半径 \(\varepsilon\) 的球覆盖至少 \(n-k\) 个点，即 smallest \(k\)-enclosing ball 的判定版；维数作为输入时强 NP-hard。[^6] 这证明即使线性、单 view、单模型，离线多数节点综合也一般困难。

这还是抽象线性估计归约；把 hardness 完整嵌入“本地可见 \(H,R\) 标签的物理 WLS”尚未证明，不能越级宣称 WLS-specific NP-hardness。

## 6. 一个越过纯 quantile 的谱结构 iff 候选

本节是本报告唯一达到“结构条件 + 构造算法 + 匹配必要性”形态的结果。

### 6.1 假设与精确行误差

对每个 \(\theta\in\Theta\)，令 \(J_\theta\in\mathbb R^{V_\theta\times V_\theta}\) 对称正定，并统一满足

\[
0<J_\theta\le\Lambda I.
\tag{6.1}
\]

不假定统一正下界。若要解释成通信算法，再假定存在固定 \(s\)，使

\[
(J_\theta)_{ij}=0\quad\text{当 }\operatorname{dist}_{G_\theta}(i,j)>s.
\tag{6.2}
\]

节点 \(i\) 的谱测度为

\[
\mu_{\theta,i}(B)
:=\langle e_i,{\bf1}_B(J_\theta)e_i\rangle.
\tag{6.3}
\]

它是 \((0,\Lambda]\) 上的概率测度。对多项式 \(p\)，用 \(p(J_\theta)b\) 估计 \(J_\theta^{-1}b\)。在全局输入 \(\|b\|_2\le1\) 下，节点最坏误差记为

\[
R_{\theta,i}(p)
:=\sup_{\|b\|_2\le1}
|e_i^\top(J_\theta^{-1}-p(J_\theta))b|.
\tag{6.4}
\]

谱定理给出标准恒等式

\[
R_{\theta,i}(p)^2
=\int_{(0,\Lambda]}
|\lambda^{-1}-p(\lambda)|^2\,d\mu_{\theta,i}(\lambda).
\tag{6.5}
\]

定义零附近 inverse-square 节点质量

\[
X_{\theta,i}(\eta)
:=\int_{(0,\eta)}\lambda^{-2}\,d\mu_{\theta,i}(\lambda).
\tag{6.6}
\]

为同时覆盖法方程右端和白化 WLS，固定 \(q>0\)，再定义

\[
\mathcal R_{\theta,i}^{(q)}(p)^2
:=\int_{(0,\Lambda]}
|1-\lambda p(\lambda)|^2\lambda^{-q}\,
d\mu_{\theta,i}(\lambda),
\qquad
X_{\theta,i}^{(q)}(\eta)
:=\int_{(0,\eta)}\lambda^{-q}\,d\mu_{\theta,i}(\lambda).
\tag{6.6a}
\]

于是 (6.4)--(6.6) 正是 \(q=2\)：

\[
\mathcal R_{\theta,i}^{(2)}(p)=R_{\theta,i}(p),
\qquad
X_{\theta,i}^{(2)}(\eta)=X_{\theta,i}(\eta).
\tag{6.6b}
\]

### 6.2 主候选定理

定义 Richardson/Neumann 多项式

\[
p_t(\lambda)
:=\frac{1-(1-\lambda/\Lambda)^t}{\lambda}
=\frac1\Lambda\sum_{k=0}^{t-1}(1-\lambda/\Lambda)^k,
\qquad t\ge1.
\tag{6.6c}
\]

右端说明它确为次数 \(t-1\) 的多项式。

**定理 6.1（多数节点共同多项式逆逼近 iff；证明闭合，原创性待核）。** 固定 \(q>0\)。在 (6.1) 下，以下三条等价：

**(S) 零附近质量在节点测度上一致依概率消失：** 对每个 \(a>0\)，

\[
\lim_{\eta\downarrow0}
\sup_{\theta\in\Theta}
m_\theta\{i:X_{\theta,i}^{(q)}(\eta)>a\}=0.
\tag{6.7}
\]

**(P) 任意精度和任意失败比例都有一个共同有限次数多项式：** 对每个 \(\varepsilon>0,\delta>0\)，存在实多项式 \(p\)，使

\[
\sup_{\theta\in\Theta}
m_\theta\{i:\mathcal R_{\theta,i}^{(q)}(p)>\varepsilon\}
\le\delta.
\tag{6.8}
\]

**(N) 甚至 Richardson/Neumann 这一条固定多项式序列已经足够：** 对每个 \(\varepsilon,\delta>0\)，存在整数 \(t\)，使

\[
\sup_{\theta\in\Theta}
m_\theta\{i:\mathcal R_{\theta,i}^{(q)}(p_t)>\varepsilon\}
\le\delta.
\tag{6.8a}
\]

若另有有限传播 (6.2)，次数 \(d\) 的一般 \(p\) 给出至多 \(ds\) hop 的确定性共同局部算法；\(p_t\) 则由标准 Richardson 递推

\[
x_0=0,\qquad
x_{k+1}=x_k+\Lambda^{-1}(b-Jx_k)
\tag{6.8b}
\]

实现。各节点只需知道共同上界 \(\Lambda\)、本地 \(J\) 条目和逐轮收到的向量值，不需要全局拓扑。

**充分性及显式算法，(S) \(\Rightarrow\) (N)。** 恒等式

\[
1-\lambda p_t(\lambda)=(1-\lambda/\Lambda)^t
\tag{6.9}
\]

给出

\[
\begin{aligned}
\mathcal R_{\theta,i}^{(q)}(p_t)^2
&=\int(1-\lambda/\Lambda)^{2t}
\lambda^{-q}\,d\mu_{\theta,i}(\lambda)\\
&\le
X_{\theta,i}^{(q)}(\eta)
+\frac{(1-\eta/\Lambda)^{2t}}{\eta^q}.
\end{aligned}
\tag{6.10}
\]

固定 \(\varepsilon,\delta\)，先由 (S) 选 \(0<\eta<\Lambda\)，使除质量至多 \(\delta\) 的节点外都有
\(X_{\theta,i}^{(q)}(\eta)\le\varepsilon^2/2\)，再取

\[
t\ge
\frac{\log\!\left(\sqrt2/(\varepsilon\eta^{q/2})\right)}
{-\log(1-\eta/\Lambda)}
\tag{6.11}
\]

（若对数分子非正，取 \(t=1\) 即可），便使 (6.10) 的第二项不超过 \(\varepsilon^2/2\)。于是这些节点的误差不超过 \(\varepsilon\)，得到 (N)；而 (N) 显然蕴含 (P)。有限传播结论来自 \(p_t\) 的次数。∎

用 \(-\log(1-\eta/\Lambda)\ge\eta/\Lambda\)，一个更易读但较粗的显式选择是

\[
t=
O\!\left(
1+\frac{\Lambda}{\eta}
\log_+\frac{1}{\varepsilon\eta^{q/2}}
\right).
\tag{6.11a}
\]

其中 \(\log_+x=\max\{0,\log x\}\)。每步只做一次有限传播矩阵-向量乘法，故轮数为 \(O(ts)\)。这里的 universality 是**定性的**：它说 Richardson 最终能达到任何可由共同多项式达到的多数节点精度；它不说 Richardson 的次数最优。若在 \([\eta,\Lambda]\) 上使用 Chebyshev 半迭代，通常可把有效条件数依赖从线性改善到平方根量级，但零附近尾项和具体谱空隙仍决定真正最优轮数。

**必要性证明，(P) \(\Rightarrow\) (S)。** 固定 \(a>0\) 和任意 \(\delta>0\)。由 (P)，对 \(\varepsilon=\sqrt a/3\) 取一个共同多项式 \(p\)。令

\[
M_p=\max_{\lambda\in[0,\Lambda]}|p(\lambda)|<\infty
\]

并选 \(\eta>0\) 使 \(\eta M_p\le1/2\)。于是对 \(0<\lambda<\eta\)，

\[
|1-\lambda p(\lambda)|
\ge1-\lambda M_p
\ge\frac12.
\tag{6.12}
\]

结合 (6.6a)，

\[
X_{\theta,i}^{(q)}(\eta)
\le4\mathcal R_{\theta,i}^{(q)}(p)^2.
\tag{6.13}
\]

所以在所有满足 \(\mathcal R_{\theta,i}^{(q)}(p)\le\sqrt a/3\) 的节点上都有 \(X_{\theta,i}^{(q)}(\eta)\le4a/9<a\)。从而

\[
\sup_\theta m_\theta\{X_{\theta,i}^{(q)}(\eta)>a\}\le\delta.
\]

又因 \(X(\eta)\) 随 \(\eta\downarrow0\) 单调下降，而 \(\delta>0\) 任意，得到 (6.7)。∎

若研究的是单个序列 \((J_n)\) 且只要求“大 \(n\)”成立，可把定理中每个 \(\sup_{\theta\in\Theta}\) 一致替换为 \(\limsup_{n\to\infty}\)；上述证明逐字不变。前者是对整个模型目录的 uniform 保证，后者是图极限语境下的 asymptotic 保证，二者不应混写。

### 6.3 固定允许失败比例的稳定版本

精确边界“至多恰好 \(\delta_0\)”可能遇到 infimum 不取到。自然且闭合的说法允许任意小失败率松弛。

**推论 6.2（固定 \(\delta_0\)）。** 给定 \(0\le\delta_0<1\)，以下等价：

1. 对每个 \(\varepsilon>0\) 和 \(\tau>0\)，存在共同多项式 \(p\)，使所有模型中 \(\mathcal R^{(q)}(p)>\varepsilon\) 的节点质量至多 \(\delta_0+\tau\)；等价地，可以把 \(p\) 限制为某个 \(p_t\)；
2. 对每个 \(a>0\)，
   \[
   \lim_{\eta\downarrow0}\sup_\theta
   m_\theta\{X_{\theta,i}^{(q)}(\eta)>a\}\le\delta_0.
   \tag{6.14}
   \]

证明沿定理 6.1，充分性先把 (6.14) 放宽到 \(\delta_0+\tau\)，必要性再令 \(\tau\downarrow0\)。因此 \(\delta_0\) 真正表示可以永久牺牲的局部病态节点比例，而非事后对一个固定误差数组取分位数。

全节点一致版本也立即得到：

\[
\forall\varepsilon>0\ \exists p:
\sup_{\theta,i}\mathcal R_{\theta,i}^{(q)}(p)\le\varepsilon
\quad\Longleftrightarrow\quad
\forall a>0\ \exists\eta>0:
\sup_{\theta,i}X_{\theta,i}^{(q)}(\eta)\le a.
\tag{6.15}
\]

这里左侧同样可以把任意多项式替换为 Richardson 序列中的 \(p_t\)。这是一个额外的 universality 结论：**只要任何共同多项式滤波器定性可行，最简单的固定步长 Richardson 也定性可行。**

### 6.4 匹配反例与“全局条件数不是答案”

以下先取 primitive-RHS 的 \(q=2\)。

**反例 A：一个低秩但完全离域的近零模态可让所有节点失败。** 取固定度 expander \(G_n\)，其归一化 Laplacian \(L_n\) 除常数模态外有统一谱隙，令

\[
J_n=L_n+\alpha_n I.
\]

常数特征向量在每个节点的谱质量都是 \(1/n\)。若 \(\alpha_n=n^{-1}\)，则对每个节点

\[
X_{n,i}(\eta)\ge\frac1{n\alpha_n^2}=n
\]

只要 \(\alpha_n<\eta\)。故 (6.7) 失败，而且任何固定多项式在充分大 \(n\) 上都有行误差至少约 \(1/(2\alpha_n\sqrt n)=\sqrt n/2\)。近零特征空间只有一维、占总维数 \(1/n\)，仍会污染**所有节点**。

**例 B：条件数发散却仍可局部逼近。** 同一例中若 \(\alpha_n=n^{-1/4}\)，常数模态对 (6.6) 的贡献是 \(n^{-1/2}\to0\)，其余谱被统一谱隙隔开。因此 (6.7) 成立，虽然 \(\kappa(J_n)\to\infty\)。这直接说明统一条件数只是强充分条件，不是多数节点局部估计的必要条件。

**例 C：少量高度局域病态模态可以被合法牺牲。** 令

\[
J_n=\operatorname{diag}(\underbrace{\alpha_n,\ldots,\alpha_n}_{m_n},1,\ldots,1),
\qquad \alpha_n\downarrow0,\quad m_n/n\to0,
\]

即使 \(m_n\alpha_n^{-2}/n\to\infty\)，(6.7) 仍成立，因为巨大尾质量只落在 \(o(n)\) 个坐标节点上。它展示多数节点标准严格弱于平均 Frobenius/\(L^2\)-trace 误差标准。

### 6.5 与有限 von Neumann 代数“依测度收敛”的区别

Nelson 的可测算子/measure topology 允许删掉小 trace 的谱投影，是处理无界逆的经典框架。[^7] 但它控制的是小谱**子空间的维数比例**，本定理控制的是该子空间经 \(1/\lambda\) 放大后在**节点基上的 leverage/行能量分布**。两者不能互换：

- 例 C 中 trace-\(L^2\) 尾可能发散，而多数节点行误差仍可消失；
- 反例 A 中近零谱投影秩比 \(1/n\to0\)，measure topology 可把它当作小投影丢掉，但它均匀铺在每个节点上，节点行误差反而处处很大。

因此定理 6.1 不是标准 trace 谱分布或非交换 measure topology 的直接重述。它也比只看经验谱分布更细：需要根谱测度 \(\mu_{\theta,i}\) 的节点间分布。顶点谱测度和局部弱收敛本身是经典对象；Bordenave--Lelarge 等证明了局部弱收敛与经验谱/resolvent 的联系。[^8][^9] 本次检索未发现这些工作给出 (6.7)--(6.8) 的节点多数 iff。

### 6.6 能否直接称为 WLS 的充要条件

primitive RHS \(b\) 对应 \(q=2\)。对白化原始观测，定理其实有一个同样精确、而且更贴近 WLS 的 \(q=1\) 版本。

考虑无先验 WLS，并令

\[
\xi=R^{-1/2}z,\qquad
G=H^*R^{-1/2},\qquad
J=H^*R^{-1}H=GG^*.
\tag{6.16}
\]

中心估计是 \(J^{-1}G\xi\)，局部多项式估计是 \(p(J)G\xi\)。在 \(\|\xi\|_2\le1\) 下，节点行误差平方恰为

\[
\begin{aligned}
\|e_i^\top(J^{-1}-p(J))G\|_2^2
&=e_i^\top(J^{-1}-p(J))J(J^{-1}-p(J))e_i\\
&=\int |1-\lambda p(\lambda)|^2\lambda^{-1}\,
d\mu_i(\lambda)
=\mathcal R_i^{(1)}(p)^2.
\end{aligned}
\tag{6.17}
\]

所以定理 6.1 逐字给出：

> 白化 WLS 在多数节点上存在共同多项式局部估计，当且仅当
> \(\int_{(0,\eta)}\lambda^{-1}d\mu_i\) 在节点测度上一致依概率消失；而且 Richardson 序列已经普适。

这比 primitive-RHS 的 \(\lambda^{-2}\) 条件更弱，原因是测量映射 \(G\) 在近零信息方向上同时缩小了输入能量。

例如 §6.4 的 expander 加 grounding 模型在 \(q=1\) 下，常数模态贡献为 \(1/(n\alpha_n)\)：\(\alpha_n=n^{-1}\) 仍使所有节点保持常数量级误差，而只要 \(n\alpha_n\to\infty\)，该单一离域模态不再阻止多数节点逼近。

测量局部性仍需单独检查：若每行 \(H\) 只涉及固定半径内的状态，且 \(R\) 按传感器块对角或 \(R^{-1/2}\) 本身有限传播，则各节点可局部形成 \(G\xi\)。若 \(R\) 的精度/白化是稠密的，(6.17) 虽代数正确，但算法不再是局部通信。

有先验 \(J=H^*R^{-1}H+J_0\) 时，若 \(J_0=B_0B_0^*\) 具有局部因子，并把 \(B_0\) 对应的虚拟白化观测并入输入，则增广 \(G=[H^*R^{-1/2},B_0]\) 仍满足 \(GG^*=J\)，故 \(q=1\) 版本保持精确。若只允许物理观测且先验项固定，则 \(GG^*=J-J_0\)，不能直接声称等号。

对一般原始观测 \(z\) 而不白化，令 \(G=H^*R^{-1}\)，多项式方案的误差为

\[
e_i^\top q(J)GG^*q(J)e_i,
\qquad q(\lambda)=\lambda^{-1}-p(\lambda).
\tag{6.18}
\]

不再只由 \(\mu_i\) 决定。若 \(cI\preceq GG^*\preceq CI\)，它与 (6.5) 双向可比；若只有统一上界，只得到充分性；若 \(GG^*\) 与 \(J\) 对易并等于谱函数 \(w(J)\)，可把 (6.6) 换成带权尾 \(\int_{(0,\eta)}\lambda^{-2}w(\lambda)d\mu_i\)。一般稠密 \(G\) 还会破坏运行时局部性。

同理，Gaussian Bayes MSE 要把 \(GG^*\) 换成输入协方差；除非对易或有 Loewner 谱窗，不能把 raw-WLS/Bayes 版本冒充为定理 6.1 的直接等价。

节点为固定维数 block 时可用矩阵值谱测度

\[
M_i(B)=E_i{\bf1}_B(J)E_i^*
\]

及尾矩阵 \(\left\|\int_{(0,\eta)}\lambda^{-q}dM_i(\lambda)\right\|\) 重复证明；标量不等式通过函数演算变成 Loewner 不等式。但 block 维数必须统一有界，且这仍是 \(p(J)\) 类算法。

### 6.7 原创性判定

定理的三块原料都已有：谱定理；Richardson/Landweber 的谱滤波表示；多项式图滤波的有限传播/分布式实现。Landweber 迭代作为 spectral regularization 及其滤波残差早已是标准逆问题理论。[^17] Shuman 等已经系统使用 Chebyshev 多项式实现分布式 inverse filtering；Emirov--Cheng--Jiang--Sun 的 Theorem 3.1、4.2、4.4 又对有限传播图滤波给出谱半径条件及指数收敛，后续 spectrum-adapted 与 inverse graph filter 工作继续发展该路线。[^10][^11][^12][^19] 现有工程文献通常在已知谱区间上做 uniform/average approximation，或默认逆乘子远离奇点；本次定向搜索没有找到以 (6.6a) 的节点间尾分布为必要充分条件、允许近零谱仅污染少数节点的结果。

所以准确说法是：

- (6.5) 是**标准恒等式**；
- Richardson 的误差核 (6.9) 和固定多项式的必要性都是**可直接证明的经典谱滤波步骤**；
- 把二者组织成“统一模型族、节点测度、多数节点、无统一谱隙”的 iff，并证明“任意多项式可行 iff 固定 Richardson 序列可行”，是**可能的新组合**；
- 是否已有 operator-algebra、graph-filter 或 inverse-problem 文献以等价形式覆盖，仍需专项逐式查重；目前不能写“首次”。

### 6.8 “可检验”应怎样理解

对一个给定有限矩阵，\(X_i^{(q)}(\eta)\) 可由特征分解或谱投影直接计算，(6.10) 又给出显式轮数证书。并且

\[
\int\lambda^k\,d\mu_i(\lambda)=(J^k)_{ii},
\]

每个固定阶矩只依赖有限邻域，所以完整根谱测度由不断增长的局部 moments 决定。

但 \(\lambda^{-q}{\bf1}_{(0,\eta)}\) 在零点无界。因此在没有额外尾 modulus 时，**不存在由某个预先固定半径的有限个 moments 自动认证 (6.7) 的结论**。换言之，定理 6.1 是可由全矩阵谱数据检验、可由局部迭代构造的结构 iff；它还不是“单个节点看一次固定半径球就能判定整个模型族是否满足条件”的局部 property tester。若投稿把它称为 locally checkable，必须再补一个有限样本谱尾估计定理。

### 6.9 为什么它不是所有 local-model-aware 算法的 iff

取

\[
J_n=\alpha_n I,\qquad \alpha_n\downarrow0.
\tag{6.19}
\]

每个节点从自己的 0-hop 系数就知道 \(\alpha_n\)，共同的局部规则“输出 \(b_i/(J_n)_{ii}\)”零轮精确恢复。因此最一般的 local-model-aware 解码器完全没有困难。

但 \(q=2\) 时每个节点的 \(X_{n,i}^{(2)}(\eta)=\alpha_n^{-2}\)（只要 \(\alpha_n<\eta\)），(6.7) 彻底失败；任何**同一个标量多项式** \(p\) 也确实不能同时逼近所有 \(1/\alpha_n\)。对白化 WLS 取 \(G_n=\sqrt{\alpha_n}I\) 时，同样有 0-hop 精确规则 \(\xi_i/\sqrt{\alpha_n}\)，而 \(q=1\) 谱尾仍失败。

所以定理 6.1 的逻辑边界不是技术细节：

> 它刻画 coefficient-oblivious 的共同谱滤波序列；它不刻画可以读取本地系数并做除法、局部消元或 view-dependent 预条件的所有局部算法。

一个自然修补是先做节点块 Jacobi 归一化。若

\[
D_\theta=\operatorname{blockdiag}_i (J_\theta)_{ii},
\qquad
K_\theta=D_\theta^{-1/2}J_\theta D_\theta^{-1/2},
\tag{6.20}
\]

各块可在 0-hop 求逆，便可对 \(K_\theta\) 应用定理 6.1，并以

\[
\widehat x
=D_\theta^{-1/2}p_t(K_\theta)D_\theta^{-1/2}b
\tag{6.21}
\]

恢复原坐标。精确误差等价需采用相应的加权输入/输出范数；若要回到普通 Euclidean 范数，还要对 \(D_\theta\) 给双向界。这个预条件版本能消除 (6.19) 的纯本地尺度病态，但仍不能涵盖任意局部 Schur 消元。图逆滤波中基于本地对角信息的预条件梯度法已有直接先例。[^20] 因此它是标准变量变换式扩展，不是“所有局部算法”的缺口闭合。

## 7. 其他跨领域语言能做什么、不能做什么

### 7.1 Sheaf 与 agreement testing

Sheaf 可把局部估计写成局部 section、把重叠一致性写成 restriction compatibility；Robinson 的 consistency radius 和 Hansen--Ghrist 的 harmonic extension 是直接相关语言。[^13] 但“把共同解码器称作 global section”本身只是重命名，没有自动处理 completion 的 minimax 量词或每模型 \(\delta\)-预算。

Agreement testing 有真正的结构性充分定理：在具有 agreement expansion 的重叠系统上，若大多数局部函数在随机重叠上相符，就存在一个全局函数与大多数局部函数一致。Dinur--Kaufman 及 Dinur--Filmus--Harsha 给出了此类结果。[^14] 若能证明 WLS completion 的局部最优系数在某个 incidence complex 上满足 agreement-expander 条件，可能产生算法；目前这个桥梁没有建立，而且已知结论是特殊有限字母/重叠系统下的充分恢复，不是一般 completion-majority iff。

### 7.2 measured coarse geometry 与 graph limits

Local-global graph limits 精确记录均匀随机根的有限 view 分布及所有局部着色输出分布；measured coarse geometry/strong almost finiteness 能表达删去小测度节点后把图分成受控块。[^15] 它们适合定义“允许 \(\delta\) 比例失败”，但只控制几何分解，不控制 WLS 的 Schur 边界影响或近零谱 leverage。定理 6.1 的 \(X_{\theta,i}\) 正说明：同样小秩的小谱空间可以局域在少数节点，也可以离域污染所有节点；纯 measured topology 不能区分二者。

### 7.3 LOCAL 不可区分性

Angluin 的 covering/不可区分思想及 Naor--Stockmeyer 的 LOCAL 模型已经覆盖“同一半径 view 的确定性算法输出相同”。[^1][^16] completion-pasting 下界可以借用这个框架，但不会自动给出连续 Gaussian/WLS 的 operator-norm 阈值或多数节点构造。

## 8. 硬门槛下的最终判决

### 8.1 抽象 completion-majority 路线：不够写论文

最强闭合结论是：

> 有限目录下，共同 \(r\)-local 规则存在，当且仅当可以在每个模型的失败预算内击中所有最小不可行约束；有限维输出/线性系数时，每条最小冲突由 Helly 定理压到维数加一。

它确有 iff、离线构造和 hardness 下界，但仍未越过论文主定理门槛，原因是：

1. 横截 + Helly 是 tolerant Helly/MaxFS 的直接应用；
2. 离散特例就是 Closest String；
3. 它需要整个 completion 目录，不能由单个节点的局部 view 检验；
4. 没有利用 Gaussian/WLS、Schur/DtN 或 M-matrix 的特有结构；
5. “逐节点半径取分位数”已被同边际反例彻底否定，只有强拼接假设才能恢复。

因此，这一部分适合放在项目的 **no-go / problem formulation / lower-bound** 章节，不应单独投稿。

### 8.2 唯一可能达到主定理形态的核

定理 6.1 同时具备：

- 非拼接、非 Chebyshev-radius 的谱结构 iff；
- 明确的共同有限轮构造 \(p(J)b\)；
- 更强的 Richardson/Neumann universality 及显式轮数界；
- 固定多项式在零附近的匹配必要性；
- 离域近零模态的失败反例与局域少数病态模态的正例；
- 严格说明为何全局条件数过强，以及为何经验谱/小秩删除仍不够。

但它尚有三个硬限制：

1. **算法类限制：** 只刻画共同多项式函数 \(p(J)\)（或预先指定的局部预条件后多项式），不是所有 view-dependent 有限传播矩阵或非线性局部算法；(6.19) 是严格分离反例；
2. **WLS 输入限制：** 对 primitive RHS 为 \(q=2\)，对白化且 \(GG^*=J\) 的 WLS 为 \(q=1\)；一般非白化输入仍需 \(GG^*\) 的对易/双向谱窗；
3. **原创性限制：** 证明短且原料经典，尚未完成以“diagonal spectral tail / local spectral measure / convergence in measure / inverse graph filtering”为关键词的系统逐式查重。

按用户“没有稍强、非拼接的主定理就不写论文”的硬门槛，本报告的诚实判决是：

- **一般 completion-majority 方向：NO-GO；**
- **原始“所有 local-model-only 算法”的目标：仍是 NO-GO；** 对角反例表明谱尾 iff 不是该最宽问题的答案；
- **若把论文明确收窄为 SPD 多项式/预条件多项式滤波：达到候选主定理门槛，但原创性 PENDING。** 在专项查重和一般相关噪声扩展完成前，不应宣称已足以投稿。

若后续要把候选提升为可靠论文核心，最小工作包是：

1. 对定理 6.1 做专项先例检索，尤其有限 von Neumann 代数中带对角条件期望的 measure convergence、graphings 上 affiliated inverse、以及 nodewise graph-filter approximation；
2. 把已闭合的 whitened-WLS \(q=1\) 版扩展到一般非白化/相关噪声，或在局部 \(G\) 且 \(cJ^\beta\preceq GG^*\preceq CJ^\beta\) 下给出统一 iff；
3. 给出可计算的有限样本证书/估计器，用局部 moments 或随机 probing 估计 \(X_i(\eta)\) 的节点分布；
4. 证明多项式类相对更广的 \(r\)-local linear rules 的 gap，或找出使多项式无损的对称性/同质性假设；
5. 做局域 tentacle 与离域低模态的数值相变实验，分别验证 RHS 的 inverse-square leverage 与 whitened-WLS 的 inverse-first-power leverage，而非条件数，预测多数节点误差。

## Sources

[^1]: Moni Naor, Larry Stockmeyer, “What Can Be Computed Locally?”, *SIAM Journal on Computing* 24(6) (1995), 1259–1277, especially §2 and Lemma 3.2. DOI [10.1137/S0097539793254571](https://doi.org/10.1137/S0097539793254571).

[^2]: Luis Montejano, Déborah Oliveros, “Tolerance in Helly-Type Theorems”, *Discrete & Computational Geometry* 45 (2011), 348–357; Definition 1.1, Theorems 1.1, 2.1 and 3.1. DOI [10.1007/s00454-010-9296-6](https://doi.org/10.1007/s00454-010-9296-6).

[^3]: Meir Katchalski, A. Liu, “A Problem of Geometry in \(\mathbb R^n\)”, *Proceedings of the AMS* 75(2) (1979), 284–288. DOI [10.1090/S0002-9939-1979-0532151-4](https://doi.org/10.1090/S0002-9939-1979-0532151-4). Gil Kalai, “Intersection Patterns of Convex Sets”, *Israel Journal of Mathematics* 48 (1984), 161–174. DOI [10.1007/BF02761162](https://doi.org/10.1007/BF02761162).

[^4]: Laurent El Ghaoui, “Inversion Error, Condition Number, and Approximate Inverses of Uncertain Matrices”, *Linear Algebra and its Applications* 343–344 (2002), 171–193; Eq. (1.4), Theorem 6.2/Eq. (6.3), Eq. (6.5). DOI [10.1016/S0024-3795(01)00273-7](https://doi.org/10.1016/S0024-3795(01)00273-7); [author PDF](https://people.eecs.berkeley.edu/~elghaoui/Pubs/InvErr_LAA02.pdf).

[^5]: Ming Li, Bin Ma, Lusheng Wang, “On the Closest String and Substring Problems”, *Journal of the ACM* 49(2) (2002), 157–171. DOI [10.1145/506147.506150](https://doi.org/10.1145/506147.506150). Jens Gramm, Rolf Niedermeier, Peter Rossmanith, “Fixed-Parameter Algorithms for Closest String and Related Problems”, *Algorithmica* 37 (2003), 25–42. DOI [10.1007/s00453-003-1028-3](https://doi.org/10.1007/s00453-003-1028-3).

[^6]: Vladimir Shenmaier, “Complexity and Approximation of the Smallest \(k\)-Enclosing Ball Problem”, *European Journal of Combinatorics* 48 (2015), 81–87. DOI [10.1016/j.ejc.2015.02.011](https://doi.org/10.1016/j.ejc.2015.02.011).

[^7]: Edward Nelson, “Notes on Non-Commutative Integration”, *Journal of Functional Analysis* 15 (1974), 103–116. DOI [10.1016/0022-1236(74)90014-7](https://doi.org/10.1016/0022-1236(74)90014-7).

[^8]: C. D. Godsil, B. Mohar, “Walk Generating Functions and Spectral Measures of Infinite Graphs”, *Linear Algebra and its Applications* 107 (1988), 191–206. DOI [10.1016/0024-3795(88)90245-5](https://doi.org/10.1016/0024-3795(88)90245-5).

[^9]: Charles Bordenave, Marc Lelarge, “Resolvent of Large Random Graphs”, *Random Structures & Algorithms* 37 (2010), 332–352. DOI [10.1002/rsa.20313](https://doi.org/10.1002/rsa.20313); [arXiv:0801.0155](https://arxiv.org/abs/0801.0155).

[^10]: David I. Shuman, Pierre Vandergheynst, Daniel Kressner, Pascal Frossard, “Distributed Signal Processing via Chebyshev Polynomial Approximation”, *IEEE Transactions on Signal and Information Processing over Networks* 4(4) (2018), 736–751. DOI [10.1109/TSIPN.2018.2824239](https://doi.org/10.1109/TSIPN.2018.2824239); [arXiv:1111.5239](https://arxiv.org/abs/1111.5239).

[^11]: Tiffany Fan, David I. Shuman, Shashanka Ubaru, Yousef Saad, “Spectrum-Adapted Polynomial Approximation for Matrix Functions with Applications in Graph Signal Processing”, *Algorithms* 13(11) (2020), 295. DOI [10.3390/a13110295](https://doi.org/10.3390/a13110295).

[^12]: Cheng Cheng, Qiyu Sun, Cong Zheng, “Iterative Polynomial Approximation Algorithms for Inverse Graph Filters”, SAMPTA 2025. DOI [10.1109/SAMPTA64769.2025.11133514](https://doi.org/10.1109/SAMPTA64769.2025.11133514); [arXiv:2504.14341](https://arxiv.org/abs/2504.14341). 另见 Junzheng Jiang, David B. Tay, “Decentralised Signal Processing on Graphs via Matrix Inverse Approximation”, *Signal Processing* 167 (2020), 107273. DOI [10.1016/j.sigpro.2019.07.010](https://doi.org/10.1016/j.sigpro.2019.07.010).

[^13]: Michael Robinson, “Sheaves Are the Canonical Data Structure for Sensor Integration”, *Information Fusion* 36 (2017), 208–224; Definition 20, Proposition 21, Theorem 29. DOI [10.1016/j.inffus.2016.12.002](https://doi.org/10.1016/j.inffus.2016.12.002). 本地：literature/07_sheaves_topology/2016_robinson_sheaves_sensor_integration.pdf。 Jakob Hansen, Robert Ghrist, “Toward a Spectral Theory of Cellular Sheaves”, *Journal of Applied and Computational Topology* 3 (2019), 315–358. DOI [10.1007/s41468-019-00038-7](https://doi.org/10.1007/s41468-019-00038-7). 本地：literature/07_sheaves_topology/2018_hansen_ghrist_spectral_theory_cellular_sheaves.pdf。

[^14]: Irit Dinur, Tali Kaufman, “High Dimensional Expanders Imply Agreement Expanders”, FOCS 2017. DOI [10.1109/FOCS.2017.94](https://doi.org/10.1109/FOCS.2017.94). Irit Dinur, Yuval Filmus, Prahladh Harsha, “Agreement Tests on Graphs and Hypergraphs”, *SIAM Journal on Computing* 54(2) (2025). DOI [10.1137/21M1397684](https://doi.org/10.1137/21M1397684), especially Theorems 1.1–1.2.

[^15]: Hamed Hatami, László Lovász, Balázs Szegedy, “Limits of Locally-Globally Convergent Graph Sequences”, *Geometric and Functional Analysis* 24 (2014), 269–296; Theorem 3.2. DOI [10.1007/s00039-014-0258-7](https://doi.org/10.1007/s00039-014-0258-7). 本地：literature/02_graphical_models_local_inference/2012_hatami_lovasz_szegedy_local_global_graph_limits.pdf。Gábor Elek, Ádám Timár, “Strong Almost Finiteness”, *Journal of Functional Analysis* 289 (2025), 111116. DOI [10.1016/j.jfa.2025.111116](https://doi.org/10.1016/j.jfa.2025.111116).

[^16]: Dana Angluin, “Local and Global Properties in Networks of Processors”, STOC 1980, 82–93. DOI [10.1145/800141.804655](https://doi.org/10.1145/800141.804655).

[^17]: Heinz W. Engl, Martin Hanke, Andreas Neubauer, *Regularization of Inverse Problems*, Kluwer, 1996, especially the spectral-filter treatment of Landweber iteration. DOI [10.1007/978-94-009-1740-8](https://doi.org/10.1007/978-94-009-1740-8).

[^18]: K. Yu. Osipenko, “Generalized Adaptive versus Nonadaptive Recovery from Noisy Information”, *Journal of Complexity* 53 (2019), 162–172; Lemma 2 and Eq. (7). DOI [10.1016/j.jco.2018.12.001](https://doi.org/10.1016/j.jco.2018.12.001).

[^19]: Nazar Emirov, Cheng Cheng, Junzheng Jiang, Qiyu Sun, “Polynomial Graph Filters of Multiple Shifts and Distributed Implementation of Inverse Filtering”, *Sampling Theory, Signal Processing, and Data Analysis* 20 (2022), article 2; Theorem 3.1, Remark 3.2, Theorems 4.2 and 4.4. DOI [10.1007/s43670-021-00019-x](https://doi.org/10.1007/s43670-021-00019-x). 本地：literature/03_sparse_inverse_graph_filters/2020_emirov_cheng_jiang_sun_polynomial_inverse_graph_filter.pdf。

[^20]: Cheng Cheng, Nazar Emirov, Qiyu Sun, “Preconditioned Gradient Descent Algorithm for Inverse Filtering on Spatially Distributed Networks”, *IEEE Signal Processing Letters* 27 (2020), 1834–1838. DOI [10.1109/LSP.2020.3029699](https://doi.org/10.1109/LSP.2020.3029699).
