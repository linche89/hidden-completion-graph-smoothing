# 从共同多项式到所有局部线性算子：可成立的主定理、边界与撞车审计

**日期：** 2026-09-13  
**研究对象：** 有限传播 SPD 正规矩阵的局部求逆，以及 Gaussian/WLS 中由白化量测诱导的逐节点误差。  
**结论口径：** 下文把标准恒等式、直接推论、可能的新组合和未闭合问题逐项标明；“未检索到同句定理”不等于原创。

## 0. 结论先行

这条线最后形成了一个肯定结果和三个否定边界。

1. **全局知道完整矩阵后离线定制局部系数（A 类）没有新结构可挖。** 对任意目标线性算子 \(T\)，节点 \(i\) 的最佳 \(r\)-hop 行算子就是把 \(T_i\) 在 \(r\)-球外的列删掉。最佳误差恰为该行的球外 \(\ell_2\) 尾；允许 \(\delta\) 比例节点失败时，充要条件和最小轮数就是这些行尾的 \((1-\delta)\)-分位数。这是坐标 Hilbert 投影，不是新的图论定理，也通常不是算子代数中的条件期望。

2. **对白化 WLS 的一个自然且很宽的算法类，确有比“共同多项式”更强的 iff。** 算法先在本地形成充分统计量 \(b=G\xi\)，再用任意有限传播线性算子 \(L_\theta b\) 求解。只要求每个被选中的规则在整个模型族上有有限的统一增益

   \[
   M=\sup_{\theta}\|L_\theta\|<\infty,
   \]

   并不要求 \(L_\theta\) 是多项式，也不要求它与 \(J_\theta\) 对易。此时下列三件事等价：

   - 对任意 \(\varepsilon,\delta>0\)，某个稳定的有限轮局部线性规则能在每个模型的至少 \(1-\delta\) 节点达到白化输入最坏误差 \(\varepsilon\)；
   - 小特征值的节点逆一次矩

     \[
     X_{\theta,i}^{(1)}(\eta)
     =\int_{(0,\eta)}\lambda^{-1}\,d\mu_{\theta,i}(\lambda)
     \]

     在整个模型族上对随机根一致依概率趋于零；
   - 一个完全共同、无需训练的 Richardson 多项式族已经能做到同样的多数节点逼近。

   核心必要性只有一行但不是定义式：

   \[
   \sqrt{X_{\theta,i}^{(1)}(\eta)}
   \le
   E_{\theta,i}^{(1)}(L_\theta)+M\sqrt{\eta}.
   \]

   因而“任何稳定局部线性求解器可行”会自动推出“最简单的共同多项式可行”。这是本报告中唯一达到**结构条件 + 构造算法 + 匹配必要性**形态的候选主定理。

3. **对任意 RHS 的 \(q=2\) 行误差，上述无对易结论不能逐节点照搬。** 若 \(L_\theta\) 与 \(J_\theta\) 对易，则同样的 iff 成立，尾改为 \(\lambda^{-2}\)。距离传递网络和阿贝尔平移不变网络是自然的充分结构。若不对易，本文给出一个固定 \(\delta\) 的局部、稳定、精确反例：好节点误差为零，但其 \(\lambda^{-2}\) 小谱尾保持为常数。因此对易性不是装饰。它是否能换成明显更弱且仍可检验的条件，目前没有闭合。

4. **“所有原始量测上的局部线性算子”仍未被完全覆盖。** 若算法直接输出 \(K_\theta\xi\)，未必能在不破坏传播半径的前提下写成 \(L_\theta G_\theta\xi\)。所以第 2 点严格覆盖的是“先本地汇总 \(b=G\xi\)，再局部求解正规方程”的完整线性类，而不是所有可能利用量测零空间、边变量和方向变量的线性电路。

5. **格点临界维数不是原创核。** 在 \(\mathbb Z^d\) 上，primitive RHS 的临界维数是 \(d=4\)，白化 WLS/GFF 的临界维数是 \(d=2\)；最佳半径 \(r\) 误差分别为

   \[
   \Theta\!\left(r^{-(d-4)/2}\right),
   \qquad
   \Theta\!\left(r^{-(d-2)/2}\right).
   \]

   最优次数 \(r\) 多项式通过 Christoffel 函数也达到同阶。但 Backhausz–Virág 已从 linear factor-of-i.i.d. 的谱测度角度覆盖顶点传递图上的定性存在性，Green 函数渐近和 Christoffel 端点渐近也都是经典结果。因此这一节是漂亮例子和 sharpness 证书，不宜包装为主原创。

**总判定：** 原始口号“所有局部线性算法的最宽 iff”仍是 **NO-GO**；但在“白化 WLS + 本地充分统计量 + 稳定有限传播求解器”这一工程上自然的完整类中，下面的定理 4.1 是一个非平凡、可审稿的主定理候选。必须立即降调的是：若模型族退化为一个无限顶点传递 Gaussian 模型，则 linear-factor 的定性存在性已经落入 Backhausz–Virág Theorem 2，GFF 的 \(q=1\) 情形更由其 Proposition 26 直接覆盖。尚待专项查重的仅是不齐次有限模型族、逐节点多数语义以及“任意稳定、可不对易的正规方程规则可行 \(\Rightarrow\) 一个共同 Richardson 规则可行”这一反向压缩结论。

---

## 1. 统一模型与三种算法权限

### 1.1 网络、数据位置与有限传播

令 \(\Theta\) 是模型族。每个 \(\theta\in\Theta\) 给出：

- 有限节点集 \(V_\theta\)，节点概率测度 \(m_\theta\)；
- 图距离 \(\rho_\theta\)；
- 每个节点一个固定维数状态块；为简洁可先读成标量；
- 自伴正定矩阵 \(J_\theta\)，满足共同谱上界

  \[
  0<J_\theta\le \Lambda I;
  \tag{1.1}
  \]

- \(J_\theta\) 的传播宽度至多 \(s\)，即

  \[
  (J_\theta)_{ij}=0
  \quad\text{若}\quad
  \rho_\theta(i,j)>s.
  \tag{1.2}
  \]

节点块行抽取记为 \(E_{\theta,i}\)，其伴随 \(E_{\theta,i}^{*}\) 把节点块嵌回全局空间。标量情形就是 \(E_i=e_i^{*}\)。

若输入坐标不是节点而是边/量测，则给每个输入坐标一个位置，并令 \(S_{\theta,i,r}\) 表示节点 \(i\) 在 \(r\) 轮内可见的输入坐标集合。相应正交坐标投影记为 \(P_{\theta,i,r}\)，球外投影记为 \(Q_{\theta,i,r}=I-P_{\theta,i,r}\)。

一个线性算子 \(L_\theta\) 有传播宽度 \(r\)，指其第 \(i\) 块行只用 \(S_{\theta,i,r}\) 中的输入。无界消息的同步 LOCAL 模型里，节点收集完半径 \(r\) 的运行时数据后，恰能计算任何这样的已知行线性形式。

### 1.2 A、B、C 三类不能混写

**A：完整模型可见、离线定制。** 对每个完整 \(\theta\)，设计者可以看见全部 \(J_\theta\)，并为每个节点预装任意 \(r\)-局部系数。运行时通信仍只有 \(r\) 轮，但系数可能包含全局 advice。

**B：系数只能由 rooted \(r\)-view 产生。** 存在一个共同映射

\[
\Phi_r:\{\text{带参数的 rooted }r\text{-views}\}
\longrightarrow
\{\text{局部系数行}\},
\tag{1.3}
\]

拥有相同规范化 view 的发生项必须得到相同系数。若唯一 ID 被放进 view，B 可能退化成查表 advice；因此必须说明结论是否对所有 ID 指派成立。

**C：同一个内禀等变规则。** 除 B 外，规则还尊重所有保留模型标记的重标号：

\[
L_{\pi\theta}
=\Pi_{\mathrm{out}}L_\theta\Pi_{\mathrm{in}}^{*}.
\tag{1.4}
\]

若 B 的 view 一开始就取 rooted 标记同构类且禁止 ID advice，则 B 与 C 很接近；若 view 保留名字，两者不同。后文始终明确使用哪一个版本。

### 1.3 两种误差，不要把 \(q=1\) 与 \(q=2\) 混淆

对任意 \(q>0\)，定义正规方程求解器 \(L_\theta\) 在节点 \(i\) 的加权误差

\[
E_{\theta,i}^{(q)}(L_\theta)
:=
\left\|
E_{\theta,i}(I-L_\theta J_\theta)J_\theta^{-q/2}
\right\|_{2\to2}.
\tag{1.5}
\]

其伴随形式是

\[
E_{\theta,i}^{(q)}(L_\theta)
=
\left\|
J_\theta^{-q/2}(I-J_\theta L_\theta^{*})
E_{\theta,i}^{*}
\right\|_{2\to2}.
\tag{1.6}
\]

- \(q=2\) 时，

  \[
  E_{\theta,i}^{(2)}(L_\theta)
  =
  \|E_{\theta,i}(J_\theta^{-1}-L_\theta)\|_{2\to2},
  \tag{1.7}
  \]

  正是任意 \(\|b\|_2\le1\) RHS 的逐节点 robust 误差。

- \(q=1\) 时，它是白化 WLS 的逐节点误差。若

  \[
  G_\theta=H_\theta^{*}R_\theta^{-1/2},
  \qquad
  J_\theta=G_\theta G_\theta^{*},
  \tag{1.8}
  \]

  中心估计为 \(x^{\star}=J_\theta^{-1}G_\theta\xi\)，局部正规方程算法为
  \(\widehat x=L_\theta G_\theta\xi\)，则

  \[
  \begin{aligned}
  \|E_i(J^{-1}-L)G\|_{2\to2}^{2}
  &=
  \left\|
  E_i(J^{-1}-L)J(J^{-1}-L)^{*}E_i^{*}
  \right\|\\
  &=
  \bigl(E_i^{(1)}(L)\bigr)^2.
  \end{aligned}
  \tag{1.9}
  \]

先验项可作为虚拟白化量测并入 \(G\)。但若 \(R^{-1/2}\) 稠密，白化本身不是局部操作，(1.8) 的网络解释便失效；应要求 \(G\) 有共同有限传播宽度，或直接从局部独立量测建模。

---

## 2. A 类的精确答案：最佳行就是坐标截断

### 命题 2.1（任意目标算子的 oracle 行最优）

令 \(T_\theta\) 是从全局输入 Hilbert 空间到节点状态的中心线性算子。固定 \((\theta,i,r)\)，所有 \(r\)-局部块行 \(K_i\) 满足 \(K_i=K_iP_{i,r}\)。则

\[
\inf_{K_i=K_iP_{i,r}}
\|E_iT_\theta-K_i\|_{2\to2}
=
\|E_iT_\theta Q_{i,r}\|_{2\to2},
\tag{2.1}
\]

且截断

\[
K_i^{\star}=E_iT_\theta P_{i,r}
\tag{2.2}
\]

达到下确界。标量输出或 Frobenius 范数下它还是唯一最优解；块输出的谱范数下可能有别的并列最优解。

**证明。** 按可见/不可见输入正交分解，

\[
E_iT_\theta-K_i
=
\bigl[
E_iT_\theta P_{i,r}-K_i,\ 
E_iT_\theta Q_{i,r}
\bigr].
\]

对块行 \(A,B\)，

\[
\|[A,B]\|_{2\to2}^{2}
=\lambda_{\max}(AA^{*}+BB^{*})
\ge \lambda_{\max}(BB^{*})=\|B\|^2.
\]

取第一块为零即达到该下界。证毕。

**状态：标准 Hilbert 投影/直接恒等式。**

对 primitive RHS，取 \(T_\theta=J_\theta^{-1}\)。对白化原始量测，取
\(T_\theta=J_\theta^{-1}G_\theta\)，并以量测坐标的位置定义球。

### 推论 2.2（A 类多数节点 iff 与最小轮数）

定义精确行尾

\[
\tau_{\theta,i}(r)
:=
\|E_iT_\theta Q_{i,r}\|_{2\to2}.
\tag{2.3}
\]

允许每个完整模型分别离线定制系数时，存在一个 \(r\)-传播线性估计器，使每个 \(\theta\) 至少 \(1-\delta\) 的节点误差不超过 \(\varepsilon\)，当且仅当

\[
\sup_{\theta\in\Theta}
m_\theta\{i:\tau_{\theta,i}(r)>\varepsilon\}
\le\delta.
\tag{2.4}
\]

因此最小轮数精确为

\[
r_A^{\star}(\varepsilon,\delta;\Theta)
=
\inf\left\{
r:
\sup_\theta m_\theta\{i:\tau_{\theta,i}(r)>\varepsilon\}
\le\delta
\right\}.
\tag{2.5}
\]

**证明。** 必要性逐节点用命题 2.1；充分性同时取每一行的截断。证毕。

**状态：分位数恒等式；作为论文主定理是 NO-GO。**

### 2.3 为什么这不是一般的“条件期望定理”

全矩阵的半径掩膜

\[
\Pi_r(M)_{ij}
=
\mathbf 1_{\{\rho(i,j)\le r\}}M_{ij}
\tag{2.6}
\]

确是 Frobenius 内积下到给定稀疏线性子空间的正交投影。但半径 \(r>0\) 的有限传播矩阵通常不构成代数：两个 \(r\)-传播矩阵相乘可有 \(2r\) 传播。因此 \(\Pi_r\) 通常不是 \(C^{*}\)-或 von Neumann 代数意义的条件期望，也未必在谱范数下正收缩。只有对角代数或固定块分割等真正的子代数情形才可这样称呼。

稀疏 approximate inverse 文献通常求

\[
\min_{K\in\mathcal S}\|I-JK\|_F
\quad\text{或}\quad
\min_{K\in\mathcal S}\|I-KJ\|_F,
\tag{2.7}
\]

以得到可由小型 least-squares 构造的预条件器；Grote–Huckle 的 SPAI 是经典代表。[^1] 这与已知 \(J^{-1}\) 后最小化
\(\|J^{-1}-K\|_F\) 不同。后一个问题的答案虽然就是删列，但计算这些系数本身已经要求全局逆。

---

## 3. A 不等于 B/C：completion 与 advice 反例

### 3.1 Schur 补说明 A 的系数含全局信息

把可见区域 \(I=B_r(i)\) 与外部 \(O\) 分块：

\[
J=
\begin{bmatrix}
A&E\\ E^{*}&D
\end{bmatrix}.
\tag{3.1}
\]

中心逆的局部块是

\[
(J^{-1})_{II}
=
(A-ED^{-1}E^{*})^{-1}.
\tag{3.2}
\]

即使 \(A,E\) 完全相同，未知 completion \(D\) 也会改变局部最优系数。这是离散 Dirichlet-to-Neumann/Schur 边界影响的标准恒等式。

一个数值清楚的 2-by-2 族是

\[
J_d=
\begin{bmatrix}
2&M\\ M&d
\end{bmatrix},
\qquad
d\in\{M^2,2M^2\}.
\tag{3.3}
\]

根只看见 \(J_{00}=2\) 和跨边权 \(M\)，看不见外部对角 \(d\)。两 completion 的根部自系数分别是

\[
(J_{M^2}^{-1})_{00}=1,
\qquad
(J_{2M^2}^{-1})_{00}=2/3.
\tag{3.4}
\]

A 可以分别预装二者；B/C 的相同 view 必须共用一个系数，因而至少有一个 completion 的可见列误差不小于 \(1/6\)。这不依赖无环性。
这里采用“根知道 incident coupling、但不知道边外节点自参数”的 halo 口径；若 view 定义为不显示 cut edge 的纯 induced ball，只需在根与该 2-by-2 边界块之间加一段相同 buffer，并把半径整体平移一格。

### 3.2 B 的定义若允许名字，会退化成 A

若每个节点 ID 都进入 view，且模型目录有限，\(\Phi_r\) 可以把 ID 当作索引，查出预先编码的全局逆行。此时所谓“只由局部 view 生成”并没有排除 advice。要研究拓扑局部性，必须至少采用以下一种口径：

- 匿名 rooted 标记同构；
- 对所有 ID 指派都正确；
- order-invariant / ID-oblivious 规则；
- 对跨规模模型族使用同一个有限描述规则。

### 3.3 局部缩放反例：共同多项式绝不等于一般 B/C

令

\[
J_\alpha=\alpha I,\qquad \alpha\downarrow0.
\tag{3.5}
\]

每个节点从 0-hop view 看见 \(\alpha\)，故 B/C 规则 \(L_\alpha=\alpha^{-1}I\) 零轮精确。但不存在一个固定有限次数多项式在所有趋零 \(\alpha\) 上逼近 \(1/\alpha\)。这说明：

- 不加稳定性约束时，谱尾不是所有 local-model-aware 算法的必要条件；
- “多项式普适”必须来自额外结论，而不能靠定义把所有局部算法叫作图滤波器；
- 绝对误差对局部单位缩放敏感，归一化模型类必须说清。

---

## 4. 核心定理：白化 WLS 下所有稳定局部正规方程算子

### 4.1 节点小谱尾

令

\[
P_{\theta,\eta}
:=
\mathbf 1_{(0,\eta)}(J_\theta).
\tag{4.1}
\]

对任意 \(q>0\)，定义节点块的小谱逆矩

\[
X_{\theta,i}^{(q)}(\eta)
:=
\left\|
J_\theta^{-q/2}P_{\theta,\eta}E_{\theta,i}^{*}
\right\|_{2\to2}^{2}.
\tag{4.2}
\]

标量情形若 \(\mu_{\theta,i}\) 是 \(J_\theta\) 在 \(e_i\) 处的谱测度，则

\[
X_{\theta,i}^{(q)}(\eta)
=
\int_{(0,\eta)}
\lambda^{-q}\,d\mu_{\theta,i}(\lambda).
\tag{4.3}
\]

所需结构条件是

\[
\boxed{
\forall a>0:\quad
\lim_{\eta\downarrow0}
\sup_{\theta\in\Theta}
m_\theta\{i:X_{\theta,i}^{(1)}(\eta)>a\}=0.
}
\tag{ST1}
\]

它允许：

- \(\inf_\theta\lambda_{\min}(J_\theta)=0\)；
- 少量节点的尾任意大；
- 不同模型有不同坏区；
- 不要求全局条件数有界。

### 定理 4.1（稳定的任意局部线性求解器 iff 谱尾 iff Richardson）

假设 (1.1)–(1.2)，并假设每个节点块维数一致有界。考虑任何算法类 \(\mathcal C\)，满足：

1. \(\mathcal C\) 至少包含所有共同 Richardson 算子 \(p_t(J_\theta)\)；
2. \(\mathcal C\) 中每个被选择的跨模型规则 \(\theta\mapsto L_\theta\) 都有某个有限增益

   \[
   M(L):=\sup_{\theta}\|L_\theta\|<\infty.
   \tag{4.4}
   \]

这里 \(M(L)\) 可以随所选精度、失败率和轮数改变；不要求所有算法共用一个固定 \(M\)。

则下列三项等价。

**(F1) 任意稳定局部规则可行。** 对每个 \(\varepsilon,\delta>0\)，存在某个有限传播规则
\(\{L_\theta\}\in\mathcal C\)，使

\[
\sup_\theta
m_\theta\{i:E_{\theta,i}^{(1)}(L_\theta)>\varepsilon\}
\le\delta.
\tag{4.5}
\]

**(ST1) 节点小谱逆一次矩一致依概率消失。**

**(R1) Richardson 本身可行。** 对每个 \(\varepsilon,\delta>0\)，存在同一个整数 \(t\)，使对所有 \(\theta\)，共同多项式

\[
p_t(\lambda)
=
\frac{1-(1-\lambda/\Lambda)^t}{\lambda}
=
\frac1\Lambda\sum_{k=0}^{t-1}
(1-\lambda/\Lambda)^k
\tag{4.6}
\]

满足

\[
\sup_\theta
m_\theta\{i:E_{\theta,i}^{(1)}(p_t(J_\theta))>\varepsilon\}
\le\delta.
\tag{4.7}
\]

其传播宽度不超过 \((t-1)s\)。

#### 证明：F1 推出 ST1

固定 \(a>0\) 和任意 \(\rho>0\)。在 F1 中取

\[
\varepsilon=\sqrt a/2,\qquad \delta=\rho,
\]

得到规则 \(L_\theta\) 及有限 \(M=\sup_\theta\|L_\theta\|\)。

对每个节点块，由三角不等式和谱投影收缩性，

\[
\begin{aligned}
\sqrt{X_{\theta,i}^{(1)}(\eta)}
&=
\|J_\theta^{-1/2}P_{\theta,\eta}E_i^{*}\|\\
&\le
\|P_{\theta,\eta}J_\theta^{-1/2}
(I-J_\theta L_\theta^{*})E_i^{*}\|
+
\|P_{\theta,\eta}J_\theta^{1/2}
L_\theta^{*}E_i^{*}\|\\
&\le
E_{\theta,i}^{(1)}(L_\theta)
+
M\sqrt\eta .
\end{aligned}
\tag{4.8}
\]

这里没有交换 \(P_{\eta}\) 与 \(L_\theta\)；只用了
\(\|P_\eta J^{1/2}\|\le\sqrt\eta\)。

取 \(M\sqrt\eta\le\sqrt a/2\)。所有误差不超过 \(\sqrt a/2\) 的节点都有
\(X_{\theta,i}^{(1)}(\eta)\le a\)，故

\[
\sup_\theta m_\theta\{i:X_{\theta,i}^{(1)}(\eta)>a\}
\le\rho.
\]

\(\rho\) 任意即得 ST1。

#### 证明：ST1 推出 R1

Richardson 残差恒等式为

\[
I-p_t(J_\theta)J_\theta
=
(I-J_\theta/\Lambda)^t.
\tag{4.9}
\]

谱区间按 \((0,\eta)\cup[\eta,\Lambda]\) 分开，得到每个节点块的显式证书

\[
\boxed{
\bigl(E_{\theta,i}^{(1)}(p_t(J_\theta))\bigr)^2
\le
X_{\theta,i}^{(1)}(\eta)
+
\frac{(1-\eta/\Lambda)^{2t}}{\eta}.
}
\tag{4.10}
\]

给定 \(\varepsilon,\delta\)，先用 ST1 选 \(\eta\)，使除至多 \(\delta\) 节点外
\(X_i^{(1)}(\eta)\le\varepsilon^2/2\)；再选共同 \(t\)，使第二项不超过
\(\varepsilon^2/2\)。于是得 R1。又因
\(\sup_{0<\lambda\le\Lambda}|p_t(\lambda)|\le t/\Lambda\)，该共同规则满足
(4.4)，故 R1 推出 F1。证毕。

**状态：可能的新组合。** 证明成分全是谱投影、三角不等式和 Richardson；尚未找到把“任意不对易、view-dependent、稳定有限传播 \(L_\theta\)”直接压缩成同一 Richardson，并以**逐节点对角尾的 uniform in-probability**为 iff 的原始论文。它不是 radius-of-information 的定义重述，因为 ST1 是只依赖 \(J\) 的结构证书，且给出一个共同构造。

### 4.2 固定 \(\delta\) 版本的精确边界

若固定 \(\delta\) 且对每个 \(\varepsilon>0\) 都有 (4.5)，则对每个 \(a>0\)

\[
\lim_{\eta\downarrow0}
\sup_\theta
m_\theta\{i:X_{\theta,i}^{(1)}(\eta)>a\}
\le\delta.
\tag{4.11}
\]

反向若右侧仅知 \(\le\delta\)，一般只能保证任意 \(\zeta>0\) 的
\(\delta+\zeta\) 失败率；要得到恰好 \(\delta\)，还需分位数达到性或对
\(\eta\) 的统一右连续控制。这是多数节点量词中不可忽略的边界，不应把
\(\limsup\le\delta\) 写成自动达到 \(\delta\)。

### 4.3 all-node 版本

把节点概率条件换成最大值，就得到

\[
\forall a>0:\quad
\lim_{\eta\downarrow0}
\sup_{\theta,i}
X_{\theta,i}^{(1)}(\eta)=0
\tag{4.12}
\]

当且仅当任意精度的 all-node 稳定局部逼近存在，当且仅当 Richardson
all-node 可行。这里不再允许病态区躲进少数节点。

### 4.4 为什么稳定性不能删

在 (3.5) 的 \(J_\alpha=\alpha I\) 族中，

\[
X_{\alpha,i}^{(1)}(\eta)=\alpha^{-1}
\quad(\alpha<\eta),
\]

故 ST1 完全失败；但 \(L_\alpha=\alpha^{-1}I\) 零轮精确。唯一被违反的正是

\[
\sup_{\alpha}\|L_\alpha\|<\infty.
\]

这表明稳定性不是技术方便，而是排除“从本地读到任意小单位后无限放大”的最小工程语义。它也与能量/噪声放大要求吻合。

### 4.5 该定理到底覆盖多少“所有局部算法”

定理 4.1 的必要性甚至没有用有限传播、B 或 C；因此它覆盖 A/B/C 中所有满足 (4.4) 的正规方程算子。充分性则给出 C 类共同规则。

但对直接量测算法

\[
\widehat x=K_\theta\xi
\tag{4.13}
\]

不能无条件声称存在同半径 \(L_\theta\) 使 \(K_\theta=L_\theta G_\theta\)。全局上可先投影掉
\(\ker G_\theta\)，但所需投影和伪逆一般稠密，会破坏局部性。因此：

- **已证明：** 所有稳定的“先算 \(b=G\xi\)，再用局部 \(L b\)”线性算法；
- **未证明：** 所有直接读取原始边量测、允许额外边状态或利用 \(\ker G\) 的局部线性算法。

### 4.6 应用推论：轮数上下界、energy/error tradeoff 与复杂度

先把 STq 的节点分位数写成可计算参数。给定 \(a,\delta>0\)，令

\[
\eta_q(a,\delta)
:=
\sup\left\{
\eta\in(0,\Lambda):
\sup_\theta
m_\theta\{i:X_{\theta,i}^{(q)}(\eta)>a\}
\le\delta
\right\}.
\tag{4.14}
\]

若集合为空就令它为零。取任意
\(0<\eta<\eta_q(\varepsilon^2/2,\delta)\)，由
\((1-\eta/\Lambda)^{2t}\le \exp(-2t\eta/\Lambda)\)，只要

\[
t
\ge
\left\lceil
\frac{\Lambda}{2\eta}
\log_+\!\left(\frac{2}{\varepsilon^2\eta^q}\right)
\right\rceil,
\qquad
\log_+(x):=\max\{0,\log x\},
\tag{4.15}
\]

Richardson 就在至少 \(1-\delta\) 节点达到误差 \(\varepsilon\)。对白化 WLS 取
\(q=1\)；对满足定理 5.1 对易假设的 primitive RHS 取 \(q=2\)。

设 \(G\) 的量测到节点传播宽度为 \(h\)。则 sufficient-statistic-first 实现

\[
\xi\ \xrightarrow[\le h\text{ 轮}]{}\ b=G\xi
\ \xrightarrow[\le (t-1)s\text{ 轮}]{}\ 
p_t(J)b
\tag{4.16}
\]

总通信半径至多

\[
r_{\rm Rich}\le h+(t-1)s.
\tag{4.17}
\]

另一方面，直接以原始白化量测为输入，令

\[
T_\theta^{\rm raw}=J_\theta^{-1}G_\theta,\qquad
\tau_{\theta,i}^{\rm raw}(r)
=
\|E_iT_\theta^{\rm raw}Q_{i,r}^{\rm meas}\|.
\tag{4.18}
\]

命题 2.1 给任何 \(r\)-轮线性映射 \(K_\theta\xi\)——包括不因子化成
\(L_\theta G_\theta\xi\) 的更大类——以下不可绕过的下界：

\[
r_{\rm any}^{\star}(\varepsilon,\delta)
\ge
\inf\left\{
r:
\sup_\theta
m_\theta\{i:\tau_{\theta,i}^{\rm raw}(r)>\varepsilon\}
\le\delta
\right\}.
\tag{4.19}
\]

因此已证明的通信轮数夹逼是

\[
\boxed{
r_A^\star(\varepsilon,\delta;T^{\rm raw})
\le
r_{\rm any}^\star(\varepsilon,\delta)
\le
r_{K=LG}^\star(\varepsilon,\delta)
\le
h+(t-1)s.
}
\tag{4.20}
\]

第一项允许全局离线 advice，所以只能作信息论下界；最后一项是无需全局系数表的显式 C 类上界。若两边同阶，才可宣称通信最优；一般网络上目前没有自动同阶。

在常见的“每传一个标量耗能 \(c_{\rm tx}\)，每条 \(J\)-边每轮交换一个
\(d_x\)-维状态块”模型中，Richardson 的量级为

\[
\begin{array}{ll}
\text{标量传输/能量:}
&
O\!\left(t\,{\rm nnz}_{\rm blk}(J)d_x\right),
\\[2mm]
\text{全网算术:}
&
O\!\left(t\,{\rm nnz}_{\rm blk}(J)d_x^2\right),
\\[2mm]
\text{节点 }i\text{ 工作内存:}
&
O\!\left((\deg_J(i)+1)d_x\right)
\ \text{个动态标量，另存本地块系数。}
\end{array}
\tag{4.21}
\]

量测汇总 \(G\xi\) 的一次性代价再加
\(O({\rm nnz}_{\rm blk}(G))\)。递推

\[
x_{k+1}=(I-J/\Lambda)x_k+b/\Lambda,\qquad x_0=0,
\tag{4.22}
\]

只保留当前 \(x_k\) 和本地 \(b\)，不需存整个 \(r\)-球。相比之下，直接实现 A 类最优截断行通常要先收集并存储球内全部输入，节点内存至少随
\(|B_r(i)|\) 增长，而且系数表来自全局逆。

若对固定 \(\varepsilon,\delta\)，
\(\eta_1(\varepsilon^2/2,\delta)\) 在网络规模上有正下界，则 (4.15) 给出与规模无关的轮数和有界度网络上的每节点常数量级通信；全局原始数据泛洪至少受网络直径约束，并需增长的 payload/存储。反之若谱尾分位数趋零，(4.15) 会明确显示延迟和能耗怎样恶化。这个比较依赖固定大小消息和每跳能耗模型；不能脱离物理层直接声称“局部通信一定更省电”。

**状态：定理 4.1 与 oracle 投影的直接应用推论；给出了可量化收益，但不是新的独立数学定理。**

---

## 5. primitive RHS 的 \(q=2\)：对易时成立，不对易时会坏

### 定理 5.1（对易稳定规则的 \(q>0\) 版本）

把 ST1 中的 \(X^{(1)}\) 换成

\[
\boxed{
\forall a>0:\quad
\lim_{\eta\downarrow0}
\sup_\theta
m_\theta\{i:X_{\theta,i}^{(q)}(\eta)>a\}=0.
}
\tag{STq}
\]

假设每个候选规则除 (4.4) 外还满足

\[
L_\theta J_\theta=J_\theta L_\theta.
\tag{5.1}
\]

则对任意 \(q>0\)，以下等价：

- 任意 \(\varepsilon,\delta\) 下存在稳定、对易、有限传播的多数节点局部求解器；
- STq；
- 共同 Richardson 多项式多数节点可行。

尤其 \(q=2\) 精确覆盖 primitive RHS 误差。

**证明。** 充分性同 (4.10)，一般形式是

\[
\bigl(E_{\theta,i}^{(q)}(p_t(J_\theta))\bigr)^2
\le
X_{\theta,i}^{(q)}(\eta)
+
\frac{(1-\eta/\Lambda)^{2t}}{\eta^q}.
\tag{5.2}
\]

必要性中，对易性使 \(P_{\theta,\eta}\) 约化 \(L_\theta\)。令
\(v=P_{\theta,\eta}E_i^{*}y\)，则

\[
\begin{aligned}
\|J^{-q/2}(I-JL^{*})v\|
&=
\|J^{-q/2}v-L^{*}J^{1-q/2}v\|\\
&\ge
\|J^{-q/2}v\|-M\|J^{1-q/2}v\|\\
&\ge
(1-M\eta)\|J^{-q/2}v\|.
\end{aligned}
\tag{5.3}
\]

低谱与其正交补不混合，故总误差至少为上式。若 \(M\eta\le1/2\)，则

\[
E_{\theta,i}^{(q)}(L_\theta)
\ge \tfrac12\sqrt{X_{\theta,i}^{(q)}(\eta)}.
\tag{5.4}
\]

随后量词与定理 4.1 相同。证毕。

**状态：谱演算的直接推论；“多数节点 + 跨模型稳定规则 + Richardson 普适”是可能的新组合。**

完全对易不是逻辑最弱假设。上述证明实际只需 \(L_\theta^{*}\) 约化所用的小谱子空间，即
\(L_\theta^{*}P_{\theta,\eta}=P_{\theta,\eta}L_\theta^{*}\)。但该条件依赖未知全局谱投影，不具局部可检验性；全对易是最干净的结构充分条件。

### 5.2 为什么 \(q=1\) 可不对易而 \(q=2\) 的同一证明不行

不对易时仍有一般估计

\[
\sqrt{X_{\theta,i}^{(q)}(\eta)}
\le
E_{\theta,i}^{(q)}(L_\theta)
+
M\eta^{1-q/2},
\qquad 0<q<2.
\tag{5.5}
\]

所以定理 4.1 实际可扩到所有 \(0<q<2\)。在 \(q=2\) 时第二项退化成 \(M\)，不再随 \(\eta\downarrow0\) 消失；\(q>2\) 更坏。

### 5.3 一个固定失败率下的稳定局部反例

固定整数 \(k\ge2\)，把每个 \(k\)-节点块视为一个 clique，其中一个节点由模型系数内禀地区分为 anchor。令

\[
u_\lambda
=
\frac{(1,\lambda,\ldots,\lambda)^\top}
{\sqrt{1+(k-1)\lambda^2}},
\qquad
J_\lambda
=
\lambda u_\lambda u_\lambda^{*}
+
(I-u_\lambda u_\lambda^{*}),
\tag{5.6}
\]

\(\lambda\downarrow0\)。\(J_\lambda\) 的谱是 \(\{\lambda,1,\ldots,1\}\)，在 clique 上是一跳矩阵。

对任一非 anchor 节点 \(i\)，只要 \(\lambda<\eta<1\)，

\[
X_{\lambda,i}^{(2)}(\eta)
=
\frac{|u_\lambda(i)|^2}{\lambda^2}
=
\frac1{1+(k-1)\lambda^2}
\longrightarrow1.
\tag{5.7}
\]

现在令 \(D_{\rm good}\) 删除 anchor 输出行，并取

\[
L_\lambda=D_{\rm good}J_\lambda^{-1}.
\tag{5.8}
\]

则：

- \(L_\lambda\) 对 \(k-1\) 个好节点与 \(J_\lambda^{-1}\) 完全相同，所以这些节点的 \(q=2\) 误差为零；
- \(\sup_{\lambda}\|L_\lambda\|\le 1+\sqrt{k-1}\)，故它稳定；
- clique 的 1-hop view 含整个标记块，anchor 由系数可辨认，所以这是 B/C 可实现的一跳规则；
- \(L_\lambda\) 不与 \(J_\lambda\) 对易。

因此在固定失败率 \(\delta=1/k\) 下，可以精确估计至少 \(1-\delta\) 节点，而同一批好节点的 \(X^{(2)}\) 尾不消失。任取更大的固定 \(k\)，反例失败率可预先做得任意小。

这严格否定了“把定理 4.1 的不对易证明原样改成 \(q=2\)”以及固定 \(\delta\) 的无条件谱尾必要性。它**尚未**否定一个量词为“同一模型族对每个 \(\delta>0\) 都可行”的更复杂无对易定理；该全量词版本目前未闭合。

---

## 6. 哪些自然结构让 C 类自动落入对易/多移位类

### 6.1 距离传递图：完整 C 类就是单移位多项式

设 \(G\) 是有限连通距离传递图，邻接矩阵为 \(A\)，距离-\(h\) 矩阵为 \(A_h\)。任意标量、线性、\(\operatorname{Aut}(G)\)-等变算子 \(L\) 的元素只依赖有序点对的距离，因此

\[
L=\sum_{h=0}^{D}c_hA_h.
\tag{6.1}
\]

若 \(L\) 的传播宽度至多 \(r\)，则 \(c_h=0\) 对所有 \(h>r\)。

距离传递推出距离正则；Bose–Mesner 三项递推又给出

\[
A_h=v_h(A),\qquad \deg v_h=h.
\tag{6.2}
\]

所以

\[
\boxed{
\{\text{等变的 }r\text{-传播线性算子}\}
=
\{p(A):\deg p\le r\}.
}
\tag{6.3}
\]

这是 association scheme 的经典事实，见 Brouwer–Cohen–Neumaier 第 2、4、7 章。[^2] Backhausz–Virág 也在正则树情形明确写出：每个有限支撑 radial 函数都是 \(p(A)\delta_o\)。[^3]

若 \(J=\alpha I+\beta A\) 且 \(\beta\ne0\)，(6.3) 等价于 \(J\) 的多项式。更一般地，即使 \(J=f(A)\) 未生成整个 Bose–Mesner 代数，所有 C 算子仍与 \(J\) 对易，故定理 5.1 仍覆盖完整 C 类。

**限制：**

- 单个传递图上所有节点同质，\(\delta\)-多数语义退化成 all-or-none；
- 标记、边权或节点势打破对称后，结论不再自动成立；
- 仅顶点传递远远不够，置换群的 centralizer 可非交换；
- 唯一 ID 会破坏 (6.1) 的内禀等变含义。

### 6.2 阿贝尔 Cayley/格点：完整 C 类是多移位 Laurent 多项式

在有限阿贝尔群或 \(\mathbb Z^d\) 上，平移等变线性算子恰是卷积：

\[
(Lx)(v)=\sum_{g}k(g)x(v-g).
\tag{6.4}
\]

有限传播等价于 \(k\) 有限支撑，故 \(L\) 是坐标平移的 Laurent 多项式。阿贝尔卷积彼此对易，所以若 \(J\) 也是有限支撑卷积，定理 5.1 覆盖所有平移等变有限传播线性算子。

在 \(d\ge2\) 时，这一类一般严格大于单个 Laplacian \(J\) 的多项式：各向异性的有限支撑卷积不必只依赖 Laplacian 符号。因此正确术语是“多移位滤波”，不是硬说所有算子都是 \(p(J)\)。Emirov–Cheng–Jiang–Sun 系统研究了多个可交换 shift 的多项式滤波及分布式逆滤波。[^4]

### 6.3 一般图上的反例

在带节点势的图上令

\[
J=A+D.
\]

0-hop 等变规则 \(L=f(D)\) 是局部对角算子；只要相邻节点的 \(D\) 不同，通常
\([L,J]\ne0\)。即便无标记， irregular 图上的 degree-based 对角规则也通常不与邻接或 Laplacian 对易。

图滤波文献本身也明确区分 node-invariant polynomial filters 与更宽的 node-variant filters；Segarra–Marques–Ribeiro 的 Proposition 4 给出 node-variant 精确实现条件。[^5] 所以在一般网络上把“所有局部线性算子”直接等同于单移位多项式，会与既有定义冲突。

---

## 7. 格点上的 sharp phase transition：完整局部类与多项式同阶

这一节给出一个能检验定理 2.1、4.1、5.1 是否尺度正确的模型，但不作为原创主核。

### 7.1 统一 \(q\) 模型

令 \(\Delta_d\) 是 \(\mathbb Z^d\) 上的正离散 Laplacian，其 Fourier 符号

\[
\lambda(\omega)
=
2\sum_{j=1}^{d}(1-\cos\omega_j)
\asymp|\omega|^2
\quad(\omega\to0).
\tag{7.1}
\]

考虑平稳目标

\[
T_q=\Delta_d^{-q/2}.
\tag{7.2}
\]

严格说，无限格点上的 \(\Delta_d\) 没有有界逆；以下把
\(\Delta_d^{-q/2}\) 当作谱演算得到的稠密定义算子，并询问其根行是否属于
\(\ell_2\)。有限版本可取 torus 上常数模态正交补的 Moore–Penrose 逆，或先加质量
\(\kappa I\) 再令 \(\kappa\downarrow0\)。所以本节是 (1.1) 的临界极限模型，不是把奇异 Laplacian 偷写成 SPD。

- \(q=2\)：\(T_2=\Delta^{-1}\)，对应 primitive white RHS；
- \(q=1\)：\(T_1=\Delta^{-1/2}\)，其输出协方差是 \(\Delta^{-1}\)，对应 GFF，也与白化 edge-WLS 的
  \(\Delta^{-1}B^{*}\) 有相同节点误差谱权。

根谱测度在零附近满足

\[
d\mu_o(\lambda)
\asymp
\lambda^{d/2-1}\,d\lambda.
\tag{7.3}
\]

因此

\[
\int_0 \lambda^{-q}\,d\mu_o(\lambda)<\infty
\quad\Longleftrightarrow\quad
d>2q.
\tag{7.4}
\]

### 命题 7.1（所有有限支撑卷积的最佳误差）

当 \(d>2q\) 时，\(T_q\) 的卷积核满足

\[
|t_q(x)|\asymp |x|^{q-d}
\tag{7.5}
\]

（对 \(q=2\) 是 Green 核；\(q=1\) 可读成 fractional Green，edge-WLS 则是普通 Green 的梯度）。于是命题 2.1 给出

\[
\inf_{\operatorname{prop}(K)\le r}
\|e_0^{*}(T_q-K)\|_2^2
=
\sum_{|x|>r}|t_q(x)|^2
=
\Theta(r^{2q-d}),
\tag{7.6}
\]

从而最优误差

\[
\boxed{
\Theta\!\left(r^{-(d-2q)/2}\right).
}
\tag{7.7}
\]

这里下确界可取遍所有平移等变有限传播线性算子，最优者就是核的空间截断。Green 函数
\(|x|^{2-d}\) 的经典渐近可见 Lawler–Limic 第 4 章或 Uchiyama。[^6]

在临界 \(d=2q\) 时，无限格点行不在 \(\ell_2\)；边长 \(N\) 的均值零 torus 上固定 \(r\) 的平方误差按
\(\Theta(\log(N/r))\) 发散。若 \(d<2q\)，主阶为 \(\Theta(N^{2q-d})\)。因此：

| 输入/目标 | \(q\) | 可统一局部逼近的维数 | 临界维数 | 最优误差 |
|---|---:|---:|---:|---:|
| primitive \(b\mapsto\Delta^{-1}b\) | 2 | \(d>4\) | 4 | \(\Theta(r^{-(d-4)/2})\) |
| whitened edge-WLS / GFF | 1 | \(d>2\) | 2 | \(\Theta(r^{-(d-2)/2})\) |

### 命题 7.2（最佳单移位多项式达到同阶）

次数至多 \(r-1\) 的 \(p\) 对应残差多项式

\[
Q(\lambda)=1-\lambda p(\lambda),
\qquad
\deg Q\le r,\quad Q(0)=1.
\tag{7.8}
\]

其平方误差是

\[
\int
|Q(\lambda)|^2\lambda^{-q}\,d\mu_o(\lambda).
\tag{7.9}
\]

所以最小值正是权

\[
w(\lambda)\asymp
\lambda^{\beta},
\qquad
\beta=d/2-1-q
\tag{7.10}
\]

在硬端点 \(0\) 的 Christoffel 函数。若 \(d>2q\)，即 \(\beta>-1\)，经典 Jacobi 型端点渐近给出

\[
\inf_{\deg p\le r-1}
\int|1-\lambda p(\lambda)|^2\lambda^{-q}\,d\mu_o
=
\Theta(r^{-2(\beta+1)})
=
\Theta(r^{-(d-2q)}).
\tag{7.11}
\]

因此单移位最优多项式虽然只是所有有限支撑卷积的子类，却达到 (7.7) 的相同阶；最小轮数为

\[
r^\star(\varepsilon)
=
\Theta\!\left(
\varepsilon^{-2/(d-2q)}
\right).
\tag{7.12}
\]

Christoffel 函数的定义本身就是带 \(Q(0)=1\) 的最小 \(L_2\) 多项式问题；广义 Jacobi 权的端点阶见 Lubinsky 等工作。[^7] 简单 Richardson 由 (5.2) 可给定性充分性，但直接优化该上界通常多一个对数；要得到 (7.11) 的 sharp 阶应使用 Christoffel/Jacobi 最优残差。

**状态：经典 Green 渐近 + 经典 Christoffel 渐近的直接组合。**

### 7.3 与 factor-of-i.i.d. 的直接撞车

Backhausz–Virág 的 Theorem 2 在一个固定的无限顶点传递图上证明（过程均为零均值且有有限二阶矩）：一个有限 Borel 测度是某个 linear factor-of-i.i.d. 过程的谱测度，当且仅当它相对于图谱测度绝对连续；该定理甚至把这与一般 factor-of-i.i.d. 及其 \(\bar d_2\) 极限的谱测度刻画列为等价。Definition 14 和 Proposition 24 还处理多项式生成的 spherical linear factors 及其
\(\bar d_2\) 闭包。[^3]

对 \(Y=J^{-q/2}Z\)，目标谱测度密度就是 \(\lambda^{-q}\)。有限方差条件正是

\[
\int\lambda^{-q}\,d\mu_o<\infty.
\tag{7.13}
\]

其 Proposition 26 直接证明 transient 顶点传递图上的 GFF 是 linear factor-of-i.i.d.；这就是单一齐次无限模型中 \(q=1\) 的定性侧。在 \(\mathbb Z^d\) 上 transience 临界为 \(d=2\)。\(q=2\) 虽未以“primitive RHS state estimation”命名，也立即落在 Theorem 2 的同一谱密度框架中。因而定理 4.1 不能以“首次给出 \(q=1\) 谱判据”表述；候选增量只可能是它的 uniform finite-family、majority-node、非对易稳定算法反推共同 Richardson 这组量词。

所以：

- \(d>2\) / \(d>4\) 的定性阈值不是新发现；
- 空间截断最优是 Hilbert 投影；
- polynomial sharp rate 是 Christoffel 端点问题；
- 可保留的价值是把这些结果翻译成 WLS 通信轮数并与一般非齐次多数节点定理 4.1 对照。

---

## 8. 与现有领域的精确关系

### 8.1 图滤波与 distributed inverse

Segarra–Marques–Ribeiro 把 node-invariant graph filter 定义成 shift 的矩阵多项式，研究任意线性网络算子的精确/近似实现，并进一步引入 node-variant filters；论文还明确把 distributed state estimation 列为动机。[^5] 因而：

- “多项式可在有限轮实现”已知；
- “node-variant 比共同多项式更宽”已知；
- 本文不能把局部线性算子按定义缩成 polynomial；
- 定理 4.1 的新意若成立，在于反向结论：任意稳定 node-variant 正规方程求解器的多数节点可行性，强迫一个共同 Richardson 也可行。

Emirov 等对多个可交换 graph shifts、分布式 polynomial inverse 和迭代 Chebyshev/Richardson 型实现给出系统算法。[^4] 这覆盖定理 5.1 的构造侧，不覆盖本文的逐节点多数必要性。

### 8.2 shift-enabled 与重复特征值

“与 shift 对易”不自动等于“是 shift 的多项式”。shift-enabled 文献专门研究何时二者相等；重复特征值时，commutant 内可有在同一特征空间中非标量作用的算子。[^8] 因此定理 5.1 有意只要求对易，充分算法才选 \(p(J)\)。

### 8.3 非交换 \(L_p\)、affiliated inverse 与 measure topology

Nelson 的非交换积分和 Fack–Kosaki 的广义奇异数理论包含以下标准框架：无界 affiliated operator 可由谱截断逼近；trace/measure topology 和非交换 \(L_p\) 用谱投影的 trace 控制尾。[^9][^10]

这与本文相邻，但不是同一个量：

- 非交换 \(L_p\) 的典型条件是 trace 平均

  \[
  \tau(J^{-q}\mathbf1_{(0,\eta)}(J))\to0;
  \]

- measure topology 主要看 \(|J^{-1}|\) 的大值谱投影 trace；
- 本文看的是固定节点基底中的对角块

  \[
  X_i^{(q)}(\eta)
  =
  \|J^{-q/2}P_\eta E_i^{*}\|^2
  \]

  对随机节点是否依概率消失。

trace 尾通过 Markov 不等式推出节点 in-probability 尾；反向不成立，因为少量节点可承载任意大的能量。本文条件依赖网络节点 masa/坐标，不是酉不变的非交换 \(L_p\) 条件。因此 affiliated inverse 的 bounded spectral truncation 是重要先例和语言来源，但目前没有发现它直接包含定理 4.1 的逐节点、允许少量病态区、finite-propagation Richardson 结论。

### 8.4 局部弱收敛与谱测度

Benjamini–Schramm 收敛对有界连续谱函数或多项式矩给出自然收敛；负次矩在零处无界，必须额外要求 uniform integrability。本文 ST1 正是一种**节点级**而非 trace 平均级的负矩尾条件。已有谱收敛文献提供背景，但仅有局部弱收敛不推出 ST1。

Berthier–Bach–Gaillard 已用逐顶点 Laplacian 谱测度定义定量 spectral dimension，并据此给 gossip/averaging 的规模无关多项式收敛率；其 Proposition 1 专门验证 torus，Corollary 1 给算法率。[^11] 这与第 7 节的“谱维数控制迭代复杂度”高度相邻，但目标是平均过程而非不稳定逆，条件是上界
\(\mu_i((0,E])\lesssim E^{d/2}\)，没有给定理 4.1 的负矩节点尾 iff 或“任意稳定求解器推出 Richardson”反向结论。

---

## 9. 反例清单：每个假设在挡什么

### 9.1 去掉稳定性

\[
J_\alpha=\alpha I,\quad
L_\alpha=\alpha^{-1}I.
\]

零轮精确，谱尾失败。说明定理 4.1 必须控制跨模型增益。

### 9.2 把 A 当成 B

Schur 族 (3.3) 的可见块相同、外部 completion 不同，局部逆自系数相差 \(1/3\)。A 可预装不同值，B/C 不可。

### 9.3 把 C 当成 polynomial

一般带权图的 \(L=f(D)\) 是 0-hop 等变规则，却通常不与 \(J=A+D\) 对易，更不必是 \(p(J)\)。

### 9.4 把 \(q=1\) 证明复制到 \(q=2\)

clique-anchor 族 (5.6)–(5.8) 在 \(1-1/k\) 节点精确、稳定、一跳，但这些好节点的
\(X^{(2)}\to1\)。这是逐节点 sharp 反例。

### 9.5 把白化量测算法等同于任意 RHS

白化 WLS 有 \(GG^{*}=J\)，故误差只含 \(\lambda^{-1}\)；primitive white RHS 含
\(\lambda^{-2}\)。两者临界维数差 2，混用会直接得出错误拓扑要求。

### 9.6 忽略量测局部性

若 \(R^{-1/2}\) 稠密，即使 \(H\) 稀疏，\(G=H^{*}R^{-1/2}\) 也可能全局；此时先形成
\(b=G\xi\) 已经需要全局通信，定理 4.1 不能被解释为原网络上的局部算法。

---

## 10. 原创性判决与最小可发表版本

### 10.1 已经撞车或明显不足

以下内容不能单独做论文主定理：

1. A 类最优行是 \(J^{-1}\) 的 \(r\)-hop 截断；
2. A 类 majority iff 是行尾分位数；
3. 距离传递图的等变有限支撑核是邻接多项式；
4. 阿贝尔平移等变有限支撑算子是多移位卷积；
5. GFF/primitive forcing 在格点的 \(d=2/d=4\) 阈值；
6. Green 尾与 Christoffel 多项式给出的 sharp 幂率；
7. 共同 polynomial graph filter 的分布式实现。

它们分别属于 Hilbert 投影、Bose–Mesner 代数、Fourier 卷积、factor-of-i.i.d.、随机游走 Green 函数和正交多项式的经典范围。

### 10.2 最小候选主定理

建议只考虑如下精确标题范围：

> 对统一上界、有限传播 SPD 正规矩阵族，白化 WLS 的本地充分统计量由任意跨模型稳定的有限传播线性求解器在几乎所有节点逼近，当且仅当节点小谱逆一次矩一致依概率消失；若成立，一个共同 Richardson 规则已经普适。

它比若干已知部件多出的实质是：

- 算法侧量化遍历所有稳定、view-dependent、甚至不对易的有限传播 \(L_\theta\)；
- 结构侧只看每个节点的低谱逆一次矩尾，允许条件数发散和少数任意病态点；
- 结论把任意可行规则压成同一个显式多项式，而不是仅证明存在某个局部近似；
- 必要性与充分性逐节点匹配，并自然产生 fixed-\(\delta\) 的 slack 边界。

**当前标签：可能的新组合，不能宣称原创。**

### 10.3 是否已跨过“不是拼接恒等式”的门槛

数学上，定理 4.1 已跨过单纯分位数恒等式门槛：它从任意稳定算法的存在推出一个可计算的小谱条件，再推出一个固定算法族。其必要性不依赖候选算法是多项式，这一点是真正的压缩结论。

但要支撑整篇论文，还差至少一项：

- 把直接原始量测算子 \(K\xi\) 也纳入，而不要求先因子化为 \(LG\xi\)；
- 给出 ST1 可由局部样本/有限球可靠认证的附加结构；
- 在非齐次随机网络或缺陷网络上给出非平凡 \(\delta\) phase diagram 与 matching communication lower bound；
- 找到该精确 iff 的现有先例或完成更系统的 operator-algebra / local weak convergence 查重。

若这些都做不到，最诚实结论是：**它适合作为强 lemma 或短理论 note，尚不足以支撑“彻底终结所有局部状态估计算法”的大标题。**

### 10.4 未闭合问题

1. \(q=2\) 在不对易规则下、量词为“同一模型族对所有 \(\delta>0\) 可行”时，是否仍可推出某个比 ST2 弱的可检验条件？
2. 对任意局部 \(K_\theta\xi\)，能否在保持多数节点误差和传播半径至常数倍的同时投影成
   \(L_\theta G_\theta\xi\)？
3. ST1 是节点基底依赖的 in-probability 条件；它与 measured graph operator algebra 中哪一种带 masa 的局部 \(L_2\) 可积性完全等价？
4. 在有界度局部弱收敛模型上，哪些可检验条件保证 ST1，而不退回全局统一谱隙？

---

## Sources

[^1]: M. J. Grote and T. Huckle, “Parallel Preconditioning with Sparse Approximate Inverses,” *SIAM Journal on Scientific Computing* 18(3), 838–853 (1997), DOI [10.1137/S1064827594276552](https://doi.org/10.1137/S1064827594276552). 该文是 SPAI 的经典来源；其目标是使 \(AM\) 或 \(MA\) 接近单位阵，而非假定已知精确逆后做坐标截断。

[^2]: A. E. Brouwer, A. M. Cohen, and A. Neumaier, *Distance-Regular Graphs*, Springer (1989), especially Chapters 2, 4, and 7, DOI [10.1007/978-3-642-74341-2](https://doi.org/10.1007/978-3-642-74341-2). 距离矩阵形成 Bose–Mesner 代数以及 \(A_h=v_h(A)\) 是标准背景。

[^3]: Á. Backhausz and B. Virág, “Spectral measures of factor of i.i.d. processes on vertex-transitive graphs,” *Annales de l’Institut Henri Poincaré Probabilités et Statistiques* 53(4), 2260–2278 (2017), DOI [10.1214/16-AIHP790](https://doi.org/10.1214/16-AIHP790), arXiv [1505.07412](https://arxiv.org/abs/1505.07412). 关键位置：Theorem 2（固定无限顶点传递图、有限二阶矩过程的谱测度 iff）；Definition 13–14（linear/spherical factors，论文 pp. 7–8）；Proposition 24（\(\bar d_2\) 闭包，p. 19）；Definition 25 与 Proposition 26（transient 图上的 GFF，pp. 20–21）。本地 PDF：`literature/02_graphical_models_local_inference/2017_backhausz_virag_spectral_measures_factor_iid.pdf`；抽取文本：`research/corpus/02_graphical_models_local_inference/2017_backhausz_virag_spectral_measures_factor_iid.txt`。

[^4]: N. Emirov, C. Cheng, J. Jiang, and Q. Sun, “Polynomial graph filters of multiple shifts and distributed implementation of inverse filtering,” *Sampling Theory, Signal Processing, and Data Analysis* 20, Article 2 (2022), DOI [10.1007/s43670-021-00019-x](https://doi.org/10.1007/s43670-021-00019-x). 关键位置：Eq. (1.1)–(1.2) 定义可交换 shifts 与多变量多项式；Algorithm 2.2 给分布式实现；Theorems 3.1, 4.2, 4.4 给逆滤波迭代收敛。

[^5]: S. Segarra, A. G. Marques, and A. Ribeiro, “Optimal Graph-Filter Design and Applications to Distributed Linear Network Operators,” *IEEE Transactions on Signal Processing* 65(15), 4117–4131 (2017), DOI [10.1109/TSP.2017.2703660](https://doi.org/10.1109/TSP.2017.2703660). 关键位置：Section II 定义 polynomial graph filters；Sections III–IV 分别处理 node-invariant 与 node-variant filters；Proposition 4 给 node-variant 精确实现条件。

[^6]: G. F. Lawler and V. Limic, “The Green’s function,” Chapter 4 in *Random Walk: A Modern Introduction*, Cambridge University Press (2010), DOI [10.1017/CBO9780511750854.005](https://doi.org/10.1017/CBO9780511750854.005), especially Theorem 4.3.1; K. Uchiyama, “Green’s Functions for Random Walks on \(\mathbb Z^N\),” *Proceedings of the London Mathematical Society* 77(1), 215–240 (1998), DOI [10.1112/S0024611598000458](https://doi.org/10.1112/S0024611598000458).

[^7]: D. S. Lubinsky, “Pointwise Asymptotics for Orthonormal Polynomials at the Endpoints of the Interval via Universality,” *International Mathematics Research Notices* 2020(4), 961–982, DOI [10.1093/imrn/rny042](https://doi.org/10.1093/imrn/rny042); T. Danka, “Universality limits for generalized Jacobi measures,” *Advances in Mathematics* 316, 613–666 (2017), DOI [10.1016/j.aim.2017.06.026](https://doi.org/10.1016/j.aim.2017.06.026). 两者给幂型端点权的 Christoffel/核渐近；本文只用其阶。

[^8]: L. Chen, S. Cheng, V. Stanković, and L. Stanković, “Shift-Enabled Graphs: Graphs Where Shift-Invariant Filters are Representable as Polynomials of Shift Operations,” *IEEE Signal Processing Letters* 25(9), 1305–1309 (2018), DOI [10.1109/LSP.2018.2849685](https://doi.org/10.1109/LSP.2018.2849685).

[^9]: E. Nelson, “Notes on Non-Commutative Integration,” *Journal of Functional Analysis* 15(2), 103–116 (1974), DOI [10.1016/0022-1236(74)90014-7](https://doi.org/10.1016/0022-1236(74)90014-7).

[^10]: T. Fack and H. Kosaki, “Generalized \(s\)-Numbers of \(\tau\)-Measurable Operators,” *Pacific Journal of Mathematics* 123(2), 269–300 (1986), DOI [10.2140/PJM.1986.123.269](https://doi.org/10.2140/PJM.1986.123.269).

[^11]: R. Berthier, F. Bach, and P. Gaillard, “Tight Nonparametric Convergence Rates for Stochastic Gradient Descent under the Noiseless Linear Model,” *Advances in Neural Information Processing Systems 33* (2020), [official proceedings paper](https://proceedings.neurips.cc/paper/2020/hash/1b33d16fc562464579b7199ca3114982-Abstract.html), arXiv [2006.08212](https://arxiv.org/abs/2006.08212). 相关位置：Section 3.2、Proposition 1、Corollary 1。

## 本地材料

- `literature/02_graphical_models_local_inference/2017_backhausz_virag_spectral_measures_factor_iid.pdf`
- `research/corpus/02_graphical_models_local_inference/2017_backhausz_virag_spectral_measures_factor_iid.txt`
- `literature/02_graphical_models_local_inference/SOURCE_URLS.md`
- `literature/03_sparse_inverse_graph_filters/2015_segarra_marques_ribeiro_distributed_linear_network_operators.pdf`
- `research/corpus/03_sparse_inverse_graph_filters/2015_segarra_marques_ribeiro_distributed_linear_network_operators.txt`
- `literature/03_sparse_inverse_graph_filters/2020_emirov_cheng_jiang_sun_polynomial_inverse_graph_filter.pdf`
- `research/corpus/03_sparse_inverse_graph_filters/2020_emirov_cheng_jiang_sun_polynomial_inverse_graph_filter.txt`
- `literature/03_sparse_inverse_graph_filters/1984_demko_moss_smith_decay_inverse_band_matrices.pdf`
- `research/agent_reports/majority_local_iff_audit.md`
- `research/agent_reports/theory_phase2.md`
