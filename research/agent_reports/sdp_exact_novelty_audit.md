# Theorem 4 / `(SDP-cont)` 专项原创性审计

**审计日期：** 2026-09-13  
**被审计对象：** `research/agent_reports/block_robust_boundary_attack.md` 的引理 3、定理 4 与 `(SDP-cont)`  
**判定目的：** 判断该 exact SDP 能否作为论文的强主定理，而不只是一个正确的应用性推论。

## 0. 执行结论

结论很明确：**定理 4 数学上是正确的，但就目前表述而言，不宜作为强主定理。** 它被单个 norm-bounded full block 的非严格 Petersen lemma（等价地，lossless full-block S-procedure）直接包含；所需工作不是新的 multiplier 理论，而是一次变量重写和 Schur complement。

更具体地，令

\[
S=\Sigma^{-1},\qquad
\alpha=M^{-1},\quad \beta=m^{-1},\quad
h=\frac{\alpha+\beta}{2},\quad d=\frac{\beta-\alpha}{2},
\]

并定义

\[
B_Q=
\begin{bmatrix}
hF^TC^T-Q^T\\ hC^T
\end{bmatrix},\qquad
G=d\begin{bmatrix}F^T\\I_n\end{bmatrix},\qquad
H=C^T.
\tag{A}
\]

这里对原矩阵的第二个输出块施加了一个无害的正交符号翻转
\(\operatorname{diag}(I_k,-I_n)\)；它不改变谱范数，所以可把下块统一写成正号。

则原问题精确等价于

\[
R(C,F;m,M)
=\inf_Q\sup_{\Delta^T\Delta\preceq I_n}
\|B_Q+G\Delta H\|_2.
\tag{B}
\]

对候选上界 \(\sqrt\tau\)，把谱范数写成标准 Schur LMI 后，(B) 正是 Petersen lemma 的教科书型输入

\[
\mathcal C+\mathcal E\Delta\mathcal G
+\mathcal G^T\Delta^T\mathcal E^T\preceq0
\quad\forall\Delta^T\Delta\preceq I.
\]

van Waarde–Camlibel–Eising–Trentelman 2023 的 Proposition 4.16(b) 于是逐字给出一个标量 multiplier 的必要充分 LMI。[^1] 该论文还明确说明，这个命题就是非严格 Petersen lemma，并把它置于 Scherer 的 single-full-block lossless S-procedure 谱系中；其 Proposition 4.14 进一步说明 Scherer 2001, Theorem 5.3 对单个 full block 的 relaxation 是 exact。[^1][^2]

因此本次审计的判定是：

- **正确性：高置信通过。**
- **“首次得到有限 exact SDP”：不通过。** 公式可能没有以完全相同的 \((C,F,Q,m,M)\) 符号发表过，但作为数学定理，它是一行现成 Petersen lemma 的直接代入。
- **locality-constrained \([Q,0]\) 的作用：** 它令应用问题有意义，但 \(Q\) 只仿射进入 nominal block \(B_Q\)；再加仿射稀疏/局部性约束不会改变 full-block 消元。因此它尚未产生新的 lossless result。
- **任意 self-adjoint matrix interval 的自然推广：也不是新主定理。** 其椭球像是既有 QMI image / contraction completion 的直接特例，随后仍是一块 Petersen LMI。
- **可保留的论文价值：** 作为“网络边界近似问题如何精确落入 single-full-block robust model matching”的应用性命题、算法接口和反例工具，价值仍然不错；但不能承担全文最强数学贡献。

## 1. 从原问题到一个标准 full-block 问题

原 robust radius 是

\[
R=\inf_Q\sup_{\alpha I\preceq S\preceq\beta I}
\|[CSF-Q,-CS]\|_2.
\tag{1}
\]

对任意矩阵 \(A\in\mathbb R^{p\times r}\)，

\[
\|A\|_2=\sup_{\|u\|=1}\|A^Tu\|_2.
\]

固定左测试向量 \(u\)，置 \(v=C^Tu\) 和 \(z=Sv\)。报告的引理 3 给出

\[
\{Sv:\alpha I\preceq S\preceq\beta I\}
=\{hv+dw:\|w\|\le\|v\|\}.
\tag{2}
\]

任意满足 \(\|w\|\le\|v\|\) 的 \(w\) 都可写成 \(w=\Delta v\)，其中 \(\Delta^T\Delta\preceq I\)。反过来任意 contraction 当然满足该范数界。虽然原 \(S\) 是 self-adjoint，而这里的 \(\Delta\) 未要求 self-adjoint，但二者对**一个给定向量**的像完全相同；Householder 构造能把每个可行 \(w\) 由一个 self-adjoint contraction 实现。因此这个替换没有扩大联合 supremum。

代入 \(z=hC^Tu+d\Delta C^Tu\)，有

\[
\begin{aligned}
&\|F^Tz-Q^Tu\|^2+\|z\|^2\\
&\quad=\left\|
\left(
\begin{bmatrix}hF^TC^T-Q^T\\hC^T\end{bmatrix}
+d\begin{bmatrix}F^T\\I_n\end{bmatrix}\Delta C^T
\right)u\right\|^2.
\end{aligned}
\]

交换两个 supremum 后即得到 (B)。这一等价是整个审计的关键：**引理 3 已经把 Loewner interval 变成一个单 full block；定理 4 的剩余部分就是标准 robust-performance elimination。**

### 1.1 与 Petersen lemma 的逐变量映射

记 \(r=k+n\)，并按 (A) 定义 \(B_Q\in\mathbb R^{r\times p}\)、\(G\in\mathbb R^{r\times n}\)、\(H\in\mathbb R^{n\times p}\)。对固定 \((Q,\tau)\)，令

\[
\mathcal C=
\begin{bmatrix}
-\tau I_p&B_Q^T\\
B_Q&-I_r
\end{bmatrix},\quad
\mathcal E=\begin{bmatrix}0_{p\times n}\\G\end{bmatrix},\quad
\mathcal G=\begin{bmatrix}H&0_{n\times r}\end{bmatrix}.
\tag{3}
\]

则各对象与 van Waarde et al. 2023, Proposition 4.16 的符号一一对应如下。

| 本问题 | Proposition 4.16 中的对象 | 含义 |
|---|---|---|
| \(\Delta\in\mathbb R^{n\times n}\) | \(F\) | norm-bounded full block |
| \(\Delta^T\Delta\preceq I_n\) | \(F^TF\preceq\bar F\), \(\bar F=I_n\) | uncertainty set |
| \(\mathcal C\) | \(C=C^T\) | nominal symmetric LMI |
| \(\mathcal E\) | \(E\) | uncertainty left factor |
| \(\mathcal G\) | \(G\) | uncertainty right factor |
| Petersen multiplier \(\varepsilon>0\) | \(\lambda>0\) | scalar lossless multiplier |

Schur complement 给出

\[
\|B_Q+G\Delta H\|_2^2\le\tau
\iff
\mathcal C+\mathcal E\Delta\mathcal G
+\mathcal G^T\Delta^T\mathcal E^T\preceq0.
\tag{4}
\]

Proposition 4.16(b) 在 \(\mathcal E\ne0\)、\(I_n\succ0\)、\(\mathcal G\ne0\) 下断言，(4) 对所有 contraction 成立，当且仅当存在 \(\varepsilon>0\) 使

\[
\boxed{
\mathcal C+\varepsilon\mathcal E\mathcal E^T
+\varepsilon^{-1}\mathcal G^T\mathcal G\preceq0.}
\tag{5}
\]

本报告定理 4 的非退化假设 \(m<M\) 和 \(C\ne0\) 正好保证这些条件：\(G\ne0\) 是因为其下块含 \(dI_n\)，而 \(H=C^T\ne0\)。再对 \(\varepsilon^{-1}\mathcal G^T\mathcal G\) 作一次 Schur complement，(5) 变为完全仿射的 LMI

\[
\boxed{
\begin{bmatrix}
-\tau I_p&B_Q^T&H^T\\
B_Q&-I_r+\varepsilon GG^T&0\\
H&0&-\varepsilon I_n
\end{bmatrix}\preceq0.}
\tag{6}
\]

最小化 \(\tau\) 并添加任意仿射 \(Q\)-约束，仍是 exact SDP。

### 1.2 (6) 与报告 `(SDP-cont)` 的关系

报告使用变量 \(w=[z;u]\)，先以球约束

\[
\|z-hC^Tu\|^2\le d^2\|C^Tu\|^2
\]

构成一个标量二次 implication，再调用普通 S-lemma。式 (6) 则先写

\[
z=hC^Tu+d\Delta C^Tu,qquad\Delta^T\Delta\preceq I,
\]

再调用 Petersen lemma。两条路线是同一个单二次约束消元的两种坐标表示。对 \(d>0\)，变量替换 \(z\leftrightarrow\Delta C^Tu\) 保持可行像，S-lemma 证书在相应 congruence 与 multiplier 重标度下给出同一 epigraph。

换言之，`(SDP-cont)` 并非偶然“长得像”已有 LMI；它与现成 Petersen LMI 证明的是同一个 robust inequality。报告中允许 \(\lambda\ge0\) 不构成区别：在非退化问题中 \(\lambda=0\) 不可能满足其 LMI，因为 \(L_Q^TL_Q\) 的 \(z\)-块包含 \(FF^T+I_n\succ0\)，所以可行 multiplier 实际必为正。

## 2. 最近的现成定理及碰撞强度

### 2.1 最直接：Petersen lemma / single full block

最接近且足以终止 novelty claim 的来源是：

> H. J. van Waarde, M. Kanat Camlibel, J. Eising, H. L. Trentelman, “Quadratic Matrix Inequalities with Applications to Data-Based Control,” *SIAM Journal on Control and Optimization* 61(4), 2251–2281 (2023), DOI 10.1137/22M1486807.[^1]

关键位置：

- **Theorem 3.4，journal p. 2258 左右（作者稿 pp. 8–9）：** QMI solution set 在右乘线性映射下的 exact image；当映射满列秩或 QMI 的右下块可逆时，包含关系升级为等号。
- **Proposition 4.14，journal p. 2265 左右（作者稿 pp. 15–16）：** 把 single full-block S-procedure 的 exactness 写成 matrix S-lemma。正文明确说明它由 Scherer 2001, Eq. (1.2) 与 Theorem 5.3 直接得到。
- **Proposition 4.16(b)，journal p. 2267 左右（作者稿 p. 17）：** 非严格 Petersen lemma；在 \(E\ne0,\bar F\succ0,G\ne0\) 时是必要充分条件。这正是式 (3)–(5) 所用的一行定理。

更早的直接谱系包括 Petersen 1987 的 uncertain-system stabilization lemma、Scherer 1997 的 full-block S-procedure，以及 Scherer 2001 对单 full block 的 losslessness。[^2][^3][^4] Bisoffi–De Persis–Tesi 2022 的 Facts 1–2 又以现代、便于直接使用的形式列出 strict/nonstrict Petersen lemma。[^5]

**碰撞等级：A（直接包含）。** 即使没有找到完全相同的网络术语和 \([CSF-Q,-CS]\) 公式，审稿人只需做 (A) 的三项定义，即可从 Proposition 4.16 复现主 LMI。这不足以支持“新的 exact robust SDP theorem”。

### 2.2 Matrix S-lemma / full-block S-procedure

van Waarde–Camlibel–Mesbahi 2022 已系统给出 matrix S-lemma 的必要充分形式，并用于把矩阵不确定集上的二次要求化成 LMI。[^6] van Waarde et al. 2023 又把 matrix S-lemma、full-block S-procedure 与 Petersen lemma 的关系明确接通。El Ghaoui–Oustry–Lebret 1998 的 uncertain SDP 工作也早已区分一般结构化不确定性下的充分 relaxation 与 full uncertainty 下可成为必要充分条件的情形。[^7]

因此，本问题“只有一个 quadratic uncertainty，所以 scalar multiplier lossless”的证明路线也不是一个独立的新机制。它是 classical S-lemma 的正确应用；从 full-block 角度看，则是 Petersen lemma 的正确应用。

### 2.3 Robust approximate inverse / model matching

El Ghaoui 2002 的 Eq. (1.4) 已把“选择一个共同矩阵，使一族不确定逆的最大 operator-norm 误差最小”作为 common approximate inverse 问题；其 Theorem 6.2（p. 184）通过 LFR、elimination 与 SDP multiplier 处理 structured inversion error，并对 unstructured full perturbation 给出 exact 情形。[^8]

本问题不与该文的具体 Eq. (6.5) 完全相同：这里仅观察 \(CS\)，被近似对象是耦合行块 \([CSF,-CS]\)，中心还被限制成 \([Q,0]\)。然而这些差异只说明**应用实例不同**，不能恢复定理 4 的基础 novelty，因为 (B) 已经把本实例送进 Petersen lemma。

用 robust model matching 的语言，(B) 是

\[
\min_{Q\in\mathcal Q_{\rm loc}}
\sup_{\|\Delta\|\le1}\|B_Q+G\Delta H\|,
\]

其中设计变量只进入 nominal term。单 full block 的 exact LMI 与 affine design constraints 的组合是标准 robust convex synthesis 操作。

### 2.4 2025 QMI Chebyshev center 结果

Shakouri–van Waarde–Camlibel 2025 研究由 QMI 诱导的矩阵集合在 unitarily invariant norms 下的 Chebyshev center/radius；其 Proposition 2(c) 给出 matrix-ball 参数化，Theorem 7 给出紧 QMI 集的 unconstrained center/radius 闭式。[^9]

当前问题可看成矩阵族

\[
\mathcal A=\{[CSF,-CS]:\alpha I\preceq S\preceq\beta I\}
\]

到受限 affine center set \(\{[Q,0]:Q\in\mathcal Q_{\rm loc}\}\) 的距离。受限 center 并不是该文 Theorem 7 的原样结论，这是一个真实区别；但这个区别已经由 (5) 的 affine \(Q\) 联合优化直接解决，故尚不能把它包装成新的强定理。

## 3. 任意 self-adjoint matrix interval：像公式正确，但也是经典特例

考虑增强后的不确定集

\[
S_-\preceq S\preceq S_+,qquad
H_0=\frac{S_++S_-}{2},\qquad
\Delta_0=\frac{S_+-S_-}{2}\succ0.
\]

令 \(X=S-H_0\)。则

\[
-\Delta_0\preceq X\preceq\Delta_0
\iff
X=\Delta_0^{1/2}K\Delta_0^{1/2},
\quad K=K^T,quad -I\preceq K\preceq I.
\tag{7}
\]

对固定 \(v\)，self-adjoint contraction 的单向量像是整个球，因此

\[
\boxed{
\{Sv:S_-\preceq S\preceq S_+\}
=\left\{z:
(z-H_0v)^T\Delta_0^{-1}(z-H_0v)
\le v^T\Delta_0v
\right\}.}
\tag{8}
\]

所以父任务提出的椭球公式是正确的。

但 (8) 不是可主张原创的 matrix-interval theorem。去掉 self-adjoint 记号后，(7) 是标准 contraction/matrix-ball 参数化；van Waarde et al. 2023, Theorem 3.3 给出一般 QMI matrix-ball 参数化，而 Theorem 3.4 给出其在线性映射下的 exact image。具体取 QMI 变量 \(X\) 以及

\[
\Pi=
\begin{bmatrix}
\Delta_0&0\\0&-\Delta_0^{-1}
\end{bmatrix},qquad
\begin{bmatrix}I\\X\end{bmatrix}^{T}
\Pi
\begin{bmatrix}I\\X\end{bmatrix}
=\Delta_0-X^T\Delta_0^{-1}X\succeq0.
\tag{9}
\]

右乘 \(v\) 后，Theorem 3.4（因 \(-\Delta_0^{-1}\) 可逆）给出像 QMI

\[
v^T\Delta_0v-y^T\Delta_0^{-1}y\ge0,qquad y=Xv.
\]

该 QMI 暂时允许一般 contraction；self-adjoint 限制不会缩小单向量像，可由一个二维 self-adjoint contraction/Householder completion 实现。Davis–Kahan–Weinberger 1982 和 Parrott 1978 的 operator/block-contraction completion 理论提供更一般的经典背景。[^10][^11]

在 (8) 之后，robust approximation 又可写成单 full block。令

\[
B_Q=
\begin{bmatrix}
F^TH_0C^T-Q^T\\H_0C^T
\end{bmatrix},\quad
G=
\begin{bmatrix}
F^T\Delta_0^{1/2}\\\Delta_0^{1/2}
\end{bmatrix},\quad
H=\Delta_0^{1/2}C^T,
\]

即有同样的 \(\sup_{\Delta^T\Delta\preceq I}\|B_Q+G\Delta H\|\) 表示和 Petersen LMI。

**结论：** 从标量谱窗推广到一般 self-adjoint interval 扩大了应用 scope，却不改变理论难度等级；“matrix interval image + lossless SDP”两步均已有现成结果。可以作为实用 generalization/corollary，但不能单独挽救主定理 novelty。

## 4. DKW、Parrott 与 operator-valued Nevanlinna–Pick 的位置

Davis–Kahan–Weinberger 1982 的主 completion theorem（pp. 447 起）刻画所有能把给定块补成范数不超过 \(\mu\) 的未知块；Parrott 1978 则从 quotient norm 与 Sz.-Nagy–Foiaş lifting 的角度给出相关 operator block completion。[^10][^11] 这些结果比“指定一个 contraction 在单向量上的作用”更一般，因此能解释引理 3/式 (8) 的 completion 几何。

但它们不是本问题 outer minimax LMI 的最直接来源：completion theorem 负责回答“某个向量像能否由 contraction 实现”，Petersen lemma 负责回答“该 full block 下的 robust quadratic inequality 是否成立”。两者组合仍然是经典 machinery。

matrix/operator-valued Nevanlinna–Pick、commutant lifting 和 strong Parrott theorems 处理的是带解析函数、算子代数或多点插值约束的动态问题；例如 Delsarte–Genin–Kamp 1979 系统处理 matrix-valued Pick interpolation，Foiaş–Tannenbaum 1989 则把 strong Parrott completion 与 commutant lifting/插值联系起来。[^13][^14] 本问题只有一个静态有限维 full block，在进入那些更重的理论以前就已被 Petersen lemma 完全消元。除非后续模型真正加入频率变量、因果性、插值节点或 commutation constraints，否则引用 operator-valued Nevanlinna–Pick 不会增加原创性，只会使叙述显得过重。

## 5. Robust least squares 的相邻性

El Ghaoui–Lebret 1997 对 norm-bounded uncertain least squares 给出 exact robust counterpart；在 unstructured data uncertainty 下，minimax residual 可化为显式 convex conic problem。[^12] 它与当前工作共享以下思想：

- 一个共同决策必须同时适应一族 norm-bounded 数据；
- 目标是 worst-case Euclidean/operator norm；
- full unstructured uncertainty 使 robust counterpart 可 exact convexify。

它不直接给出 `(SDP-cont)`，因为当前目标是 uniform operator mapping，并且 uncertainty 在两块中耦合出现。不过从 novelty 评判看，这一文献再次说明“full norm ball + worst-case norm + convex design”本身是成熟范式。真正可发表的贡献必须来自网络结构、耦合多块、闭式分类或 sharp lower bound，而不能只来自 robustification。

## 6. Locality-constrained Schur structure 到底新增了什么

本问题最有应用识别度的结构是：

1. \(S\) 是外部网络消元后的 inverse Schur/completion block；
2. \(C\) 只选择边界可见行；
3. \([CSF,-CS]\) 由同一个 \(S\) 耦合；
4. 实现器被迫是 \([Q,0]\)，而且 \(Q\) 可带通信半径、稀疏或局部 opcode 限制。

这些结构足以形成一篇有意义的网络估计论文，但在定理 4 中，它们只决定 Petersen 模板里的常数 \(B_Q,G,H\) 和 affine feasible set \(\mathcal Q_{\rm loc}\)。对于每个 admissible \(Q\)，uncertainty elimination 完全相同；再联合最小化 \(Q\) 仍是一个 SDP。

因此要区分两种 claim：

- **可以主张（经更完整应用查重后）：** 首次把某类局部 state-estimation/Schur-boundary approximation 精确重写为这一 robust model-matching SDP；首次给出它与通信轮次/能耗/边界宽度的可计算联系；首次在相应网络数据上验证规模规律。
- **不宜主张：** 新的 lossless S-procedure、新的 exact full-block robust-performance theorem、首次把 self-adjoint interval 的单向量像写成球/椭球，或首次得到 norm-bounded robust approximation 的有限 LMI。

换言之，locality 目前贡献的是**问题建模和决策集**，不是新的 robust-elimination 数学。

## 7. 关于 \(\operatorname{span}(\operatorname{range}C^T,\operatorname{range}F)+1\) 压缩

对标量谱窗 \(hI\pm dI\)，确有一个 exact 的低维 representer reduction。令

\[
\mathcal U=operatorname{span}(\operatorname{range}C^T+\operatorname{range}F),
\qquad s=\dim\mathcal U.
\tag{10}
\]

固定 \(u\)，写

\[
z=hC^Tu+dy,qquad \|y\|\le\|C^Tu\|,qquad
y=a+b, a\in\mathcal U, b\perp\mathcal U.
\]

因为 \(C^Tu\in\mathcal U\)、\(\operatorname{range}F\subseteq\mathcal U\)，有

\[
F^Tb=0,qquad
\|z\|^2=\|hC^Tu+da\|^2+d^2\|b\|^2.
\]

所以目标只通过标量 \(\|b\|\) 看见整个 \(\mathcal U^\perp\)。给定 \(a\) 后，worst case 会使用剩余半径，而方向无关；可把 \(b\) 精确替换成某个额外坐标 \(\rho e_\perp\)。于是 uncertainty/output 空间可从 \(n\) 维压到

\[
\boxed{s+\mathbf 1_{\{s<n\}}}
\]

维，而不是总写成 \(s+1\)。这会相应缩小 S-lemma/Petersen LMI。证明只是正交分解和旋转对称性，不需要 relaxation。

更具体地，取 \(U\in\mathbb R^{n\times s}\) 为 \(\mathcal U\) 的正交基，置

\[
\widehat C^T=U^TC^T,\qquad \widehat F=U^TF.
\]

若 \(s<n\)，定义

\[
\bar C^T=\begin{bmatrix}\widehat C^T\\0_{1\times p}\end{bmatrix},
\qquad
\bar F=\begin{bmatrix}\widehat F\\0_{1\times k}\end{bmatrix}.
\tag{11}
\]

则原 S-lemma 中的 \(z\in\mathbb R^n\) 可精确替成
\(\bar z=[z_{\mathcal U};z_\perp]\in\mathbb R^{s+1}\)，并且

\[
\|\bar z-h\bar C^Tu\|^2\le d^2\|\bar C^Tu\|^2,qquad
\left\|\begin{bmatrix}
\bar F^T\bar z-Q^Tu\\\bar z
\end{bmatrix}\right\|^2
\]

分别与原问题的可行像和 worst-case objective 完全相同。故可直接在 `(SDP-cont)` 中以
\((\bar C,\bar F,s+1)\) 替代 \((C,F,n)\)。若 \(s=n\)，则直接用
\((U^TC^T,U^TF,n)\)，无需添加虚拟坐标。

这条压缩有计算价值，但应谨慎评价其数学 novelty：

- 它是标准“目标只看一个子空间、球在正交补上旋转不变”论证的直接结果；目前未找到完全相同 \((C,F)\) 符号下的 primary theorem，不等于可以声称高原创性。
- 单独作为主定理强度仍不够，更适合作为 complexity proposition。
- 对一般 matrix interval \(\Delta_0\) 不能无条件沿用。椭球度量破坏了 \(\mathcal U^\perp\) 的旋转对称性；除非相关子空间对 \(H_0\)、\(\Delta_0^{1/2}\)（或相应 generalized metric）不变，否则一个额外标量不再能代表整个正交补。

一个有希望升级的方向是：给出**一般 anisotropic interval 下的最小 exact compression space**，证明其维数必要且充分，并给出可由 block Krylov/invariant closure 计算的形式；若还能把维数绑定到 graph cut 而非 exterior size，这会比标量窗的 \(s+1\) 观察强得多。

## 8. 如果要形成强主定理，应怎样增强

因为 single-full-block exactness 已撞车，简单扩大 \(C,F\) 尺寸、添加 affine locality constraints、或把 \(\alpha I\preceq S\preceq\beta I\) 换成 \(S_-\preceq S\preceq S_+\) 都不够。以下方向才可能把工作提升为主定理。

### 8.1 首选：结构化多块仍 lossless 的图论条件

实际网络往往产生多个耦合 Schur completions、多个 cut 或多个局部 uncertainty blocks。一般多个 quadratic/full blocks 下，标量 multiplier S-procedure 不再 lossless。可尝试证明：对由 chordal/tree decomposition、running-intersection property 或某种 separator nesting 诱导的 blocks，某类 sparse multiplier 是必要充分的。

强主定理应包含：

- 明确的图条件；
- exactness 的 iff 证明，而不只是 sufficient SDP；
- 违反条件时的反例或保守 gap；
- LMI 尺寸随 separator/treewidth 而不是全局节点数增长。

这是真正超出 single-full-block Petersen lemma 的地方，也是最贴合“局部通信”的数学升级。

### 8.2 Constrained Chebyshev center 的闭式/分类定理

不要只给 SDP；刻画何时 \(Q^*\) 有闭式、何时唯一、何时简单 Dirichlet/local inverse 已最优，以及最坏 completion 的秩和方向。一个合格形式可以是：

\[
Q^*=\text{某个 GSVD/canonical-angle 公式}
\quad\Longleftrightarrow\quad
\text{关于 }(C,F,\mathcal Q_{\rm loc})\text{ 的可检验条件}.
\]

若能给 necessary-and-sufficient classification，而不是若干 sufficient special cases，就可能成为主定理。

### 8.3 网络拼接定理与匹配 converse

把单 cut 的 robust radius 嵌入全网估计：证明给定局部失败概率、通信半径和 completion model 时，全局误差/成功概率如何由各 separator radius 精确或近精确地拼接，并给出任何局部算法都必须付出的匹配下界。仅有上界不够；upper/lower bound 同阶甚至常数 sharp，才有“终结问题”的力度。

### 8.4 非仿射实现约束下的 sharp result

若局部 opcode 还必须满足因果性、量化、正性、稳定性或有限字长，\(Q\) 的可行集不再只是 affine subspace。可以追求 exact reformulation、可证明 approximation ratio，或复杂性二分定理。这样的结果不再被 Petersen lemma 一行解决。

### 8.5 Minimal compression theorem

把第 7 节提升为一般 theorem：对 anisotropic interval、多个 blocks 和 graph-local observation，找出最小 sufficient subspace，并证明任何 exact algorithm 至少需要该维数。若维数由 local cut rank/treewidth 控制，同时给出构造算法和 tight lower bound，这可成为主线贡献，而不只是 SDP 尺寸优化。

## 9. 保守 novelty 置信度与 go/no-go

| 判断 | 置信度 | 说明 |
|---|---:|---|
| `(SDP-cont)` 被现成 Petersen/full-block theorem 直接包含 | **0.97** | 已找到 primary theorem、非严格条件和逐变量映射 |
| 一般 self-adjoint interval 的椭球像不是新定理 | **0.95** | QMI image + contraction completion 直接覆盖 |
| locality affine constraints 本身不产生新的 exactness theorem | **0.98** | 只改变 \(B_Q\) 的 affine design set |
| 文献中未出现完全相同的网络符号和应用叙述 | **0.55** | 不能靠关键词阴性检索确认，且对核心数学 novelty 无决定性意义 |
| 标量窗的 \(s+1\) 压缩可作为未发表新命题 | **0.35** | 精确但很可能属于常见 symmetry/representer reduction；需单独查复杂性与 robust optimization 文献 |

**Go/no-go：**

- 把定理 4 作为工具性 proposition/corollary：**GO**。
- 把定理 4 或一般 matrix interval 版本作为强主定理并宣称新的 exact LMI：**NO-GO**。
- 以“structured multi-cut losslessness / sharp graph composition / constrained center 闭式分类 / minimal anisotropic compression”之一重建主定理：**GO，值得继续。**

## 10. 检索边界与仍需查的库

本审计优先查了能直接判定公式包含关系的 primary sources：SIAM、IEEE/Automatica、Springer、ScienceDirect 上的 QMI/matrix S-lemma、Petersen lemma、full-block S-procedure、approximate inverse、robust least squares、DKW/Parrott completion 与 2025 QMI Chebyshev center。核心 no-go 不依赖关键词阴性结果，而依赖 Proposition 4.16 的正面逐式包含。

若未来要对“完全相同网络应用是否首见”或第 7 节 compression 的 novelty 作投稿级声明，仍应补：

1. MathSciNet 与 zbMATH 的 cited-by/MSC 交叉检索；
2. IEEE Xplore 与 IFAC PapersOnLine 的 distributed estimation、localized control、network inverse/Schur complement 全文检索；
3. ProQuest theses、HAL、DART-Europe 与博士论文中的 unpublished variants；
4. operator theory 专著中 Parrott/Douglas lemma 的 self-adjoint completion 变体；
5. robust optimization 与 numerical linear algebra 中 symmetry reduction、effective uncertainty dimension、facial reduction/low-rank robust SDP 文献。

这些补查可能影响“应用首见”和“compression 首见”，但不会改变本审计的核心结论：**单 full-block exact SDP 已被现成 theorem 直接覆盖。**

## Primary sources

[^1]: H. J. van Waarde, M. Kanat Camlibel, J. Eising, and H. L. Trentelman, “[Quadratic Matrix Inequalities with Applications to Data-Based Control](https://doi.org/10.1137/22M1486807),” *SIAM Journal on Control and Optimization* 61(4) (2023), 2251–2281. See Theorem 3.4, Proposition 4.14, and especially Proposition 4.16(b). [Official SIAM page](https://epubs.siam.org/doi/10.1137/22M1486807).

[^2]: C. W. Scherer, “[LPV Control and Full Block Multipliers](https://doi.org/10.1016/S0005-1098(00)00176-X),” *Automatica* 37(3) (2001), 361–375. See p. 367 and Theorem 5.3 for exactness of the relevant single-full-block relaxation.

[^3]: C. W. Scherer, “[A Full Block S-Procedure with Applications](https://doi.org/10.1109/CDC.1997.657686),” in *Proceedings of the 36th IEEE Conference on Decision and Control* (1997), 2602–2607.

[^4]: Ian R. Petersen, “[A Stabilization Algorithm for a Class of Uncertain Linear Systems](https://doi.org/10.1016/0167-6911(87)90102-2),” *Systems & Control Letters* 8(4) (1987), 351–357.

[^5]: Andrea Bisoffi, Claudio De Persis, and Pietro Tesi, “[Data-Driven Control via Petersen’s Lemma](https://doi.org/10.1016/j.automatica.2022.110537),” *Automatica* 145 (2022), 110537. See Facts 1–2 on pp. 3–4 for strict and nonstrict forms.

[^6]: Henk J. van Waarde, M. Kanat Camlibel, and Mehran Mesbahi, “[From Noisy Data to Feedback Controllers: Nonconservative Design via a Matrix S-Lemma](https://doi.org/10.1109/TAC.2020.3047577),” *IEEE Transactions on Automatic Control* 67(1) (2022), 162–175. See the matrix S-lemma results, in particular Theorems 9 and 13 in the accepted manuscript numbering.

[^7]: Laurent El Ghaoui, François Oustry, and Hervé Lebret, “[Robust Solutions to Uncertain Semidefinite Programs](https://doi.org/10.1137/S1052623496305717),” *SIAM Journal on Optimization* 9(1) (1998), 33–52.

[^8]: Laurent El Ghaoui, “[Inversion Error, Condition Number, and Approximate Inverses of Uncertain Matrices](https://doi.org/10.1016/S0024-3795(01)00273-7),” *Linear Algebra and its Applications* 343–344 (2002), 171–193. See Eq. (1.4), Lemma 6.1, and Theorem 6.2 on p. 184. [Author’s primary manuscript](https://people.eecs.berkeley.edu/~elghaoui/Pubs/InvErr_LAA02.pdf).

[^9]: Amir Shakouri, Henk van Waarde, and M. Kanat Camlibel, “[Chebyshev Centers and Radii for Sets Induced by Quadratic Matrix Inequalities](https://doi.org/10.1007/s00498-025-00424-w),” *Mathematics of Control, Signals, and Systems* (2025). See Proposition 2(c) and Theorem 7. [Official full text](https://link.springer.com/article/10.1007/s00498-025-00424-w).

[^10]: Chandler Davis, W. M. Kahan, and H. F. Weinberger, “[Norm-Preserving Dilations and Their Applications to Optimal Error Bounds](https://doi.org/10.1137/0719029),” *SIAM Journal on Numerical Analysis* 19(3) (1982), 445–469. See the main block-completion theorem beginning around p. 447.

[^11]: Stephen Parrott, “[On a Quotient Norm and the Sz.-Nagy–Foiaş Lifting Theorem](https://doi.org/10.1016/0022-1236(78)90060-5),” *Journal of Functional Analysis* 30(3) (1978), 311–328.

[^12]: Laurent El Ghaoui and Hervé Lebret, “[Robust Solutions to Least-Squares Problems with Uncertain Data](https://doi.org/10.1137/S0895479896298130),” *SIAM Journal on Matrix Analysis and Applications* 18(4) (1997), 1035–1064. [Author’s primary manuscript](https://people.eecs.berkeley.edu/~elghaoui/Pubs/rob-ls.pdf).

[^13]: Philippe Delsarte, Yves Genin, and Yves Kamp, “[The Nevanlinna–Pick Problem for Matrix-Valued Functions](https://doi.org/10.1137/0136005),” *SIAM Journal on Applied Mathematics* 36(1) (1979), 47–61.

[^14]: Ciprian Foiaş and Allen Tannenbaum, “[A Strong Parrott Theorem](https://doi.org/10.1090/S0002-9939-1989-0972228-9),” *Proceedings of the American Mathematical Society* 106(3) (1989), 777–784. See Theorems 1–2 and the introduction’s connection to commutant lifting and interpolation.
