# 独立权重低秩信息更新的逆矩阵顶点结论：证明、边界与原创性审计

## 结论先行

考虑

\[
J(w)=H+\sum_{e=1}^{m}w_ea_ea_e^T,\qquad
T(w)=H^{1/2}J(w)^{-1}H^{1/2},\qquad H\succ0.
\]

这条路线得到四个可以严格闭合的结论，但它们的原创性强弱差别很大。

1. **有限盒上的顶点结论成立，而且有矩阵值凸包恒等式。** 对
   \(w\in\prod_e[\ell_e,u_e]\)，每个 \(T(w)\) 都是 \(2^m\) 个角点矩阵
   \(T(v)\) 的显式凸组合。因此任意连续凸函数 \(\Phi\) 的最坏值都在角点取得。
   但这一部分是 rational-multiaffine 顶点方法的直接矩阵值强化，不应作为原创主定理：
   Barmish 的 multiaffine 顶点引理和 Blanchini--Colaneri--Giordano--Zorzan
   2022 的 totally-multiaffine 灵敏度顶点定理已经非常接近。[^2][^3]
2. **非负无界权重的闭凸包可以完全写出。** 令
   \(B=H^{-1/2}[a_1\ \cdots\ a_m]\)，则
   \[
   \overline{\operatorname{conv}}\{T(w):w\ge0\}
   =\operatorname{conv}\{P_{\ker B_S^T}:S\subseteq[m]\}.
   \]
   每个互不相同的投影都是这个多面体的 exposed vertex。由此任意连续凸损失在
   \(w\ge0\) 上的 supremum 精确等于这些零空间投影上的最大值。
3. **决定“端点是否足够”的本质不是字面上的 rank one。** 对单参数
   \((A+tU)^{-1}\)，所有连续凸损失都只需检查两端点，当且仅当
   \(A^{-1/2}UA^{-1/2}=cP\)，其中 \(P\) 是正交投影。rank one 自动满足，
   但较高秩的 scaled projection 也满足。一般 rank-two 共享权重会失败。
4. **对 common local Schur/WLS 规则，这给出精确有限场景 SDP。** 它能严格避免
   把真实结构粗暴放大成 \(0\preceq T\preceq I\) 的保守性；一个二维 WLS 例子中，
   这种 full-block 外逼近的保守比可以随参数 \(M\) 无界增长。

截至 2026-09-14 的审计结论是：有限盒顶点公式、DPP 均值投影、图随机森林平均均有
明确先例；“无界闭凸包 + 全部不同投影均 exposed + common Schur 规则精确 LMI +
full-block 无界保守比”这一整套组合尚未检索到同样陈述，但各证明都较短，不能仅凭
阴性检索宣称首次发现。它适合作为更大 completion/local-WLS 主定理的核心工具；若单独
投稿应用数学论文，当前数学厚度仍偏薄。

## 1. 模型、符号与量词

记 \(\mathbb S^n\) 为实对称矩阵空间，内积为
\(\langle X,Y\rangle=\operatorname{tr}(XY)\)。有限盒部分假设

\[
\mathcal W=\prod_{e=1}^{m}[\ell_e,u_e],\qquad
0\le \ell_e<u_e<\infty.
\]

更一般地，证明只需 \(J(w)\) 在整个盒上正定；非负权重只是保证这一点的自然 WLS
情形。零宽度坐标可直接删除。对 \(\sigma\in\{0,1\}^m\)，定义角点

\[
v_e^\sigma=\ell_e+\sigma_e(u_e-\ell_e),
\quad
x_e=\frac{w_e-\ell_e}{u_e-\ell_e},
\quad
\lambda_\sigma(x)=\prod_{e=1}^{m}x_e^{\sigma_e}(1-x_e)^{1-\sigma_e}.
\]

无界部分固定下界为零。若原问题有 \(\ell_e>0\)，应先把
\(\sum_e\ell_ea_ea_e^T\) 吸收到新基准
\(H_\ell=H+\sum_e\ell_ea_ea_e^T\) 中，再用 \(H_\ell\) 白化；若仍用原来的
\(H\) 白化，极限一般是投影的合同变换，而不是正交投影本身。

## 2. 有限盒：精确角点凸组合

### 定理 2.1（有限盒的矩阵值凸包恒等式）

对每个 \(w\in\mathcal W\)，有

\[
T(w)=\sum_{\sigma\in\{0,1\}^m}\pi_\sigma(w)T(v^\sigma),
\tag{2.1}
\]

其中

\[
\pi_\sigma(w)=
\frac{\lambda_\sigma(x)\det J(v^\sigma)}{\det J(w)},
\qquad
\pi_\sigma(w)\ge0,
\qquad
\sum_\sigma\pi_\sigma(w)=1.
\tag{2.2}
\]

因而

\[
\operatorname{conv}\{T(w):w\in\mathcal W\}
=\operatorname{conv}\{T(v^\sigma):\sigma\in\{0,1\}^m\}.
\tag{2.3}
\]

特别地，对任意连续凸函数 \(\Phi:\mathbb S^n\to\mathbb R\)，

\[
\max_{w\in\mathcal W}\Phi(T(w))
=\max_{\sigma\in\{0,1\}^m}\Phi(T(v^\sigma)).
\tag{2.4}
\]

#### 证明

固定除 \(w_e\) 外的所有参数。因为第 \(e\) 个系数矩阵是 rank one，\(J(w)\)
的每个方形子矩阵关于 \(w_e\) 都是“固定矩阵 + rank-at-most-one 矩阵乘
\(w_e\)”；其行列式关于 \(w_e\) 是 affine 的。因此 \(\det J(w)\) 和
\(\operatorname{adj}J(w)\) 的每个元素都分别关于每个 \(w_e\) affine，即都是
multiaffine 多项式。

标准张量积 Bernstein 插值于是给出

\[
\det J(w)=\sum_\sigma\lambda_\sigma(x)\det J(v^\sigma),
\qquad
\operatorname{adj}J(w)=
\sum_\sigma\lambda_\sigma(x)\operatorname{adj}J(v^\sigma).
\]

用 \(J^{-1}=\operatorname{adj}J/\det J\)，并在每一项中乘除
\(\det J(v^\sigma)\)，得到 (2.1)--(2.2)。整个盒上 \(J\succ0\)，所以角点和
内部的行列式均为正；权重非负。权重和为一则由行列式插值式直接得到。固定合同变换
\(X\mapsto H^{1/2}XH^{1/2}\) 保留同一组凸组合权重。

(2.3) 的一个方向来自 (2.1)，另一个方向来自角点本身属于原像。
最后由凸性

\[
\Phi(T(w))\le\sum_\sigma\pi_\sigma(w)\Phi(T(v^\sigma))
\le\max_\sigma\Phi(T(v^\sigma)),
\]

而角点可行，故 (2.4) 成立。证毕。

### Sherman--Morrison 的逐坐标解释

上述显式全局权重也可由逐坐标线段性质得到。固定其他参数，令
\(R_\ell=(M+\ell aa^T)^{-1}\)、\(\delta=t-\ell\)、
\(\Delta=u-\ell\)、\(c=a^TR_\ell a\)。Sherman--Morrison 公式[^1]给出

\[
R(t)=R_\ell-\frac{\delta}{1+\delta c}R_\ell aa^TR_\ell
=(1-\theta(t))R_\ell+\theta(t)R_u,
\tag{2.5}
\]

其中

\[
\theta(t)=\frac{\delta(1+\Delta c)}{\Delta(1+\delta c)}\in[0,1].
\tag{2.6}
\]

依次消去每个内部坐标，就得到角点凸组合。这个证明直观地说明了“每个不确定标量只
控制一个 rank-one 信息因子”的作用。

### 已知结果碰撞

这一节不能作为论文的原创核。

- Barmish 1994 的 Lemma 14.5.5 已给出 hyperrectangle 上 multiaffine 标量函数的
  极值顶点性质。[^2]
- Blanchini et al. 2022 的 Definition 1 把“所有 minors 都 multiaffine”定义为
  totally multiaffine；其 BDC 分解明确包含由独立标量乘 rank-one 矩阵的情形。
  Theorem 2 在 robust nonsingularity 下证明标量
  \(-H_\mathrm{out}J(\delta)^{-1}E_\mathrm{in}\) 的紧上下界由参数盒角点取得。
  本地归档 PDF 的 Definition 1 在 PDF p.6，Theorem 2 在 PDF p.12。[^3]
- Alamo--Tempo--Rodríguez--Camacho 2008 也系统研究 interval/multiaffine robust
  constraints 的顶点归约和 LMI 应用。[^4]

(2.1) 比单个标量灵敏度界更强：它给所有矩阵元素同一组凸组合权重，从而一次覆盖任意
矩阵凸损失。不过其证明正是 adjugate/determinant 的 rational-multiaffine 机制，属于
自然强化，而不是足以单独支撑主定理的全新原理。

## 3. 单参数的最宽端点条件

rank one 是一个简单而稳健的充分条件，但并非点态必要条件。下面的 iff 给出单参数
情形真正的边界。

### 定理 3.1（单参数逆曲线的端点 iff）

令 \(A\succ0\)、\(U\succeq0\)、\(U\ne0\)，并令 \(0<\tau<\infty\)。下列条件等价：

1. 对每个 \(t\in[0,\tau]\)，
   \((A+tU)^{-1}\in\operatorname{conv}\{A^{-1},(A+\tau U)^{-1}\}\)；
2. 对每个连续凸函数 \(\Phi:\mathbb S^n\to\mathbb R\)，
   \[
   \max_{0\le t\le\tau}\Phi((A+tU)^{-1})
   =\max\{\Phi(A^{-1}),\Phi((A+\tau U)^{-1})\};
   \]
3. \(X=A^{-1/2}UA^{-1/2}\) 的所有非零特征值都等于同一个 \(c>0\)，即
   \(X=cP\)，其中 \(P\) 是正交投影；
4. \(UA^{-1}U=cU\) 对某个 \(c>0\) 成立。

#### 证明

固定合同变换是线性双射，所以只需考察
\(Y(t)=(I+tX)^{-1}\)。设 \(X=V\operatorname{diag}(\lambda_j)V^T\)，
\(\lambda_j\ge0\)。若某个内点满足

\[
Y(t)=(1-\theta)I+\theta Y(\tau),
\]

则对每个 \(\lambda_j>0\)，

\[
\theta=\frac{t(1+\tau\lambda_j)}{\tau(1+t\lambda_j)}.
\tag{3.1}
\]

当 \(0<t<\tau\) 时，(3.1) 关于 \(\lambda_j\) 严格递增，所以同一个 \(\theta\)
能适用于所有正特征值，当且仅当这些正特征值全相等。反过来，若 \(X=cP\)，则

\[
(I+tX)^{-1}=I-\frac{tc}{1+tc}P

\]

显然位于两个端点之间。这证明 1 与 3 等价。条件 3 等价于 \(X^2=cX\)，合同变换后
恰为 \(UA^{-1}U=cU\)，所以 3 与 4 等价。

1 推出 2 是凸性的直接结果。若 1 不成立，则某个逆矩阵点落在两个端点的闭线段之外；
有限维空间中的严格分离定理给出一个线性泛函，在该内部点上的值严格大于两个端点，
而线性泛函本身就是连续凸函数，故 2 失败。证毕。

若参数区间是 \([\ell,u]\)，把 \(A+\ell U\) 当作新基准即可。rank-one
\(U=aa^T\) 总满足

\[
UA^{-1}U=(a^TA^{-1}a)U.
\]

但 \(U\) 也可以是任意秩的“白化后 scaled projection”。因此严谨措辞应是：

> 独立 rank-one 因子在不需要额外代数结构时保证顶点结论；一般 shared rank-two
> 因子不保证它。rank one 不是单参数逐点必要条件。

### 多参数的可检查充分条件

更一般地，令

\[
J(w)=H+\sum_e w_eU_e,\qquad U_e\succeq0.
\]

固定 \(w_{-e}\)，记 \(A_{-e}(w)=H+\sum_{f\ne e}w_fU_f\)。若对每个 \(e\) 和
每个允许的 \(w_{-e}\)，都存在 \(c_e(w_{-e})\ge0\) 使

\[
U_eA_{-e}(w)^{-1}U_e=c_e(w_{-e})U_e,
\tag{3.2}
\]

则每条坐标逆曲线都位于自己的端点线段中，逐坐标递归仍推出有限盒角点凸包结论。
(3.2) 是充分条件；多参数整体顶点性质可能因其他角点的凸包几何而成立，所以这里不把
它误写成多参数全局必要条件。

### 一般 shared rank two 的显式反例

令

\[
J(t)=I+t\operatorname{diag}(1,4),\qquad 0\le t\le1,
\]

并取线性（因而凸）函数 \(\Phi(X)=X_{11}-X_{22}\)。则

\[
\Phi(J(0)^{-1})=0,\qquad
\Phi(J(1)^{-1})=\frac3{10},\qquad
\Phi(J(1/2)^{-1})=\frac13>\frac3{10}.
\]

因此一般 rank-two 共享标量权重时，内部点可以严格超过两个端点。失败原因正是白化后
两个正特征值 \(1,4\) 不相同。

定向检索没有发现定理 3.1 以完全相同的“所有凸函数端点 iff”形式发表；但它是短的
谱分解与分离定理推论。最稳妥的论文定位是 sharp boundary lemma，原创性仍标记
**待核实**，不能仅因没有搜到就声称首创。

## 4. 无界非负权重：零空间投影的精确闭凸包

令

\[
A_0=[a_1\ \cdots\ a_m],\qquad B=H^{-1/2}A_0,
\qquad W=\operatorname{diag}(w),
\]

则

\[
T(w)=(I+BWB^T)^{-1}.
\tag{4.1}
\]

对 \(S\subseteq[m]\)，定义

\[
P_S=P_{\ker B_S^T}
=I-B_S(B_S^TB_S)^\dagger B_S^T.
\tag{4.2}
\]

### 定理 4.1（无界族的精确闭凸包与 exposed 场景）

有

\[
\overline{\operatorname{conv}}\{T(w):w\in\mathbb R_+^m\}
=\operatorname{conv}\{P_S:S\subseteq[m]\}.
\tag{4.3}
\]

对任意连续凸函数 \(\Phi\)，

\[
\sup_{w\ge0}\Phi(T(w))=\max_{S\subseteq[m]}\Phi(P_S).
\tag{4.4}
\]

此外，每个互不相同的 \(P_S\) 都是右侧多面体的 exposed vertex。

#### 确定性证明

先固定任意有限 \(w\ge0\)，取 \(U\ge\max_e w_e\)。定理 2.1 应用于盒
\([0,U]^m\)，说明 \(T(w)\) 属于 \(2^m\) 个矩阵
\(T(U1_S)\) 的凸包。对固定 \(S\)，谱分解或 Woodbury 公式给出

\[
\lim_{U\to\infty}(I+UB_SB_S^T)^{-1}=P_{\ker B_S^T}=P_S.
\tag{4.5}
\]

把 \(U\to\infty\)。相应凸组合系数位于紧的有限维概率单纯形中，故可取一个收敛
子列；(4.5) 表明子列极限把固定的 \(T(w)\) 写成 \(P_S\) 的凸组合。因此所有
\(T(w)\) 都属于 \(\operatorname{conv}\{P_S\}\)。反过来，(4.5) 表明每个
\(P_S\) 都是原族的极限，故得到 (4.3)。由凸性和连续性立即得到 (4.4)。注意若某些
子集张成相同列空间，它们给出同一个投影，应去重。

最后固定一个投影 \(P\)，用线性泛函

\[
L_P(X)=\operatorname{tr}((2P-I)X).
\]

对任意另一个正交投影 \(Q\)，

\[
L_P(P)-L_P(Q)
=\operatorname{tr}(P)+\operatorname{tr}(Q)-2\operatorname{tr}(PQ)
=\|P-Q\|_F^2.
\tag{4.6}
\]

若 \(P\ne Q\)，差严格为正，所以 \(L_P\) 唯一暴露 \(P\)。证毕。

式 (4.6) 有一个重要的负面含义：对“任意凸损失”这一宽量词，不存在一个统一的更小
场景表可以删掉某个不同投影而仍保持精确，因为总有线性损失把它单独选为最坏点。
实际问题可以利用特定 \(C,F,Q\) 的结构删场景，但那不再是目标无关结论。

## 5. 行列式抽样 / DPP 均值恒等式

定理 4.1 还有一个精确概率表示。定义

\[
\mathbb P_w(S)=
\frac{\det(B_S^TB_S)\prod_{e\in S}w_e}
{\det(I+BWB^T)},
\tag{5.1}
\]

空集的行列式和乘积均定义为 1；线性相关子集的概率为零。Cauchy--Binet 保证分母
正好是所有分子之和。若
\(L=W^{1/2}B^TBW^{1/2}\)，(5.1) 就是标准 L-ensemble DPP：
\(\mathbb P(S)=\det(L_{S,S})/\det(I+L)\)。[^7]

### 命题 5.1（均值零空间投影）

若 \(S\sim\mathbb P_w\)，则

\[
T(w)=\mathbb E_w[P_{\ker B_S^T}].
\tag{5.2}
\]

#### 证明与准确归属

令 \(G=BW^{1/2}\)。Dereziński--Khanna--Mahoney 2020 的 Supplemental
Lemma 5（本地合并 PDF pp.15--16）证明：若
\(S\sim\operatorname{DPP}(G^TG)\)，则所选列空间投影满足[^5]

\[
\mathbb E[P_{\operatorname{span}G_S}]
=G(I+G^TG)^{-1}G^T
=I-(I+GG^T)^{-1}.
\]

概率为正的子集不会包含零权列，故
\(\operatorname{span}G_S=\operatorname{span}B_S\)。取正交补便得到 (5.2)。同一公式
后来在 Dereziński--LeJeune--Needell--Rebrova 2025 的 Lemma 15 中明确重述。[^6]
Kassel--Lévy 2023 还梳理了 mean projection theorem 在电网络、DPP 和随机线性代数中的
更早谱系。[^8]

所以 (5.2) **不是新结果**。而且准确说，(5.1) 是可变样本大小的 L-ensemble；给定
\(|S|=k\) 后才是 fixed-size 的加权 volume sampling。不能把二者不加限定地混称为同一个
fixed-\(k\) 恒等式，因为固定大小下的期望一般没有 (5.2) 这么简单。

### 图上的已知特例

若 \(H=\mu I\)，且 \(a_e\) 是无向图的边关联向量，则 \(\ker B_S^T\) 是“在子图
\(([n],S)\) 的每个连通分量上为常数”的空间，\(P_S\) 就是在每个分量内取平均的块投影。
Pilavcı--Amblard--Barthelmé--Tremblay 2021 的 Proposition 2（本地 PDF p.6，
Eqs. (19)--(21)）已经证明：随机 spanning forest 的分区平均算子期望为
\((L+Q)^{-1}Q\)；当 \(Q=\mu I\) 时，条件算子正是这种分量平均投影。[^9]

因此图上的“逆矩阵是随机森林分区投影的平均”也不能作为原创点。可能保留的内容是
对所有非负权重取 supremum 的 sharp converse、投影场景的 exposed 性，以及它们和
common local robust design 的结合。

## 6. common local Schur/WLS 规则的精确 SDP

考虑一个分块正规矩阵

\[
\mathcal J(w)=
\begin{bmatrix}A&E\\E^T&D(w)\end{bmatrix}\succ0,
\qquad
\Sigma(w)=D(w)-E^TA^{-1}E.
\]

令

\[
F=E^TA^{-1},\qquad C=R_iA^{-1}E,
\]

其中 \(R_i\) 选取本地目标。分块逆公式表明中央解的目标行是

\[
\left[
R_iA^{-1}+C\Sigma(w)^{-1}F,
\ -C\Sigma(w)^{-1}
\right].
\tag{6.1}
\]

若算法必须在看到远端权重前选定同一个本地修正 \(Q\)，并只使用本地区域的 primitive
normal-equation RHS，则规则 \([R_iA^{-1}+Q,0]\) 的误差矩阵为

\[
\mathcal E_Q(w)=
\left[C\Sigma(w)^{-1}F-Q,\ -C\Sigma(w)^{-1}\right]
\in\mathbb R^{p\times(k+n)}.
\tag{6.2}
\]

这里 \(Q\in\mathcal Q_{\rm loc}\subseteq\mathbb R^{p\times k}\)，
\(C\in\mathbb R^{p\times n}\)，\(F\in\mathbb R^{n\times k}\)。若输入是原始观测
\(z\) 而不是 primitive RHS，还必须把 \(b=H_{\rm meas}^TR^{-1}z\) 的右因子保留在
误差中；下面的 operator norm 不能无条件冒充 raw-measurement WLS 范数。

### 命题 6.1（有限盒的精确 robust-local SDP）

若

\[
\Sigma(w)=H+\sum_ew_ea_ea_e^T,qquad w\in\mathcal W,
\]

则

\[
\inf_{Q\in\mathcal Q_{\rm loc}}
\sup_{w\in\mathcal W}\|\mathcal E_Q(w)\|_2
=
\inf_{Q\in\mathcal Q_{\rm loc}}
\max_{\sigma\in\{0,1\}^m}\|\mathcal E_Q(v^\sigma)\|_2.
\tag{6.3}
\]

若 \(\mathcal Q_{\rm loc}\) 是 affine/convex 可表示集合，右侧有精确半正定 epigraph：

\[
\begin{array}{ll}
\operatorname{minimize}_{Q,\varepsilon}&\varepsilon\\[1mm]
\operatorname{subject\ to}&Q\in\mathcal Q_{\rm loc},\\[1mm]
&
\begin{bmatrix}
\varepsilon I_p&\mathcal E_Q(v^\sigma)\\
\mathcal E_Q(v^\sigma)^T&\varepsilon I_{k+n}
\end{bmatrix}\succeq0,
\quad\forall\sigma\in\{0,1\}^m.
\end{array}
\tag{6.4}
\]

证明只需注意，对固定 \(Q\)，\(\mathcal E_Q\) 是
\(T=H^{1/2}\Sigma^{-1}H^{1/2}\) 的 affine 函数，谱范数是连续凸函数；应用定理
2.1 即得 (6.3)。(6.4) 是 \(\|X\|_2\le\varepsilon\) 的标准 Schur LMI，两个对角
块的维数分别是 \(p\) 和 \(k+n\)。

### 命题 6.2（无界权重的精确投影场景 SDP）

若 \(w\ge0\) 无上界，令

\[
S_S=H^{-1/2}P_SH^{-1/2}.
\]

则逐个固定 \(Q\) 有

\[
\sup_{w\ge0}\|\mathcal E_Q(w)\|_2
=\max_{S\subseteq[m]}
\left\|[CS_SF-Q,-CS_S]\right\|_2,
\tag{6.5}
\]

因此再对 \(Q\) 取 infimum 仍保持等号，并可用与 (6.4) 相同的有限组 LMI。这里左侧
一般是 supremum 而不是由有限 \(w\) 达到的 maximum。

### 实际收益与不能回避的代价

- 连续的可靠度/精度区间不再需要 gridding；(6.4) 是确定性、无漏点的证书。
- 场景只与不确定的边界因子数 \(m\) 有关。当一个节点的边界端口很小，离线设计对
  \(m\) 是 exact fixed-parameter tractable；角点可并行计算，也可按 Gray code 用
  rank-one 更新把稠密逆更新成本降到每场景 \(O(n^2)\)。
- 最坏情况下仍有 \(2^m\) 个不同投影。定理 4.1 的 exposed 结论说明，对任意凸损失
  的统一理论不能保证继续删去某个不同场景。大 \(m\) 时需要利用具体损失、拟阵 flats
  去重、cut generation 或近似抽样；DPP 抽样可做近似评估，但不是 exact robust certificate。

## 7. 忽略方向结构会产生无界保守比

下面给出一个可直接放入论文的命题，说明保留投影方向不仅是理论美化。

### 命题 7.1（full-block Loewner 外逼近的无界保守比）

真实结构化无界族取

\[
\mathcal T_{\rm str}
=\{\operatorname{diag}(t,1):0\le t\le1\},
\]

它由 \(H=I_2\)、\(a=e_1\)、\(w\ge0\) 产生，其中
\(t=(1+w)^{-1}\)。令

\[
C=[1,0],\qquad F=\operatorname{diag}(1,M),\qquad Q=[q_1,q_2],
\]

并使用损失

\[
\left\|[CTF-Q,-CT]\right\|_2.
\tag{7.1}
\]

则真实结构化 minimax 值恰为 1；若把集合外逼近为
\(\mathcal T_{\rm full}=\{T:0\preceq T\preceq I\}\)，则相同 common \(Q\) 的
minimax 值至少为 \(M/2\)。所以二者比值至少为 \(M/2\)，随 \(M\to\infty\) 无界。

#### 证明

对任意 \(Q\)，取 \(T=I\) 时第二块 \(-CT=[-1,0]\)，所以真实 minimax 至少为 1。
取 \(Q_*=[1,0]\) 后，对 \(T=\operatorname{diag}(t,1)\)，(7.1) 等于

\[
\sqrt{(1-t)^2+t^2}\le1,
\]

故真实值正好为 1。

另一方面，\(\mathcal T_{\rm full}\) 包含

\[
P_\pm=\frac12
\begin{bmatrix}1&\pm1\\\pm1&1\end{bmatrix}.
\]

对任意共同 \(Q=[q_1,q_2]\)，\(CP_\pm F-Q\) 的第二个坐标分别是
\(\pm M/2-q_2\)。因此

\[
\max_{s\in\{+,-\}}
\|[CP_sF-Q,-CP_s]\|_2
\ge
\max\{|M/2-q_2|,|-M/2-q_2|\}
\ge M/2.
\]

证毕。

这个例子是合法的 WLS/Schur 实例，而不是随意拼出的矩阵函数。取

\[
K(w)=I_2+we_1e_1^T,
\qquad
\mathcal J(w)=
\begin{bmatrix}
I_2&F^T\\
F&K(w)+FF^T
\end{bmatrix}.
\tag{7.2}
\]

上左块的 Schur complement 恰为 \(K(w)\)。取 \(R_i=[1,0]\)，则
\(E=F^T\)、\(F=E^TI_2^{-1}\)、\(C=R_iI_2^{-1}E=[1,0]\)。而且

\[
\mathcal J(w)=
\begin{bmatrix}I_2&F\\0&K(w)^{1/2}\end{bmatrix}^{T}
\begin{bmatrix}I_2&F\\0&K(w)^{1/2}\end{bmatrix},
\]

所以它直接是单位噪声 WLS 的 Gram 信息矩阵；每行最多连接一个 local coordinate 和
一个 exterior coordinate，剩余行为 unary factors。

## 8. 原创性与投稿风险审计

| 结论 | 最近的已知结果 | 当前判断 | 论文中建议定位 |
|---|---|---|---|
| 有限盒角点最大值 | Barmish 1994 multiaffine 顶点引理；Blanchini et al. 2022 Theorem 2；Alamo et al. 2008 | **已明显撞车** | 基础引理，不宣称原创 |
| 显式矩阵凸组合 (2.1) | positive-denominator rational multiaffine / Bernstein 机制；Blanchini 的标量 ratio 证明非常接近 | 可能是更强表述，但增量小 | 可写“matrix-valued strengthening”，不可当主核 |
| 单参数端点 iff | 未检索到同措辞；谱证明很短 | **novelty pending，数学正确性高** | sharp boundary lemma |
| DPP 均值投影 (5.2) | Dereziński et al. 2020 Lemma 5；2025 Lemma 15；Kassel--Lévy 2023 | **已知** | 概率解释与证明捷径 |
| 图分区平均表示 | Pilavcı et al. 2021 Proposition 2；matrix-forest 文献 | **已知** | 图特例背景 |
| 无界闭凸包 (4.3) | 可由已知 DPP 均值 + 权重极限短推；未见相同 sharp-converse 包装 | **可能新组合，置信度中低** | 与应用定理合并，不单独夸大 |
| 每个不同投影均 exposed | 基础投影几何，一行 Frobenius 恒等式；未见在此族中强调 | 正确但单独较薄 | 证明场景不可统一删减 |
| common Schur/WLS exact LMI | robust control 有广泛 vertex/LMI 与 uncertain inverse 文献；未见本结构同式 | **novelty pending** | 最有应用价值的命题之一 |
| full-block 保守比无界 | 未检索到该二维 common-rule 构造 | **novelty pending** | 用来证明结构化 hull 的实质收益 |

必须避免的表述：

- “首次发现 rank-one inverse 的顶点原理”；
- “rank one 是端点性质的必要条件”；
- “DPP/随机森林表示是本文新发现”；
- “无论 \(m\) 多大都得到高效算法”；
- “该结果解决一般相关噪声 raw-WLS”，除非把观测到 primitive RHS 的右因子明确加入。

较稳妥的贡献表述是：

> 对独立标量信息因子的 bounded/unbounded 可靠度不确定性，给出矩阵值精确凸包及其
> 零空间投影极限，并把它转化为 common local Schur 规则的有限、无保守 SDP；同时
> 给出 shared higher-rank 权重的精确单参数边界和一个证明方向无关 full-block
> 外逼近可任意保守的 WLS 实例。

即使采用这段措辞，投稿前仍应针对 “convex hull of inverses of rank-one affine matrix
families”“regularized volume sampling inverse representation”“robust common approximate
inverse under diagonal uncertainty” 做数据库级逐定理查重。当前 web/本地语料的阴性结果
只能支持 novelty pending，不能证明不存在先例。

## 9. 可复现检查与归档

数值脚本：

```text
python experiments/theory_search/rank_one_inverse_vertices.py
```

截至 2026-09-14，以下五项均输出 `PASS`：

1. finite-box recursive corner decomposition；
2. determinant/DPP projection mixture；
3. unbounded projection limits；
4. shared rank-two endpoint failure；
5. structured-vs-full-block unbounded gap。

Blanchini et al. 2022 已归档为：

```text
literature/08_robust_control_boundary/2022_blanchini_colaneri_giordano_zorzan_vertex_results.pdf
```

文件为 43 页、1,721,303 bytes，SHA-256 为
`0cb75f87351596418aee15913ff314f333f469be37c300d180ff897b45bf59db`；对应文本已写入
`research/corpus/08_robust_control_boundary/`，并已更新 category 的 `SOURCE_URLS.md`、
`SHA256SUMS` 和全局 `research/corpus/manifest.tsv`。

## 10. 最终建议

这条线值得保留，但不建议把“有限盒 \(2^m\) 角点”本身升格为主定理。更合理的结构是：

1. 用定理 2.1 作已知 machinery 的自包含矩阵引理；
2. 用定理 3.1 精确交代 rank one 与 shared block 的数学边界；
3. 把定理 4.1、命题 6.2 和命题 7.1 合成一个应用主结果：真实 completion 的极端点
   是哪些、common local rule 如何精确优化、为什么保留这些方向能产生不可替代的收益；
4. 再增加一个非枚举算法结果、场景分离复杂度结果，或 hidden-node/partial-partition 的
   更强 converse，才更接近“强应用数学主定理”的标准。

当前判断：**数学正确性 GO；有限盒原创性 NO-GO；整套无界结构化 robust-WLS 组合
NOVELTY-PENDING；作为更大主定理的引擎 GO。**

## 注释

[^1]: Sherman, J.; Morrison, W. J. “Adjustment of an Inverse Matrix Corresponding to a Change in One Element of a Given Matrix.” *Annals of Mathematical Statistics* 21(1), 1950, 124–127. https://doi.org/10.1214/aoms/1177729893
[^2]: Barmish, B. R. *New Tools for Robustness of Linear Systems*. Macmillan, 1994, Lemma 14.5.5.
[^3]: Blanchini, F.; Colaneri, P.; Giordano, G.; Zorzan, I. “Vertex results for the robust analysis of uncertain biochemical systems.” *Journal of Mathematical Biology* 85, Article 35, 2022. https://doi.org/10.1007/s00285-022-01799-z; open full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC9485104/
[^4]: Alamo, T.; Tempo, R.; Rodríguez Ramírez, D.; Camacho, E. F. “A new vertex result for robustness problems with interval matrix uncertainty.” *Systems & Control Letters* 57(6), 2008, 474–481. https://doi.org/10.1016/j.sysconle.2007.11.003
[^5]: Dereziński, M.; Khanna, R.; Mahoney, M. W. “Improved guarantees and a multiple-descent curve for Column Subset Selection and the Nyström method.” *NeurIPS 33*, 2020, 4953–4964; Supplemental Lemma 5 and Appendix C. https://proceedings.neurips.cc/paper/2020/hash/342c472b95d00421be10e9512b532866-Abstract.html
[^6]: Dereziński, M.; LeJeune, D.; Needell, D.; Rebrova, E. “Fine-grained Analysis and Faster Algorithms for Iteratively Solving Linear Systems.” *Journal of Machine Learning Research* 26(144), 2025, 1–49, Lemma 15. https://www.jmlr.org/papers/volume26/24-1906/24-1906.pdf
[^7]: Kulesza, A.; Taskar, B. *Determinantal Point Processes for Machine Learning*. Foundations and Trends in Machine Learning 5(2–3), 2012, 123–286. https://doi.org/10.1561/2200000044
[^8]: Kassel, A.; Lévy, T. “On the mean projection theorem for determinantal point processes.” *ALEA* 20, 2023, 497–504. https://doi.org/10.30757/ALEA.v20-17
[^9]: Pilavcı, Y.; Amblard, P.-O.; Barthelmé, S.; Tremblay, N. “Graph Tikhonov Regularization and Interpolation via Random Spanning Forests.” *IEEE Transactions on Signal and Information Processing over Networks* 7, 2021, 359–374, Proposition 2. https://doi.org/10.1109/TSIPN.2021.3084879

## Sources

1. Sherman, J.; Morrison, W. J. “[Adjustment of an Inverse Matrix Corresponding to a Change in One Element of a Given Matrix](https://doi.org/10.1214/aoms/1177729893).” 1950.
2. Barmish, B. R. *New Tools for Robustness of Linear Systems*. 1994, Lemma 14.5.5.
3. Blanchini, F.; Colaneri, P.; Giordano, G.; Zorzan, I. “[Vertex results for the robust analysis of uncertain biochemical systems](https://doi.org/10.1007/s00285-022-01799-z).” 2022. [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9485104/).
4. Alamo, T.; Tempo, R.; Rodríguez Ramírez, D.; Camacho, E. F. “[A new vertex result for robustness problems with interval matrix uncertainty](https://doi.org/10.1016/j.sysconle.2007.11.003).” 2008.
5. Dereziński, M.; Khanna, R.; Mahoney, M. W. “[Improved guarantees and a multiple-descent curve for Column Subset Selection and the Nyström method](https://proceedings.neurips.cc/paper/2020/hash/342c472b95d00421be10e9512b532866-Abstract.html).” NeurIPS 2020, Supplemental Lemma 5.
6. Dereziński, M.; LeJeune, D.; Needell, D.; Rebrova, E. “[Fine-grained Analysis and Faster Algorithms for Iteratively Solving Linear Systems](https://jmlr.org/papers/v26/24-1906.html).” JMLR 26(144), 2025, Lemma 15.
7. Kulesza, A.; Taskar, B. “[Determinantal Point Processes for Machine Learning](https://doi.org/10.1561/2200000044).” 2012.
8. Kassel, A.; Lévy, T. “[On the mean projection theorem for determinantal point processes](https://doi.org/10.30757/ALEA.v20-17).” 2023.
9. Pilavcı, Y.; Amblard, P.-O.; Barthelmé, S.; Tremblay, N. “[Graph Tikhonov Regularization and Interpolation via Random Spanning Forests](https://doi.org/10.1109/TSIPN.2021.3084879).” 2021, Proposition 2.
