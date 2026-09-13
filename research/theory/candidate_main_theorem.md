# Local Schur-interval minimax：候选主定理

## 0. 状态与结论

本文只记录数学闭合性，不作原创性声明。候选主定理达到以下形式：

- 对任意实对称 Loewner 区间给出 robust local estimator 的 exact if-and-only-if LMI；
- 从 SDP 的最优 \(Q\) 直接构造共同局部估计器；
- 给出 LMI 阶数、退化区间处理和 scalar spectral window 下的 cut-dimension 压缩；
- 给出一个优化过 \(Q\) 后仍击穿“只查两个各向同性端点”的解析反例。

数学链条完全闭合，但独立 novelty 审计已经定位到直接上位结果：单个 full-block 的 Petersen lemma/lossless full-block S-procedure 已包含这里的 exact LMI；van Waarde–Camlibel–Eising–Trentelman (2023) Proposition 4.16 可逐变量专门化到本问题，给 \(Q\) 加局部性仿射约束并不会逃离该结论。[^1] 因此本文的 SDP 定理是网络估计语境下的正确 specialization，不是原创主定理。一般矩阵区间只是在 full block 前后加入 \(\Delta^{1/2}\) whitening；scalar cut 压缩则是旋转对称性的初等推论，也不足以单独恢复原创性。

全文在有限维实 Hilbert 空间中陈述。复数情形需要把转置改为伴随并单独引用复 S-lemma，不在此处自动外推。

## 1. 一般矩阵区间的单向量像

令

\[
S_-=S_-^T\preceq S_+=S_+^T,\qquad
\bar S:=\frac{S_-+S_+}{2},\qquad
\Delta:=\frac{S_+-S_-}{2}\succeq0.
\]

记

\[
\mathcal I(S_-,S_+):={S=S^T:S_-\preceq S\preceq S_+\}.
\]

### 引理 1（正定半宽时的精确椭球像）

若 \(\Delta\succ0\)，则对任意 \(v\in\mathbb R^n\)，

\[
\boxed{
\{Sv:S\in\mathcal I(S_-,S_+)\}
=
\left\{z:
(z-\bar Sv)^T\Delta^{-1}(z-\bar Sv)
\le v^T\Delta v
\right\}.}
\tag{1.1}
\]

#### 证明

任意区间元素唯一写成

\[
S=\bar S+\Delta^{1/2}K\Delta^{1/2},\qquad
K=K^T,\quad -I\preceq K\preceq I.
\tag{1.2}
\]

令 \(a=\Delta^{1/2}v\)。则

\[
z-\bar Sv=\Delta^{1/2}Ka,
\]

所以左乘 \(\Delta^{-1/2}\) 后的范数不超过 \(\|a\|\)，得到“\(\subseteq\)”。

反过来，给定满足右侧的 \(z\)，令

\[
b=\Delta^{-1/2}(z-\bar Sv),\qquad \|b\|\le\|a\|.
\]

若 \(a=0\)，则 \(v=0,z=0\)，任取 \(K\) 即可。若 \(b=0\)，取 \(K=0\)。其余情形令

\[
\rho=\frac{\|b\|}{\|a\|}\le1,
\]

并选择一个对称正交 Householder 矩阵 \(H\)，使
\(H(a/\|a\|)=b/\|b\|\)。取 \(K=\rho H\)，则 \(K\) 是自伴收缩且
\(Ka=b\)。代回 (1.2) 即构造出满足 \(Sv=z\) 的区间元素。证毕。

### 引理 2（半宽奇异时的精确退化椭球）

设 \(\operatorname{rank}\Delta=r\)，取任意满列秩
\(B\in\mathbb R^{n\times r}\) 使 \(\Delta=BB^T\)。则

\[
\boxed{
\{Sv:S\in\mathcal I(S_-,S_+)\}
=\{\bar Sv+By:\|y\|^2\le v^T\Delta v\}.}
\tag{1.3}
\]

等价地，位移 \(w=z-\bar Sv\) 必须满足

\[
w\in\operatorname{ran}\Delta,\qquad
w^T\Delta^\dagger w\le v^T\Delta v.
\]

#### 证明

Loewner 序迫使 \(S-\bar S\) 在 \(\ker\Delta\) 上为零，并给出表示

\[
S-\bar S=BK B^T,\qquad K=K^T,\quad\|K\|_2\le1.
\]

于是位移为 \(B K B^Tv\)。集合
\(\{KB^Tv:K=K^T,\|K\|\le1\}\) 是半径 \(\|B^Tv\|\) 的欧氏球，证明与引理 1 的 Householder 构造相同。又因
\(\|B^Tv\|^2=v^T\Delta v\)，得 (1.3)。证毕。

特别地，\(\Delta=0\) 时区间只有一个矩阵 \(S=\bar S\)。

当 \(\Delta\succ0\) 时，(1.2) 还是仿射双射，因此区间极点恰为

\[
S=\bar S+\Delta^{1/2}K\Delta^{1/2},\qquad
K=K^T,\quad K^2=I.
\tag{1.4}
\]

固定 \(Q\) 的目标关于 \(S\) 凸，所以可选一个上述极点作为最坏矩阵。
注意在非交换情形，这只表示 whitened contraction 的特征值为
\(\{\pm1\}\)，并不表示 \(S\) 与 \(S_\pm\) 共享特征向量或其普通特征值逐个等于端点特征值。

## 2. 一般 Loewner 区间的 lossless finite SDP

给定

\[
C\in\mathbb R^{p\times n},\qquad
F\in\mathbb R^{n\times k},\qquad
Q\in\mathbb R^{p\times k},
\]

定义

\[
\mathcal R(C,F;S_-,S_+)
:=\inf_Q\sup_{S_-\preceq S\preceq S_+}
\|[CSF-Q,-CS]\|_2.
\tag{2.1}
\]

### 定理 3（正定半宽：exact if-and-only-if SDP）

假设 \(\Delta\succ0\) 且 \(C\ne0\)。令

\[
\mathcal D=
\begin{bmatrix}
-\Delta^{-1}&\Delta^{-1}\bar S C^T\\
C\bar S\Delta^{-1}&
C(\Delta-\bar S\Delta^{-1}\bar S)C^T
\end{bmatrix},
\qquad
\mathcal E=
\begin{bmatrix}0_{n\times n}&0\\0&I_p\end{bmatrix},
\tag{2.2}
\]

并定义

\[
\mathcal L_Q=
\begin{bmatrix}
F^T&-Q^T\\
I_n&0
\end{bmatrix}
\in\mathbb R^{(k+n)\times(n+p)}.
\tag{2.3}
\]

则对每个 \(\varepsilon\ge0\)，以下两件事等价：

1. 存在一个 \(Q\) 使
   \[
   \sup_{S_-\preceq S\preceq S_+}
   \|[CSF-Q,-CS]\|_2\le\varepsilon;
   \tag{2.4}
   \]
2. 存在 \(Q\) 和 \(\lambda\ge0\) 使
   \[
   \boxed{
   \begin{bmatrix}
   -\varepsilon^2\mathcal E+\lambda\mathcal D&\mathcal L_Q^T\\
   \mathcal L_Q&-I_{k+n}
   \end{bmatrix}\preceq0.}
   \tag{2.5}
   \]

因此

\[
\boxed{
\mathcal R(C,F;S_-,S_+)^2
=\min_{\tau,Q,\lambda}\ \tau}
\tag{2.6}
\]

满足

\[
\tau\ge0,qquad\lambda\ge0,qquad
\begin{bmatrix}
-\tau\mathcal E+\lambda\mathcal D&\mathcal L_Q^T\\
\mathcal L_Q&-I_{k+n}
\end{bmatrix}\preceq0.
\tag{2.7}
\]

#### 证明

固定 \(S,Q\)，由谱范数的左奇异向量变分式，

\[
\|[CSF-Q,-CS]\|_2^2
=\sup_{\|u\|=1}
\left(\|F^TSC^Tu-Q^Tu\|^2+\|SC^Tu\|^2\right).
\tag{2.8}
\]

交换两个 supremum 不改变值。对固定 \(u\)，令

\[
v=C^Tu,\qquad z=Sv,\qquad w=[z;u].
\]

由引理 1，允许的 \(z\) 满足

\[
(z-\bar SC^Tu)^T\Delta^{-1}(z-\bar SC^Tu)
\le u^TC\Delta C^Tu.
\tag{2.9}
\]

展开 (2.9)，恰为

\[
w^T\mathcal Dw\ge0.
\tag{2.10}
\]

平方损失则为 \(w^T\mathcal L_Q^T\mathcal L_Qw\)，而
\(\|u\|^2=w^T\mathcal Ew\)。约束齐次，并且 \(u=0\) 时 (2.9) 迫使
\(z=0\)，所以 (2.4) 等价于

\[
w^T\mathcal Dw\ge0
\Longrightarrow
w^T(\mathcal L_Q^T\mathcal L_Q-\varepsilon^2\mathcal E)w\le0.
\tag{2.11}
\]

因为 \(\Delta\succ0,C\ne0\)，可取 \(C^Tu\ne0\) 和
\(z=\bar SC^Tu\)，使 (2.9) 严格成立。因此单二次约束的 S-lemma 在此 lossless，(2.11) 当且仅当存在 \(\lambda\ge0\) 使

\[
\mathcal L_Q^T\mathcal L_Q-arepsilon^2\mathcal E
+\lambda\mathcal D\preceq0.
\]

对固定负定块 \(-I_{k+n}\) 取 Schur complement 即得 (2.5)。把
\(\varepsilon^2\) 替换为变量 \(\tau\) 并最小化，得到 (2.6)–(2.7)。证毕。

若 \(C=0\)，取 \(Q=0\) 即有 \(\mathcal R=0\)，不应为了形式统一而调用一个没有严格可行点的 S-lemma。

### 定理 4（奇异半宽的 exact SDP）

设 \(\Delta=BB^T\)、\(B\in\mathbb R^{n\times r}\) 满列秩，并令

\[
\widetilde{\mathcal D}=
\begin{bmatrix}
-I_r&0\\0&C\Delta C^T
\end{bmatrix},qquad
\widetilde{\mathcal E}=
\begin{bmatrix}0_{r\times r}&0\\0&I_p\end{bmatrix},
\tag{2.12}
\]

\[
\widetilde{\mathcal L}_Q=
\begin{bmatrix}
F^TB&F^T\bar SC^T-Q^T\\
B&\bar SC^T
\end{bmatrix}
\in\mathbb R^{(k+n)\times(r+p)}.
\tag{2.13}
\]

若 \(C\Delta C^T\ne0\)，则

\[
\mathcal R(C,F;S_-,S_+)^2
=\min_{\tau,Q,\lambda\ge0}\tau
\]

满足 exact LMI

\[
\boxed{
\begin{bmatrix}
-\tau\widetilde{\mathcal E}
+\lambda\widetilde{\mathcal D}
&\widetilde{\mathcal L}_Q^T\\
\widetilde{\mathcal L}_Q&-I_{k+n}
\end{bmatrix}\preceq0.}
\tag{2.14}
\]

该退化表示的 LMI 阶数为 \(r+p+k+n\)，而不是把伪逆公式放回一个
\(2n+p+k\) 阶 LMI；当 \(r\ll n\) 时它同时给出精确的 uncertainty-rank 压缩。

#### 证明

由引理 2 写

\[
z=\bar SC^Tu+By,\qquad
\|y\|^2\le u^TC\Delta C^Tu.
\]

令 \(\widetilde w=[y;u]\)。可行性是
\(\widetilde w^T\widetilde{\mathcal D}\widetilde w\ge0\)，损失是
\(\|\widetilde{\mathcal L}_Q\widetilde w\|^2\)，归一化仍由
\(\widetilde{\mathcal E}\) 给出。条件
\(C\Delta C^T\ne0\) 恰好提供严格可行点，故重复定理 3 的 lossless S-lemma 和 Schur complement 即得。证毕。

若 \(C\Delta C^T=0\)，则 \(\Delta^{1/2}C^T=0\)，所有允许扰动均满足
\((S-\bar S)C^T=0\)。不确定性对目标完全不可见，问题退化为单场景

\[
\inf_Q\|[C\bar SF-Q,-C\bar S]\|_2,
\]

用一个标准 operator-norm epigraph LMI 精确求解。该分支也覆盖
\(\Delta=0\)。

## 3. Cut 子空间：何时真能压到 \(\dim U+1\)

设 \(U\subseteq\mathbb R^n\) 是 cut 子空间，且

\[
\operatorname{ran}C^T\subseteq U,\qquad
\operatorname{ran}F\subseteq U.
\tag{3.1}
\]

### 命题 5（标量谱窗的精确固定维数压缩）

若

\[
S_-=\alpha I,\qquad S_+=\beta I,\qquad \alpha<\beta,
\]

令 \(q=\dim U\)，并定义

\[
n_{\rm eff}=
\begin{cases}
q,&U^\perp=\{0\},\\
q+1,&U^\perp\ne\{0\}.
\end{cases}
\tag{3.2}
\]

则 (2.1) 的值可以在任意固定的 \(n_{\rm eff}\) 维子空间

\[
W=U\quad\text{或}\quad W=U\oplus\operatorname{span}\{e\},
\qquad e\in U^\perp,\ \|e\|=1,
\]

中精确计算。相应 exact SDP 的 LMI 阶数从 \(2n+p+k\) 降为

\[
\boxed{2n_{\rm eff}+p+k.}
\tag{3.3}
\]

#### 证明

此时引理 1 的椭球是球：对 \(v=C^Tu\in U\)，

\[
z\in B(hv,d\|v\|),\qquad
h=\frac{\alpha+\beta}{2},\quad d=\frac{\beta-\alpha}{2}.
\]

分解 \(z=z_U+z_\perp\)。因 \(F\) 支撑于 \(U\)，损失

\[
\|F^Tz-Q^Tu\|^2+\|z\|^2
\]

只通过 \(z_U\) 和标量 \(\|z_\perp\|\) 依赖正交部分；球约束也只通过同一个范数依赖它。若 \(U^\perp\ne\{0\}\)，映射

\[
z\longmapsto z_U+\|z_\perp\|e
\]

同时保持可行性和损失。反向包含显然，因为右端本就是原空间中的向量。
故任取一个固定 \(e\) 即可同时对所有 \(u,Q\) 压缩。若
\(U^\perp=\{0\}\)，没有额外方向，直接取 \(W=U\)。证毕。

### 命题 6（任意矩阵区间下，support 条件本身不足）

对一般 \(\bar S,\Delta\)，必须区分两种“压缩”：

1. **单个见证的集合论压缩总成立。** 任意 \(z\) 都可写成
   \(z_U+z_\perp\)，而目标只依赖 \(z_U,\|z_\perp\|\)。所以对一个已经选定的最坏
   \((u,z)\)，确实只需 \(U\) 加方向
   \(e=z_\perp/\|z_\perp\|\)。
2. **预先固定的低维线性/矩阵区间模型一般不存在。** 椭球约束包含
   \(\Delta^{-1}\) 和中心 \(\bar SC^Tu\)，通常不在 \(U^\perp\) 上旋转不变。
   把 \(z_\perp\) 替换为它的范数虽保留目标，却不保留原椭球的方向相关可行性；消去方向后得到的投影集合也未必仍是一个矩阵区间椭球。因此仅有 (3.1) 不能像命题 5 那样，把整个 robust SDP 预先替换成固定的
   \(U\oplus\operatorname{span}\{e\}\) 问题。

一个直接秩障碍如下。取 \(n=4\)、
\(U=\operatorname{span}\{e_1,e_2\}\)、\(C=[I_2\ 0]\)、\(F=0\)，并令

\[
\bar S=3I+
\begin{bmatrix}0&I_2\\I_2&0\end{bmatrix},\qquad
\Delta=0.1I.
\]

则 \(S_-=\bar S-\Delta\succ0\)。但

\[
P_{U^\perp}\bar S|_U=I_2
\]

的秩为 2；任何只保留一个固定正交方向的**线性矩阵表示**，其对应 Gram 项秩至多为 1，不能对所有 \(v\in U\) 同时保持
\(\|P_{U^\perp}\bar Sv\|^2\)。因此不存在仅由 support 假设保证的统一
\(q+1\) 线性矩阵区间压缩。这个反例不否认把每个已知见证非线性地编码成
\((z_U,\|z_\perp\|)\)；它否认的是由 support 条件自动得到一个较小、同类型、可直接代入定理 3 的 SDP 数据集。

对每一个固定 \(Q\)，当然总能选一个最坏二元组 \((u,z)\)，而单个
\(z\) 属于 \(U\oplus\operatorname{span}\{P_{U^\perp}z\}\)。这是“最坏见证至多多一个自适应方向”，不是可在求解前固定的 SDP 降维；两者不可混淆。

若额外假设 \(\bar S,\Delta\) 在 \(U^\perp\) 上具有相应旋转对称性，或给出一个已知低维 reducing subspace 包含
\(\bar S U+\operatorname{ran}\Delta\)，则可按该子空间另行精确压缩。一般情形使用定理 3/4 的全维或低秩 \(\Delta\) 版本。

## 4. 嵌回 r-hop 局部状态估计

固定根节点 \(i\) 和收集球 \(B=B_r(i)\)，令 \(O=V\setminus B\)。对 primitive normal equation

\[
Jx=b,\qquad
J=\begin{bmatrix}A&E\\E^T&D\end{bmatrix}\succ0,
\]

定义

\[
\Sigma=D-E^TA^{-1}E,\qquad
S=\Sigma^{-1},\qquad
C=R_iA^{-1}E,\qquad
F=E^TA^{-1}.
\tag{4.1}
\]

Schur 逆公式给出根的精确全局行

\[
\left[R_iA^{-1}+CSF,\ -CS\right].
\tag{4.2}
\]

任何只使用 \(b_B\) 的线性规则都可唯一写成

\[
\widehat x_i=(R_iA^{-1}+Q)b_B,
\tag{4.3}
\]

其单位 \(\ell_2\) primitive-input 最坏误差恰为

\[
\|[CSF-Q,-CS]\|_2.
\tag{4.4}
\]

### 哪些量来自 local view

采用 coefficient-closed view：根在 \(r\) 轮加固定 halo 后获得

\[
A=J_{BB},\qquad E=J_{BO}
\]

及 cut-port 坐标。于是 \(R_i,A,E\) 和通过局部 principal solves 得到的
\(C,F\) 都由 view 确定。又有

\[
\operatorname{ran}C^T,\operatorname{ran}F
\subseteq U:=\operatorname{ran}E^T,
\]

所以 scalar spectral window 可使用命题 5 的 cut 压缩。

根并不知道实际 \(D,\Sigma,S\)，也看不到 \(b_O\)。一般矩阵端点
\(S_-,S_+\) 或标量 \(\alpha,\beta\) 必须作为 completion 类的统一先验/证书给出，不能从 \(A,E\) 自动推出。若先验以 Schur 补区间

\[
0\prec\Sigma_-\preceq\Sigma\preceq\Sigma_+
\]

给出，则逆序性质给出完全等价的

\[
\Sigma_+^{-1}\preceq S\preceq\Sigma_-^{-1};
\]

这里不要求两个端点交换。一般非标量端点还要求 completion 类共享已识别的 exterior 坐标；若不同 completion 的远端维数不同，必须先给出共同端口压缩或只使用 scalar window，不能假装全矩阵端点是 local data。

## 5. 一体化候选主定理

### 定理 7（Schur-interval optimal local estimator）

固定一个 coefficient-closed \(r\)-view，因而固定有限维实矩阵
\(A\succ0,E,R_i\)。在共同 exterior 空间上，令 completion 类为

\[
J_S=
\begin{bmatrix}
A&E\\
E^T&E^TA^{-1}E+S^{-1}
\end{bmatrix},qquad
S_-\preceq S\preceq S_+,\qquad S_-\succ0.
\tag{5.1}
\]

令 \(C,F\) 如 (4.1)，并假设 \(\Delta=(S_+-S_-)/2\succ0\)。则对每个
\(\varepsilon\ge0\)，以下三件事等价：

1. 存在一个对所有 completions 共用、只依赖 \(b_B\) 的线性
   \(r\)-hop 规则，其单位 \(\ell_2\) 输入误差至多 \(\varepsilon\)；
2. 存在 \(Q\) 使
   \[
   \sup_{S_-\preceq S\preceq S_+}
   \|[CSF-Q,-CS]\|_2\le\varepsilon;
   \]
3. 存在 \(Q,\lambda\ge0\) 满足定理 3 的 LMI (2.5)。

最优共同局部误差等于 SDP (2.6) 最优值的平方根。若
\((Q^*,\tau^*,\lambda^*)\) 是 SDP 最优解，则

\[
\boxed{
\widehat x_i=(R_iA^{-1}+Q^*)b_B}
\tag{5.2}
\]

是一个达到该 minimax 值的共同 estimator。反之，任意共同线性 estimator 都对应某个 \(Q\)，所以该值也是 matching lower bound，而不只是构造性上界。

一般 LMI 阶数为

\[
2\dim\mathsf X_O+\dim\mathsf X_i+\dim\mathsf X_B.
\]

其决策变量为 \(Q\) 的
\(\dim\mathsf X_i\times\dim\mathsf X_B\) 个系数、\(\tau\) 和一个 S-lemma multiplier。若
\(S_-=\alpha I,S_+=\beta I\)，令
\(q=\dim\operatorname{ran}E^T\)，则命题 5 把 LMI 阶数精确降为

\[
2(q+\mathbf 1_{U^\perp\ne0})
+\dim\mathsf X_i+\dim\mathsf X_B.
\tag{5.3}
\]

这给出 SDP 表示的多项式规模，而不是对特定求解器 wall-clock 的承诺。在线阶段只需收集 \(b_B\) 并应用 (5.2)；求解 SDP 和对 \(A\) 的 principal solves 可离线按 view 复用。

#### 证明

Schur 公式给出 (4.2)。任意局部线性系数 \(L:b_B\mapsto\widehat x_i\) 都可写成
\(L=R_iA^{-1}+Q\)，故 1 与 2 等价，且没有遗漏其他线性规则。定理 3 给出 2 与 3 的 exact iff。取 SDP 最优 \(Q^*\) 得到 (5.2)；反向对应说明任何规则都不能低于 SDP 值。维数结论来自 (2.5) 和命题 5。证毕。

若 \(\Delta\) 奇异，使用定理 4 替换定理 3；若不确定性对
\(C^T\) 不可见，则退化为单场景。该分支不应塞进 \(\Delta^{-1}\) 公式。

## 6. Matching 反例：端点谱不等于两个标量端点

取

\[
n=2,\quad p=k=1,\quad
S_-=\tfrac12I,\quad S_+=I,\quad
C=[1\ 0],\quad F=[0\ 4]^T.
\]

只检查 \(S=S_-,S_+\) 时，\(CSF=0\)，优化 \(q\) 后得到值 1。完整区间关于
\(q\) 的最坏损失是偶凸函数，故仍在 \(q=0\) 最小。取区间极点

\[
S=\frac12I+\frac12uu^T,\qquad
u=(\cos\theta,\sin\theta)^T,\qquad
y=\cos^2\theta,
\]

平方损失为

\[
\frac14+\frac{19}{4}y-4y^2.
\]

它在 \(y=19/32\) 达到 \(425/256\)，所以真实 minimax 为

\[
\boxed{\frac{5\sqrt{17}}{16}\approx1.288470508>1.}
\]

最坏 \(S\) 的特征值仍全在 \(\{1/2,1\}\)，对应 Schur 补
\(\Sigma=S^{-1}\) 的特征值全在 \(\{1,2\}\)。因此：

- “最坏矩阵可选 endpoint spectrum”成立；
- “只检查 \(mI,MI\)”严格错误；
- 连续投影方向正是 exact SDP 必须处理的非交换内容。

`experiments/theory_search/block_boundary_counterexample.py` 数值复现上述解析值，并验证
\(\tau=425/256,\lambda=17,q=0\) 时 (2.7) 的最大特征值约为机器零。

## 7. 完全闭合与未闭合事项

### 已完全闭合

- 任意正定半宽矩阵区间的精确椭球像；
- 半宽奇异时的 support/伪逆形式；
- 连续矩阵区间 robust minimax 的 lossless finite SDP；
- 退化 Slater 分支；
- estimator 的直接构造和 matching lower bound；
- scalar spectral window 下 \(\dim U\) 或 \(\dim U+1\) 的统一固定压缩；
- 一般各向异性区间只靠 cut support 无法统一压到 \(\dim U+1\) 的秩障碍；
- local-view 数据与 class-level interval prior 的量词划分；
- 优化后仍成立的 2×2 非各向同性端点反例。

### 尚未闭合，但不影响定理正确性

- 一般 \(C,F,S_\pm\) 下 \(Q^*\) 的解析闭式；
- 非标量区间的更强 cut-only 降维所需最弱 reducing-subspace 条件；
- 复 Hilbert 空间版本；
- SDP 的专门化一阶/大规模算法和 bit-complexity。

因此数学包达到了“exact iff + 构造 + matching lower bound/counterexample + 可计算维数”的形式强度，但 novelty 审计已经判定其核心 exact SDP 被经典 machinery 直接包含。它可以作为严谨工具定理或方法节，不能作为论文的原创主定理。

## 8. 已确认的撞车与真正可能逃离的位置

### 8.1 为什么这组增强仍未逃离 full-block machinery

定理 3 的椭球换元把

\[
S=\bar S+\Delta^{1/2}K\Delta^{1/2},\qquad K=K^T,\quad\|K\|\le1
\]

代入一个二次性能不等式。它仍然只有一个不受结构限制的 contraction block。Petersen lemma/full-block S-procedure 对这种“对所有 \(\|K\|\le1\)”的二次矩阵不等式给出 lossless multiplier；本文件的 \(\lambda\) 正是该 multiplier。[^1] 因此：

- 从球推广到各向异性椭球只是可逆 whitening，不产生新的 uncertainty coupling；
- \(Q\) 的 sparsity/locality 若只是一个独立仿射可行集，仍可直接与已知 LMI 联立；
- cut 的 \(q+1\) 压缩只使用 isotropic ball 对 \(U^\perp\) 的旋转不变性，是计算简化而非新的 robust iff；
- 2×2 反例纠正“只查 \(mI,MI\)”的误解，但并不超出 full-block 定理。

### 8.2 至少还要加入什么，才可能形成原创主定理

后续工作必须引入一个不能被单 full-block 独立处理的网络结构耦合，并同时给出 exact converse。可信方向至少包括以下之一：

1. **多个共享或重复 completion blocks。** 多个 cut、多个根或多轮协议共同依赖同一组稀疏/重复不确定块，且同一 completion 必须同时满足所有约束。此时不能为每个根独立选择一个 full block；一般 block-diagonal multiplier 只是上界。需要证明针对某类网络交叠图的 lossless structured multiplier 或给出 sharp gap。
2. **bounded-degree sparse realization。** 不把 \(S\) 允许为整个 Loewner interval，而要求它必须来自某个有界度、有限 block size、指定端口和拓扑规则的 exterior Schur complement。需要刻画哪些投影/椭球点真能由这类网络实现，并给 matching construction/no-go；这会切断“任意 full contraction”假设。
3. **跨节点同一协议。** 不允许每个 view 独立选择自由 \(Q_v\)，而要求所有 \(Q_v\) 由同一个有限轮、固定 message dimension 的 local opcode/图滤波器/更新规则生成。共享参数和一致性约束把多个 robust approximation 问题耦合起来。
4. **全局 completion compatibility 与允许失败。** 对节点分位数不能逐节点各选最坏 completion；必须要求一个全局稀疏网络同时实现所有坏局部 Schur 块。若能给出 completion-pasting 的必要充分条件，再与 \(\delta\)-节点失败量词匹配，才超出独立 full-block SDP。
5. **资源 matching theorem。** 在同一个通信模型中同时限制轮数、每边 bit/message dimension、内存或能量，并证明 estimator construction 与不可改进的 lower bound。仅把 SDP 离线算快或给 \(Q\) 加零模式不够。

最低发表门槛应是上述至少一个方向中的“网络结构 iff + 构造 + matching 反例/下界”，而不是继续对单 Loewner interval 换坐标或加更多端点公式。

## References

[^1]: Henk J. van Waarde, M. Kanat Camlibel, Jaap Eising, and Harry L. Trentelman, “[Quadratic Matrix Inequalities with Applications to Data-Based Control](https://doi.org/10.1137/22M1486807),” *SIAM Journal on Control and Optimization* 61(4) (2023), 2251–2281, especially Proposition 4.16 and its comparison with Petersen’s lemma and the full-block S-procedure. The collision claim here follows the independent variable-by-variable specialization audit; this file does not claim that the paper discusses local state estimation.
