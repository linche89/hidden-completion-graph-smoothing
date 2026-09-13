# 局部状态估计的第二阶段研究方案

## 1. 总体判断

数学研究、大规模模拟和高性能实现都值得做，但三者的作用不同。数学部分决定论文是否有新的定理；模拟负责寻找反例、检查条件是否过强并展示误差随通信半径的变化；高性能实现则证明方法在数万乃至更大网络上仍可运行。计算规模本身不能弥补定理不足，因此实验应从一开始就围绕定理中的量设计，而不是只比较几种算法的运行时间。

建议把论文的主要问题固定为：

> 对一族有环的块 Gaussian/WLS 网络，当算法参数只能由有限邻域内的模型确定时，什么条件保证至少 \(1-\delta\) 比例节点能在与网络规模无关的通信轮数内逼近中心估计？

实验的首要任务是测量三件事：远端数据对每个节点中心估计的实际影响、只使用局部模型的算法与理论最优局部规则之间的差距，以及达到给定误差所需的通信轮数。通信量、计算时间、内存和能量是第二层指标。

## 2. 理论与实验必须共同回答的问题

令中心 WLS 写成

\[
Jx^\star=b,
\qquad
J=H^\top R^{-1}H+Q,
\]

其中 \(Q\succeq0\) 可以表示先验。节点 \(i\) 只输出 \(C_i x^\star\)。对半径 \(r\)，至少需要比较以下四个量：

本节关于唯一中心解和 Schur 补的讨论先假设 \(J\succ0\)。于是每个主子矩阵 \(J_{OO}\) 也正定，下面的消元公式才可直接使用。若 \(J\) 秩亏但指定目标仍可估，需要改在可辨识商空间上求解、固定规范条件或使用广义逆；不能把下面的逆矩阵公式原样套用。

1. **信息本身是否局部**：固定全局模型后，中心估计系数在 \(r\) 邻域外的块行范数有多大；
2. **局部模型是否足够**：保持节点的 \(r\) 邻域不变、改变远端网络时，中心答案最多能变化多少；
3. **算法损失有多大**：可实际运行的局部算法与最佳 \(r\) 轮规则之间差多少；
4. **多少节点成功**：误差不超过 \(\epsilon\) 的节点比例是否至少为 \(1-\delta\)。

第 1 项可以作为“知道全局系数时”的理论下限。第 2 项才真正反映算法只知道局部模型的困难。若把二者混在一起，实验可能看似成功，实际却把全局逆矩阵预先装进了每个节点。

## 3. 已得到的数学结论

### 3.1 固定模型：最宽的精确充要条件

必须把网络模型与本次观测数据分开。令 \(\theta\) 表示通信图、系数、数据归属和目标坐标，令

\[
T_i(\theta):z\longmapsto x_i^\star
\]

为中心估计在节点 \(i\) 的线性算子，\(S_{i,r}^\theta z\) 抽取 \(r\) 跳内可获得的数据。对固定的 \((\theta,i)\) 和归一化输入 \(\|z\|_2\le1\)，允许任意确定性解码器时有精确等式

\[
\inf_f\sup_{\|z\|_2\le1}
\|T_i z-f(S_{i,r}z)\|_2
=\|T_i(I-S_{i,r}^*S_{i,r})\|_{2\to2}. \tag{3.1}
\]

因此，固定模型上节点 \(i\) 能被 \(r\) 轮局部通信精确恢复，当且仅当中心估计算子在球外的数据块为零；能以误差 \(\varepsilon\) 恢复，当且仅当该块行范数不超过 \(\varepsilon\)。这就是给定通信模型与单位球输入后的最宽条件，既不要求图无环，也不要求每个节点的本地观测矩阵满秩。

但式 (3.1) 允许算法预先知道全局算子 \(T_i\)。它刻画的是“数据是否局部”，还不是“算法只知道本地网络系数”时的答案。

### 3.2 只知道局部模型：共同规则的精确条件

把半径 \(r\) 的带标记局部视图记为 \(v\)。所有拥有同一视图的全局扩展必须使用同一个局部线性系数 \(\ell_v\)，故其精确 minimax 值是

\[
a_r^{\rm lin}(v)=
\inf_{\ell_v}\sup_{(\theta,i):\mathcal V_r(i,\theta)=v}
\|T_i(\theta)-\ell_v S_{i,r}^{\theta}\|_{2\to2}. \tag{3.2}
\]

在有限维且该模型类的算子范数一致有界时，存在一个仅由局部模型决定的线性 \(r\) 轮规则、对所有相容的远端扩展和所有单位输入误差不超过 \(\varepsilon\)，当且仅当 \(a_r^{\rm lin}(v)\le\varepsilon\)。若允许任意非线性规则，精确对象不是 (3.2)，而是对每个本地数据值，所有相容全局输出所成集合的最小包围球半径；这一区别不能省略。

“允许 \(\delta\) 节点失败”也有不同量词。本文建议采用最强且最清楚的一种：先为每种局部视图固定一个规则，再要求每个全局模型中至少 \(1-\delta\) 比例的节点，对所有单位输入都满足误差界。令 \(Q_{1-\delta}\) 表示节点误差的 \((1-\delta)\)-分位数，则相应精确值为

\[
A^{\rm loc}_{r,\delta}=
\inf_{(\ell_v)}\sup_{\theta}
Q_{1-\delta}^{\mu_\theta}
\!\left(
\|T_i(\theta)-\ell_{\mathcal V_r(i,\theta)}S_{i,r}^{\theta}\|
\right). \tag{3.3}
\]

条件 \(A^{\rm loc}_{r,\delta}\le\varepsilon\) 就是这个版本的充要条件。它不能与“对每个输入临时允许不同的失败节点”或“只对随机数据平均的 Bayes 失败率”混用。

### 3.3 SPD 线性方程：Schur 补给出的可计算双向界

先研究正规方程的原始充分统计量 \(Jx=b\)，令 \(J=J^*\succ0\)。把根节点的邻域和外部分别记为 \(B,O\)，写成

\[
J=\begin{bmatrix}A&E\\E^*&D\end{bmatrix},
\qquad T=D-E^*A^{-1}E.
\]

固定模型、只保留 \(b_B\) 时的最佳球外误差与局部零边界解的误差分别为

\[
g_i(B)=\|P_iJ^{-1}P_O\|,
\qquad
d_i(B)=\|P_iJ^{-1}-P_iA^{-1}P_B\|.
\]

令 \(Y=P_iJ^{-1}P_O=-P_iA^{-1}ET^{-1}\)、\(F=E^*A^{-1}\)，则有精确分解

\[
P_iJ^{-1}-P_iA^{-1}P_B=Y[-F,I],
\]

从而

\[
g_i(B)\le d_i(B)
\le\sqrt{1+\|E^*A^{-1}\|^2}\,g_i(B). \tag{3.4}
\]

右侧的局部解和倍数只使用 \(B\) 内的 principal block 及割边系数。进一步，对所有共享同一闭合局部视图的远端扩展，(3.4) 给出

\[
\sup_{\theta\succ v}g_i^\theta(B)
\le a_r^{\rm lin}(v)
\le \gamma_i(B)\sup_{\theta\succ v}g_i^\theta(B),
\quad
\gamma_i(B)=\sqrt{1+\|E^*A^{-1}\|^2}. \tag{3.5}
\]

若 \(\gamma_i\le\Gamma\)，对节点分位数同样有

\[
G^{\rm off}_{r,\delta}
\le A^{\rm loc}_{r,\delta}
\le\Gamma G^{\rm off}_{r,\delta}, \tag{3.6}
\]

其中 \(G^{\rm off}_{r,\delta}\) 是固定模型最佳球外误差的最坏节点分位数。这说明：在局部边界放大因子受控的模型类中，中心算子的行尾对多数节点消失，当且仅当同一个仅看局部模型的 Dirichlet 规则对多数节点成功，误差只差一个显式倍数。也不必要求所有节点都好：\(g_i\le a\) 的失败比例为 \(\delta_1\)、\(\gamma_i\le b\) 的失败比例为 \(\delta_2\) 时，至少 \(1-\delta_1-\delta_2\) 的节点误差不超过 \(ab\)。

统一谱窗 \(mI\preceq J\preceq MI\) 会自动给出仅依赖条件数的 \(\Gamma\)，并通过经典的稀疏逆矩阵衰减导出指数充分条件[^10]；但这条充分条件自 1984 年起已经成熟，不能作为论文的新意。更接近本项目的已有 finite-section 工作也要求全局稳定性和规则图增长条件，并未直接给出任意远端扩展和 \(\delta\)-节点版本[^12]。式 (3.4)–(3.6) 的线性代数成分本身也是标准工具；论文价值要落在严格的局部模型量词、completion 下界、多数节点分位数及资源含义的组合，而不能宣称“发明了 Schur 补”。

### 3.4 Grounded Laplacian：图论与概率论版本

对标量相对测量和少量绝对锚点，信息矩阵为

\[
J=L_W+\operatorname{diag}(\kappa)=S(I-P),
\]

其中 \(P\) 是带吸收的次随机转移矩阵。把输入归一成 \(c=S^{-1}b\)，并令从节点 \(i\) 出发的随机游走在吸收前寿命为 \(\zeta\)、首次走出 \(B_r(i)\) 的时间为 \(\tau_B\)。在 \(\|c\|_\infty\le1\) 下，三个误差有精确含义：

\[
\begin{aligned}
o_i(r)&=\mathbb E_i[\text{吸收前在球外的总访问次数}],\\
d_i(r)&=\mathbb E_i[(\zeta-\tau_B)\mathbf1\{\tau_B<\zeta\}],\\
n_i(r)&=\mathbb E_i[(\zeta-r-1)_+],
\end{aligned}
\qquad
o_i(r)\le d_i(r)\le n_i(r). \tag{3.7}
\]

它们依次对应：知道全局系数后的最佳球内截断、只解局部零边界子问题，以及做 \(r\) 步 Neumann/message-passing。Green 函数、walk-sum 和随机游走占用表示本身都是经典结果[^11][^13]；这里有用之处是把三种误差和多数节点量词放在同一比较中。

在走出球后的期望剩余寿命至多为 \(L\) 时，本地可计算的出球概率 \(h_i(r)=\Pr_i(\tau_B<\zeta)\) 满足

\[
h_i(r)\le o_i(r)\le d_i(r)\le Lh_i(r). \tag{3.8}
\]

因此，对这个模型类，多数节点的出球概率趋零是局部估计成功的常数因子充要条件。它比“全网条件数统一有界”更贴近实际局部影响，但剩余寿命约束不能删掉。

### 3.5 不可能性边界

如果同一个局部视图允许任意规模的远端扩展，又没有外部 Schur 补、resolvent 或剩余寿命的统一控制，那么只要根到边界的传递非零，就可以让外部 Schur 补趋于奇异，使 \(\|P_iJ^{-1}P_O\|\) 任意大。在 grounded Laplacian 中，即使度数和非零权重都有界，只需在球外接一条任意长且没有锚点的悬挂路径，根的球外行尾便可按 \(\sqrt N\) 或 \(N\) 增长。

所以一般模型类中不存在一个只靠有限局部拓扑、同时对任意远端 completion 有效的非平凡统一误差证书。最宽的正确答案必然包含某种远端影响的紧性条件，例如 Schur 补稳定性、吸收后剩余寿命、completion 的分布或规模约束。允许 \(\delta\) 个节点失败只能排除真正受远端病态区域影响的节点，不能仅凭“病态节点数量少”推出其余节点安全；单锚点路径就是直接反例。

这也准确回答了原始两个假设：**无环不是必要条件，每个 \(A_i\) 满列秩也不是必要条件；真正必要的是中心目标可辨识，以及所选通信半径之外的影响在指定范数、模型类和失败量词下足够小。** 不存在一个脱离这些量词的纯拓扑万能条件。

### 3.6 必须优先构造的反例

理论成立前至少要测试以下情形：

- 只有 \(o(n)\) 个节点病态，但接近零的全局模态影响正比例节点；
- 两个网络在所有半径 \(r\) 的根邻域内相同，远端闭合方式不同，中心 WLS 答案相差常数；
- 信息矩阵逆矩阵的块行很长程，但乘上测量矩阵后发生抵消，使最终估计系数仍然局部；
- 测量噪声本身存在长程相关，导致 \(R^{-1}\) 稠密；
- 每个节点的平均误差小，但存在一个全网一致的输入方向使联合误差不消失；
- 少数坏点的数量很小，边界却很大，从而影响大量正常节点。

只要其中一个反例击穿候选等价，就应收紧定理，而不是用更多假设掩盖。

## 4. 测试集的选择

在本次核查的 PGLib、MATPOWER、ACTIVSg、PEGASE、SimBench、JuliaGrid、PowerModelsDistributionStateEstimation、DPLib、SuiteSparse 和 pose-graph 官方生态中，没有发现一个被广泛采用、同时固定“网络参数、真值、测量位置、噪声/协方差与种子、通信图或分区”的完整测试集。这是对已核查范围的判断，不是声称世界上不存在零散作者数据。最合适的做法是使用公开网络参数，按公开脚本生成状态、测量位置、测量噪声、测量的初始持有者和通信图。这样既保留真实拓扑和物理参数，又能控制可观性、误差和病态程度。

| 数据来源 | 可用规模与内容 | 与本题的关系 | 建议用途 |
|---|---|---|---|
| PGLib-OPF | MATPOWER 格式的输电网，从小型 IEEE case 到 78,484 节点；数据采用 CC BY 4.0，软件采用 MIT[^1] | 有真实或合成的母线、支路和电气参数，但没有固定的 WLS 测量方案 | 核心物理基准；由脚本生成 DC/AC 测量和 Gaussian 噪声 |
| Texas A&M ACTIVSg | 完全合成的 2,000、10,000、25,000、70,000 及组合 82,000 节点等输电网；官方允许免费商用或非商用，并要求保留引用说明[^2] | 大规模拓扑和物理参数非常匹配；部分下载需要填写表单 | 最终大规模验证；优先采用 PGLib/MATPOWER 已收录且可稳定下载的版本 |
| SimBench | 从低压到超高压的网络，并带全年负荷、发电和储能曲线；数据库采用 ODbL/DbCL，代码采用 BSD-3-Clause[^3] | 适合研究配电网、异质节点、时间变化工作点 | 泛化与时序压力测试，不作为第一个理论实验 |
| pandapower 自带网络 | IEEE、CIGRE、低压和配电网模型；其状态估计模块支持母线/线路/变压器的电压、功率和电流测量[^4] | 能直接生成和核对 AC WLS 结果 | 小中型正确性基线与 AC 扩展；需逐个核对随附 case 的数据许可 |
| SuiteSparse Matrix Collection | 大量真实稀疏矩阵；矩阵本身总体采用 CC BY 4.0，但还应保留每个矩阵头部的原始引用[^5] | 可检验稀疏求解和逆矩阵衰减，却不天然带有测量持有者和通信图 | 矩阵级压力测试；不能替代状态估计实验 |
| RTS-GMLC | 开放 CSV 及多种电网格式，并提供时间序列；数据可免费使用、复制和分发，但必须保留官方声明并在论文中致谢 DOE/NREL/ALLIANCE[^6] | 小型但资料完整，适合检查时序和测量生成流程 | 端到端可复现实例和单元测试 |
| JuliaGrid | MIT 许可的电力系统分析软件，支持 AC/DC 状态估计、观测性、坏数据分析和测量生成[^7] | 不是主要数据集，但可作为独立中心求解器核对结果 | 小中型交叉验证；本机当前没有 Julia，可后置 |
| DPLib | 提供多区域 PGLib/MATPOWER 网络、区域映射和联络线，是最接近分区通信实验的公开工作[^9] | 没有固定状态估计测量、噪声、真值或中心参考解；仓库也没有标准 LICENSE 文件 | 只引用其分区设计；当前不归档、不再分发代码或派生数据 |
| g2o/GTSAM pose graph | 图优化与 SLAM 因子图；g2o 和 GTSAM 主体采用 BSD 许可[^8] | 在线性化以后也是块稀疏最小二乘，但变量在旋转/位姿流形上 | 作为跨领域外部验证，不能与线性 WLS 主定理混写 |

### 4.1 推荐的三层基准

**第一层：可完全控制的合成网络。** 使用路径、环、二维网格、随机几何图、随机正则图、社区图、哑铃图以及成对的图覆盖实例。合成测量直接写成节点测量与边差分测量，从而可以独立控制锚点数量、每个节点测量矩阵的秩、噪声、条件数、坏节点比例和坏区边界大小。这一层用于验证定理和制造反例。

**第二层：电力网络。** 使用 PGLib 的 14、118、300、1,354、2,869、9,241、13,659、30,000 和 78,484 节点 case，形成从单元测试到单机极限的规模梯度。先做线性的 DC 状态估计，再在少量中型 case 上做 AC Gauss–Newton 线性化。

当前已经按不可变的 v23.07 标签归档 14、118、200、300、1,354、2,000、2,869 和 10,000 节点八个 case；30,000 和 78,484 节点 case 等稀疏生成器与求解流程稳定后再按需下载，避免为尚不能运行的规模提前堆数据。

**第三层：外部泛化。** 使用 SimBench 的配电网与时间序列，另选少量 pose-graph 数据，检查结论是否依赖输电网结构。第三层若不符合主定理的假设，应明确写成经验性扩展。

## 5. 如何从网络参数生成可复现的 WLS 实例

### 5.1 线性图模型

对图的加权关联矩阵 \(B\)，可以生成

\[
z_e=W^{1/2}Bx+\eta_e,
\qquad
z_a=A x+\eta_a,
\]

其中 \(A\) 只在锚点或具有本地传感器的节点上非零。必须同时保存每个测量依赖的状态变量和该测量在第 0 轮由哪些通信节点持有；这两个集合一般不同。信息矩阵为

\[
J=B^\top W B+A^\top R_a^{-1}A+Q.
\]

这正好复现“边上的相对测量 + 少量节点绝对测量”的模型，并允许许多局部 \(A_i\) 不满秩。块状态版本把每条边的标量权重替换为小型矩阵即可。每个实例都应报告局部 \(A_i\) 的秩分布、全局可估计性以及目标功能是否可估，不能只报告信息矩阵条件数。

### 5.2 电力系统 DC 状态估计

状态取母线相角；测量由支路有功潮流、母线有功注入和少量 PMU 相角组成。根据 case 的支路电抗生成相应 Jacobian，在一个参考母线固定相角。测量位置、噪声标准差和随机种子必须写入实例文件。至少设置以下四种传感器方案：

1. 高冗余且均匀覆盖；
2. 只有少量锚点，其余主要是支路测量；
3. 固定比例的局部低质量或缺失测量；
4. 小规模病态区域，其噪声方差或边权呈数量级变化。

### 5.3 AC 扩展

AC WLS 用中心 Gauss–Newton 迭代得到每一步的稀疏正规方程。理论若只针对线性模型，实验应把 AC 结果准确称为“在每个线性化点应用局部求解器”，不能据此声称已经证明非线性版本。

## 6. 算法和基线

实验至少需要四个层次，缺一层就难以解释结果：

1. **中心解**：稀疏直接法或高精度预条件共轭梯度，作为数值真值；
2. **知道全局系数时的最佳局部误差**：计算或估计中心增益块行在半径外的范数；在明确的单位球输入和无限带宽模型下，这既是任何 \(r\) 轮算法的下限，也可由截断正确的中心系数达到；
3. **只知道局部模型的算法**：局部子问题、边界截断、局部 Schur 补近似和低次多项式方法；若 \(J\) 相对通信图的作用距离为 \(\rho\)，\(k\) 次多项式需要的通信半径至多为 \(k\rho\)，不能一概把次数等同于轮数；
4. **已有分布式基线**：Richardson/共轭梯度式迭代、Gaussian BP，以及适用时的重叠 Schwarz 方法。

不能显式形成完整的 \(J^{-1}\)。若需要误差分位数，应随机抽取根节点，并对每个抽到的节点做转置稀疏求解以得到相应块行。Hutchinson 一类随机探针适合估计平均迹或 Frobenius 量，不能未经额外证明就代替逐节点的 \(2\to2\) 范数及其分位数。报告 \(1-\delta\) 成功比例时，应给二项置信区间或 Dvoretzky–Kiefer–Wolfowitz 界，说明抽样节点数量带来的统计误差。

多项式方法所需的谱区间必须注明来源。`eigsh` 给出的近似极端特征值可以用于探索性实验，却不是自动包含完整谱的严格界；需要认证时，应使用构造时已知的界、带残差的外包界或保守扩大的区间。如果谱界来自全局预处理和广播，其代价也必须计入，而不能把多项式迭代写成完全由局部模型产生。

## 7. 评价指标

### 7.1 理论指标

- 每个节点相对中心解的块误差；
- 原始测量单位球下，中心估计系数的球外范数不超过 \(\epsilon\) 的节点比例；
- 达到 \((\epsilon,\delta)\) 所需的最小半径或通信轮数；
- 中心增益块行的邻域外范数；
- 局部误差证书与实际误差的比值；
- 局部算法相对于“知道全局系数的最佳局部下限”的倍数差距。

### 7.2 资源指标

- 每轮、每边和全网消息数；
- 传输标量数、比特数及其与跳数的乘积；
- 每节点峰值内存和本地浮点运算量；
- 端到端时间，以及构建模型、预处理和迭代时间的分项；
- 单机 CPU/GPU 的吞吐量和内存带宽；
- 多机时的强扩展与弱扩展效率。

能耗只能在明确模型下报告。如果硬件计数器、权限和采样频率经过核验，单机可以测量“这台机器上的实现能耗”；当前尚未完成这种核验，因此只能先报告时间和内存。无线网络能耗则必须另选路径损耗、发送/接收成本和介质接入模型，不能把服务器功耗直接解释成无线节点能耗。

## 8. 高性能实现路线

当前工作站为 AMD Ryzen 9 9950X3D（16 核/32 线程）、约 126 GiB 内存和 NVIDIA RTX 5080 16 GiB，足以承担大部分单机实验。现有 Python 环境包含 NumPy、SciPy、NetworkX、Pandas 和 Numba；CUDA 12.8 已安装，但当前 PyTorch 是 CPU 版本，尚无 CuPy、MPI、PETSc 或 pandapower。

一个秒级 CPU smoke benchmark 已经跑通，见 `experiments/locality_benchmark.py`。它在 72 节点有环块 WLS 上通过了中心残差、Schur 消元结果与中心解一致、抽样中心估计系数尾部随半径不增三项检查。修订后的普通配置只让 20% 节点拥有满秩本地锚点；另一个更严格的配置只给 72 个节点中的 1 个节点配置满秩锚点，信息矩阵仍正定，三项检查也全部通过。代码还把“测量依赖的状态变量”和“测量的初始持有者”分开保存。这个结果只证明局部满秩不是数值可解性的必要条件，并验证了接口和公式；当前生成器仍有稠密临时矩阵，不能据此声称数万节点或 GPU/MPI 已经跑通。

建议分三步实现：

1. **正确性版本**：Python + SciPy 稀疏矩阵，完成实例生成、中心解、半径邻域和误差统计；当前 smoke 生成器仍先建立稠密因子行，只能用于小规模接口检查，正式扩展前必须改为直接生成 COO/CSR；
2. **单机加速**：将图遍历和小块运算用 Numba/C++ 并行，矩阵乘法使用优化 CSR 核；需要 GPU 时在项目独立环境中安装与 CUDA 12.8 兼容的 CuPy 或 GPU 稀疏库；
3. **多机版本**：在 Linux 集群上用 PETSc/MPI 实现分区稀疏矩阵与邻域消息，METIS/ParMETIS 只负责数据分区，算法通信仍按理论中的边或分隔集统计。

GPU 适合大量相同结构的稀疏矩阵—向量运算或许多小块批处理；不规则的逐根 BFS 和大小差异很大的局部直接分解可能更适合 CPU。是否使用 GPU 应由性能分析决定，而不是预先把 GPU 写成贡献。

## 9. 实验矩阵

### 阶段一：反例和小规模精确计算

- 在内存和分解开销允许的小实例上对所有节点计算精确中心增益行；规模升高后立即改为均匀或分层抽样节点；
- 扫描半径、锚点比例、噪声和条件数；
- 找出候选定理失败的最小实例。

### 阶段二：真实拓扑的统计验证

- PGLib/SimBench 中 \(n=10^3\) 到约 \(10^5\)；先以实际非零元、分解填充和内存决定上限，不预先保证所有 case 都能直接分解；
- 随机抽取根节点估计误差分位数；
- 固定 \((\epsilon,\delta)\)，比较最小轮数、算法误差和理论上下界；
- 分别报告正常区域和人工植入病态区域。

### 阶段三：高性能扩展

- 单机强扩展：1、2、4、8、16 核；
- 单机规模扩展：逐步增加节点与边直到内存或时间上限；
- GPU 与 CPU 只比较相同数值容差和相同算法；
- 集群弱扩展：每个进程保持固定节点/边数，增加进程数；
- 集群强扩展：固定大实例，增加进程数并记录通信占比。

## 10. 复现要求

每个实验实例都应保存：数据源与版本、原文件校验和、转换脚本版本、随机种子、测量位置、测量依赖的状态变量、测量的初始持有者、噪声模型、通信图、数值精度、停止条件和硬件信息。公开数据保持原样存放在 `datasets/raw/`，所有修改写入 `datasets/processed/`，并在文件名和元数据中注明变换。不得把手工修改后的 case 冒充原始数据。

结果表不直接依赖交互式 notebook。核心实验应由命令行配置文件启动，并产生机器可读的 CSV/Parquet 指标、日志和图表数据。小规模单元测试应验证局部结果在半径覆盖全图时回到中心解。成对图覆盖反例必须保存覆盖映射、选定根节点和局部带标记邻域的显式同构证书，并逐项验证局部系数与本地数据相同；不能把两个独立随机提升图直接称为不可区分反例。

从固定图中有放回均匀抽取节点时，DKW 界可以控制该图上整条误差分布的经验误差。若同一批节点用于多个半径并要声称所有半径同时覆盖，需要做联合校正；不同随机图、测量噪声和系数带来的不确定性则必须通过独立实例重复另行报告。

## 11. 继续或停止的判断标准

下列任一结果都说明需要收缩或调整论文：

- 只能重新得到“统一条件数有界时指数衰减”的成熟充分条件；
- 局部可检验量对中心误差没有必要性，且无法用一个清楚的附加假设补足；
- 只知道局部模型的算法与预装全局系数的下限存在随规模发散的差距；
- 所谓多数节点成功完全来自宽松平均指标，在固定误差分位数下消失；
- HPC 加速只减少常数时间，却不能让可测试网络规模或通信轮数明显扩大。

反之，如果能得到一个允许少数病态区域的双向结构条件、一个本地误差证书、一个达到同阶半径的算法，并在合成反例和真实电网上同时成立，这个项目就有足够坚实的理论与实验价值。

## 12. 近期执行顺序

1. 独立复核 (3.1)–(3.8) 的证明、量词和常数，并对最接近的 finite-section、local linear solver 与 potential-theory 文献做针对定理表述的二次查重；
2. 暂时并行保留两条数学线：一般 block-SPD 的 Schur 双向界作为通用框架，grounded Laplacian 的随机游走刻画与远端扩展反例作为更直观的主定理候选；
3. 把当前小规模生成器改成全程 COO/CSR，加入路径、网格、随机几何图、随机正则图、社区图、哑铃图以及严格验证的成对局部不可区分实例；
4. 实现 PGLib/MATPOWER case 到 DC-WLS 的确定性转换；公开数据的固定版本归档已经完成；
5. 在 14、118、300 和 1,354 节点 case 上完成中心解、局部解、理论界和统计抽样的全部正确性检查；
6. 再扩展到 10,000–78,484 节点，使用随机抽根估计“原始测量单位球下，中心估计系数的球外范数超过阈值”的节点比例；
7. CPU 稀疏版本稳定后再投入 GPU，确认算法确为批量稀疏运算且受算力或带宽限制后才做 CUDA 优化；
8. 单机结果明确受通信或内存限制后，再在 Linux 集群上做 MPI/PETSc 强弱扩展；最后根据反例与结果决定是否加入 AC、SimBench 和 pose graph。

## 13. 当前本地成果

- `datasets/raw/pglib_opf_v23.07/`：8 个原始 PGLib case、上游 LICENSE/README/CHANGELOG、逐文件来源和 SHA-256；
- `datasets/raw/suitesparse_ccby4/`：`usroads`、`wathen100`、`bcsstk18` 原始压缩包、CC BY 4.0 法律文本、来源和 SHA-256；
- `experiments/locality_benchmark.py`：CPU 小规模正确性骨架；
- `experiments/configs/smoke.json`：20% 节点有本地满秩锚点；
- `experiments/configs/single_anchor_smoke.json`：72 个节点中仅 1 个具有本地满秩锚点；
- `experiments/smoke_result.json` 与 `experiments/single_anchor_smoke_result.json`：两个本机正确性运行的机器可读结果；
- `research/agent_reports/theory_phase2.md`：精确 minimax、Schur 双向界、带吸收随机游走表示、反例及文献碰撞审计的完整数学底稿；
- `research/agent_reports/datasets_phase2.md`：逐数据源的完整许可、格式和适用性审计；
- `research/agent_reports/hpc_experiments_phase2.md`：本机能力、CPU/GPU/MPI 分层实现与扩展实验设计；
- `research/agent_reports/phase2_adversarial_notes.md`：对数学量词、测试代码、统计结论和 HPC 声称的反方审计。

两个数据目录中的 `SHA256SUMS` 已在本机重新计算并逐项比对通过。当前共归档约 10.24 MiB；原始文件保持不变，后续生成的 WLS 实例将另存到 `datasets/processed/`。

## Sources

[^1]: Power Grid Lib, “[PGLib-OPF](https://github.com/power-grid-lib/pglib-opf),” v23.07; [license](https://github.com/power-grid-lib/pglib-opf/blob/master/LICENSE). Data: CC BY 4.0; software: MIT.
[^2]: Texas A&M University, “[Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/)”; “[ACTIVSg2000](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg2000/).” The repository describes synthetic cases from 9 to 82,000 buses and states the ACTIVSg data are free for commercial or non-commercial use subject to the displayed use and citation terms.
[^3]: SimBench, “[Datasets](https://simbench.de/en/download/datasets/)”; e2nIEE, “[SimBench License](https://github.com/e2nIEE/simbench/blob/develop/LICENSE).” Database: ODbL/DbCL; code: BSD-3-Clause.
[^4]: pandapower, “[Networks](https://pandapower.readthedocs.io/en/latest/networks.html)”; “[State Estimation](https://pandapower.readthedocs.io/en/docs/estimation.html).” The documentation lists supported network families and measurement types. MATPOWER notes that its software license does not automatically cover every bundled case file; see the [MATPOWER manual](https://matpower.org/docs/MATPOWER-manual-8.0b1.pdf).
[^5]: SuiteSparse Matrix Collection, “[About and License](https://suitesparse-collection-website.herokuapp.com/about).” Matrices are distributed under CC BY 4.0 while matrix-specific metadata and citations must be retained.
[^6]: Grid Modernization Laboratory Consortium, “[RTS-GMLC](https://github.com/GridMod/RTS-GMLC).” The repository contains source CSV and tool-specific formats together with its Data Use Disclaimer Agreement.
[^7]: M. Cosovic et al., “[JuliaGrid](https://mcosovic.github.io/JuliaGrid.jl/stable/),” MIT-licensed power-system analysis framework; see the [state-estimation API](https://mcosovic.github.io/JuliaGrid.jl/dev/api/stateEstimation/).
[^8]: R. Kümmerle et al., “[g2o](https://github.com/RainerKuemmerle/g2o)”; BorgLab, “[GTSAM](https://github.com/borglab/gtsam).” Both projects describe their main libraries as BSD-licensed, with component-level exceptions documented in their repositories.
[^9]: M. Hasanzadeh and A. Kargarian, “[DPLib: A Standard Benchmark Library for Distributed Power System Analysis and Optimization](https://arxiv.org/abs/2506.20819),” 2025; [project repository](https://github.com/LSU-RAISE-LAB/DPLib). The repository describes academic/research use but, as checked on 2026-09-13, does not include a standard root LICENSE file granting general redistribution terms.
[^10]: S. Demko, W. F. Moss, and P. W. Smith, “Decay Rates for Inverses of Band Matrices,” *Mathematics of Computation* 43 (1984), 491–499. [DOI: 10.1090/S0025-5718-1984-0758197-9](https://doi.org/10.1090/S0025-5718-1984-0758197-9).
[^11]: D. M. Malioutov, J. K. Johnson, and A. S. Willsky, “Walk-Sums and Belief Propagation in Gaussian Graphical Models,” *JMLR* 7 (2006), 2031–2064. [Article and PDF](https://jmlr.org/papers/v7/malioutov06a.html).
[^12]: C. Cheng, Y. Jiang, and Q. Sun, “Spatially Distributed Sampling and Reconstruction,” *Applied and Computational Harmonic Analysis* 47 (2019), 109–148. [DOI: 10.1016/j.acha.2017.07.007](https://doi.org/10.1016/j.acha.2017.07.007).
[^13]: Á. Carmona, A. M. Encinas, M. J. Jiménez, and À. Martín, “Random Walks Associated with Symmetric M-matrices,” *Linear Algebra and its Applications* 693 (2024), 324–338. [DOI: 10.1016/j.laa.2023.10.009](https://doi.org/10.1016/j.laa.2023.10.009).
