# Phase 2：局部通信 Gaussian/WLS 状态估计公开测试集尽调

> 尽调截止：2026-09-13（Asia/Shanghai）  
> 目标：为“每个节点只经局部通信，估计自身目标分量”的静态 Gaussian/WLS 项目建立许可清楚、可复现、可扩展的公开 benchmark。  
> 证据标记：**A** = 官方数据页、官方文档或仓库原文件；**B** = 作者论文/预印本；**I** = 本报告基于前述材料作出的建模判断。  
> 归档原则：只保存体量合理且数据许可明确的原始文件；原文件不改名、不改内容，另存来源、许可证和 SHA-256。

## 0. 执行结论

**结论一：在本次审计覆盖的主流官方生态内，没有发现一个被广泛采用、同时固定下列全部对象的 power-system state-estimation benchmark：**

\[
(\text{网络参数},\ \text{真值},\ \text{测量位置/类型},\ \text{噪声与协方差/随机种子},\ \text{通信图与分区}).
\]

这是一个有明确证据边界的“未发现”结论，不是对全球所有论文和私有工业数据不存在性的形式证明。审计范围包括 PGLib-OPF、MATPOWER 及其 state-estimator 示例、IEEE PES test feeders、UW PSTCA、TAMU ACTIVSg、PEGASE、SimBench、JuliaGrid、PowerModelsDistributionStateEstimation.jl（PMDSE）和 DPLib；还抽查了 SuiteSparse、g2o/SE-Sync 与图生成器。最接近的是 2025/2026 版 DPLib：它固定网络分区、区域数据、边界节点和联络线，但仍不固定 SE 测量、噪声、随机实现和中央 SE 参考解；仓库又没有标准 LICENSE 文件，只在 README 写 “Academic and Research Use Only”，所以本次只引用论文/仓库，**不归档其代码或派生数据**。[DPLib 论文，证据 B](https://arxiv.org/abs/2506.20819)；[DPLib 仓库，证据 A](https://github.com/LSU-RAISE-LAB/DPLib)

**结论二：首选物理网络底座是固定版本的 PGLib-OPF v23.07，而不是笼统地说“MATPOWER cases”。** PGLib 数据文件明确采用 CC BY 4.0，软件部分采用 MIT，文件头保留各案例来源；相反，MATPOWER 软件自 5.1 起采用 BSD-3-Clause，但其手册明确警告：随软件分发的 case files 不由该 BSD 许可证统一覆盖。[PGLib v23.07 LICENSE，证据 A](https://raw.githubusercontent.com/power-grid-lib/pglib-opf/v23.07/LICENSE)；[MATPOWER 8.1 手册，第 11--12 页，证据 A](https://matpower.org/docs/MATPOWER-manual-8.1.pdf)

**结论三：PGLib 是网络/OPF benchmark，不是现成的状态估计 benchmark。** 它给出拓扑、支路、母线、机组和运行约束，但没有统一的 SE 测量部署、协方差、噪声种子或通信覆盖。必须由本项目公开冻结一个生成协议；DC 线性 WLS 应是主基准，AC WLS/雅可比线性化是物理外推层。MATPOWER 的旧 `state_est` 示例和 JuliaGrid/PMDSE 可以帮助生成测量和交叉验证，但它们本身没有替我们固定唯一的 benchmark 配置。[PGLib README，证据 A](https://raw.githubusercontent.com/power-grid-lib/pglib-opf/v23.07/README.md)；[MATPOWER state_est 示例，证据 A](https://matpower.org/docs/ref/matpower6.0/extras/state_estimator/state_est.html)；[JuliaGrid measurement model，证据 A](https://mcosovic.github.io/JuliaGrid.jl/dev/tutorials/measurementModel/)

**结论四：推荐三层核心套件。**

| 层 | 数据 | 研究作用 | 是否已归档 |
|---|---|---|---|
| L1 物理电网 | PGLib 14/118/200/300/1354/2000/2869/10000 bus；自行冻结 DC/AC 测量与通信配置 | 现实稀疏拓扑、物理权重、跨规模验证 | 是，v23.07 原文件 |
| L2 通用稀疏/图局部性 | SuiteSparse usroads 的 anchored-Laplacian WLS；wathen100 与 bcsstk18 作 SPD/病态数值压力测试 | 将“拓扑局部性”与“电力语义”解耦；检验条件数和求解器鲁棒性 | 是，原始 Matrix Market tarball |
| L3 可控真值合成 | NetworkX 的树、格、环/团环、几何图、小世界、正则/扩展图、SBM、barbell、BA | 精确控制直径、环、分区割、anchor 密度、秩、条件数、噪声和故障概率 | 按配置生成，不预存大包 |

g2o/SE-Sync pose graph 和真实无线定位数据建议作为 **L4 非线性扩展门**，不纳入核心三层：模型具有群结构、规范自由度和非凸性，且常见仓库的软件许可证不能自动证明其中第三方数据集可再分发。

## 1. 什么才算本项目可复现的 benchmark

仅有一张图或一个稀疏矩阵不够。建议把一个实例冻结为

\[
\mathcal B =
(G_{\rm phys},p,x^\star,\mathcal M,h/H,R,\xi,
G_{\rm comm},\Pi,\mathcal T,x_{\rm ref}),
\]

其中：

- \(G_{\rm phys}\)、\(p\)：物理网络和参数；
- \(x^\star\)：生成无噪声量测的真值；
- \(\mathcal M\)、\(h/H\)：测量位置、类型和非线性/线性算子；
- \(R\)、\(\xi\)：协方差、噪声实现和产生它的种子；
- \(G_{\rm comm}\)：实际允许消息通过的通信图，不能默认与物理图相同；
- \(\Pi\)：区域/代理/控制中心的分区；
- \(\mathcal T\)：每个节点要求输出的目标坐标或线性功能；
- \(x_{\rm ref}\)：中央高精度 WLS/MAP 参考解和求解容差。

对于研究“有限半径局部估计的充要条件”，还必须记录每个节点初始知道哪些模型参数、消息精度/维度、同步模型和失败机制。这些是算法实验协议，不应被悄悄塞进数据文件。

本次采用五条入库准则：

1. 官方、版本化且 URL 稳定；
2. **数据**许可证明确，不能只看到求解器代码开源就推断数据可再分发；
3. 原始文件规模合理，现阶段不为“以后也许用”下载几十 GB；
4. 能直接或经有记录的变换构成稀疏局部因子 WLS；
5. 来源、许可证、哈希和变换链都可审计。

## 2. 电力系统网络数据候选总表

| 来源 | 典型规模/格式 | 数据许可与下载 | 能否构成 node/edge WLS | 优点 | 主要缺口与本次处理 |
|---|---|---|---|---|---|
| PGLib-OPF v23.07 | 3--78,484 buses；MATPOWER v2 `.m`；另有 API/SAD 运行变体 | 数据 CC BY 4.0；软件 MIT；Git tag 可固定 | 可以。支路潮流是端点因子，电压/相角是 unary 因子，注入量测是 1-hop star/hyperedge 因子 | 官方 IEEE PES task-force 库、规模跨度大、源头与版本清楚 | 不含统一 SE 测量/噪声/通信分区；归档精选 8 例 |
| MATPOWER | 大量 `.m` cases；`makeYbus`、`makeBdc`、PF/OPF 工具 | 软件 BSD-3；**case 文件许可不统一** | 可以，并可生成 AC/DC truth 和 Jacobian | 工具成熟、文档完善、中央交叉验证方便 | 不应整包当成统一 BSD 数据；本次不另归档 cases |
| IEEE PES distribution feeders | 13/34/37/123、8500、9500 node 等；常见 OpenDSS/Excel/专用格式 | 官方旧开发仓库的 `LICENSE.md` 内容仍是 “TBD”；IEEE data-sharing 条目标记 license N/A | 适合三相配电 SE，边/节点测量自然 | 行业辨识度高、配电不平衡场景真实 | 许可与格式不够统一；不归档，只有在逐 case 核清权利后再纳入 |
| [UW PSTCA](https://labs.ece.uw.edu/pstca/) | 9--300 bus 经典 transmission cases 为主；MATPOWER/文本 | 公开教学研究服务，但不是现代统一开放数据许可；页面带使用限制/免责声明 | 可生成小型 DC/AC SE | 经典、论文对照多 | 规模较小、权利不如 PGLib 清楚；不归档重复副本 |
| TAMU ACTIVSg | 200/500/2000/10,000/25,000 bus 合成电网；PowerWorld、MATPOWER、PSS/E、PSLF 等 | 官方页面称商业/非商业均可免费使用，但需下载表单并鼓励引用 | 很适合物理规模与地理结构测试 | 完全合成、格式多、具有现实统计特征 | 表单访问；本次不绕过表单。使用 PGLib 中许可清楚的 200-bus ACTIVSg 副本 |
| PEGASE | PGLib 常见 89/1354/2869/9241/13659 bus；MATPOWER | 采用 PGLib 版本时受 PGLib CC BY 4.0 与文件头署名约束 | 适合大规模 transmission SE | 欧洲高压网形态、环路丰富 | 原工程目标是 OPF/潮流；本次归档 1354 与 2869 |
| SimBench | 低压到超高压、多场景和时序 profiles；CSV/pandapower 等 | 数据库 ODbL 1.0、单项内容 DbCL 1.0；代码 BSD-3 | 可生成配电/输电 AC SE 和时序伪量测 | 13 个底层网族、运行场景和 profiles 丰富 | share-alike 与数据结构更复杂；先作为可选扩展，不在最小核心归档 |
| DPLib 2025/2026 | 40 个 5--20,758 bus 多区域 case；`.mat/.csv/.m`，bus-region map 与 tie lines | 仓库无标准 LICENSE；README 仅称 Academic and Research Use Only | 分区层非常适合；SE 测量仍需自行添加 | 目前最贴近“局部区域通信”的公开工作 | 只有 DC/AC OPF 验证器；不固定 SE truth/noise/placement；**不归档** |

IEEE PES feeder 的许可判断来自其[官方原型仓库](https://github.com/ieee-pes-amps/dtf-dev)和内容为 [“TBD” 的 LICENSE.md](https://raw.githubusercontent.com/ieee-pes-amps/dtf-dev/master/LICENSE.md)。这不等于所有 IEEE feeder 永远不可用，而是说不能在没有逐项授权核对时把整套镜像进本项目。

ACTIVSg2000 官方页面明确说明它是完全合成的 2000-bus Texas-footprint 系统，提供多种格式，并称可免费用于商业或非商业用途，但下载要提交表单并引用来源。[ACTIVSg2000 官方页，证据 A](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg2000/) PGLib v20.07 的变更日志又说明 500/2000/10000 bus 的旧 ACTIVSg case 已被相近的 GOC cases 取代，因此本次规模档采用 PGLib 的 2000/10000 GOC，而保留 case200_activ 作为 ACTIVSg 小型代表。[PGLib changelog，证据 A](https://raw.githubusercontent.com/power-grid-lib/pglib-opf/v23.07/CHANGELOG.md)

SimBench 的许可证必须分开理解：数据库 ODbL、数据库内单项内容 DbCL、代码 BSD-3，不能只说“SimBench 是 BSD”。[SimBench LICENSE，证据 A](https://github.com/e2nIEE/simbench/blob/develop/LICENSE) 官方下载页显示网络和 time series 可通过 CSV/pandapower 等方式获得，适合以后增加日内多工况，但不是当前静态线性 WLS 的最短路径。[SimBench datasets，证据 A](https://simbench.de/en/download/datasets/)

## 3. PGLib：首选物理层及规模设计

PGLib 官方 README 明确说它由 IEEE PES 的 emerging-algorithm benchmark task force 维护、case 采用 MATPOWER 格式、原目标是 AC-OPF，并提醒迁移到其他问题时要审慎。因此这里把它当作“物理网络底座”，不把 OPF baseline 误称为 SE ground truth。[证据 A](https://raw.githubusercontent.com/power-grid-lib/pglib-opf/v23.07/README.md)

### 3.1 已归档的八个档位

| case | buses | branches/edges | 上游谱系 | 用途 |
|---|---:|---:|---|---|
| case14_ieee | 14 | 20 | IEEE | 单元测试、手工检查 |
| case118_ieee | 118 | 186 | IEEE | 小型标准对照 |
| case200_activ | 200 | 245 | ACTIVSg | 合成现实拓扑、小规模分区 |
| case300_ieee | 300 | 411 | IEEE | 小中规模回归 |
| case1354_pegase | 1,354 | 1,991 | PEGASE | 真实形态环网 |
| case2000_goc | 2,000 | 3,639 | ARPA-E GOC synthetic | 中型合成 transmission |
| case2869_pegase | 2,869 | 4,582 | PEGASE | 中型物理压力测试 |
| case10000_goc | 10,000 | 13,193 | ARPA-E GOC synthetic | 首轮规模上限 |

节点/边数取自 PGLib v23.07 的[官方 BASELINE](https://github.com/power-grid-lib/pglib-opf/blob/v23.07/BASELINE.md)。库还覆盖 8,387/9,241-bus PEGASE，10,192/20,758/78,484-bus EPIGRIDS，10,480/19,402/24,464/30,000-bus GOC 和 13,659-bus PEGASE；最大 78,484 个节点、126,146 条边。当前不必把这些全部下载：10k 已能验证稀疏线性代数和通信规模趋势；更大档位等实现稳定后按需拉取即可。

### 3.2 MATPOWER 许可边界

MATPOWER 8.1 手册写明软件采用三条款 BSD，但紧接着说明：随 MATPOWER 分发的 case files **不受该 BSD 许可证覆盖**，数据通常是经许可纳入或由公开来源转换；ACTIVSg、PEGASE、RTE 等还要求额外引用文件头中的文献。[MATPOWER manual 8.1，第 11--12 页，证据 A](https://matpower.org/docs/MATPOWER-manual-8.1.pdf) 因而：

- 可以使用 MATPOWER 的开源代码生成 \(Y_{\rm bus},Y_f,Y_t,B_{\rm bus},B_f\) 和中央 PF/OPF；
- 不应写“所有 MATPOWER cases 都是 BSD”；
- 对已进入 PGLib 且有 CC BY 4.0 声明的副本，以固定 PGLib tag 和该文件头为准；
- 新增 MATPOWER-only case 时必须逐文件审核来源和再分发许可。

## 4. 为什么仍需我们冻结一套 SE 测量协议

### 4.1 现有工具各自固定了什么

| 生态 | 网络 | 可复现真值 | 测量位置/类型 | 噪声/协方差 | 通信分区 | 最终判断 |
|---|---|---|---|---|---|---|
| PGLib | 固定 | 需重新求 PF/OPF 并保存 | 无 | 无 | 无 | 最佳网络底座，不是 SE benchmark |
| MATPOWER `state_est` 示例 | 特定 30-bus 示例 | 由脚本求得 | 示例中硬编码 | 给出 sigma 规则并调用 `normrnd`，但未冻结通用 seed/suite | 无 | 很好的测量配方原型，不是大规模标准 |
| JuliaGrid | 用户选择/可导入 MATPOWER、PSS/E | 可由 PF 产生 | API 让用户添加 voltmeter/wattmeter/varmeter/PMU | `noise=true` 可加白 Gaussian；PMU 可处理相关误差 | 无 | 首选生成器和中央解交叉验证工具 |
| PMDSE | OpenDSS/MATPOWER/JSON | 可由 distribution PF 产生 | CSV 明确记录 placement/accuracy | `add_measurements!` 支持 seed；可从 PF 写测量 | 无 | 三相/配电与非线性扩展优秀 |
| DPLib | 固定 | 只有网络/OPF 工作流中的运行数据，不是统一 SE truth | 无统一 SE 配置 | 无统一 SE 配置 | **固定** | 可借鉴分区格式，不能直接充当完整 SE benchmark |

MATPOWER 旧示例的量测向量包括两端有功/无功支路潮流、母线有功/无功注入、相角和电压幅值；标准差采用约 2% 幅值加 floor 的规则，随后调用 `normrnd` 加 Gaussian error。[示例源码，证据 A](https://matpower.org/docs/ref/matpower6.0/extras/state_estimator/state_est.html) 但源码也明确标注后续量测索引选择是“specific to the 30-bus system”，所以应借鉴其统计形式，而不是声称它定义了通用 benchmark。

JuliaGrid 的 2026 文档支持 AC/DC/PMU WLS、观测性、坏数据分析，并允许人工生成

\[
\epsilon_i\sim\mathcal N(0,v_i),\qquad z_i=e_i+\epsilon_i.
\]

其代码采用 MIT 许可证；它适合将我们生成的数据与另一套实现交叉验证，但 meter placement 和 variance 仍由用户提供。[JuliaGrid 首页](https://mcosovic.github.io/JuliaGrid.jl/stable/)；[measurement model](https://mcosovic.github.io/JuliaGrid.jl/dev/tutorials/measurementModel/)；[state-estimation API](https://mcosovic.github.io/JuliaGrid.jl/dev/api/stateEstimation/)

PMDSE 的官方文档把网络、测量、SE settings 分开；测量 CSV 显式包含 component、quantity、distribution 和参数，`add_measurements!` 有 seed，`write_measurements!` 可从 PF 结果生成 voltage/injection measurements。代码为 BSD-3-Clause。它非常适合未来的三相不平衡 distribution-SE 层，但仍没有通信图/区域协议。[PMDSE 仓库](https://github.com/Electa-Git/PowerModelsDistributionStateEstimation.jl)；[输入与测量文档](https://electa-git.github.io/PowerModelsDistributionStateEstimation.jl/v0.4/input_data_format/)

### 4.2 DPLib 是“几乎撞车”，但没有真正撞车

DPLib v4 预印本宣称 40 个 5--20,758 bus 多区域 case，输出区域 MATPOWER 结构、local/global bus mapping、tie-line、boundary 和 bus-to-region map；验证器是 ADMM DC/AC OPF。[论文摘要，证据 B](https://arxiv.org/abs/2506.20819) 仓库 README 也把 distributed state estimation 列为潜在用途，但同时说明用户需在网络/区域结构上添加 application-specific data，仓库展示的验证求解器仍是 OPF。[仓库 README，证据 A](https://github.com/LSU-RAISE-LAB/DPLib)

它没有固定：

- 哪些母线/支路安装何种测量；
- \(R\)、噪声随机种子和坏数据；
- 用哪个 AC/DC 解作 \(x^\star\)；
- 每个区域/节点要求输出的 target；
- distributed SE 的中央参考解和验收容差。

许可上，截止审计日仓库根目录没有标准 LICENSE 文件，GitHub 许可证识别为空；README 只写 “Academic and Research Use Only”。这既不是常见 OSI license，也没有给出完整的复制、修改、再分发条款。论文称 “open-source” 不能替代法律许可文本。因此本项目目前只引用，不 clone、不镜像；若作者补充明确许可证，再重新评估。即便原始 PGLib 数据是 CC BY 4.0，DPLib 自己的工具代码和组织后的全部文件也不能据此被自动推定具有相同许可。

### 4.3 对“完整标准 benchmark 不存在”的证据边界

本报告能诚实支持的命题是：

> 截止 2026-09-13，在上述官方且广泛使用/高度相关的生态中，未发现同时冻结网络、SE truth、measurement placement、noise realization/covariance 和 communication partition 的大规模通用 benchmark。

不能从中推出：

- 没有某一篇论文的作者仓库保存过其特定实验的全部数组；
- 没有电网公司内部数据或受限数据集；
- 世界上任何语言/地区的仓库都不存在类似产物。

真正可发表的贡献不应写成“我们发明了首个电力 SE 数据集”，而应更精确地写成：**在开放许可、版本化网络之上，发布首个面向 fixed-radius/local-information Gaussian/WLS 的完整、资源可审计 protocol 和实例 manifest**。这是否能写“首个”仍应在投稿前再做一次针对性检索。

## 5. 建议冻结的 PGLib 测量生成协议

### 5.1 DC 线性 WLS：核心主任务

对每个 PGLib case：

1. 固定 slack/reference bus，运行 DC PF 或由 AC solution 投影得到 \(\theta^\star\)；
2. 用 MATPOWER `makeBdc` 或 JuliaGrid 构造

\[
P_{\rm bus}=B_{\rm bus}\theta+P_{\rm businj},\qquad
P_f=B_f\theta+P_{\rm finj};
\]

3. 选择量测行并生成

\[
z=H\theta^\star+c+\varepsilon,\qquad
\varepsilon\sim\mathcal N(0,R);
\]

4. 去除 reference angle 后检查 \(H^\top R^{-1}H\) 的秩、最小特征值与条件数；
5. 用高精度 sparse QR/Cholesky/LU 得到中央 \(x_{\rm ref}\)，同时保存 normal-equation residual 和原始 rectangular residual。

因子语义必须标清：

- branch flow \(b_e(\theta_u-\theta_v)\)：真正的 pairwise edge factor；
- bus angle/PMU：unary anchor factor；
- bus injection：涉及该母线及全部一跳邻居，是 star/hyperedge factor，不应伪称 pairwise；
- 若理论只允许 node/edge factor，可用“只选 branch flows + anchors”版本；若允许 1-hop local factor，再另开 injection-inclusive 版本。

### 5.2 AC WLS：物理外推与线性化层

1. 运行 AC PF/OPF，保存收敛的 \(V^\star\)、支路潮流和 solver/version/tolerance；
2. 生成 \(|V|,\angle V,P/Q\) bus injection、两端 \(P/Q\) branch flow 和可选 PMU voltage/current；
3. 噪声用独立 Gaussian 基线，同时保留可选相关 PMU block covariance；
4. 做两条任务：

   - 完整非线性 AC WLS；
   - 在 \(x^\star\) 或统一 flat-start 附近固定 Jacobian \(J\)，形成 \(J^\top R^{-1}J\) 的线性 Gaussian 子问题。

后一条与本项目的静态线性定理最贴近，也能把 power-grid 物理结构和通用 SPD/Green-function 分析接起来。AC Gauss--Newton 是否收敛是另一个问题，不能和“固定线性 WLS 是否可被局部近似”混成一个成功率。

### 5.3 必须至少有四种 measurement placement

| 配置 | 目的 |
|---|---|
| Dense-observable | 全/高比例 branch flows + injections + 少量 anchors，验证正确性与规模 |
| Sparse-observable | 在保证全局可观/满列秩的前提下减少量测，测试接近奇异和远程信息依赖 |
| Locally-rank-deficient | 每个节点/区域单独不满秩，但聚合后可辨识，直接击中“每个 \(A_i\) 满列秩是否过强” |
| Boundary-starved/adversarial | 刻意减少跨区或目标附近量测，制造 separator、远程 anchor 和 local indistinguishability 压力 |

每一类都应存 placement file，而不是在运行时随机挑完就丢掉。量测密度、anchor 密度与 \(R\) 同时影响条件数，必须分别扫描并报告。

### 5.4 通信结构不可默认等于电气结构

每个物理实例至少冻结三种 communication views：

1. **bus-level physical overlay**：通信边等于在运支路，最适合理论与小消息测试；
2. **regional quotient graph**：先将母线分成 \(k\) 区，每区一个 agent，跨区 tie-line 形成邻接；保存 bus-to-region map；
3. **degraded overlay**：在保持连通或明确失联事件的前提下删边、限带宽、加丢包/节点失败，用于 \((\epsilon,\delta)\) 与能耗实验。

分区可用 deterministic spectral partition、METIS（若环境和许可证允许）或本项目自己的 seeded partitioner，但必须保存最终 map；只记录“用了 spectral clustering”不足以复现。DPLib 可作为字段设计参考，不能直接复制其未明确许可代码。

### 5.5 建议派生实例目录模式

~~~text
benchmark_id/
  network_source.json       # upstream tag, URL, source hash, attribution
  physical_network.*        # unmodified or referenced source
  truth.npz                 # x*, PF/OPF status and residual
  measurements.npz          # z, H/J or factor list, R, placement labels
  communication.graphml     # G_comm
  partition.csv             # bus -> agent
  targets.csv               # node/agent -> desired state coordinates
  protocol.json             # seeds, versions, tolerances, failure model
  central_reference.npz     # x_ref and verification residuals
  LICENSES/
~~~

原始 PGLib 文件仍留在 `datasets/raw`；所有噪声、测量和分区放 `datasets/derived`，并用新的 benchmark version 管理，避免误把派生数据说成未经修改的上游副本。

## 6. SuiteSparse：怎样用、怎样不能用

SuiteSparse Matrix Collection 官方说明矩阵数据采用 CC BY 4.0，并要求保留矩阵内 metadata/引用；它鼓励为了重复性不要修改原矩阵，若修改则改名并说明。[官方 About/License，证据 A](https://sparse.tamu.edu/about) 因此本次保留原始 Matrix Market tarball。

| 矩阵 | 规模 | 官方属性 | 本项目映射 | 限制 |
|---|---:|---|---|---|
| Gleich/usroads | \(129{,}164^2\)，330,870 expanded nnz，56 components；另含 xy | binary undirected US road graph | 定向后构造 incidence \(B\)，做相对差量测与 anchored graph-Laplacian WLS | 原矩阵不是 SPD；需每个 component 一个 anchor，或只取最大连通分量 |
| GHS_psdef/wathen100 | \(30{,}401^2\)，471,601 nnz，SPD，\(\kappa_2\approx5.82\times10^3\) | random 2D/3D Wathen finite-element matrix | 作中等条件数 sparse precision/normal system | 仅有 \(Q\) 不等于提供了局部 measurement factorization |
| HB/bcsstk18 | \(11{,}948^2\)，149,090 nnz，SPD，\(\kappa_2\approx3.46\times10^{11}\) | nuclear-station structural stiffness matrix | 病态数值与收敛压力测试 | 无传感器语义；不能拿它支持 power-SE 物理结论 |

数据来自 [usroads](https://sparse.tamu.edu/Gleich/usroads)、[wathen100](https://sparse.tamu.edu/GHS_psdef/wathen100) 和 [bcsstk18](https://sparse.tamu.edu/HB/bcsstk18) 的官方页。网页 nnz 对称展开计数；tar 内 Matrix Market header 只存一个三角部分，所以 stored entries 分别为 165,435、251,001 和 80,519，不是下载损坏。

### 6.1 usroads 的严格 node/edge WLS 构造

给每条无向边任选方向，令 \(B\) 为 incidence matrix：

\[
z_e=x_u-x_v+\epsilon_e,\qquad
z_a=x_a+\epsilon_a .
\]

若 \(C\) 选择 anchor 节点，则

\[
Q=B^\top W_eB+C^\top W_aC,\qquad
b=B^\top W_ez_e+C^\top W_az_a.
\]

当每个连通分量至少一个正权 anchor 时，\(Q\succ0\)。这是一个完全可解释的 pairwise Gaussian WLS，通信图可以直接取 road graph。它最适合检验大直径、空间局部性与远程 anchor 影响。原 usroads 有 56 个 components；核心实验建议同时保留：

- 最大连通分量 + 单 anchor；
- 全图 + 每分量至少一个 anchor；
- 多 anchor 密度扫描。

### 6.2 为什么两个 SPD 矩阵只是数值副基准

任意 SPD \(Q\) 都可写为 \(Q=L^\top L\)，但 Cholesky \(L\) 可能 fill-in，不能据此声称每行测量只依赖图上邻居。wathen100 背后虽有局部有限元装配，下载包却没有给本项目可直接采用的 sensor-factor ownership；bcsstk18 同理。因此二者可比较 GaBP/CG/稀疏直接法、谱近似和病态鲁棒性，却不应和具有明确局部测量语义的 PGLib/usroads 混在一张“物理 SE 精度”主表里。

## 7. Sensor localization、SLAM、g2o 与 SE-Sync

### 7.1 模型上为什么有价值

pose graph 的边测量是两个 pose 的相对变换，天然局部；问题有全局 gauge，自然需要 anchor；回环边控制环结构，和本项目关心的“树上精确、含环局部近似、远程信息衰减”高度相似。g2o 的文本格式用 VERTEX/EDGE 记录状态、相对测量和 information matrix，适合解析成 factor graph。[g2o file format，证据 A](https://github.com/RainerKuemmerle/g2o/wiki/File-Format)

SE-Sync 明确求解 \(SE(2)/SE(3)\) 上的 pairwise relative-pose maximum likelihood，并提供全局最优性的可认证方法；仓库示例包含常用 `sphere2500.g2o` 等数据。[SE-Sync 仓库，证据 A](https://github.com/david-m-rosen/SE-Sync)；[论文，证据 B](https://arxiv.org/abs/1612.07386)，DOI [10.1177/0278364918784361](https://doi.org/10.1177/0278364918784361)

### 7.2 为什么不能当主线性 benchmark

- state 在 \(SE(d)\) 流形上，残差非线性且有旋转；
- 有 gauge ambiguity，需固定 pose 或在 quotient 上比较；
- 求全局非凸 ML 与解一个固定 SPD 线性 WLS 不是同一问题；
- pose-graph 数据常只有测量与初始值，不一定有可信 ground truth；
- 多机器人通信/区域 ownership 通常没有随单机器人 g2o 文件一起定义。

因此它适合作为“我们的局部性思想能否迁移到非线性 factor graph”的次级挑战，不能用它替代线性定理边界的主实验。

### 7.3 许可为何暂不够

g2o 核心代码主要是 BSD，但仓库明确列出 LGPL/GPL 子部件；更关键的是，**软件许可证不自动覆盖从其他项目收集的 benchmark 数据**。[g2o 仓库许可说明，证据 A](https://github.com/RainerKuemmerle/g2o) SE-Sync 仓库的 C++/MATLAB 实现是 LGPL-3-or-later，但其 `data/` 下汇集的若干第三方 pose graphs 没有逐文件数据许可证清单。故本次不下载/镜像这些数据；日后若使用，需回到每个原始数据作者页面逐项核权。

真实 WiFi/RSSI/UWB/视觉惯导定位数据通常解决“移动设备相对固定 anchors 的位置”或端到端感知问题，需要校准、同步、特征提取，并不直接给一个多节点静态 Gaussian pairwise WLS。对 sensor-network localization 的核心数学压力测试，带坐标真值的 random geometric graph + noisy ranges/bearings 更干净；真实数据只在算法具备相应 observation model 后加入。

## 8. 可控真值图生成器

NetworkX 3.6.1 官方文档覆盖经典图、格图、随机图、几何图、community/SBM 和 expander 等生成器，代码采用 BSD-3-Clause。[生成器目录，证据 A](https://networkx.org/documentation/stable/reference/generators.html)；[LICENSE，证据 A](https://github.com/networkx/networkx/blob/main/LICENSE.txt)

| 图族 | 要控制的参数 | 检验的数学现象 |
|---|---|---|
| path / balanced tree | \(n\)、分支数、高度 | 无环基线、直径和远程 anchor；有限传播下界最清楚 |
| cycle / circulant | \(n\)、chord offsets | 最小含环反例、girth、双路径信息 |
| grid / hex lattice | 长宽、periodic | 空间网络、边界/体积比、\(D=\Theta(\sqrt n)\) |
| ring of cliques / caveman | 团数、团大小、桥数 | 强局部相关 + 弱全局 cut，Schur 边界影响 |
| random geometric | \(n\)、维数、连接半径 | 无线/传感器局部通信，坐标和链路距离可用于能耗代理 |
| connected Watts--Strogatz | \(n,k,p\) | 从规则大直径到少量 shortcut 的连续过渡 |
| random \(d\)-regular / expander | \(n,d\)、seed | 有界度、谱隙、快速 mixing；与格图形成对照 |
| SBM / planted partition | block sizes、\(p_{\rm in},p_{\rm out}\) | 区域结构、cut size、分区通信 |
| barbell / lollipop | 团和连接路径大小 | 极端 bottleneck 与慢混合 |
| Barabási--Albert | \(n,m\) | hub、度异质和消息拥塞；不属于 bounded-degree 主定理族 |

每个随机实例保存 generator 名称、NetworkX 版本、完整参数、seed 和最终 edge list；若生成图不连通，不能悄悄重抽，应记录 retry count 或明确取最大连通分量。

### 8.1 三类合成测量模板

1. **相对差 + anchor（严格 pairwise）**

\[
z_{uv}=C_{uv}x_u-D_{uv}x_v+\epsilon_{uv},\qquad
z_i=A_ix_i+\epsilon_i.
\]

可从 scalar difference 到 \(d\)-维 block factors，专门构造每个 \(A_i\) 不满列秩但全局信息矩阵可逆的实例。

2. **一般 1-hop local factor**

\[
z_i=H_i x_{\{i\}\cup N(i)}+\epsilon_i,
\]

用于模拟 power injection 和局部 cluster sensor；需把 factor graph 与通信 graph 分开保存。

3. **非线性 range/bearing**

\[
r_{uv}=\|p_u-p_v\|+\epsilon_{uv},
\]

固定 anchors 消除刚体 gauge；只作为 nonlinear extension，不混入线性精度主表。

### 8.2 可控轴

- \(n\in\{10^2,10^3,10^4,10^5\}\)，平均度和最大度；
- diameter、girth、谱隙、separator/cut size、community strength；
- state block dimension \(d_x\) 和 measurement dimension；
- anchor density、measurement density、局部秩缺失模式；
- edge/unary weights 的动态范围与 \(\kappa(Q)\)；
- SNR、相关噪声 block、outlier fraction；
- 通信边丢包、永久 edge failure、node failure 和消息量化。

测量缺失、通信丢包、节点失败是三个不同概率空间，必须分别报告。所谓 \((\epsilon,\delta)\) 成功也要明确 \(\delta\) 是对随机噪声、随机根节点、随机图还是通信故障取概率。

## 9. 推荐的三层 benchmark suite

### L1：Power-WLS physical core

**Correctness 档：** case14、case118、case300。  
**Topology/provenance 档：** case200_activ、case1354_pegase、case2869_pegase。  
**Scaling 档：** case2000_goc、case10000_goc。

每例冻结：

- DC branch-flow + angle-anchor 严格 edge/unary 版；
- DC injection-inclusive 1-hop-factor 版；
- AC Jacobian-at-truth 版；
- 可选完整 AC WLS；
- dense、sparse-observable、locally-rank-deficient、boundary-starved 四种 placement；
- bus-level、\(k\)-region、degraded 三种 communication view；
- 至少 10 个公开 noise seeds，另有 noise-free correctness seed。

### L2：Sparse locality and conditioning

- usroads 最大连通分量 + 1 anchor：大直径 pairwise WLS；
- usroads 全 56 components + 每分量 anchor：多分量/目标覆盖；
- wathen100：中等条件数 SPD；
- bcsstk18：\(\kappa\approx 3.46\times10^{11}\) 的病态压力。

报告时把 usroads 放“语义合法的 graph WLS”，后两者放“generic precision/solver”；不要把三者不加区分地平均。

### L3：Controlled theorem boundary

用 path/tree/cycle/grid/RGG/WS/random-regular/SBM/barbell/BA，逐轴扫描：

\[
(r,\epsilon,\delta,\text{rounds},\text{message bytes},
\text{anchor density},\kappa,\text{cut size}).
\]

这一层最能检验原创数学结论，因为可以构造一对半径 \(r\) 邻域完全相同、但远端数据改变中央 WLS 的 indistinguishable instances；也能检验谱隙或 inverse-decay 条件何时足以给出近似。

### 可选 L4：Nonlinear transfer gate

只在逐数据集许可核清后加入 g2o/SE-Sync；PMDSE/SimBench 可提供三相 distribution AC SE。L4 的成败不能反向决定线性 iff 定理是否成立。

## 10. 实验指标与发布规则

### 10.1 估计质量

- \(\|P_i\hat x-P_ix_{\rm ref}\|\) 的绝对、相对和 \(R/Q\)-加权误差；
- 节点误差分布、最坏节点、median、95/99 percentile；
- \(F_n(r,\epsilon)\)：半径 \(r\) 内达到误差阈值的节点比例；
- 在 random-root 口径下报告 \(\Pr(e_i(r)\le\epsilon)\)，同时保留 per-node 数据；
- AC 情况另报 objective gap、measurement residual 和是否落入同一局部解。

### 10.2 通信与能耗代理

- rounds / radius；
- 每条消息维度、bytes 和峰值；
- total messages、byte-hops、max-link load；
- 若有几何距离，报告 \(\sum b_e d_e^\alpha\) 的 radio-energy proxy；
- packet loss/retransmission 单列，不用 wall-clock 替代通信开销。

### 10.3 计算

- 每 agent 的时间、峰值内存、factorization/iteration count；
- 中央 sparse direct、CG/MINRES 和选定分布式基线；
- 记录 hardware、线程数、BLAS、solver、容差和 stopping rule。

### 10.4 数据发布

- raw 永不覆盖；derived 每次生成使用 semantic version；
- 上游 CC BY attribution、原论文和修改说明随包发布；
- seed 不只是一个整数：同时保存 RNG family/library/version；
- 所有中央解保存 primal/normal residual，防止把求解器失败当局部算法误差；
- benchmark leaderboard 不把物理层、generic SPD 和非线性层汇成一个无解释总分。

## 11. 本次本地归档

总归档体量 10,740,304 bytes（约 10.24 MiB），所有 `SHA256SUMS` 条目已在 2026-09-13 本机重新计算并逐项比对通过。

### 11.1 PGLib-OPF v23.07

本地目录：`E:\Code\MAS\datasets\raw\pglib_opf_v23.07`

固定下载基址为 `https://raw.githubusercontent.com/power-grid-lib/pglib-opf/v23.07/`，在其后拼接表中的文件名即可得到每个 case 的 immutable-tag URL；完整逐文件 URL 已存于同目录 `SOURCE_URLS.md`。

| 文件 | bytes |
|---|---:|
| pglib_opf_case14_ieee.m | 13,781 |
| pglib_opf_case118_ieee.m | 80,336 |
| pglib_opf_case200_activ.m | 89,620 |
| pglib_opf_case300_ieee.m | 154,407 |
| pglib_opf_case1354_pegase.m | 666,250 |
| pglib_opf_case2000_goc.m | 781,594 |
| pglib_opf_case2869_pegase.m | 1,464,409 |
| pglib_opf_case10000_goc.m | 3,421,442 |

同时保存上游 `LICENSE`、`README.md`、`CHANGELOG.md`，以及本项目的 `SOURCE_URLS.md` 和 `SHA256SUMS`。case 与三个上游 metadata 文件均来自 immutable `v23.07` tag，未修改。

### 11.2 SuiteSparse CC BY 4.0 subset

本地目录：`E:\Code\MAS\datasets\raw\suitesparse_ccby4`

稳定官方 URL 分别为：

- `https://sparse-files.engr.tamu.edu/MM/Gleich/usroads.tar.gz`；
- `https://sparse-files.engr.tamu.edu/MM/GHS_psdef/wathen100.tar.gz`；
- `https://sparse-files.engr.tamu.edu/MM/HB/bcsstk18.tar.gz`。

| 文件 | bytes | SHA-256 前 12 位 |
|---|---:|---|
| usroads.tar.gz | 1,013,706 | 3ee6913382a9 |
| wathen100.tar.gz | 2,316,717 | 2498caf696db |
| bcsstk18.tar.gz | 709,409 | ebd14190b838 |
| LICENSE_CC-BY-4.0.txt | 18,657 | 9ba9550ad484 |

另有 `README.md`、`SOURCE_URLS.md`、`SHA256SUMS`。tarball 未解包重写；Matrix Market header 中的矩阵来源和作者信息完整保留。官方 archive server 的 HTTPS 在本机发生 TLS transport failure，实际下载回退到同一官方服务器的 HTTP endpoint；稳定 HTTPS URL 与最终内容哈希均已记录，避免把传输方式当作来源变化。

### 11.3 明确没有归档的对象

- DPLib：无标准 LICENSE 文件；
- g2o/SE-Sync sample datasets：逐数据集再分发权不清；
- IEEE PES feeder 镜像：旧官方仓库许可证为 TBD；
- TAMU ACTIVSg 直接包：需通过官方表单，本次不绕过；
- SimBench：许可明确但 ODbL/DbCL obligations 和体量/结构超出最小核心；
- 78k PGLib 等大文件：当前无必要，固定 tag 可随时按需获取。

## 12. 下一步建议

1. 先实现一个版本化 `benchmark-builder`，输入固定 PGLib file + YAML/JSON protocol，输出第 5.5 节目录；不要先写分布式算法。
2. 用 case14 做三方一致性：本项目中央解、MATPOWER、JuliaGrid；DC 应达到数值精度一致，AC 对齐量测定义和 reference convention。
3. 冻结 benchmark v0.1：case14/118/300、usroads 小子图、path/cycle/grid/RGG，各含 locally-rank-deficient 与 boundary-starved 配置。
4. 在 v0.1 的 schema、哈希和残差全部稳定后，再扩到 1354/2869/10000 与完整 usroads。
5. 若要采用 DPLib 分区，先向作者请求明确 LICENSE；否则从 CC BY 4.0 的 PGLib 原文件用我们自己的 seeded partitioner 重新生成并公开配置。
6. 投稿前再做一次专门检索，核查 2026 年后是否出现完整 distributed-SE benchmark；当前报告的“不存在”措辞必须保留审计范围限定。

## 13. 主要官方来源

- [PGLib-OPF 官方仓库](https://github.com/power-grid-lib/pglib-opf)
- [PGLib v23.07 README](https://raw.githubusercontent.com/power-grid-lib/pglib-opf/v23.07/README.md)
- [PGLib v23.07 LICENSE](https://raw.githubusercontent.com/power-grid-lib/pglib-opf/v23.07/LICENSE)
- [PGLib v23.07 BASELINE](https://github.com/power-grid-lib/pglib-opf/blob/v23.07/BASELINE.md)
- [PGLib paper / archive report](https://arxiv.org/abs/1908.02788)
- [MATPOWER license](https://matpower.org/license/)
- [MATPOWER 8.1 manual](https://matpower.org/docs/MATPOWER-manual-8.1.pdf)
- [MATPOWER state_est example](https://matpower.org/docs/ref/matpower6.0/extras/state_estimator/state_est.html)
- [TAMU ACTIVSg2000](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg2000/)
- [IEEE PES distribution-feeder prototype](https://github.com/ieee-pes-amps/dtf-dev)
- [SimBench datasets](https://simbench.de/en/download/datasets/)
- [SimBench license](https://github.com/e2nIEE/simbench/blob/develop/LICENSE)
- [DPLib paper](https://arxiv.org/abs/2506.20819)
- [DPLib repository](https://github.com/LSU-RAISE-LAB/DPLib)
- [JuliaGrid documentation](https://mcosovic.github.io/JuliaGrid.jl/stable/)
- [JuliaGrid measurement model](https://mcosovic.github.io/JuliaGrid.jl/dev/tutorials/measurementModel/)
- [PMDSE repository](https://github.com/Electa-Git/PowerModelsDistributionStateEstimation.jl)
- [PMDSE input/measurement documentation](https://electa-git.github.io/PowerModelsDistributionStateEstimation.jl/v0.4/input_data_format/)
- [SuiteSparse Matrix Collection license/about](https://sparse.tamu.edu/about)
- [SuiteSparse usroads](https://sparse.tamu.edu/Gleich/usroads)
- [SuiteSparse wathen100](https://sparse.tamu.edu/GHS_psdef/wathen100)
- [SuiteSparse bcsstk18](https://sparse.tamu.edu/HB/bcsstk18)
- [g2o repository](https://github.com/RainerKuemmerle/g2o)
- [g2o file format](https://github.com/RainerKuemmerle/g2o/wiki/File-Format)
- [SE-Sync repository](https://github.com/david-m-rosen/SE-Sync)
- [NetworkX graph generators](https://networkx.org/documentation/stable/reference/generators.html)

## 14. 最终判决

**Benchmark availability：PARTIAL。** 网络、分区工具、测量生成器和数值矩阵分别存在，但本次审计范围内没有一个广泛采用的包同时冻结完整 distributed Gaussian/WLS tuple。置信度：**中高（0.85）**；受限于无法形式排除零散论文仓库或非英语/非公开数据。

**Dataset readiness：SOLVED for Phase 2 core.** PGLib v23.07 精选案例和 SuiteSparse 精选矩阵已按许可、来源与哈希规范归档；三层 suite 和测量/通信生成协议已明确。置信度：**高（0.95）**。

**DPLib reuse：NO-GO pending license clarification.** 可引用其论文并借鉴字段，不应在本项目中再分发当前仓库内容。置信度：**高（0.95）**。
