# 局部通信网络状态估计：研究尽调工作区

## 研究目标

判断下述问题的哪些部分已经被解决、哪些只是算法特例、哪些仍可能构成真正的理论缺口：

> 对网络线性/高斯状态估计，在节点仅与邻居通信、只输出自身状态且允许误差与少量节点失败时，存在低轮数、低通信、低能耗、低本地计算估计器的最宽充要条件是什么？

本项目首先做可证伪的文献尽调，不预设结果具有原创性。“本次未检索到”不会写成“首次提出”。

## 主要交付物

- `../paper/main.tex`：scope-frozen 正式论文稿；一个几何主定理、三个应用推论，含完整证明、
  sharp Schatten 边界、exact SDP 与保守原创性表述。
- `../hidden_completion/`：finite-\(h\)/unbounded 精确场景枚举与 spectral SDP 求解器；支持
  zero pattern、\(r\)-hop、shared coefficient，并输出 active worst-case scenarios。
- `p1_p2_delivery.md`：本轮论文落地与精确求解器的逐项验收、复现命令和测试快照。
- `main_theorem_decision_v2.md`：**已冻结的论文主线**；finite-\(h\) 精确 hull、完整输入
  all-Schatten 归约、exact SDP、\(q/(q+h)\) 收敛率，以及无界 partial-partition sharp limit。
- `agent_reports/finite_hidden_budget_theorem.md`：有限隐藏节点预算主定理的完整证明、
  exposed-vertex characterization、路径 converse 与解析收益例。
- `agent_reports/direct_smoother_novelty_audit.md`：截至 2026-09-14 的原始文献专项查重与
  审稿风险边界。
- `survey_2011_2026.md`：2011–2026 主尽调报告，含问题分层、最宽基础充要条件、直系与跨领域撞车、资源边界及拟议研究缺口。
- `formal_problem_spec.md`：后续写定理与论文时使用的统一问题定义。
- `next_research_program.md`：下一阶段定理、反例、算法和下界的研究顺序；使用普通术语，不为候选条件另造名称。
- `phase2_theory_experiments.md`：数学、公开测试集、大规模模拟与 CPU/GPU/MPI 实现的一体化方案及当前执行状态。
- `phase3_execution.md`：数学原创性审计与本地稀疏 CPU/GPU 实验的执行和验收规范；近期不承诺集群结果。
- `phase3_results.md`：第三阶段数学判决、一般块 exact SDP 的 Petersen-lemma 碰撞、PGLib CPU 与百万状态 GPU 实跑结果，以及下一主定理的硬门槛。
- `corpus/manifest.tsv`：本地全文库清单、页数、SHA-256 与文本索引。
- `agent_reports/theory_phase2.md`：精确局部 minimax、Schur 双向界、带吸收随机游走表示与不可能性反例的数学底稿。
- `agent_reports/datasets_phase2.md`：公开测试集、许可、格式、下载与测量生成方案的逐项审计。
- `agent_reports/hpc_experiments_phase2.md`：本机能力、稀疏 CPU/GPU/MPI 路线和大规模实验矩阵。
- `agent_reports/phase2_adversarial_notes.md`：对数学量词、实验实现、统计解释和 HPC 声称的独立反方审计。
- `agent_reports/theorem_novelty_audit.md`：逐定理、逐量词的原创性碰撞审计；旧定理包判 NO-GO，并把一般块 exact SDP 标为待专项核验；最终判决见下一项。
- `agent_reports/block_robust_boundary_attack.md`：一般实 block 鲁棒局部估计的 lossless SDP、有限 completion 对偶/KKT 与二维端点反例。
- `agent_reports/sdp_exact_novelty_audit.md`：将 exact SDP 逐变量归约到 single-full-block Petersen lemma 的专项审计。
- `agent_reports/structured_completion_main_theorem_search.md`：一般 fixed pairwise factors 的 matroid-flat reduction、uniform/normalized-heterogeneous hidden graph 的 exact partial-partition theorem、cut-size FPT LMI、degree-2 sharp converse 与 Petersen strict gap；当前最强候选主定理。
- `theory/candidate_main_theorem.md`：任意矩阵区间、退化区间与 cut 维数压缩的完整正确性底稿；明确作为已知 full-block 工具而非原创主定理。
- `agent_reports/scalable_experiments_implementation.md`：PGLib 14--10,000 节点稀疏 DC-WLS 的实现语义、CPU/GPU 结果和复现边界。
- `agent_reports/local_gpu_scale.md`：约百万状态的 CSR SpMV/SpMM 与固定轮局部迭代 CPU/GPU 实跑。
- `../datasets/raw/`：已核验许可与 SHA-256 的 PGLib v23.07 和 SuiteSparse 原始数据子集。
- `../experiments/`：不显式形成完整逆矩阵的 CPU 正确性骨架、配置与运行结果。

## 必须区分的四种“局部”

1. **逐轮局部**：每轮只与一跳邻居通信，但轮数可达网络直径；这仍可实现全网泛洪。
2. **半径局部**：轮数 `r` 不随网络规模增长，节点输出只依赖 `r`-邻域。
3. **资源局部**：每条消息、节点内存和节点计算量只依赖局部维数/度数，而不随全网规模增长。
4. **知识局部**：算法只知道局部模型与局部拓扑，不预先持有全局矩阵、谱界或全局分解。

只有第 2–4 项同时受到约束时，“局部算法”才对应我们真正关心的可扩展性问题。

## 核心判定轴

每篇论文按以下字段记录：

- 模型：静态/动态、线性/非线性、Gaussian/WLS、节点变量/全局公共变量；
- 目标：每节点全局重构或仅重构自身分量；
- 信息：本地测量、边测量、先验、全局模型知识；
- 拓扑：树、有环图、连通/强连通、树宽、谱隙、几何/平移不变；
- 保证：精确有限时、渐近精确、有限半径近似、均方/高概率；
- 条件性质：充分、必要、充要，针对“问题本身”还是针对“某个算法收敛”；
- 资源：轮数、bit 数/bit-hop、能耗、消息维度、时间和空间复杂度；
- 失败语义：所有节点、随机失败、恶意节点，或至少 `1-δ` 比例节点成功；
- 反例/边界：结论不覆盖的实例与隐含全局代价。

## 工作流

1. **全文本地化**：65 篇公开全文已归档，并转换为可由 `rg` 检索的 UTF-8 文本；`08_robust_control_boundary` 保存 Petersen lemma/QMI 碰撞来源，`09_network_completion` 保存 DPP/随机森林、Kron reduction 与 partition-polytope 近邻；见 `corpus/manifest.tsv`。
2. **直系谱系审计**：先审读原始论文及 2015、2018、2020、2022、2023、2026 后续工作，确定“无环”和“每个局部矩阵满列秩”分别被放宽到哪里。
3. **跨领域撞车审计**：Gaussian BP/空间混合、稀疏逆衰减、图滤波、有限时线性变换、functional observability、localized control、通信复杂度和率失真。
4. **统一数学翻译**：把各领域结论翻译成同一对象——中心估计算子（在线性 Gaussian 模型中由信息矩阵逆矩阵构成）能否由有限半径数据逼近。
5. **缺口判决**：分别说明已经解决、只解决一部分、看起来仍开放或问题定义不完整，并给出证据和置信度。
6. **研究立项门槛**：只有当候选定理包含可检验结构条件、构造算法和匹配下界，并且不只是“最优局部风险的定义式”时，才进入原创定理阶段。

## 当前初步警报

- “树上有限步精确”已经是经典结果，不应作为新贡献。
- “每个节点本地可观/满秩”早已可放宽为整体可估计性；标量边测量还有更弱的单锚点条件。
- “有环图也能算出精确解”在允许 `O(n)` 轮传播或全局消元时已经被解决。
- “局部近似等价于相关性衰减”在有限字母 Gibbs 模型中已有形式化充要结论。
- “稀疏正定矩阵的逆可被局部多项式近似”是数值线性代数与图滤波的成熟结论。
- 可能仍未统一解决的是：任意有环 Gaussian/WLS 图族、仅局部模型、`1-δ` 节点成功，以及轮数—bit-hop—能耗—计算的匹配上下界。
