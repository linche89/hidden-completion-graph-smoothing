# Phase 2 二次反方审计：必须修正项

> 只读审计对象：`research/phase2_theory_experiments.md` 与 `experiments/locality_benchmark.py`。未修改二者。  
> 总结：当前 smoke 的 Schur 代数和 `R=I,Q=0` 下的转置求解均正确；但有 9 个会造成指标偷换、错误外推或实验无法验证核心主张的问题。

## M1. Schur 公式必须显式假设 `J≻0`（至少 `J_OO` 可逆）

主报告先写 `Q⪰0`，但没有明确 `J` 正定，随后直接使用 `J_OO^{-1}`。若 WLS 仅半正定、目标虽可估但全状态不唯一，该逆和唯一中心解都未必存在。代码的生成器因每点有正权锚点而使 `J≻0`，所以 smoke 没暴露此问题。

必须改成：

> 本节先假设 `J≻0`，于是每个主子矩阵 `J_OO≻0`，Schur 公式成立。秩亏但目标可估的情形需另用商空间/约束规范化或广义逆处理，不能直接套此公式。

代码中的 `schur_oracle` 公式本身正确：

`S=J_BB-J_BO J_OO^{-1}J_OB`，`g=b_B-J_BOJ_OO^{-1}b_O`。

它形成 dense `J_OB` 解和 dense Schur，只能是小规模 oracle，不能被描述成大规模局部实现。

## M2. 默认 smoke 仍满足“每个本地 `A_i` 满列秩”，因此尚未验证核心放宽

代码为每个节点添加

`sqrt(anchor_weight*scale) I_d`。

smoke 中 `anchor_weight=.2`，坏区 `scale=1e-4>0`，所以 **所有节点的局部锚点块仍为满列秩**；缩小奇异值不等于秩亏。当前实验只能验证有环和病态性，不能支持“局部满秩不必要”的经验结论。

必须在正式生成器中独立加入 `anchor_fraction/anchor_nodes` 和 `anchor_rank∈{0,…,d}`，至少包含：只有一个全秩锚点、其余 `A_i=0`；每点只测不同低维方向但全局目标可估；全局状态秩亏但指定 `C_i x` 可估。每个实例应报告 `rank(A_i)` 分布与全局/功能可估计性，而不是只报告 `κ(J)`。

把 `bad_anchor_scale` 手动设为 0 可让坏区 `A_i=0`，但现有默认 smoke 没有这样做，且非坏区仍全部满秩。

## M3. `measurement_support` 被错误地兼任 measurement-holder map

代码用因子依赖的状态节点 `[u,v]` 判定边测量只要与 `B_r(i)` 相交就可见。这在当前特殊约定“同一条边测量在两个端点初始时都完整可用”下正确；但一般 `supp(H_k)` 与“数据 `z_k` 由哪个通信节点持有”是两个不同对象。主报告自己要求保存 holder map，代码尚未落实。

必须把二者拆成：

- `measurement_variable_support[k]`：第 `k` 个测量依赖哪些状态块；
- `measurement_holders[k]`：哪些通信节点在第 0 轮持有该数值。

球内可见性应按 `holders[k]∩B_r(i)≠∅`，而不是按因子 support。若边测量只存于一个端点、独立 edge agent 或中心 PMU，当前 mask 会错误低估 raw-`z` 尾部。当前结果必须注明“双端复制 holder”语义。

## M4. raw-`z` tail、`b` tail、实例误差与 Bayes MSE 不能共用“失败比例”名称

代码正确计算了当前模型中的两个不同量：

- `state_rhs_tail = ||E_iJ^{-1}Q_out^b||`，输入域是节点信息向量 `b` 的欧氏单位球；
- `raw_wls_tail = ||E_iJ^{-1}H^TQ_out^z||`，输入域是原始 `z` 的欧氏单位球，因为当前 `R=I`。

一般模型应为 `(R^{-1}HY_i)^T`，且先验 `Q` 对 `J` 和 `b` 的贡献也需写清。两个 tail 的单位、输入归一化和通信 holder 都不同，不能互换。

代码的 `failure_curves/sample_failure_fraction` 实际是“抽样节点中 raw-`z` worst-case operator tail 超过阈值的比例”，不是某次状态估计失败率、随机数据失败概率或 Bayes MSE 超阈值概率。必须重命名或在输出 schema 中写全 `raw_z_operator_tail_exceedance_fraction`，并在主报告为每张图声明输入范数和 `ε` 的单位。

此外，tail 在“全局系数预装、无限带宽收集球内 raw data”模型下是最佳局部 worst-case 误差，既是 lower bound 也可由截断中心系数达到；主报告只称“理论下限”会掩盖这个量词。

## M5. Chebyshev 的传播声明只在当前 `prop(J)=1` 时成立，且当前谱界不是认证界

代码直接令 `degree=radius`。这只因当前 pairwise factor 与通信边一致，使 `prop(J)≤1`；一般有限作用距离 `ρ` 时，degree `k` 需要 `kρ` 通信半径。若 factor graph 与 communication graph 不同，必须按相对于通信度量的 `prop(J)` 换算。

`eigsh(...,tol=1e-5)` 返回的是近似极端特征值，不保证给出包住整个谱的下/上界。将其直接映射到 `[-1,1]` 可以做经验 polynomial，但不能作为带理论误差保证的 Chebyshev 证书；病态时轻微低估谱区间就可能导致外插和不稳定。认证实验应采用构造时已知界、带 Ritz residual 的外包界或保守 padding。

当前 `chebyshev_global_relative_error` 只是在一个随机 RHS 上相对中心解的实例误差，不是多项式的 operator-norm certificate，也不是节点 `δ` 风险。主报告不能用该曲线验证谱充分条件。

最后，代码的 bit-hop proxy 没计入全局极端特征值计算/谱界广播；若谱界不是模型先验，必须把它列为 global preprocessing，而不能把 `r·2|E|·d·bits` 称为总通信。

## M6. DKW 使用本身正确，但目前只给“固定图、固定半径”的条件置信带

代码对固定实例从节点集合均匀、有放回抽根，重复根也按重复样本计数；此时标准 DKW 对该实例节点 tail 的经验 CDF 有效，尾分布 `P(e>ε)` 也可用同一半宽。

但同一组根被用于多个半径。每个半径各自标 `α=.05` 只给 marginal 95% 覆盖；若正文把整条 `(r,ε)` 曲面称为同时 95% 置信区域，必须对半径数作 Bonferroni/union-bound 调整，或证明联合经验过程界。该 DKW 还条件于已经生成的单张图，不覆盖 graph seed、噪声参数或随机系数的 ensemble uncertainty；图间重复需另报。

## M7. 当前 `cover` 不是 paired-cover converse 生成器，且补连通步骤会破坏 covering property

代码只生成一个随机 lift，没有生成两个 rooted `r`-ball 明确同构而中央答案不同的成对实例；若随机 lift 不连通，代码又任意给分量加边，这些新边不来自 base lift，会破坏 covering map。注释“local balls can agree”不能作为不可区分性验证。

正式 paired-cover 测试必须保存 covering map、选定 roots、rooted colored/coefficient ball 的显式同构证书，并保证本地数据也相同；两个 cover 的远端闭合或锚点再有控制地不同。无法连通的 lift应重采样或按理论允许保留分量，不能事后任意加边仍称 cover。

## M8. 当前生成器是 `O(n²d²)` dense 临时构造，主报告的 `10^4–10^5` 路线不能由它支持

每个因子先建立宽度 `nd` 的 dense block，随后 `vstack`，最后才转 CSR。即使最终 `H` 稀疏，临时内存和写入量仍是二次的。当前 skeleton 只能用于 smoke，不能用其成功推断“数万节点已可运行”。必须先改为 COO/CSR triplet 直接生成。

主报告阶段一还写 `n=10²…10⁴` 时“对所有节点计算精确中心增益行”。这需要约 `nd` 个转置 RHS；即便不保存完整逆，总计算量也可能不可接受，若保存块行则重新接近显式逆的存储。必须改成：低端小规模全节点；规模升高后随机抽根/分层抽根，并报告 DKW/有限总体误差。不能用“高精度替代量”掩盖 exact tail 与 probe surrogate 的差异。

## M9. 随机探针不能未经证明直接给 `2→2` 块行分位数，功率/能耗能力也被过度承诺

主报告说平均量或分位数可“随机抽根并用随机探针估计”。Hutchinson/Hutch++ 自然估计 trace/Frobenius 平均量；它不会自动给每个节点的块行 `2→2` norm 或其精确 quantile。分位数应由抽样根的逐根块转置求解估计，或为 sketch 误差给单独浓缩界。必须把 root sampling 与 matrix probing 两种误差源分开。

同样，“单机可读取 CPU/GPU 功率计数器”未经本机审计验证。RTX 可由 `nvidia-smi` 采瞬时功率，但 Windows 上 Ryzen package energy、采样频率、权限和积分误差未确认。必须改成“若硬件计数器和权限可用则采样并校准”；否则只报 wall time/bit-hop，不得声称已有实现能耗测量。

## 核查通过、无需修改的核心实现

- `schur_oracle` 在当前 `J≻0` 模型中的左右端消元符号正确。
- `J^TY_i=E_i^T` 后，当前 `R=I` 时 `(HY_i)^T=E_iJ^{-1}H^T`，raw-`z` 块行公式正确。
- 当前“双端均持有边测量”语义下，球内任一端点可见即保留测量的 mask 正确。
- `np.linalg.norm(d×m block,ord=2)` 给块行 `2→2` 范数，而非误用 Frobenius。
- 有放回均匀抽根时 DKW 半宽公式 `sqrt(log(2/α)/(2m))` 正确。
- 当前 pairwise 图模型中 `prop(J)=1`，所以 degree `r` 多项式的传播半径不超过 `r`；问题只在不能外推到一般 factor/communication graph。

修完 M1–M9 后，Phase 2 报告才可以把 smoke 定位为“理论量的接口验证”，而不是对局部秩放宽、大规模扩展或多数节点定理的实证支持。
