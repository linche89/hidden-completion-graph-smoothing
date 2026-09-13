# 本地可扩展 WLS 实现实录

> 完成日期：2026-09-13  
> 范围：直接稀疏合成、PGLib MATPOWER 解析与 DC-WLS、CPU 正确性/规模跑、项目隔离 CuPy GPU 路线。没有使用或声称 MPI。

## 1. 已交付代码

| 文件 | 作用 |
|---|---|
| `experiments/scalable/wls_core.py` | `WLSInstance`、直接 COO/CSR 的块 synthetic WLS、`J=H^TR^{-1}H`、抽样根转置求解与 `b/raw-z/whitened-z` 尾部 |
| `experiments/scalable/pglib_dc.py` | 最小 PGLib/MATPOWER v2 数字 parser、参考母线处理、稀疏 DC-WLS 构造 |
| `experiments/scalable/cpu_benchmark.py` | 六档 PGLib CPU correctness/scale、时间/RSS 分项、抽样尾部和 DKW 元数据 |
| `experiments/scalable/gpu_sparse_benchmark.py` | 同一 CSR、同容差的 CPU/CuPy SpMV 与 Jacobi-CG；预处理、H2D、warmup、算法、D2H 分开 |
| `experiments/scalable/test_correctness.py` | 14-bus parser/参考母线/holder/support/全半径零尾，以及少锚点块 WLS 单测 |
| `experiments/scalable/requirements-gpu.txt` | 项目隔离 GPU 环境的精确依赖 |

静态搜索确认新模块没有 `toarray`、`todense`、`inv` 或完整逆构造。`H` 和 `J` 直接以 triplet 形成 COO/CSR；允许的 dense 对象只有长度 `N` 的向量、单个抽样逆行以及小型 `d×d` 局部块。

原 `experiments/locality_benchmark.py` 保留作历史 smoke，没有被覆盖；新的大规模入口全部在 `experiments/scalable/`。

## 2. 模型语义

### 2.1 PGLib parser

parser 只解析已归档 PGLib 文件中的数字 `mpc.baseMVA`、`mpc.bus` 和 `mpc.branch`，不是任意 MATLAB 解释器。它会拒绝缺块、ragged matrix、重复 bus id 和不合法参数。八个本地 case 的 SHA-256 均重新核对并与 `datasets/raw/pglib_opf_v23.07/SHA256SUMS` 一致。

只采用 status 为 1 且 `|x|>10^{-12}` 的支路。tap 的 MATPOWER 零值按 1 处理，phase shift 被保存为已知 offset。通信图由 active modeled branches 形成；每个连通分量优先选 type-3 slack，否则选最小内部节点为参考母线，并从未知角状态中消去。

### 2.2 DC-WLS 测量

- **from-end branch flow**：`baseMVA/(x·tap)(θ_f-θ_t-shift)`；默认数值只由 from bus 持有；
- **bus injection**：所有相邻 DC branch flows 的代数和；由该 bus 持有；
- **bus angle**：非参考母线角度；由该 bus 持有。

观测中保留相移产生的 `known_offset`，求 WLS 前使用 `z-offset`。三种标准差分别进入对角 precision，中心系统是

`J=H^T diag(precision)H`, `b=H^T diag(precision)(z-offset)`。

每条测量分别保存：

1. `measurement_variable_support`：参考角消去后真正出现非零列的节点；
2. `measurement_physical_support`：原物理因子涉及的节点；
3. `measurement_holders`：第 0 轮拥有测量数值的通信节点。

因此不会再把“因子依赖谁”和“数据放在哪里”偷换成同一个集合。CPU 跑采用 `flow_holder=from`，不是默认假设边测量在两端复制。

PGLib-OPF 多数 case 的 `Va` 是 flat start，所以 benchmark 用固定 seed 生成 smooth synthetic DC angle truth；这不是声称 PGLib 提供了真实状态轨迹。测量拓扑和电气参数来自 PGLib，状态/测量位置/噪声由本脚本生成。

### 2.3 三种尾部

对抽样标量根解 `J^Ty=e_i`，不形成逆：

- `state_rhs_tail_l2=||e_i^TJ^{-1}Q_out^b||_2`；
- `physical_raw_z_tail_l2=||e_i^TJ^{-1}H^TR^{-1}Q_out^z||_2`；
- `whitened_z_tail_l2=||e_i^TJ^{-1}H^TR^{-1/2}Q_out^z||_2`。

潮流/注入以 MW、角度以 rad 表示，直接对 physical raw coordinates 取欧氏单位球会混合单位。因此跨测量类型的 `ε/δ` exceedance 只采用无量纲 whitened tail；physical raw tail 仍保留供明确物理扰动模型使用。measurement mask 按 holder 是否进入根的通信球计算。

## 3. CPU 结果

命令固定 `OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`，使用 Python 3.13.5、NumPy 2.3.3、SciPy 1.17.1。中心解为 SuperLU 分解加三角解；尾部用 32 个均匀、有放回根样本和半径 `{0,1,2,4}`。以下为单次 correctness/scale run，不是严谨 microbenchmark；RSS 是显式检查点最大值，不是 OS peak。

| case | buses / state | measurements | `nnz(H)` / `nnz(J)` | parse / build (s) | factor / solve (s) | 32 roots × 4 radii tails (s) | case total (s) | checkpoint RSS MiB | relative residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE 14 | 14 / 13 | 26 | 61 / 73 | .0006 / .0015 | .00007 / .00001 | .0010 | .0041 | 67.3 | `4.56e-16` |
| IEEE 118 | 118 / 117 | 226 | 518 / 707 | .0010 / .0017 | .00017 / .00001 | .0049 | .0126 | 68.0 | `1.19e-14` |
| IEEE 300 | 300 / 299 | 499 | 1,129 / 1,741 | .0018 / .0030 | .00030 / .00002 | .0104 | .0160 | 68.9 | `1.15e-13` |
| PEGASE 1354 | 1,354 / 1,353 | 2,395 | 5,277 / 7,359 | .0077 / .0139 | .00113 / .00006 | .0448 | .0685 | 73.1 | `1.08e-13` |
| PEGASE 2869 | 2,869 / 2,868 | 5,442 | 12,020 / 16,300 | .0173 / .0307 | .00269 / .00015 | .0983 | .1510 | 79.0 | `1.55e-12` |
| GOC 10000 | 10,000 / 9,999 | 16,202 | 35,685 / 49,557 | .0487 / .1113 | .01027 / .00038 | .3001 | .4745 | 96.2 | `8.72e-14` |

六个 case 全部 `status=ok`，抽样 whitened tail 对递增半径非增。14 和 118 的全直径检查中 `b`、physical raw-`z` 与 whitened-`z` 球外尾部均精确为零。机器可读结果在 `experiments/results/pglib_cpu_final.json`。

这些时间很短的原因是当前 DC-WLS 状态每 bus 只有一个标量且 `J` 极稀疏；不能外推为 AC 块状态、稠密噪声 precision 或高 fill-in 图也同样便宜。

## 4. 少锚点 synthetic 验证

新 synthetic generator 的节点 anchor 比例和 anchor rank 是独立参数。单测使用 40 节点 cycle、块维 2、5% 全秩 anchor，其余 95% 节点无本地 anchor；全局 `J` 仍可解，残差低于 `1e-10`。这才真正离开“每个 `A_i` 满列秩”，不同于旧 smoke 仅把每点 anchor 权重缩小但保持正秩。

边因子使用小型正交 block `x_u-Q_ex_v`；这些 `d×d` 局部矩阵可以 dense，但全局 `H/J` 始终 sparse。

## 5. 项目隔离 GPU 路线与实跑

官方 CuPy 14.2 安装文档列出 Windows wheel、Python 3.13 和 CUDA 12.8 支持，CUDA 12.x wheel 包名为 `cupy-cuda12x`：[CuPy installation](https://docs.cupy.dev/en/stable/install.html)。`cupyx.scipy.sparse.linalg.cg` 官方接口支持 SPD sparse matrix/LinearOperator：[CuPy CG](https://docs.cupy.dev/en/latest/reference/generated/cupyx.scipy.sparse.linalg.cg.html)。

没有改系统 Python。隔离环境位于 `experiments/.venv-gpu`，freeze 为：

```text
cuda-pathfinder==1.8.1
cupy-cuda12x==14.2.0
networkx==3.6.1
numpy==2.3.3
psutil==7.2.2
scipy==1.17.1
```

CuPy 自检识别：RTX 5080、compute capability 12.0、cuSPARSE available；CuPy wheel linked runtime 12.9，系统 local toolkit 12.8，driver API 13.2。CSR identity SpMV 自检正确。

真实 GPU benchmark 使用 PGLib GOC 10000 转换得到的同一 `9999×9999`、49,303-nnz SPD `J`。CPU/GPU 都用 Jacobi-preconditioned CG、`rtol=1e-8`、`atol=0`、相同 RHS；转移、warmup 和算法时间分开：

| 指标 | CPU | RTX 5080 |
|---|---:|---:|
| 1,000 次 SpMV algorithm time | 0.0178 s | 0.0393 s |
| Jacobi-CG algorithm time | 0.2406 s | 1.8071 s |
| iterations | 7,470 | 7,476 |
| final relative residual | `9.28e-9` | `8.68e-9` |

另有 H2D matrix+vectors 0.0578 s、首次 warmup 0.0870 s、D2H 约 0.0001 s；SpMV CPU/GPU 相对差 `1.11e-16`。结果在 `experiments/results/pglib_gpu_case10000_final.json`。

结论不是“GPU 加速成功”，而是 **GPU 数值路线已打通，但该矩阵太小、每次只有约 49k nnz，且 CuPy CG 是许多细粒度 Python 驱动迭代，所以 RTX 5080 明显慢于 CPU**。GPU 应用于更大的 bounded-degree synthetic 图、batched RHS、Chebyshev/SLQ，而非用这次结果宣称 speedup。

负结果也保留：无预条件 10k CG 在 CPU/GPU 均达到 20,000 iteration 上限，残差约 `2.5e-4/2.7e-4`，见 `pglib_gpu_case10000.json`；它不能作为同容差性能结果。

## 6. 验收与剩余边界

已满足：

- 全局 `H/J` 直接 COO/CSR，无 dense 构造、无完整逆；
- 14/118/300 必跑及 1354/2869/10000 扩展全部完成；
- DC reference、flow/injection/angle、noise、offset、support 与 holder 语义显式；
- CPU residual、转换/构造、factor/solve、tail 与 RSS 分项；
- 项目隔离官方 GPU 路线、真实 PGLib SpMV/CG、相同容差与 CPU 对照；
- 明确 `mpi=false`，没有 MPI 性能主张。

仍不能声称：

- 当前最小 parser 等价于 MATPOWER/pandapower 的完整 AC/DC 引擎；
- 单次 wall time 是稳定性能统计；
- 32 根样本给出窄的 `δ` 置信带（95% DKW 半宽仍约 0.240）；
- GPU 对该问题已经有速度优势；
- 本机实现验证了 MPI、多机、无线 joule 或 AC 非线性状态估计。

下一步最有价值的是：用新 sparse synthetic generator 产生 `N≥10^6` 的 bounded-degree SPD 系列，比较 batched Chebyshev/SLQ 的 CPU/GPU 吞吐；同时把 PGLib root sample 增加到至少 185，才获得每固定半径约 ±0.10 的 95% DKW 带。
