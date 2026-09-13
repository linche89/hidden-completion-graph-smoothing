# Phase 2：HPC 与理论验证实验专项

> 审计日期：2026-09-13，Asia/Shanghai  
> 审计原则：只读检查本机硬件和软件；没有安装、卸载或更新任何系统组件。  
> 本阶段产物：环境能力边界、理论驱动实验设计、CPU/GPU/MPI 分层方案、复现与正确性规范，以及一个秒级 CPU smoke benchmark。

> 后续修订说明：二次审计后，主线已把 smoke 配置改为只有 20% 节点具有本地满秩锚点，并把测量的变量支撑与初始持有者分开；三项正确性检查重新运行后仍通过。下文 6.1 节的具体数值来自修订前“每点有弱锚点”的首次运行，只保留为环境与接口记录，不作为修订后配置的科学结果。

## 0. 执行结论

这台工作站现在就足以开展 **CPU 版中等规模稀疏 WLS、图族扫描、抽样根节点 Green/WLS 块行尾部、局部截断/多项式算法与 `δ` 分位曲线**。它不缺 CPU 和内存，当前主要缺口是软件栈：Windows Python 的 PyTorch 是 CPU build，没有 CuPy、mpi4py、PETSc、CHOLMOD/PyAMG；本机 Windows 和 WSL 都没有可调用的 MPI。RTX 5080 与 CUDA toolkit 本身可用，但 Python GPU 稀疏计算尚未接通。

因此建议按三层推进：

1. **现在**用 SciPy/NetworkX 做正确性、小中规模定理验证和参数空间筛查；
2. **GPU 环境就绪后**把有限传播多项式、批量随机探针和 batched SpMV 放到 RTX 5080；
3. **集群阶段**用 PETSc/Trilinos + MPI 做分布式 SpMV/迭代解和强弱扩展，作业数组负责图族、seed、病态参数的外层并行。

实验不能只画 runtime。主科学对象应是

`e_i(r)`、`F_n(r,ε)`、`q_{r,δ}`、局部低谱质量、边界 Schur 影响和 cover indistinguishability gap，

然后才把 rounds、message dimension、bit-hop、时间与内存作为实现代价叠加上去。

## 1. 本机只读审计

### 1.1 硬件与操作系统

| 项目 | 实测 | 判断 |
|---|---:|---|
| OS | Windows 11 Pro 64-bit, 10.0.26100, build 26100 | 原生 SciPy/CUDA 可用；生产 MPI 更适合集群或单独配置的 WSL/Linux 环境 |
| CPU | AMD Ryzen 9 9950X3D，16 cores / 32 logical processors | 很适合实例/seed 外层并行和共享内存 SpMV；需防止进程 × BLAS 线程过度订阅 |
| RAM | 125.65 GiB，总审计时约 82.52 GiB 空闲 | 足够中等规模稀疏直接解和大规模矩阵-free 迭代；不能据此忽略 fill-in |
| Page file | 8 GiB，审计时约 179 MiB 使用 | 不应把 swap 当算法内存；生产作业应在接近 RAM 上限前停止 |
| GPU | NVIDIA GeForce RTX 5080，16,303 MiB；compute capability 12.0 | 适合多 RHS SpMV、Chebyshev/Lanczos、随机探针；16 GiB 不适合无约束稀疏直接分解 |
| NVIDIA driver | 596.36；驱动报告 CUDA 13.2 | 可向后运行 CUDA 12.8 构建；正式实验仍需记录 driver/runtime 双版本 |
| CUDA toolkit | 12.8, nvcc 12.8.61 | 原生 Windows 工具链存在；cuSPARSE DLL、header 和 import library 均已找到 |
| 主工作盘 | E: Samsung SSD 990 EVO Plus NVMe，约 3.73 TiB，约 1.57 TiB 空闲 | 适合结果、稀疏矩阵和 checkpoint；不要把海量小文件作为主结果格式 |
| WSL2 | Ubuntu，32 logical CPUs；可见内存约 61 GiB + 16 GiB swap；可见 RTX 5080 | 当前 WSL 内存上限约为宿主的一半，做集群前原型时要按 61 GiB 预算 |

Win32 `AdapterRAM` 对显存的报告发生 32-bit 截断，不能采用；这里使用 `nvidia-smi` 的 16,303 MiB。

### 1.2 语言、编译器与数值库

| 栈 | 实测状态 | 能做什么 / 缺什么 |
|---|---|---|
| Python | 3.13.5, `C:\Python313\python.exe` | 可直接运行本阶段 skeleton |
| NumPy / SciPy | NumPy 2.3.3, SciPy 1.17.1 | 有 CSR/CSC、SuperLU、`spsolve/splu/cg/minres/gmres/eigsh/lobpcg` |
| BLAS/LAPACK | SciPy OpenBLAS 0.3.30, 64-bit integer, `MAX_THREADS=24` | CPU dense/block 运算可多线程；默认线程数未显式固定，benchmark 必须记录并控制 |
| NetworkX / pandas / psutil | 3.6.1 / 3.0.0 / 7.2.2 | 图生成、结果表和 RSS 采样可用 |
| Numba / joblib | 0.62.1 / 1.5.3 | 可优化生成器与做进程级外层并行 |
| PyTorch | 2.10.0+cpu，`torch.cuda.is_available()==False` | **当前不能用 RTX 5080**；不要看到安装了 torch 就把实验标成 GPU |
| GPU Python | CuPy、JAX 均缺失 | Python GPU SpMV/稀疏解当前不可用 |
| 稀疏扩展 | 无 `scikit-sparse/CHOLMOD`、PyAMG、PETSc/petsc4py、pymetis | 缺多重网格、分布式稀疏解和高质量图划分接口 |
| 数据/绘图 | matplotlib、h5py、pyarrow、zarr 均缺失 | 当前 skeleton 输出 JSON；生产结果存储和绘图环境需另行固定 |
| Power-grid 数据 | pandapower 缺失 | 现阶段做合成 WLS；真实电网 benchmark 需以后单独准备环境/数据 |
| C/C++ | 当前普通 PowerShell 找不到 `cl`，但 VS Community 2022 17.14 完整安装；MSVC 14.44 路径存在 | 从 VS Developer Prompt 可用 MSVC；普通 shell 需初始化环境。另有 Strawberry MinGW GCC/G++ 13.2 |
| Build tools | CMake 3.29.2、Ninja 可见 | 原生 C++ 小程序可构建；CUDA+MSVC 应用 Developer Prompt |
| Julia / MATLAB | 未发现 | 本阶段不把它们列为可复现依赖 |
| R language | 未发现；PowerShell 的 `R` 名称不是 R runtime | 不可作为统计绘图后端 |
| MPI | Windows 无 `mpiexec/mpicc/mpicxx`；WSL 也未发现 `mpiexec` | 本机当前不能做真正 MPI 扩展实验 |

WSL 的系统 Python 是 3.12.3；当前 shell 中 `python3` 优先命中 Anaconda Python 3.11.7，后者有 NumPy 1.26.4、SciPy 1.11.4、NetworkX 3.1。另有 GCC/G++ 13.3，但没有确认到 CMake、Ninja、MPI、nvcc 或 Python GPU 数值栈。WSL 的 `nvidia-smi` 可见 GPU 只说明驱动透传正常，不等于 CUDA 开发环境完整。建议先在 WSL 建立项目独立 venv 做 Linux 单机版，集群环境另行部署；若以后需要利用宿主全部内存，再显式调整 `.wslconfig`，本次没有改动。

### 1.3 当前能力边界

**现在可靠可跑：**

- bounded-degree 图上 CSR/CSC 合成、稀疏乘法、SuperLU 中等规模直接解；
- CG/MINRES/GMRES 等矩阵-free 迭代；
- 抽样根节点的转置稀疏求解与 WLS/Green 块行尾部；
- Chebyshev 多项式、局部 Dirichlet 截断、小规模 Schur oracle；
- 多进程并行的图族 × 参数 × seed 扫描；
- JSON/CSV 级结果记录。

**当前不能诚实声称可跑：**

- Python/CUDA 版稀疏算法；
- MPI 强/弱扩展；
- PETSc/Hypre/Trilinos/MUMPS/CHOLMOD 基线；
- MATLAB/Julia 对照；
- pandapower 标准电网实例；
- 大规模 Parquet/HDF5/Zarr 数据流水线。

## 2. 实验要验证的理论命题

实验按 falsifiable hypotheses 组织，而不是按算法名单组织。

### H1：算子尾部确实刻画固定实例的最佳局部 worst-case 误差

对抽样根节点 `i` 和半径 `r`，计算

`e_i^b(r)=||E_iJ^{-1}Q_out||_{2→2}`

以及原始 WLS 数据版本

`e_i^z(r)=||E_iJ^{-1}H^TR^{-1}Q_out^z||_{2→2}`。

用局部截断算法的实测误差和针对尾部奇异向量构造的 adversarial input 同时验证上、下界。两类尾部不能混用：前者把本地初始量视为信息向量 `b`，后者才对应原始测量 `z`。

### H2：局部低谱质量比全局条件数更能预测多数节点局部性

构造少量病态区域，使 `κ(J_n)` 随 `n` 发散，但远离病态区的大多数节点维持小的

`∫ λ^{-2} dμ_{n,i}(λ)` 或低谱窗口质量。

比较：全局 `κ`、节点到病态区距离、局部 Lanczos 谱测度、Schur 边界影响与 `e_i(r)` 的解释力。核心检验是：固定 `(r,ε)` 时，`F_n(r,ε)` 是否只随坏区测度而不随全局最坏条件数恶化。

### H3：存在 `δ` 相变而不是 uniform locality

估计

`F_n(r,ε)=n^{-1}|{i:e_i(r)>ε}|`

和经验 `q_{r,δ}`。调节病态区比例 `θ`，测试当 `δ>θ` 时所需半径是否保持有界，而 uniform (`δ=0`) 半径随图规模/病态强度增长。该实验直接区分“多数节点可局部”与“全局 operator norm 可局部”。

### H4：cover/local-view 不可区分性给出匹配轮数下界

生成具有相同 rooted `r`-ball、但远端闭合方式/锚点/病态区不同的 graph covers。对成对实例保持根节点局部系数和数据完全一致，计算两个中央目标之差。任何仅知局部视图的算法，其一侧误差至少为目标差的一半。将此 converse 与截断/多项式 achievability 曲线比较，而不是只比较算法之间谁快。

### H5：结构性算法界能预测误差，而非只拟合 wall time

测试 Dirichlet 截断、局部 Schur 近似、固定 Chebyshev、阻尼 Jacobi/GaBP 与可选局部 Krylov。每个算法同时报告：

- 实例误差；
- 它对应的 worst-case/operator certificate；
- 与最优 sampled tail 的 gap；
- 是否用了全局谱界、全局预处理或离线系数。

### H6：通信资源前沿与局部性相变一致

固定浮点消息或明确量化器后，画出

`(rounds, message dimension, bit-hop, error, failure fraction)`

的 Pareto 集。若加入无线几何，再单独定义 weighted bit-hop / radio model energy；不能把抽象 bit-hop 直接标成 joules。

## 3. 合成图族与块 WLS 生成器

### 3.1 图族矩阵

| 族 | 作用 |
|---|---|
| path / balanced tree | 检查直径下界、长程锚点传播、树不等于常数半径 |
| cycle / torus | 隔离“有环”效应，避免把环与高 treewidth 混为一谈 |
| 2-D / 3-D grid | 测试几何 Green 衰减和 separator/Schur 代价 |
| bounded-degree random regular | expander、高混合、局部树状但全局闭合快 |
| barbell / lollipop / weak cut | 制造小谱隙和跨瓶颈远端影响 |
| stochastic block / geometric | 控制 community 与空间通信代价 |
| random lifts / paired covers | 构造相同局部视图、不同全局目标的不可区分下界 |
| chordal / bounded-treewidth | 对照精确消元与 fill-in，验证环本身不是障碍 |
| power-grid cases（以后） | 真实结构外部验证；当前缺 pandapower，不列入 smoke |

规模取 `n=2^k` 或近似等比序列，degree、block dimension 和坏区比例分别控制。不要让 graph family 随 `n` 无意改变平均度。

### 3.2 块 WLS 因子

基础模型用有限作用距离因子：

- 节点锚点：`z_i=A_ix_i+v_i`；
- 边/局部超边：`z_e=B_{eu}x_u+B_{ev}x_v+v_e`；
- `R` 先取块对角，之后再加有限范围相关噪声；
- `J=H^TR^{-1}H`，目标 `C_ix` 可取本地块或低维功能。

必须同时保存 communication graph、factor/support graph 和 measurement holder map。若边测量被两端共同持有，`r`-球外 mask 是“其 support 与球不相交”；若只有一个 holder，mask 按 holder 定义。两种语义会改变 `T_i` 尾部，不能混为一个数据集。

### 3.3 病态区域注入

至少包含四种机制：

1. 连通子区的 anchor 权重乘 `η_a≪1`；
2. 穿过坏区/边界的 coupling 权重乘 `η_e≪1`；
3. 块测量方向近线性相关，控制最小奇异值而不改变拓扑；
4. 少量高 leverage 或强相关噪声块。

控制量为坏区测度 `θ`、直径、到根距离、`η_a,η_e` 和块维 `d`。应设计两组 `J` 拥有相似全局条件数但不同坏区位置，以及相同局部谱质量但不同远端 cover，以区分相关性和因果机制。

### 3.4 当前 smoke 生成器的限制

现有 skeleton 为代码清晰使用 dense block rows 后转 CSR，只适合 smoke/小规模；正式大规模生成器必须直接写 COO/CSR triplets，避免 `O(n^2d^2)` 临时存储。现有 `cover` 是四节点 base graph 的随机 lift，足以测试接口，但 paired-cover converse 需要显式保证两实例 rooted `r`-ball 同构并同步局部系数。

## 4. 中心解、局部算法与证书

### 4.1 中心稀疏解

- 小中规模：`splu/spsolve`，分解一次、多 RHS 复用；记录 symbolic/numeric factor time 和 fill ratio。
- 大规模 SPD：PCG + 可说明来源的预条件器；比较 residual error 与 estimation error，不能只用迭代残差替代统计误差。
- 块不定 KKT/约束模型：MINRES/GMRES。
- 集群：PETSc KSP + Hypre/GAMG，必要时 MUMPS/SuperLU_DIST 作较小基线。

**禁止显式形成 `J^{-1}`。** 中央 `x̂` 只解 `Jx=b`。

### 4.2 抽样根节点的精确块行尾部

对目标块 `E_i` 解

`J^TY_i=E_i^T`。

于是 `E_iJ^{-1}=Y_i^T`，而

`T_i=E_iJ^{-1}H^TR^{-1}=(R^{-1}HY_i)^T`。

对球外坐标/测量列做小型 SVD 得到块行 `2→2` tail。一个根节点只需 `d` 个转置 RHS；相同 factorization 可供全部抽样根复用。大规模迭代解应把线性求解误差传播成 tail 的上下界，不能把近似块行当精确 ground truth。

当根数过多时使用两条路线：

- 随机均匀抽根，逐根转置求解，估计 `F_n`/quantile；
- Hutchinson/Hutch++ 或 block probes 估计平均 Frobenius tail/trace。

第二条不能冒充逐节点 `2→2` tail；它只回答平均平方风险。

### 4.3 `r`-局部算法

1. **Dirichlet truncation**：只解 `J_BB x_B=b_B`，便宜但边界偏差可非单调。
2. **Exact Schur oracle**：用全局外部消元形成 `J_BB-J_BOJ_OO^{-1}J_OB`；它只用于小实例正确性，不是局部算法。
3. **Localized Schur**：只消元厚度 `s` 的 shell，或压缩边界 Dirichlet-to-Neumann map；报告 `(r,s)`。
4. **Fixed polynomial**：Chebyshev/Jackson/Remez 逼近 `1/x`；degree `k` 对传播半径为 `kρ(J)`。全局谱界若由离线计算提供，必须计为 global knowledge。
5. **Stationary message passing**：Jacobi/Richardson/GaBP；区分算法收敛条件与问题可行性。
6. **Krylov**：CG 的全局 inner products 不属于纯邻居 LOCAL；只有改成局部/去同步版本后才能与相同通信模型比较。

### 4.4 局部谱量

不能全特特征分解。对抽样根用 block Lanczos / stochastic Lanczos quadrature 估计 `μ_{n,i}`，记录 quadrature steps 与误差；对全局谱端点用 `eigsh`、Gershgorin 或构造时已知界。检验窗口包括

- `μ_i([0,λ_0])`；
- `∫_{λ>0}λ^{-2}dμ_i(λ)` 的截断/正则版本；
- 最优次数 `k` 的局部 `L²(μ_i)` 多项式误差。

这些首先是 polynomial achievability 证书；是否对任意协议必要必须由 paired-cover 或其他不可区分下界验证。

## 5. 指标与统计协议

### 5.1 误差指标必须分层

- `operator_wc_b`：从本地 `b` 到目标的 worst-case tail；
- `operator_wc_z`：从原始 `z` 到目标的 WLS worst-case tail；
- `bayes_mse_gap`：局部与全局条件均值的风险差；
- `instance_rel_l2` / nodewise block error：具体抽样数据误差；
- adversarial paired-view gap：局部不可区分 converse。

任何图都不能用一次随机 `x,z` 的 RMSE 代替 operator/Bayes 结论。

### 5.2 `δ` 曲线的抽样误差

令根节点从 `V_n` 均匀、有放回抽样 `m` 次，对每个根精确/高精度估计 `e_i(r)`。经验 CDF 的 DKW 带满足

`P(sup_t |F̂_m(t)-F_n(t)|>η)≤2exp(-2mη²)`，

所以 `η=sqrt(log(2/α)/(2m))`。它同时覆盖全部 `ε`/quantile；若只关心预先固定的一个 `(r,ε)`，可用 Wilson 或 Clopper–Pearson 二项区间得到更窄的 pointwise 报告。例：`α=.05` 时，DKW 半宽约 `.10` 需 `m≥185`，半宽 `.05` 需 `m≥738`。

抽样必须独立于图值和坏区；如采用分层抽样，应报告权重并使用相应有限总体/分层置信区间。多个半径复用同一根样本有利于 paired comparison，但跨半径 simultaneous inference 要使用统一 DKW 带或做多重性校正。

### 5.3 通信与资源

每个结果至少记录：

- sequential rounds；
- 每条边每轮 message scalar count 和 quantization bits；
- 总 bit-hop、最大节点 bit、峰值边负载；
- 若为无线几何模型：距离加权 transmit/receive energy，明确路径损耗和 MAC 假设；
- preprocessing、factorization、solve/iteration、tail certification 分项时间；
- peak RSS、GPU peak allocated/reserved、matrix/vector bytes；
- 本地 FLOP/SpMV 次数和全局 reduction 次数。

“32-bit 浮点数 × directed edge-round”只叫 bit-hop proxy，不能叫实际能量。

## 6. 分层实现和扩展实验

### 6.1 Tier 0：当前 Python smoke

文件：

- `experiments/locality_benchmark.py`
- `experiments/configs/smoke.json`
- `experiments/README.md`

它不形成逆矩阵；抽样根通过 `J^TY=E_i^T` 得到精确块行；包含 central、Dirichlet、Schur oracle、Chebyshev、raw-WLS tail、DKW failure curve、bit-hop proxy 和正确性断言。

本次实跑：72 节点 cycle，块维 2，状态维 144，`κ≈2068.9`，总 wall time 约 0.08 s。三项检查全部通过：中央相对残差 `1.6e-16`，Schur oracle 最大根误差约 `3.2e-15`，抽样 raw-WLS tail 随半径非增。抽样中位 raw-WLS tail 从 `r=0` 的约 `.740` 降到 `r=8` 的约 `.0213`；`ε=.1` 的抽样失败比例从 `r≤4` 的 100% 降到 `r=8` 的 6.25%。由于仅 `m=16`，95% DKW 半宽约 `.340`，这正说明正式曲线必须增加根样本而不能过度解读 smoke。

Dirichlet 样本误差前几步不必单调；Chebyshev 全局误差在该病态例上仍很大且不严格单调。这不是 smoke 失败，而是提示正式实验要将“最优算子 tail”“具体局部边界近似”和“固定多项式”分开。

### 6.2 Tier 1：工作站 CPU

外层用 process pool/joblib 并行 `(graph family, n, θ, η, d, seed)`；每个进程内部把 BLAS/OpenMP 线程固定为 1–2。单个大实例才反过来使用较多 BLAS/solver 线程。必须在结果中保存 `OMP_NUM_THREADS/OPENBLAS_NUM_THREADS`，否则时间不可复现。

建议阶段：

1. `N=nd≤10^4`：精确 sampled tails、Schur oracle、paired-cover converse；
2. `10^4<N≤2×10^5`：factor reuse + 64–256 sampled roots，视 fill-in 调整；
3. `N≈10^6` 或更大：bounded-degree CSR + iterative/matrix-free，主要做 polynomial/probes，中央真值用严格 tolerance 迭代解和残差传播。

这些只是保守起点，不是硬上限。路径/树和 2-D grid 的 fill-in 完全不同；每次扩容依据实际 `nnz(L+U)` 和 RSS，而不是只看 `N`。当前 smoke 生成器在升级为直接 COO 前不能用于后两档。

### 6.3 Tier 2：本机 GPU

GPU 最适合：

- 同一 `J` 上大量向量的 batched SpMV；
- Chebyshev/Lanczos/SLQ；
- block-CG/多 RHS iterative tails；
- 参数相同的小图 batch。

不建议先做 GPU 稀疏 LU：16 GiB 和 fill-in 很快成为瓶颈，且与“局部算法”的主要计算模式不一致。WDDM 下显示/桌面当前约占数 GiB，单作业按可用显存而非标称 16 GiB 配 batch。

当前只能把 GPU 实现列为下一环境阶段：硬件、driver、nvcc、cuSPARSE 已在，但 Python 是 CPU-only。若以后采用 C++/CUDA，可从 VS Developer Prompt 调用 MSVC + nvcc；若采用 Python，应固定与 RTX 5080/compute 12.0 兼容的 CUDA runtime，并先做 SpMV 数值一致性测试。这里不安装任何组件。

### 6.4 Tier 3：MPI 集群

推荐 Linux + Slurm，软件栈为编译器/MPI、PETSc（KSP + GAMG/Hypre）、ParMETIS/PT-Scotch，可选 CUDA-aware MPI。分层并行为：

1. **作业数组**：图族、规模、病态参数、seed，几乎无通信；
2. **节点内线程/GPU**：本地 CSR SpMV、多 RHS；
3. **节点间 MPI**：按图划分矩阵，每轮 polynomial/message passing 只交换 halo；
4. **少量全局 reduction**：谱端点、CG/Lanczos；单独计数，不混入邻居通信。

示意作业数组：

```bash
#SBATCH --array=0-255
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
export OMP_NUM_THREADS=1
srun python run_manifest_item.py --manifest configs.jsonl --index "$SLURM_ARRAY_TASK_ID"
```

分布式单实例：

```bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
srun --cpu-bind=cores python mpi_locality.py --config case.json
```

GPU 版本通常设一 rank/GPU，并记录 GPU topology 与 MPI GPU-direct 是否实际启用。脚本只是部署模板；当前仓库没有 MPI driver，不能把它作为已运行结果。

### 6.5 强/弱扩展

**强扩展：** 固定同一个大图、相同 partition seed、相同 tolerance 和 sampled roots，增加 ranks/GPUs；报告 total、SpMV、halo、reduction、imbalance，直到通信主导。不能改变迭代数而只画总时间。

**弱扩展：** 固定每 rank 节点数、边数、坏区比例和根样本数；扩大 cover/grid。报告每轮 halo bytes 与跨分区 edge cut，避免图划分质量变化伪装成算法扩展性。

**算法扩展与机器扩展分离：** `r` 增大导致的 rounds/bit-hop 是理论资源；rank 增大导致的 MPI 时间是实现资源。二者分别画图。

## 7. 正确性与复现规范

### 7.1 必过的数值检查

1. `J` 对称误差和 SPD/目标可估计性检查；
2. 中央解 `||Jx-b||/||b||` 与 forward-error proxy；
3. 小实例用 dense/高精度参考交叉验证，但生产路径不形成逆；
4. Schur oracle 的根块与中央解一致；
5. 精确 sampled tail 对嵌套球非增，直径半径后为零；
6. 构造尾部顶奇异向量，验证 worst-case 下界达到；
7. paired covers 的 rooted ball、局部系数和本地数据逐项一致；
8. `R≠I` 时显式核对 `T_i=(R^{-1}HY_i)^T`，不能漏掉 `R^{-1}`；
9. float32/float64、CPU/GPU 的误差差异按 condition number 缩放解释；
10. 迭代中央解的 residual tolerance 必须明显小于被比较的 locality error。

### 7.2 数据与 provenance

每个 case 保存规范化 config、SHA-256、seed、代码 commit、OS/CPU/GPU、包版本、线程环境、solver/tolerance、图摘要与 measurement-holder 语义。当前目录不是 Git repository，所以 smoke 不能记录 commit；正式研究应在版本库中运行或保存源码 hash。

随机数用主 seed 派生 graph / coefficient / state / noise / root-sampling 子 seed，避免改变一个模块后其余随机流漂移。结果建议一行一 case 的 JSONL/Parquet，矩阵只对失败案例保存 `.npz`；不要默认保存每个大图的所有根块行。

### 7.3 Benchmark 卫生

- 生成、预处理、factorization、solve 和 certificate 分项计时；
- GPU 计时前 warm-up 和 synchronize；
- 至少 3 次独立进程重复，报告 median 和离散度；
- 固定 CPU affinity/threads，记录后台负载；
- 先检查 correctness，再纳入 timing 汇总；
- timeout、OOM 和不收敛是数据点，不能静默删除。

## 8. 推荐执行顺序

1. 把 smoke 生成器改成直接 COO，并新增 unit tests；
2. 做 `path/cycle/grid/random_regular/barbell/cover × d∈{1,2,4} × n≤10^4` 的 CPU correctness sweep；
3. 实现抽样根 `e_i^z(r)` + DKW/Wilson，先用 `m=185` 得到 ±0.10 的 simultaneous CDF 精度；
4. 实现局部 Lanczos 谱量和 paired-cover converse，验证 H2–H4；
5. 冻结理论图表后再优化 CPU 并行；
6. 单独建立 CUDA 环境并只迁移 SpMV/polynomial/probes；
7. 最后移植 PETSc/MPI，做强弱扩展和通信字节审计。

先完成第 2–4 步，才有足够证据判断“measured Green locality”是否值得写成定理。若局部谱量不能预测 tail，或 paired-cover 下界与任何谱证书不匹配，应及时修改理论，而不是靠更大机器掩盖。
