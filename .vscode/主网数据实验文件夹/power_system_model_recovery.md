# 环网下数据中心承载力评估模型

### 引言

### 研究背景与目标

随着云计算、大数据和人工智能技术的快速发展，数据中心作为数字经济的核心基础设施，其能耗和电力需求呈现快速增长趋势。数据中心的高功率密度和动态负载特性对电网的稳定运行和承载能力提出了新的挑战。特别是在环网结构下，数据中心的接入位置、容量配置和运行策略会显著影响电网的潮流分布、电压稳定性和供电可靠性。

本研究旨在建立环网结构下数据中心承载力评估的数学模型，通过优化数据中心的功率分配和电网的运行状态，评估电网在满足安全约束条件下所能承载的最大数据中心负荷。该模型综合考虑了交流潮流的非线性特性、能量存储系统的动态特性以及数据中心与电网的耦合关系，为电网规划、运行和扩容提供理论依据和决策支持。

### 数学模型

考虑一个包含$n$个节点的环网系统，时间周期为$t = 1, 2, \dots, T$。模型的目标是最大化数据中心的总负荷承载能力，同时满足电网的安全运行约束。数学模型表述如下：

#### 模型概述与SOCP松弛原理

交流最优潮流（AC-OPF）问题本质上是非凸非线性规划问题，直接求解计算代价高且难以保证全局最优。本模型采用**二阶锥规划（SOCP）松弛**方法，将原始非凸问题转化为可高效求解的凸优化问题。

**松弛的核心思路**：

交流潮流方程要求辅助变量满足等式约束（锥面约束）：

$$
T_{ij,t}^2 + Z_{ij,t}^2 = U_{i,t} \cdot U_{j,t}
$$

将上述等式约束松弛为二阶锥不等式约束（式(5)），允许解落在锥的内部或表面，从而将非凸问题转化为 SOCP 凸问题。

**松弛的几何意义**：
- 原始约束要求点 $(T_{ij,t}, Z_{ij,t}, U_{i,t}, U_{j,t})$ 位于**二阶锥的表面**（锥面约束）
- 松弛后的约束允许点位于**二阶锥的内部或表面**（锥约束）

**松弛有效性条件**：在辐射状（树形）配电网结构中，该松弛通常是**紧的（tight）**，即最优解自然满足原始等式约束，SOCP 解可直接作为原问题的最优解使用。对于环网结构，需要通过附加的**相角参考约束**（式(7)）消除潮流方程的旋转不变性，并在求解后验证松弛间隙，必要时进行**锥还原**以恢复物理可行解。

**模型变量体系**：本模型引入三类辅助变量代替原始电压相角变量：

| 辅助变量 | 物理含义 | 关系式 |
|--------|--------|------|
| $T_{ij,t}$ | 与余弦相角相关的支路变量 | $T_{ij,t} = V_{i,t}V_{j,t}\cos\theta_{ij,t}$ |
| $Z_{ij,t}$ | 与正弦相角相关的支路变量 | $Z_{ij,t} = V_{i,t}V_{j,t}\sin\theta_{ij,t}$ |
| $U_{i,t}$ | 节点电压幅值平方 | $U_{i,t} = V_{i,t}^2$ |

通过引入上述辅助变量，功率平衡方程（式(3)、式(4)）和支路潮流方程（式(15)、式(16)）均可表示为关于这些变量的线性方程，二阶锥松弛约束（式(5)）将非凸的电压相角关系转化为标准二阶锥约束，从而构成完整的 SOCP 模型。

#### 目标函数

$$
\min \quad -hc \tag{1}
$$

#### 约束条件

$T_{\text{set}} = \{0,1,2,\dots,T-1\}$，$D_{\text{set}} = \{0,1,2,\dots,D\}$。

$$
hc \leq \sum_{i \in IDC} P_{IDC, i, t}, \quad \forall t \in T_{\text{set}} \tag{2}
$$

$$
P_{G, i, t} - P_{L, i, t} - P_{IDC, i, t} - U_{i, t} G_{ii} - \sum_{\substack{j=1 \\ j \neq i}}^{n} \left( T_{ij, t} G_{ij} + Z_{ij, t} B_{ij} \right) = 0, \quad \forall i \in N, \forall t \in T_{\text{set}} \tag{3}
$$

$$
Q_{G, i, t} - Q_{L, i, t} - 0 + U_{i, t} B_{ii} - \sum_{\substack{j=1 \\ j \neq i}}^{n} \left( Z_{ij, t} G_{ij} - T_{ij, t} B_{ij} \right) = 0, \quad \forall i \in N, \forall t \in T_{\text{set}} \tag{4}
$$

$$
\left\| 2T_{ij, t}; 2Z_{ij, t}; U_{i, t} - U_{j, t} \right\|_2 - \left( U_{i, t} + U_{j, t} \right) \leq 0, \quad \forall (i,j) \in B, \forall t \in T_{\text{set}} \tag{5}
$$

$$
T_{ij, t} - T_{ji, t} = 0, \quad Z_{ij, t} + Z_{ji, t} = 0, \quad \forall (i,j) \in B, \forall t \in T_{\text{set}} \tag{6}
$$

$$
\sum_{(i,j) \in \Omega_{B_k}} \left[ \tan^{-1} \left( Z_{ij, t}^{0} / T_{ij, t}^{0} \right) + \frac{ \left( T_{ij, t}^{0} Z_{ij, t} - Z_{ij, t}^{0} T_{ij, t} \right) }{ \left( T_{ij, t}^{0} \right)^2 + \left( Z_{ij, t}^{0} \right)^2 } \right] = 0, \quad \forall k, \forall t \in T_{\text{set}} \tag{7}
$$

$$
T_{ii, t} = U_{i, t}, \quad \forall i \in N, \forall t \in T_{\text{set}} \tag{8}
$$

$$
Z_{ii, t} = 0, \quad \forall i \in N, \forall t \in T_{\text{set}} \tag{9}
$$

$$
U_{i, \min} \leq U_{i, t} \leq U_{i, \max}, \quad \forall i \in N, \forall t \in T_{\text{set}} \tag{10}
$$

$$
P_{Gi, \min} \leq P_{G, i, t} \leq P_{Gi, \max}, \quad \forall i \in N, \forall t \in T_{\text{set}} \tag{11}
$$

$$
Q_{Gi, \min} \leq Q_{G, i, t} \leq Q_{Gi, \max}, \quad \forall i \in N, \forall t \in T_{\text{set}} \tag{12}
$$

$$
U_{0,t} = 1 , \quad \forall t \in T_{\text{set}} \tag{13}
$$

$$
T_{ij, t} = V_{i, t} V_{j, t} \cos \theta_{ij, t}, \quad Z_{ij, t} = V_{i, t} V_{j, t} \sin \theta_{ij, t}, \quad U_{i, t} = (V_{i, t})^2, \quad \forall (i,j) \in B, \forall t \in T_{\text{set}} \tag{14}
$$

$$
P_{ij, t} = -U_{i, t} G_{ij, t} + \left( T_{ij, t} G_{ij, t} + Z_{ij, t} B_{ij, t} \right), \quad \forall (i,j) \in B, \forall t \in T_{\text{set}} \tag{15}
$$

$$
Q_{ij, t} = U_{i, t} B_{ij, t} + \left( -T_{ij, t} B_{ij, t} + Z_{ij, t} G_{ij, t} \right), \quad \forall (i,j) \in B, \forall t \in T_{\text{set}} \tag{16}
$$

$$
\lambda_{t,u} = \lambda_{t,u}^{\text{shape}} \cdot scale, \quad \forall t \in T_{\text{set}}, \forall u \in U \tag{17}
$$

$$
q^{(0)}_{t,u} = \lambda_{t,u} + \sum_{v:(v,u)\in\mathcal{E}} y^{(0)}_{t, v\rightarrow u}, \quad \forall u\in U,\; t \in T_{\text{set}} \tag{18}
$$

$$
q^{(d)}_{t,u} = r^{(d-1)}_{t-1,u} + \sum_{v:(v,u)\in\mathcal{E}} y^{(d)}_{t, v\rightarrow u}, \quad \forall u\in U,\; d=1,\dots,D,\; t=1,\dots,T-1 \tag{19}
$$

$$
q^{(d)}_{0,u}=0, \quad \forall u\in U,\; d=1,\dots,D \tag{20}
$$

$$
r^{(d)}_{t,u} = q^{(d)}_{t,u} - p^{(d)}_{t,u} - \sum_{v:(u,v)\in\mathcal{E}} y^{(d)}_{t, u\rightarrow v}, \quad \forall u,\; d=0,1,\dots,D,\; t \in T_{\text{set}} \tag{21}
$$

$$
r^{(d)}_{t,u}\ge 0, \quad \forall u,\; d,\; t \tag{22}
$$

$$
r^{(D)}_{t,u}=0, \quad \forall u\in U,\; t \in T_{\text{set}} \tag{23}
$$

$$
s\sum_{d=0}^{D} y^{(d)}_{t, u\rightarrow v} \le \bar B_{t,uv}, \quad \forall (u,v)\in\mathcal{E},\; t \in T_{\text{set}} \tag{24}
$$

$$
p^{(d)}_{t,u}\ge 0,\; y^{(d)}_{t, u\rightarrow v}\ge 0, \quad \forall u,(u,v),t,d \tag{25}
$$

$$
q^{(d)}_{T-1,u}=q^{(d)}_{0,u}, \quad \forall u\in U,\; d=0,\dots,D \tag{26}
$$

$$
\lambda_{t,u} \ge 0, \quad \forall u\in U,\; t \in T_{\text{set}} \tag{27}
$$

$$
\sum_{d=0}^{D} p^{(d)}_{t,u} \le \mu_{\max,u}, \quad \forall u \in U, \forall t \in T_{\text{set}} \tag{28}
$$

$$
P_{IDC, i, t} = \alpha_{u, t} \sum_{u \in U_i}\sum_{d=0}^{D} p^{(d)}_{t,u} + \beta_{u, t}, \quad \forall i \in IDC, \forall t \in T_{\text{set}} \tag{29}
$$

#### 变量及约束解释

##### 1.2.3.1 变量解释

**决策变量：** $T_{ij,t}$, $Z_{ij,t}$, $U_{i,t}$, $P_{G,i,t}$, $Q_{G,i,t}$, $P_{IDC,i,t}$, $P_{ij,t}$, $Q_{ij,t}$, $\lambda_{t,u}$, $q^{(d)}_{t,u}$, $r^{(d)}_{t,u}$, $p^{(d)}_{t,u}$, $y^{(d)}_{t, u\rightarrow v}$, $scale$, $hc$

**已知量（参数）：** $\lambda_{t,u}^{\text{shape}}$, $\bar{B}_{t,uv}$, $s$, $\mu_{\max,u}$, $\alpha_{u,t}$, $\beta_{u,t}$, $G_{ij}$, $B_{ij}$, $G_{ii}$, $B_{ii}$, $U_{i,\min}$, $U_{i,\max}$, $P_{Gi,\min}$, $P_{Gi,\max}$, $Q_{Gi,\min}$, $Q_{Gi,\max}$, $P_{L,i,t}$, $Q_{L,i,t}$, $T_{ij,t}^{0}$, $Z_{ij,t}^{0}$

**集合定义：**

- $T_{\text{set}} = \{0,1,2,\dots,T-1\}$：连续时间潮流优化的时间集合，其中$T$为时间周期总数
- $D_{\text{set}} = \{0,1,2,\dots,D\}$：可延迟计算任务的延迟时间集合，其中$D$为最大延迟级别
- $N$：电网节点集合，包含$n$个节点
- $B$：电网线路集合
- $U$：数据中心编号集合（算力节点集合）
- $U_i$：节点$i$（电力节点）下数据中心编号集合
- $\mathcal{E}$：数据中心之间的通信连接集合（算力链路集合）
- $k$：环网索引，用于标识电网中的不同环路
- $\Omega_{B_k}$：第$k$个环的支路集合，用于相角参考约束，其中支路$(i,j) \in \Omega_{B_k}$按正方向（不包括反向）构成一个环
- $IDC$：数据中心的电力节点集合

**环网潮流二阶锥松弛变量：**
式(2)-(16)为环网潮流二阶锥松弛模型。交流潮流方程是非凸的，难以直接求解。通过引入辅助变量$T_{ij,t}$、$Z_{ij,t}$和$U_{i,t}$，将原非凸问题松弛为二阶锥规划问题：

- $V_{i,t}$：节点$i$在时刻$t$的电压幅值 (p.u.)
- $\theta_{ij,t}$：支路$(i,j)$在时刻$t$的相角差 (rad)
- $T_{ij,t} = V_{i,t} V_{j,t}\cos\theta_{ij,t}$：与节点电压幅值和相角余弦相关的辅助变量
- $Z_{ij,t} = V_{i,t} V_{j,t}\sin\theta_{ij,t}$：与节点电压幅值和相角正弦相关的辅助变量
- $U_{i,t} = (V_{i,t})^2$：节点$i$电压幅值的平方
- $P_{G,i,t}$、$Q_{G,i,t}$：节点$i$在时刻$t$的有功和无功发电功率
- $P_{L,i,t}$、$Q_{L,i,t}$：节点$i$在时刻$t$的有功和无功负荷
- $P_{IDC,i,t}$：节点$i$在时刻$t$的数据中心功耗
- $P_{ij,t}$、$Q_{ij,t}$：支路$(i,j)$在时刻$t$的有功和无功功率流
- $T_{ij,t}^{0}$、$Z_{ij,t}^{0}$：辅助变量的初始值或参考值，用于相角参考约束的线性化

**数据中心任务调度变量：**
式(17)-(29)为数据中心任务调度模型，描述数据中心间的任务分配和处理过程：

- $\lambda_{t,u}$：时段$t$在算力节点$u$的任务到达量（决策变量，job/period）
- $q^{(d)}_{t,u}$：时段$t$内算力节点$u$的第$d$级缓存（队列）任务量 (job)
- $r^{(d)}_{t,u}$：时段$t$结束后算力节点$u$的第$d$级缓存剩余任务量 (job)
- $p^{(d)}_{t,u}$：时段$t$内算力节点$u$对第$d$级缓存任务的本地处理量 (job)
- $y^{(d)}_{t, u\rightarrow v}$：时段$t$开始时从算力节点$u$向$v$迁移的第$d$级缓存任务量 (job)
- $scale$：任务分布的全局缩放因子（决策变量），用于等比例缩放所有算力节点的任务到达量
- $hc$：全局承载力下界（优化目标变量），代表配网最小可承载总功率 (kW)

**数据中心任务调度参数：**

- $\lambda_{t,u}^{\text{shape}}$：算力节点$u$在时段$t$的基础任务分布参数，用于限制真实到达量$\lambda_{t,u}$的时空比例
- $\bar{B}_{t,uv}$：时段$t$算力链路$(u,v)$的带宽上限 (job/period)
- $s$：单位任务的数据量系数 (Band/job)，用于将任务迁移量换算为带宽占用
- $\mu_{\max,u}$：算力节点$u$的单时段最大处理能力 (job/period)
- $\alpha_{u,t}$、$\beta_{u,t}$：算力节点$u$的任务处理量到IDC用电功率的线性映射系数（用于构造$P^{IDC}$），包含制冷负荷映射

**电网相关参数：**

- $G_{ij}$、$B_{ij}$：支路$(i,j)$的电导和电纳（导纳参数）
- $G_{ii}$、$B_{ii}$：节点$i$的自导纳（对地导纳），包括并联电容、电抗等
- $U_{i,\min}$、$U_{i,\max}$：节点$i$电压平方的上下限，对应电压幅值范围$(V_{i,\min})^2 \le U_{i,t} \le (V_{i,\max})^2$
- $P_{Gi,\min}$、$P_{Gi,\max}$：节点$i$发电机有功出力的上下限
- $Q_{Gi,\min}$、$Q_{Gi,\max}$：节点$i$发电机无功出力的上下限


##### 1.2.3.2 约束解释

**目标函数相关约束：**

- 式(1)：定义优化目标为最大化全局承载力下界$hc$，等价于最小化$-hc$
- 式(2)：$hc$不超过所有节点数据中心在任意时刻的总功耗

**环网潮流约束：**

- 式(3)-(4)：节点有功和无功功率平衡方程，确保发电、负荷和网损之和为零
- 式(5)：二阶锥松弛约束，将非凸的电压相角关系转化为凸约束
- 式(6)：辅助变量对称性约束，$T_{ij,t}$对称，$Z_{ij,t}$反对称
- 式(7)：相角参考约束，对第$k$个环中所有支路的相角求和为0，消除潮流方程的旋转不变性。每个环的支路按正方向求和（不包括反向支路）
- 式(8)-(9)：对角线元素定义，$T_{ii,t}=U_{i,t}$，$Z_{ii,t}=0$
- 式(10)-(12)：运行安全约束，包括电压幅值、有功发电和无功发电的上下限
- 式(13)：参考节点电压设定为标幺值1
- 式(14)：辅助变量与物理量关系定义
- 式(15)-(16)：支路有功和无功功率流方程

**数据中心任务调度约束：**

- 式(17)：任务到达量定义，通过缩放因子scale对基础任务分布进行缩放
- 式(18)：第0级缓存任务量，等于新到达任务量加上从其他节点迁移来的任务量
- 式(19)：第$d$级缓存任务量状态更新方程（$d \geq 1$），等于上一时段剩余任务加上迁移来的任务
- 式(20)：初始条件，0时刻除第0级外的其他级缓存为空
- 式(21)：缓存剩余任务量更新方程，等于缓存任务量减去处理量和迁移出量
- 式(22)：缓存剩余任务量非负约束
- 式(23)：最高级缓存必须清空约束，确保所有任务最终都被处理
- 式(24)：带宽约束，所有级别任务的总迁移量不超过带宽上限
- 式(25)：任务处理量和迁移量非负约束
- 式(26)：周期性闭环约束，适用于周期性调度场景
- 式(27)：任务到达量非负约束
- 式(28)：单时段处理能力约束，所有级别任务的总处理量不超过最大处理能力
- 式(29)：数据中心功耗模型，功耗与处理的任务量成线性关系

该模型通过二阶锥松弛将非凸的交流潮流问题转化为可高效求解的凸优化问题，同时考虑了数据中心任务调度的时空特性，实现了电网运行与数据中心负载的协同优化。

### 算法设计

#### 锥还原算法

##### 问题背景

SOCP 求解得到的最优解可能落在二阶锥的**内部**而非表面，此时松弛不紧，解不对应物理可行的交流潮流解。需要通过**锥还原（Cone Recovery）**将解投影或驱动到锥面上。

**锥松弛间隙**定义为：

$$
\Delta_{t,i,j} = l^*_{t,i,j} \cdot v^*_{t,i} - \left[(P^*_{t,i,j})^2 + (Q^*_{t,i,j})^2\right]
$$

其中采用 DistFlow 支路潮流模型的记号（$l_{t,i,j}$ 为支路电流幅值平方，$v_{t,i}$ 为节点电压幅值平方）。间隙 $\Delta_{t,i,j} \approx 0$ 时松弛是紧的，无需还原；$\Delta_{t,i,j} > \epsilon$ 时需进行锥还原。

在本模型（导纳矩阵形式）中，对应的松弛间隙为：

$$
\delta_{ij,t} = \frac{(U_{i,t}+U_{j,t})^2 - [4T_{ij,t}^2 + 4Z_{ij,t}^2 + (U_{i,t}-U_{j,t})^2]}{(U_{i,t}+U_{j,t})^2}
$$

$\delta_{ij,t} \in [0,1]$，越接近 0 表示松弛越紧。判断标准：$\max \delta_{ij,t} < 10^{-4}$ 视为严格紧，$10^{-4}$ 到 $10^{-2}$ 为近似紧，超过 $10^{-2}$ 需重点关注。

##### 基于方向割的锥还原方法

锥还原方法基于 **Cauchy-Schwarz 不等式**，导出线性的**方向割约束**来收紧可行域。

对于支路 $(i,j)$，记锥表示向量 $\mathbf{x}_{ij} = (2P_{t,i,j}, 2Q_{t,i,j}, l_{t,i,j}-v_{t,i})^T$，二阶锥约束为 $\|\mathbf{x}_{ij}\|_2 \leq l_{t,i,j} + v_{t,i}$。

对任意方向向量 $\mathbf{d} = (d_P, d_Q, d_l)^T$，由 Cauchy-Schwarz 不等式得：

$$
\frac{\mathbf{d}^T \mathbf{x}_{ij}}{\|\mathbf{d}\|_2} \leq \|\mathbf{x}_{ij}\|_2 \leq l_{t,i,j} + v_{t,i}
$$

在方向 $\mathbf{d}$ 上添加**带收缩因子的方向割约束**：

$$
\frac{\mathbf{d}^T \mathbf{x}_{ij}}{\|\mathbf{d}\|_2} \geq r \cdot (l_{t,i,j} + v_{t,i}), \quad r \in [0, 1]
$$

其中 $r=1$ 时为最紧的割（仅保留方向 $\mathbf{d}$ 上的锥面点），$r<1$ 时保留方向 $\mathbf{d}$ 附近一定角度范围内的锥面点。

##### 层次化方向割锥还原算法

算法由三个阶段组成：初始化、主迭代循环和终止处理，以 SOCP 松弛解为起点，通过逐步添加方向割约束收紧可行域，驱动松弛解向锥面靠近。

**算法参数**：

| 参数 | 含义 | 默认值 |
|-----|-----|------|
| $\epsilon$ | 收敛容差 | $10^{-6}$ |
| $\epsilon_{\text{stag}}$ | 停滞检测阈值（最小改善量） | $10^{-4}$ |
| $W_{\text{stag}}$ | 停滞检测窗口大小 | 5 |
| $L_{\max}$ | 最大层数 | 3--5 |
| $K_{\max}$ | 最大迭代次数 | 50 |

**阶段一：初始化**

以初始 SOCP 解 $(P^*, Q^*, l^*, v^*)$ 为出发点，对每条支路 $(i,j)$ 提取当前锥表示向量方向 $\mathbf{d} = (2P^*_{t,i,j}, 2Q^*_{t,i,j}, l^*_{t,i,j}-v^*_{t,i})^T$，计算初始收缩因子 $r_{\min} = \|\mathbf{d}\| / (l^*_{t,i,j} + v^*_{t,i})$，将对应的方向割约束加入约束集合 $\mathcal{L}$，置层计数 $\ell=1$、迭代计数 $k=0$、间隙历史集合 $\mathcal{H}=\emptyset$。

**阶段二：主迭代循环**

每轮迭代包含以下步骤：

- **步骤 1（SOCP 求解）**：在当前割约束集合 $\mathcal{L}$ 下求解 SOCP，得到新解 $(P^k, Q^k, l^k, v^k)$
- **步骤 2（间隙计算）**：对所有支路和时段计算锥松弛间隙 $\Delta^k_{t,i,j}$，取最大间隙 $\Delta^k_{\max}$ 记录至历史集合 $\mathcal{H}$
- **步骤 3（全局收敛检验）**：若 $\Delta^k_{\max} \leq \epsilon$，则所有支路的松弛均已收紧，返回当前解并标记为收敛
- **步骤 4（层内自适应方向更新）**：对仍有间隙（$\Delta^k_{t,i,j} > \epsilon$）的支路，以当前解方向 $\hat{\mathbf{d}}$ 和当前层收缩因子 $r_\ell$ 生成新的方向割，补充进 $\mathcal{L}$。收缩因子按层线性递增：$r_\ell = r_{\min} + (\ell-1)/(L_{\max}-1) \cdot (1 - r_{\min})$
- **步骤 5--6（停滞检测与层次提升）**：若近 $W_{\text{stag}}$ 次迭代的最大间隙改善量均低于 $\epsilon_{\text{stag}}$，判定当前层已停滞。若 $\ell < L_{\max}$，则提升至下一层，增大收缩因子，生成新层全局方向集 $\mathcal{D}_\ell$ 并添加对应割约束，清空间隙历史重新监测；否则退出循环

**阶段三：终止处理**

若总迭代次数达到 $K_{\max}$ 仍未收敛，或层数达到 $L_{\max}$ 仍未收敛，返回当前最优解并标记为未完全收敛，供后续处理使用。

**算法流程图**：

```
初始化：提取初始解方向，计算 r_min，构建初始方向割集合 L
    ↓
while k < K_max:
    k ← k + 1
    求解 SOCP（含当前 L 中所有割约束）→ 得到 (P^k, Q^k, l^k, v^k)
    ↓
    计算所有支路时段的锥松弛间隙 Δ^k_{t,i,j}
    Δ^k_max ← max_{t,i,j} Δ^k_{t,i,j}；将 Δ^k_max 写入历史集合 H
    ↓
    若 Δ^k_max ≤ ε → 返回当前解，标记 Converged ✓
    ↓
    层内自适应：对有间隙的支路，以当前方向和 r_ℓ 生成新割，加入 L
    ↓
    若 |H| ≥ W_stag 且近期改善量 < ε_stag（停滞）：
        若 ℓ < L_max → ℓ += 1，生成全局方向集 D_ℓ 的割，清空 H
        否则        → break（层数耗尽）
    ↓
终止：返回当前解，标记 NotFullyConverged
```

##### 算法复杂度

**时间复杂度**：设决策变量数 $n = O(T \cdot |E|)$（辐射网），第 $k$ 次迭代时约束集合 $\mathcal{L}$ 已累积至多 $k \cdot T|E|$ 条割，对应 SOCP 求解代价为 $O(k \cdot n^3)$。全部 $K_{\max}$ 次迭代的总代价为：

$$
\sum_{k=1}^{K_{\max}} O(k \cdot n^3) = O(K_{\max}^2 \cdot n^3)
$$

由于 $K_{\max}$ 为小常数（默认 50），实践中可视为 $O(K_{\max}^2 \cdot n^3)$。

**空间复杂度**：割约束集合 $\mathcal{L}$ 以稀疏线性约束存储，每条割仅含 $O(1)$ 非零系数，$K_{\max}$ 轮后累积不超过 $K_{\max} \cdot T|E|$ 条，总空间复杂度为 $O(K_{\max} \cdot T \cdot |E|)$。

**典型迭代次数**：辐射网 SOCP 松弛通常已接近紧，实践中 10--30 次迭代即可收敛；硬上界为 $K_{\max}$（算法强制终止条件）。

---

#### 场景输入控制

本脚本（`SOCP_Mainnet_Solver_real_onetest.py`）采用**单组或少量场景**的测试模式，不需要传入完整的多场景组合，而是通过以下两个顶层参数灵活控制实际加载的场景：

```python
# ---- 场景输入控制参数（在 __main__ 段顶部设置）----
SCENARIO_INDEX = 0     # 指定运行哪一组场景（0-based 索引）；设为 None 时改用 N_SCENARIOS
N_SCENARIOS    = 3     # 当 SCENARIO_INDEX 为 None 时，从头取前 N_SCENARIOS 组场景运行
```

**逻辑说明：**

- 若 `SCENARIO_INDEX` 为非 `None` 整数，则只取 `all_scenarios[SCENARIO_INDEX]` 这一组数据，以单场景模式运行，便于调试和结果核查。
- 若 `SCENARIO_INDEX = None`，则取 `all_scenarios[:N_SCENARIOS]`，以少量场景批次模式运行，用于快速验证多场景逻辑的正确性。
- 两种模式均不需要修改场景生成逻辑，`generate_scenario_combinations` 仍完整生成全部组合，控制逻辑在主流程入口处截断。

**典型用法示例：**

| 目的 | 设置 |
|------|------|
| 只跑第 0 组（默认第一组）| `SCENARIO_INDEX = 0` |
| 只跑第 7 组（人工指定）| `SCENARIO_INDEX = 7` |
| 跑前 5 组做批量初步测试 | `SCENARIO_INDEX = None`, `N_SCENARIOS = 5` |

---

#### 输出目录

所有结果（`.npz` 变量文件、`summary.csv`、可视化图片）均保存在以下固定路径：

```
<脚本所在目录>/test_output/output_onetest/
```

若该目录不存在，脚本启动时自动调用 `os.makedirs(..., exist_ok=True)` 创建，无需手动操作。目录结构如下：

```
test_output/output_onetest/
├── summary.csv                         # 所有已求解场景的标量汇总（追加写入）
├── scenario_<name>/
│   ├── results.npz                     # 该场景所有决策变量最优值（压缩存储）
│   └── figures/
│       ├── voltage.png                 # 节点电压时序图
│       ├── capacity.png                # IDC 承载力时序图
│       ├── lambda.png                  # 任务到达量时序图
│       └── task_dc<k>.png             # 各 DC 任务处理量堆叠图
└── analysis/                           # 跨场景汇总分析（多场景模式时生成）
    ├── hc_distribution.png
    ├── hc_per_scenario.png
    ├── scale_vs_hc.png
    ├── voltage_violation.png
    └── report.csv
```

---

#### 优化求解模块

##### Gurobi 建模的约束分类

在使用 Python 调用 Gurobi 求解器构建上述 SOCP 模型时，需对各约束和变量进行工程化处理。以下从三个维度说明约束的处理方式。

**1. 应纳入 Gurobi 模型的约束（作为显式约束添加）**

以下约束不能通过变量界来表达，必须作为独立约束写入模型：

- **式(2)**：承载力下界约束 $hc \leq \sum_{i \in IDC} P_{IDC,i,t}$，对所有 $t$ 逐条添加线性约束
- **式(3)(4)**：节点功率平衡等式约束，对所有节点 $i$、时段 $t$ 添加等式线性约束
- **式(5)**：二阶锥松弛约束，使用 Gurobi 的 `addQConstr` 添加二次约束 $\|(2T_{ij,t},\; 2Z_{ij,t},\; U_{i,t}-U_{j,t})\|_2 \leq U_{i,t}+U_{j,t}$
- **式(6)**：辅助变量对称/反对称等式约束 $T_{ij,t}=T_{ji,t}$，$Z_{ij,t}=-Z_{ji,t}$，逐支路添加
- **式(7)**：相角参考约束（线性等式），对每个环 $k$、每个时段 $t$ 添加
- **式(8)(9)**：对角元素定义等式 $T_{ii,t}=U_{i,t}$，$Z_{ii,t}=0$，逐节点添加
- **式(17)**：任务到达量定义等式 $\lambda_{t,u}=\lambda_{t,u}^{\text{shape}} \cdot scale$，逐节点逐时段添加线性等式约束
- **式(18)(19)**：队列状态转移等式约束，逐节点逐时段逐延迟级添加
- **式(20)**：初始条件等式约束 $q^{(d)}_{0,u}=0$（$d\geq1$），逐节点逐级添加
- **式(21)**：缓存剩余量等式约束，逐节点逐时段逐延迟级添加
- **式(23)**：最高级缓存清空约束 $r^{(D)}_{t,u}=0$，逐节点逐时段添加等式约束
- **式(24)**：带宽约束（线性不等式），逐链路逐时段添加
- **式(26)**：周期性闭环等式约束 $q^{(d)}_{T-1,u}=q^{(d)}_{0,u}$，逐节点逐级添加
- **式(29)**：数据中心功耗等式约束，逐节点逐时段添加线性等式约束

**2. 冗余约束（不应单独添加，已被其他约束或变量界隐含）**

以下约束在 Gurobi 模型中**无需单独添加**：

- **式(14)**（变量关系定义）：$T_{ij,t}=V_{i,t}V_{j,t}\cos\theta_{ij,t}$ 等是辅助变量的物理含义说明，在 SOCP 松弛框架下 $T_{ij,t}$、$Z_{ij,t}$、$U_{i,t}$ 已作为独立决策变量处理，该等式关系由式(5)的二阶锥约束松弛替代，不再显式添加（添加会引入非线性，破坏 SOCP 结构）
- **式(15)(16)**（支路潮流定义）：$P_{ij,t}$、$Q_{ij,t}$ 若仅用于后处理（结果分析、可视化），则在优化阶段无需作为约束添加；若后续约束中不出现 $P_{ij,t}$、$Q_{ij,t}$，建议将其从决策变量中移除，在求解后通过辅助变量值计算得到
- **式(22)**（$r^{(d)}_{t,u}\geq 0$）：通过在变量声明时设置下界 $lb=0$ 即可隐含此约束，无需单独添加不等式约束
- **式(25)**（$p^{(d)}_{t,u}\geq 0$，$y^{(d)}_{t,u\to v}\geq 0$）：同上，通过变量下界 $lb=0$ 实现
- **式(27)**（$\lambda_{t,u}\geq 0$）：通过变量下界 $lb=0$ 实现

**3. 应在变量声明时通过上下界规定的范围**

以下量直接在 `model.addVar` 或 `model.addVars` 时设置 `lb`/`ub`，避免生成多余约束行：

| 变量 | 下界 `lb` | 上界 `ub` | 对应数学约束 |
|------|-----------|-----------|-------------|
| $U_{i,t}$ | $U_{i,\min}$ | $U_{i,\max}$ | 式(10) |
| $P_{G,i,t}$ | $P_{Gi,\min}$ | $P_{Gi,\max}$ | 式(11) |
| $Q_{G,i,t}$ | $Q_{Gi,\min}$ | $Q_{Gi,\max}$ | 式(12) |
| $r^{(d)}_{t,u}$ | $0$ | $+\infty$ | 式(22) |
| $p^{(d)}_{t,u}$ | $0$ | $\mu_{\max,u}$（可选） | 式(25) |
| $y^{(d)}_{t,u\to v}$ | $0$ | $+\infty$ | 式(25) |
| $\lambda_{t,u}$ | $0$ | $+\infty$ | 式(27) |
| $scale$ | $0$ | $+\infty$ | — |
| $hc$ | $0$ | $+\infty$ | — |

特别说明：式(13) $U_{0,t}=1$ 可通过将参考节点 $U_{0,t}$ 声明为 $lb=ub=1$ 的固定变量（或直接代入参数）实现，无需添加为约束。式(28) $\sum_d p^{(d)}_{t,u} \leq \mu_{\max,u}$ 为多变量线性不等式，不能仅通过单变量界实现，须作为显式约束添加。

---

#### 式(5) 松弛约束的松紧程度测试

式(5) 是本模型的核心约束，将原始非凸的交流潮流约束（$T_{ij}^2 + Z_{ij}^2 = U_i U_j$）松弛为二阶锥不等式：

$$
4T_{ij,t}^2 + 4Z_{ij,t}^2 + (U_{i,t}-U_{j,t})^2 \leq (U_{i,t}+U_{j,t})^2
$$

松弛是否"紧"（tight）决定了 SOCP 最优解是否对应真实的物理可行解。松弛紧意味着最优解处不等号几乎取等，即辅助变量满足原始等式约束，SOCP 解可直接恢复为交流潮流解。

##### 松紧度量指标

对每条支路 $(i,j)$、每个时段 $t$，定义**松弛间隙**（relaxation gap）：

$$
\delta_{ij,t} = \frac{(U_{i,t}+U_{j,t})^2 - [4T_{ij,t}^2 + 4Z_{ij,t}^2 + (U_{i,t}-U_{j,t})^2]}{(U_{i,t}+U_{j,t})^2}
$$

$\delta_{ij,t} \in [0,1]$，越接近 0 表示松弛越紧：

| 指标 | 含义 |
|------|------|
| $\max_{i,j,t} \delta_{ij,t}$ | 全网最差松弛程度 |
| $\text{mean}_{i,j,t} \delta_{ij,t}$ | 全网平均松弛程度 |
| 支路级 $\delta_{ij}^{\max} = \max_t \delta_{ij,t}$ | 各支路最差时段松弛 |

**判断标准**：若 $\max \delta_{ij,t} < 10^{-4}$，视为松弛严格紧；$10^{-4}$ 到 $10^{-2}$ 为近似紧；超过 $10^{-2}$ 需重点关注。

##### 松紧度测试流程

```
求解完成（model.Status == OPTIMAL 或 SUBOPTIMAL）
    ↓
后处理阶段：提取 T_var、Z_var、U_var 最优值
    ↓
for each branch (i,j) in E:
    for t in range(T):
        lhs = 4*T[i,j,t]**2 + 4*Z[i,j,t]**2 + (U[i,t]-U[j,t])**2
        rhs = (U[i,t]+U[j,t])**2
        gap[i,j,t] = (rhs - lhs) / rhs
    ↓
统计 gap 的最大值、均值、各支路最大值
    ↓
打印松紧度报告（含超阈值支路列表）
    ↓
生成松紧度热力图：以支路编号为横轴、时段 t 为纵轴，颜色映射 gap 大小
保存至 output_dir/scenario_<name>/figures/socp_gap.png
```

##### 松紧度输出内容

求解完成后，在控制台打印以下报告：

```
========== 式(5) SOCP 松弛紧程度报告 ==========
  最大松弛间隙 (全支路全时段): 3.21e-06
  平均松弛间隙               : 8.45e-09
  超过阈值(1e-4)的支路数量   : 0 / 38
  结论: 松弛严格紧，SOCP 解对应物理可行的交流潮流解
================================================
```

若存在超阈值支路，额外输出该支路的支路编号、节点对和最差时段信息，便于排查网络拓扑或参数异常。

---

#### 结果保存、读取与可视化模块

优化求解完成后，结果处理分为三个独立子模块：

**结果保存模块 `save_results(scenario_id, result_arrays, result_meta, output_dir)`**

从已提取的变量数组中将结果持久化存储：

- 存储格式：以场景名称为目录名，将变量数组保存至 `output_dir/scenario_{name}/results.npz`（`numpy.savez_compressed`，兼顾读写速度与压缩比）；同时将标量结果（$hc$、$scale$、求解时间、目标值）追加写入汇总文件 `output_dir/summary.csv`
- 线程安全：汇总 CSV 的写入使用文件锁（`filelock.FileLock`），避免并行场景同时写入时发生竞态

**结果读取模块 `load_results(scenario_id, output_dir, scenario_str)`**

从磁盘加载指定场景的结果：

- 从 `results.npz` 中反序列化各变量数组，重建为与求解阶段一致的索引结构（节点、时段、延迟级）
- 从 `summary.csv` 中读取标量汇总信息
- 返回 `(result_arrays, result_meta)` 二元组，供分析和可视化模块调用

**可视化模块 `visualize_results(scenario_id, result_arrays, result_meta, output_dir)`**

对单场景结果生成标准图表并保存：

- 图1 `voltage.png`：各典型节点 $V_{i,t}=\sqrt{U_{i,t}}$ 时序折线图，标注 0.9/1.1 上下限
- 图2 `capacity.png`：总 IDC 功率与 $hc$ 的时序对比（上子图）+ 各 DC 功率分解（下子图）
- 图3 `lambda.png`：各 DC 最优任务到达量 $\lambda_{k,t}$ 时序折线图
- 图4 `task_dc<k>.png`：各 DC 按延迟级 $d$ 堆叠的任务处理量柱状图
- 图5 `socp_gap.png`：式(5) 松弛间隙 $\delta_{ij,t}$ 热力图（支路 × 时段）
- 图片保存：使用 `matplotlib` 的 `Agg` 后端，调用 `fig.savefig(path, dpi=150, bbox_inches='tight')` 后立即 `plt.close(fig)`，避免内存泄漏；所有图片保存至 `output_dir/scenario_{name}/figures/`

---

#### 单场景求解流水线

```
__main__ 入口
    ↓
读取参数：SCENARIO_INDEX / N_SCENARIOS
    ↓
generate_scenario_combinations(...)        # 生成全量场景（笛卡尔积）
    ↓
按 SCENARIO_INDEX 或 N_SCENARIOS 截取      # 只保留待测场景
    ↓
for scenario in selected_scenarios:
    read_scenario_data(...)                # 从 CSV 读取单组时序数据
    build_and_solve_single(...)            # 构建 Gurobi SOCP 模型并求解
        ↓
    [求解成功]
        → 提取变量最优值（.X）
        → 计算式(5)松弛间隙（delta_arr）
        → 打印松紧度报告
        → save_results(...)               # 写 .npz + 追加 summary.csv
        → visualize_results(...)          # 生成 5 张图（含 socp_gap.png）
    [求解失败]
        → computeIIS() + 写 .ilp 文件
        → 打印失败原因，继续下一场景
    ↓
（可选）analyze_all_scenarios(output_dir)   # 多场景模式时运行跨场景汇总
```

单场景模式下不启动进程池，直接在主进程内串行执行，避免调试时的进程间通信开销。多场景模式（`N_SCENARIOS > 1`）仍支持 `ProcessPoolExecutor` 并行加速，与 `SOCP_Mainnet_Solver_real.py` 的并发策略保持一致。
