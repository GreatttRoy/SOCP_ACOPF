# 二阶锥交流潮流优化模型 (SOCP AC OPF)

## 1. 优化模型

$$
\begin{align}
\min_{P_{t,i,j}, Q_{t,i,j}, l_{t,i,j}, v_{t,i}, p_{t,i}, q_{t,i}, P_{g,t,g}, Q_{g,t,g}, P_{DR,t,dr}} \quad & \sum_{t=0}^{T-1} \sum_{(i,j) \in E} l_{t,i,j} \\
\text{s.t.} \quad \\
& v_{t,0} = 1.0, \quad \forall t \in \{0,1,\ldots,T-1\} \tag{1} \\
\\
& p_{t,j} = \sum_{g \in G_j} P_{g,t,g} + \sum_{dr \in DR_j} P_{DR,t,dr} + P_{d,t,j}, \quad \forall t, j \in N \tag{2a} \\
\\
& q_{t,j} = \sum_{g \in G_j} Q_{g,t,g} + Q_{d,t,j}, \quad \forall t, j \in N \tag{2b} \\
\\
& p_{t,j} = \sum_{(j,k) \in E} P_{t,j,k} - \sum_{(i,j) \in E} \left( P_{t,i,j} - r_{ij} \cdot l_{t,i,j} \right), \quad \forall t, j \in N \tag{3a} \\
\\
& q_{t,j} = \sum_{(j,k) \in E} Q_{t,j,k} - \sum_{(i,j) \in E} \left( Q_{t,i,j} - x_{ij} \cdot l_{t,i,j} \right), \quad \forall t, j \in N \tag{3b} \\
\\
& v_{t,j} = v_{t,i} - 2\left( r_{ij} P_{t,i,j} + x_{ij} Q_{t,i,j} \right) + \left( r_{ij}^2 + x_{ij}^2 \right) l_{t,i,j}, \quad \forall t, (i,j) \in E \tag{4} \\
\\
& \left(2P_{t,i,j}\right)^2 + \left(2Q_{t,i,j}\right)^2 + \left(l_{t,i,j} - v_{t,i}\right)^2 \leq \left(l_{t,i,j} + v_{t,i}\right)^2, \quad \forall t, (i,j) \in E \tag{5} \\
\\
& P_{g,\min}(t) \leq P_{g,t,g} \leq P_{g,\max}(t), \quad \forall t, g \in G \tag{6} \\
\\
& P_{DR,\min}(t) \leq P_{DR,t,dr} \leq P_{DR,\max}(t), \quad \forall t, dr \in DR \tag{7} \\
\\
& V_{\min}^2 \leq v_{t,i} \leq V_{\max}^2, \quad \forall t, i \in N \tag{8} \\
\\
& -2.5 \leq P_{t,i,j} \leq 2.5, \quad \forall t, (i,j) \in E \tag{9a} \\
\\
& -2.5 \leq Q_{t,i,j} \leq 2.5, \quad \forall t, (i,j) \in E \tag{9b} \\
\\
& 0 \leq l_{t,i,j} \leq 2.5, \quad \forall t, (i,j) \in E \tag{9c}
\end{align}
$$

---

## 2. 决策变量

### 2.1 支路变量
- $P_{t,i,j}$：支路 $(i,j)$ 在时段 $t$ 的有功功率 (p.u.)
- $Q_{t,i,j}$：支路 $(i,j)$ 在时段 $t$ 的无功功率 (p.u.)
- $l_{t,i,j}$：支路 $(i,j)$ 在时段 $t$ 的电流幅值平方 $|I_{ij}|^2$ (p.u.)

### 2.2 节点变量
- $v_{t,i}$：节点 $i$ 在时段 $t$ 的电压幅值平方 $|V_i|^2$ (p.u.)
- $p_{t,i}$：节点 $i$ 在时段 $t$ 的有功注入功率 (p.u.)
- $q_{t,i}$：节点 $i$ 在时段 $t$ 的无功注入功率 (p.u.)

### 2.3 发电机变量
- $P_{g,t,g}$：发电机 $g$ 在时段 $t$ 的有功出力 (p.u.)
- $Q_{g,t,g}$：发电机 $g$ 在时段 $t$ 的无功出力 (p.u.)

### 2.4 可调负荷变量
- $P_{DR,t,dr}$：可调负荷 $dr$ 在时段 $t$ 的削减量 (p.u.)（正数表示削减）

---

## 3. 模型参数

### 3.1 集合
- $N$：节点集合，$|N| = $ num_n（节点总数）
- $E$：支路集合，每条支路表示为 $(i, j, r_{ij}, x_{ij})$
- $G$：发电机集合，$|G| = $ num_gen（发电机总数）
- $DR$：可调负荷集合，$|DR| = $ num_dr（可调负荷总数）
- $G_j$：连接在节点 $j$ 的发电机集合
- $DR_j$：连接在节点 $j$ 的可调负荷集合

### 3.2 网络参数
- $r_{ij}$：支路 $(i,j)$ 的电阻 (p.u.)
- $x_{ij}$：支路 $(i,j)$ 的电抗 (p.u.)

### 3.3 时间参数
- $T$：时间断面数（优化的时段总数）
- $t$：时段索引，$t \in \{0, 1, \ldots, T-1\}$

### 3.4 负荷参数
- $P_{d,t,j}$：节点 $j$ 在时段 $t$ 的固定有功负荷 (p.u.)（负值表示消耗）
- $Q_{d,t,j}$：节点 $j$ 在时段 $t$ 的固定无功负荷 (p.u.)（负值表示消耗）

### 3.5 发电机参数
- $P_{g,\min}(t)$：发电机 $g$ 在时段 $t$ 的最小有功出力 (p.u.)
- $P_{g,\max}(t)$：发电机 $g$ 在时段 $t$ 的最大有功出力 (p.u.)

### 3.6 可调负荷参数
- $P_{DR,\min}(t)$：可调负荷 $dr$ 在时段 $t$ 的最小削减量 (p.u.)
- $P_{DR,\max}(t)$：可调负荷 $dr$ 在时段 $t$ 的最大削减量 (p.u.)

### 3.7 电压限制参数
- $V_{\min} = 0.9$ p.u.：节点电压幅值下限
- $V_{\max} = 1.1$ p.u.：节点电压幅值上限

---

## 4. 约束条件说明

**(1) 根节点电压约束**
- 节点 0 为电源节点（平衡节点），电压幅值恒定为 1.0 p.u.

**(2a-2b) 节点净注入功率约束**
- 定义节点的有功和无功净注入功率，包含发电机出力、可调负荷削减和固定负荷

**(3a-3b) 潮流平衡约束（Branch Flow Model）**
- 节点的净注入功率等于流出功率减去流入功率（考虑线路损耗）
- $\sum_{(j,k) \in E}$：从节点 $j$ 流出的所有支路
- $\sum_{(i,j) \in E}$：流入节点 $j$ 的所有支路
- $r_{ij} \cdot l_{t,i,j}$：支路 $(i,j)$ 的有功损耗
- $x_{ij} \cdot l_{t,i,j}$：支路 $(i,j)$ 的无功损耗

**(4) 电压降落约束（DistFlow Model）**
- 描述沿支路的电压降落关系
- 考虑了电阻和电抗引起的电压降以及线路损耗的影响

**(5) 二阶锥松弛约束**
- 对原始非凸约束 $l_{t,i,j} \cdot v_{t,i} = P_{t,i,j}^2 + Q_{t,i,j}^2$ 的凸松弛
- 等价形式：$\left\| \begin{bmatrix} 2P_{t,i,j} \\ 2Q_{t,i,j} \\ l_{t,i,j} - v_{t,i} \end{bmatrix} \right\|_2 \leq l_{t,i,j} + v_{t,i}$
- 在辐射网络中，该松弛通常是紧的（即可以得到原问题的精确解）

**(6) 发电机出力约束**
- 限制发电机的有功出力在其容量范围内
- 上下限可随时段 $t$ 变化（时变约束）

**(7) 可调负荷约束**
- 限制可调负荷的削减量在允许范围内
- 上下限可随时段 $t$ 变化（时变约束）

**(8) 电压上下限约束**
- 确保所有节点的电压幅值在安全运行范围内

**(9a-9c) 支路变量界约束**
- 限制支路有功、无功功率和电流平方的取值范围

---

## 5. 潮流结果计算公式

在求解优化问题得到最优解后，可以通过以下公式恢复相角信息并计算完整的潮流结果。

### 5.1 节点相角计算

对于时段 $t$，首先计算节点相角向量 $\boldsymbol{\theta}_t$（单位：弧度）：

$$
\boldsymbol{\theta}_t = (A_t^T)^{-1} \boldsymbol{\beta}_t
$$

其中：
- $A_t \in \mathbb{R}^{|N| \times |E|}$：节点-支路关联矩阵（按国内教材惯例需取转置）
- $\boldsymbol{\beta}_t \in \mathbb{R}^{|E|}$：支路相角差向量，通过相角恢复条件计算得到
- $\theta_i$：节点 $i$ 的相角（弧度）

**相角单位转换：**
- 弧度转角度：$\theta_{\text{度}} = \theta_{\text{弧度}} \times \frac{180}{\pi}$
- 角度转弧度：$\theta_{\text{弧度}} = \theta_{\text{度}} \times \frac{\pi}{180}$

### 5.2 支路功率相量

对于支路 $(i,j)$，支路复功率为：

$$
S_{ij} = P_{t,i,j}^* + j Q_{t,i,j}^*
$$

其中：
- $P_{t,i,j}^*, Q_{t,i,j}^*$：优化问题的最优解
- $j = \sqrt{-1}$：虚数单位

支路功率的相角：
$$
\angle S_{ij} = \arctan\left(\frac{Q_{t,i,j}^*}{P_{t,i,j}^*}\right) \quad \text{（弧度）}
$$

### 5.3 支路电流相量

支路电流的幅值和相角为：

$$
I_{ij} = \sqrt{l_{t,i,j}^*} \cdot e^{j(\theta_i - \angle S_{ij})}
$$

其中：
- $l_{t,i,j}^*$：优化问题得到的支路电流幅值平方
- $\theta_i$：节点 $i$ 的相角（弧度）
- $\angle S_{ij}$：支路功率的相角（弧度）

**电流的幅值和相角：**
- 幅值：$|I_{ij}| = \sqrt{l_{t,i,j}^*}$ (p.u.)
- 相角：$\angle I_{ij} = \theta_i - \angle S_{ij}$ （弧度）

### 5.4 节点电压相量

对于节点 $j$（除平衡节点外），节点电压相量为：

$$
V_j = \sqrt{v_{t,j}^*} \cdot e^{j\theta_j}
$$

其中：
- $v_{t,j}^*$：优化问题得到的节点电压幅值平方
- $\theta_j$：节点 $j$ 的相角（弧度）

**电压的幅值和相角：**
- 幅值：$|V_j| = \sqrt{v_{t,j}^*}$ (p.u.)
- 相角：$\angle V_j = \theta_j$ （弧度）

### 5.5 节点注入复功率

节点 $j$ 的注入复功率为：

$$
s_{0,j} = p_{t,j}^* + j q_{t,j}^*
$$

其中 $p_{t,j}^*, q_{t,j}^*$ 为优化问题的最优解。

### 5.6 计算流程

**步骤 1：求解优化问题**

获得所有决策变量的最优值：
- $P_{t,i,j}^*, Q_{t,i,j}^*, l_{t,i,j}^*$（支路变量）
- $v_{t,i}^*, p_{t,i}^*, q_{t,i}^*$（节点变量）
- $P_{g,t,g}^*, Q_{g,t,g}^*$（发电机变量）
- $P_{DR,t,dr}^*$（可调负荷变量）

**步骤 2：相角恢复**

通过相角恢复条件计算支路相角差向量 $\boldsymbol{\beta}_t$，进而计算节点相角：
$$
\boldsymbol{\theta}_t = (A_t^T)^{-1} \boldsymbol{\beta}_t
$$

**步骤 3：计算相量**

- **支路功率相量：**
  $$S_{ij} = P_{t,i,j}^* + j Q_{t,i,j}^*$$

- **支路电流相量：**
  $$I_{ij} = \sqrt{l_{t,i,j}^*} \cdot e^{j(\theta_i - \angle S_{ij})}$$

- **节点电压相量：**
  $$V_j = \sqrt{v_{t,j}^*} \cdot e^{j\theta_j}$$

- **节点注入复功率：**
  $$s_{0,j} = p_{t,j}^* + j q_{t,j}^*$$

**步骤 4：提取幅值和相角**

使用复数运算提取相量的幅值和相角：
- 幅值：$|z| = \sqrt{\text{Re}(z)^2 + \text{Im}(z)^2}$
- 相角：$\angle z = \arctan\left(\frac{\text{Im}(z)}{\text{Re}(z)}\right)$

---

## 6. 模型特点

1. **凸优化问题**：通过二阶锥松弛，将非凸的交流潮流问题转化为凸优化问题，可以高效求解并获得全局最优解

2. **精确性**：在辐射状配电网络中，当满足相角恢复条件时，二阶锥松弛可以得到原问题的精确解

3. **多时段优化**：考虑多个时间断面的联合优化，支持时变的负荷、发电机出力限制等

4. **适用场景**：适用于辐射状配电网络的最优潮流计算，可用于配电网规划、运行优化等场景

5. **模型依据**：基于 Branch Flow Model (BFM) 和二阶锥松弛 (SOCP relaxation) 技术

---

**注**：该模型实现基于 Gurobi 优化器求解，使用了高精度参数设置（`MIPGap=1e-8`, `OptimalityTol=1e-8`, `FeasibilityTol=1e-8`）以确保求解精度。
