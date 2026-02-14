# 二阶锥交流潮流优化模型 (SOCP-AC-OPF-Model)

## 目录

1. [优化模型](#1-优化模型)
2. [决策变量](#2-决策变量)
3. [模型参数](#3-模型参数)
4. [约束条件说明](#4-约束条件说明)
5. [潮流结果计算公式](#5-潮流结果计算公式)
6. [模型特点](#6-模型特点)
7. [基于方向割的锥还原方法](#7-基于方向割的锥还原方法)
8. [模型特点（修订）](#8-模型特点修订)

---

## 符号表

### 索引与集合

| 符号 | 描述 | 说明 |
|------|------|------|
| $t$ | 时段索引 | $t \in \{0, 1, \ldots, T-1\}$ |
| $i, j, k$ | 节点索引 | $i, j, k \in N$ |
| $g$ | 发电机索引 | $g \in G$ |
| $dr$ | 可调负荷索引 | $dr \in DR$ |
| $N$ | 节点集合 | 配电网所有节点的集合，$\|N\|$ = 节点总数 |
| $E$ | 支路集合 | 配电网所有支路的集合，每条支路表示为 $(i, j, r_{ij}, x_{ij})$ |
| $G$ | 发电机集合 | 所有发电机的集合，$\|G\|$ = 发电机总数 |
| $DR$ | 可调负荷集合 | 所有可调负荷的集合，$\|DR\|$ = 可调负荷总数 |
| $G_j$ | 节点 $j$ 的发电机集合 | 连接在节点 $j$ 的所有发电机 |
| $DR_j$ | 节点 $j$ 的可调负荷集合 | 连接在节点 $j$ 的所有可调负荷 |

### 决策变量

#### 支路变量

| 符号 | 描述 | 单位 |
|------|------|------|
| $P_{t,i,j}$ | 支路 $(i,j)$ 在时段 $t$ 的有功功率 | p.u. |
| $Q_{t,i,j}$ | 支路 $(i,j)$ 在时段 $t$ 的无功功率 | p.u. |
| $l_{t,i,j}$ | 支路 $(i,j)$ 在时段 $t$ 的电流幅值平方 $\|I_{ij}\|^2$ | p.u. |

#### 节点变量

| 符号 | 描述 | 单位 |
|------|------|------|
| $v_{t,i}$ | 节点 $i$ 在时段 $t$ 的电压幅值平方 $\|V_i\|^2$ | p.u. |
| $p_{t,i}$ | 节点 $i$ 在时段 $t$ 的有功注入功率 | p.u. |
| $q_{t,i}$ | 节点 $i$ 在时段 $t$ 的无功注入功率 | p.u. |

#### 发电机变量

| 符号 | 描述 | 单位 |
|------|------|------|
| $P_{g,t,g}$ | 发电机 $g$ 在时段 $t$ 的有功出力 | p.u. |
| $Q_{g,t,g}$ | 发电机 $g$ 在时段 $t$ 的无功出力 | p.u. |

#### 可调负荷变量

| 符号 | 描述 | 单位 |
|------|------|------|
| $P_{DR,t,dr}$ | 可调负荷 $dr$ 在时段 $t$ 的削减量（正数表示削减） | p.u. |

### 网络参数

| 符号 | 描述 | 单位 |
|------|------|------|
| $r_{ij}$ | 支路 $(i,j)$ 的电阻 | p.u. |
| $x_{ij}$ | 支路 $(i,j)$ 的电抗 | p.u. |
| $V_{\min}$ | 节点电压幅值下限 | p.u.（默认 0.9） |
| $V_{\max}$ | 节点电压幅值上限 | p.u.（默认 1.1） |

### 时间参数

| 符号 | 描述 |
|------|------|
| $T$ | 时间断面数（优化的时段总数） |
| $t$ | 时段索引，$t \in \{0, 1, \ldots, T-1\}$ |

### 负荷参数

| 符号 | 描述 | 单位 |
|------|------|------|
| $P_{d,t,j}$ | 节点 $j$ 在时段 $t$ 的固定有功负荷（负值表示消耗） | p.u. |
| $Q_{d,t,j}$ | 节点 $j$ 在时段 $t$ 的固定无功负荷（负值表示消耗） | p.u. |

### 发电机参数

| 符号 | 描述 | 单位 |
|------|------|------|
| $P_{g,\min}(t)$ | 发电机 $g$ 在时段 $t$ 的最小有功出力 | p.u. |
| $P_{g,\max}(t)$ | 发电机 $g$ 在时段 $t$ 的最大有功出力 | p.u. |

### 可调负荷参数

| 符号 | 描述 | 单位 |
|------|------|------|
| $P_{DR,\min}(t)$ | 可调负荷 $dr$ 在时段 $t$ 的最小削减量 | p.u. |
| $P_{DR,\max}(t)$ | 可调负荷 $dr$ 在时段 $t$ 的最大削减量 | p.u. |

### 相角恢复相关符号（第5节）

| 符号 | 描述 | 单位 |
|------|------|------|
| $\boldsymbol{\theta}_t$ | 时段 $t$ 的节点相角向量 | 弧度 |
| $\theta_i$ | 节点 $i$ 的相角 | 弧度 |
| $A_t$ | 节点-支路关联矩阵 | $A_t \in \mathbb{R}^{\|N\| \times \|E\|}$ |
| $\boldsymbol{\beta}_t$ | 支路相角差向量 | $\boldsymbol{\beta}_t \in \mathbb{R}^{\|E\|}$，弧度 |

### 相量相关符号（第5节）

| 符号 | 描述 |
|------|------|
| $S_{ij}$ | 支路 $(i,j)$ 的复功率相量，$S_{ij} = P_{t,i,j} + jQ_{t,i,j}$ |
| $I_{ij}$ | 支路 $(i,j)$ 的电流相量 |
| $V_j$ | 节点 $j$ 的电压相量 |
| $s_{0,j}$ | 节点 $j$ 的注入复功率，$s_{0,j} = p_{t,j} + jq_{t,j}$ |
| $j$ | 虚数单位，$j = \sqrt{-1}$ |
| $\angle z$ | 复数 $z$ 的相角（辐角） |

### 锥还原相关符号（第7节）

| 符号 | 描述 | 单位 |
|------|------|------|
| $\Delta_{t,i,j}$ | 支路 $(i,j)$ 在时段 $t$ 的锥松弛间隙 | p.u.$^2$ |
| $\mathbf{x}_{ij}$ | 支路功率向量，$\mathbf{x}_{ij} = (P_{t,i,j}, Q_{t,i,j})^T$ | p.u. |
| $\mathbf{d}$ | 方向向量，$\mathbf{d} = (d_P, d_Q)^T$ | 无量纲 |
| $r$ | 收缩因子，$r \in [0, 1]$ | 无量纲 |
| $r_{\min}$ | 初始收缩因子 | 无量纲（默认 0.5） |
| $r_k$ | 第 $k$ 层的收缩因子 | 无量纲 |
| $\mathcal{L}$ | 线性割约束集合 | - |
| $\mathcal{D}_\ell$ | 第 $\ell$ 层的方向集合 | - |
| $L_{\max}$ | 最大层数 | 正整数（默认 3-5） |
| $K_{\max}$ | 最大迭代次数 | 正整数（默认 50） |
| $\epsilon$ | 收敛容差 | p.u.$^2$（默认 $10^{-6}$） |

### 优化求解相关符号

| 符号 | 描述 |
|------|------|
| $x^*$ | 最优解（带上标 $*$ 表示优化问题的最优值） |
| $k$ | 迭代计数器 |
| $\ell$ | 层次计数器 |

### 算法参数

| 符号 | 描述 | 典型值 |
|------|------|--------|
| `MIPGap` | 混合整数规划相对间隙容差 | $10^{-8}$ |
| `OptimalityTol` | 最优性容差 | $10^{-8}$ |
| `FeasibilityTol` | 可行性容差 | $10^{-8}$ |

### 常用数学符号

| 符号 | 描述 |
|------|------|
| $\|\|\cdot\|\|_2$ | 欧几里得范数（2-范数） |
| $\mathbb{R}$ | 实数集 |
| $\mathbb{R}^n$ | $n$ 维实向量空间 |
| $\mathbb{R}^{m \times n}$ | $m \times n$ 实矩阵空间 |
| $\forall$ | 对于所有（全称量词） |
| $\in$ | 属于 |
| $\leq, \geq$ | 小于等于，大于等于 |
| $\Leftrightarrow$ | 等价于 |
| $\Rightarrow$ | 推出 |

---

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

## 7. 基于方向割的锥还原方法

### 7.1 问题背景

在二阶锥交流潮流优化模型中，二阶锥松弛约束 (5) 对原始非凸约束进行了凸松弛：

$$
l_{t,i,j} \cdot v_{t,i} = P_{t,i,j}^2 + Q_{t,i,j}^2 \quad \Rightarrow \quad \left\| \begin{bmatrix} 2P_{t,i,j} \\ 2Q_{t,i,j} \\ l_{t,i,j} - v_{t,i} \end{bmatrix} \right\|_2 \leq l_{t,i,j} + v_{t,i}
$$

**松弛的几何意义**：
- 原始约束要求点 $(P_{t,i,j}, Q_{t,i,j}, l_{t,i,j}, v_{t,i})$ 必须位于**二阶锥的表面**（锥面约束）
- 松弛后的约束允许点位于**二阶锥的内部或表面**（锥约束）

在辐射状配电网络中，该松弛通常是紧的（tight），即最优解自然满足等式约束。但在某些情况下，松弛解可能落在锥的内部，此时需要进行**锥还原**，将解投影到锥面上以获得原问题的可行解。

### 7.2 锥松弛间隙检测

设 SOCP 求解得到的最优解为 $(P^*_{t,i,j}, Q^*_{t,i,j}, l^*_{t,i,j}, v^*_{t,i})$，定义**锥松弛间隙**为：

$$
\Delta_{t,i,j} = l^*_{t,i,j} \cdot v^*_{t,i} - \left[(P^*_{t,i,j})^2 + (Q^*_{t,i,j})^2\right]
$$

**间隙判定**：
- 若 $\Delta_{t,i,j} \approx 0$（在数值误差范围内），则松弛是紧的，无需还原
- 若 $\Delta_{t,i,j} > \epsilon$（给定容差），则需要进行锥还原

### 7.3 基于方向的锥还原方法

当检测到松弛间隙时，采用基于方向的迭代方法将解还原到锥面上。该方法基于 **Cauchy-Schwarz 不等式**和**方向割约束**。

#### 7.3.1 Cauchy-Schwarz 不等式

对于支路 $(i,j)$，记 $\mathbf{x}_{ij} = (P_{t,i,j}, Q_{t,i,j})^T$，则二阶锥约束为：

$$
\|\mathbf{x}_{ij}\|_2 \leq \sqrt{l_{t,i,j} \cdot v_{t,i}}
$$

对于任意方向向量 $\mathbf{d} = (d_P, d_Q)^T$，由 Cauchy-Schwarz 不等式可得：

$$
\frac{\mathbf{d}^T \mathbf{x}_{ij}}{\|\mathbf{d}\|_2} \leq \|\mathbf{x}_{ij}\|_2 \leq \sqrt{l_{t,i,j} \cdot v_{t,i}}
$$

这表明：沿任意方向 $\mathbf{d}$ 的投影不超过 $\sqrt{l_{t,i,j} \cdot v_{t,i}}$。

#### 7.3.2 方向割约束

为了收紧可行域，在方向 $\mathbf{d}$ 上添加**带收缩因子的方向割**：

$$
\frac{\mathbf{d}^T \mathbf{x}_{ij}}{\|\mathbf{d}\|_2} \geq r \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}}
$$

其中 $r \in [0, 1]$ 为收缩因子：
- $r = 1$：最紧的割，仅保留方向 $\mathbf{d}$ 上的锥面点
- $r < 1$：较松的割，保留方向 $\mathbf{d}$ 附近一定角度范围内的锥面点

**几何意义（二维情况）**：

设方向 $\mathbf{d}$ 对应角度 $\alpha$，锥面上的点 $\mathbf{x}_{ij}$ 对应角度 $\theta$，则约束要求：

$$
\cos(\theta - \alpha) \geq r
$$

因此可行角度范围为：

$$
\theta \in [\alpha - \arccos(r), \alpha + \arccos(r)]
$$

**收缩因子的选择**：
- $r = 1$：保留范围 0°（仅方向 $\mathbf{d}$ 本身）
- $r = 1/\sqrt{2} \approx 0.707$：保留范围 90°（±45°）
- $r = 1/2$：保留范围 120°（±60°）
- $r = 0$：保留范围 180°（半圆）

#### 7.3.3 层次化方向细分策略

采用**从粗到细**的层次化策略，逐步逼近锥面。该策略的核心思想是：

**策略原则**：
1. **固定收缩因子 $r$**：在同一层内保持 $r$ 不变
2. **自适应搜索方向 $\mathbf{d}$**：通过迭代在当前解附近自适应地选择有效方向
3. **停滞判断**：当同一层内的方向搜索无法显著减小锥松弛间隙时，进入下一层
4. **层次提升**：增大 $r$ 值，提高割的约束强度，并添加更密集的全局方向集

**重要说明**：这里的"层"不是严格的顺序结构，而是**约束强度等级**。在每一层内会进行多轮迭代，每次迭代自适应地选择方向并添加割约束。只有当层内搜索停滞时，才提升到下一层。

---

**第一层（$\ell = 1$）：粗略约束阶段**

- **收缩因子**：$r_1 = r_{\min} \approx 0.5$
- **初始全局方向集** $\mathcal{D}_1$（象限等分线）：

对于二维功率空间 $(P, Q)$，使用四个基本方向：

$$
\begin{align}
\mathbf{d}^{(1)} &= (1, 1)^T & \text{(45°)} \\
\mathbf{d}^{(2)} &= (-1, 1)^T & \text{(135°)} \\
\mathbf{d}^{(3)} &= (-1, -1)^T & \text{(225°)} \\
\mathbf{d}^{(4)} &= (1, -1)^T & \text{(315°)}
\end{align}
$$

对应的约束为：

$$
\begin{align}
\frac{P_{t,i,j} + Q_{t,i,j}}{\sqrt{2}} &\geq r_1 \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}} \\
\frac{-P_{t,i,j} + Q_{t,i,j}}{\sqrt{2}} &\geq r_1 \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}} \\
\frac{-P_{t,i,j} - Q_{t,i,j}}{\sqrt{2}} &\geq r_1 \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}} \\
\frac{P_{t,i,j} - Q_{t,i,j}}{\sqrt{2}} &\geq r_1 \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}}
\end{align}
$$

- **层内迭代**：
  - 添加初始方向集 $\mathcal{D}_1$ 的所有割约束
  - 求解 SOCP 问题，得到解 $\mathbf{x}^{(k)}$
  - 识别当前解的方向 $\hat{\mathbf{d}}^{(k)}$
  - 在 $\hat{\mathbf{d}}^{(k)}$ 附近自适应生成新方向 $\mathbf{d}^{new}$
  - 添加新的方向割（使用固定的 $r_1$）
  - 重复迭代，直至收敛或停滞

- **停滞判断**：若连续多次迭代（如 3-5 次）锥松弛间隙的改进小于阈值（如 $10^{-7}$），则认为在当前层停滞

---

**第二层（$\ell = 2$）：中等约束阶段**

- **收缩因子**：$r_2 = r_{\min} + \frac{1}{L_{\max}-1} \cdot (1 - r_{\min}) \approx 0.75$
- **全局方向集** $\mathcal{D}_2$（象限细分方向）：

在相邻方向之间插入新方向，对于第一象限：

$$
\begin{align}
\mathbf{d}^{(1,1)} &= (2, 1)^T & \text{(约26.6°)} \\
\mathbf{d}^{(1,2)} &= (1, 2)^T & \text{(约63.4°)}
\end{align}
$$

类似地在其他三个象限也插入方向，共 8-12 个方向。

对应的约束示例：

$$
\begin{align}
\frac{2P_{t,i,j} + Q_{t,i,j}}{\sqrt{5}} &\geq r_2 \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}} \\
\frac{P_{t,i,j} + 2Q_{t,i,j}}{\sqrt{5}} &\geq r_2 \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}}
\end{align}
$$

- **层内迭代**：与第一层类似，在固定 $r_2$ 下进行自适应方向搜索

---

**第 $\ell$ 层：一般化策略**

- **收缩因子递增公式**：

$$
r_\ell = r_{\min} + \frac{\ell-1}{L_{\max}-1} \cdot (1 - r_{\min}), \quad \ell = 1, 2, \ldots, L_{\max}
$$

其中：
- $L_{\max} \in [3, 5]$：最大层数
- $r_{\min} \in [0.5, 0.7]$：初始收缩因子
- $r_{L_{\max}} = 1$：最后一层使用最紧的割

- **全局方向集 $\mathcal{D}_\ell$**：随着层数增加，方向集变得更加密集

- **层内迭代策略**：
  1. 添加全局方向集 $\mathcal{D}_\ell$ 的所有割约束（使用 $r_\ell$）
  2. 进行多轮自适应方向搜索（固定 $r_\ell$）
  3. 判断是否停滞
  4. 若停滞且 $\ell < L_{\max}$，则进入第 $\ell+1$ 层

---

**算法核心逻辑**：

```
初始化：层 ℓ = 1, r = r_1, 添加方向集 D_1
  ↓
【层 1 内部迭代】（固定 r = r_1）
  迭代 1: 求解 → 识别方向 → 添加新割（r = r_1）
  迭代 2: 求解 → 识别方向 → 添加新割（r = r_1）
  ...
  迭代 m: 停滞判断 → 是 → 进入层 2
  ↓
【层 2 内部迭代】（固定 r = r_2）
  添加全局方向集 D_2（r = r_2）
  迭代 m+1: 求解 → 识别方向 → 添加新割（r = r_2）
  迭代 m+2: 求解 → 识别方向 → 添加新割（r = r_2）
  ...
  迭代 n: 停滞判断 → 是 → 进入层 3 或收敛
```

**关键要点**：
- 在同一层内，$r$ **保持固定**，主要通过改变方向 $\mathbf{d}$ 来逼近锥面
- 只有当方向搜索停滞时，才增大 $r$，进入下一层
- 这种策略充分利用了当前约束强度，避免过早提升 $r$ 导致问题过度约束

### 7.4 锥还原算法流程

**算法：基于层次化方向割的锥还原方法**

---

**输入**：
- SOCP 松弛问题的最优解 $(P^*_{t,i,j}, Q^*_{t,i,j}, l^*_{t,i,j}, v^*_{t,i})$
- 收敛容差 $\epsilon > 0$（如 $10^{-6}$）
- 停滞容差 $\epsilon_{\text{stag}} > 0$（如 $10^{-7}$）
- 停滞判断窗口 $W_{\text{stag}}$（如 3-5）
- 最大层数 $L_{\max}$（如 3-5）
- 最大迭代次数 $K_{\max}$（如 50）
- 初始收缩因子 $r_{\min}$（如 0.5）

**输出**：
- 还原后的解 $(\bar{P}_{t,i,j}, \bar{Q}_{t,i,j}, \bar{l}_{t,i,j}, \bar{v}_{t,i})$
- 收敛状态（收敛/达到最大迭代次数/达到最大层数）

**步骤**：

---

**第一阶段：初始化**

1. **全局初始化**：
   - 迭代计数器 $k \leftarrow 0$
   - 当前层数 $\ell \leftarrow 1$
   - 方向割集合 $\mathcal{L} \leftarrow \emptyset$
   - 历史记录 $\mathcal{H} \leftarrow \emptyset$（用于存储每次迭代的间隙值）

2. **第一层初始化**：
   - 生成第一层全局方向集 $\mathcal{D}_1 = \{(1,1)^T, (-1,1)^T, (-1,-1)^T, (1,-1)^T\}$
   - 计算第一层收缩因子 $r_1 = r_{\min}$
   - 对每条支路 $(i,j) \in E$ 和每个时段 $t$，对每个方向 $\mathbf{d} \in \mathcal{D}_1$：
     - 添加约束 $\frac{\mathbf{d}^T \mathbf{x}_{ij}}{\|\mathbf{d}\|_2} \geq r_1 \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}}$ 到 $\mathcal{L}$

---

**第二阶段：主迭代循环**

3. **while** $k < K_{\max}$ **do**：

   **步骤 3.1：求解带割约束的 SOCP 问题**

   求解：
   $$
   \begin{align}
   \min \quad & \sum_{t=0}^{T-1} \sum_{(i,j) \in E} l_{t,i,j} \\
   \text{s.t.} \quad & \text{原 SOCP 约束 (1)-(9)} \\
   & \text{所有 } \mathcal{L} \text{ 中的方向割约束}
   \end{align}
   $$

   得到新解：$(P^{k+1}_{t,i,j}, Q^{k+1}_{t,i,j}, l^{k+1}_{t,i,j}, v^{k+1}_{t,i})$

   ---

   **步骤 3.2：计算锥松弛间隙**

   对所有支路和时段计算：
   $$
   \Delta^{k+1}_{t,i,j} = l^{k+1}_{t,i,j} \cdot v^{k+1}_{t,i} - \left[(P^{k+1}_{t,i,j})^2 + (Q^{k+1}_{t,i,j})^2\right]
   $$

   计算最大间隙：
   $$
   \Delta^{k+1}_{\max} = \max_{t,i,j} \Delta^{k+1}_{t,i,j}
   $$

   将 $\Delta^{k+1}_{\max}$ 添加到历史记录 $\mathcal{H}$

   ---

   **步骤 3.3：检查全局收敛**

   **if** $\Delta^{k+1}_{\max} \leq \epsilon$ **then**

   输出："算法收敛！"

   **return** $(P^{k+1}_{t,i,j}, Q^{k+1}_{t,i,j}, l^{k+1}_{t,i,j}, v^{k+1}_{t,i})$ 及收敛状态

   **end if**

   ---

   **步骤 3.4：层内自适应方向搜索**

   获取当前层的收缩因子：
   $$
   r_\ell = r_{\min} + \frac{\ell-1}{L_{\max}-1} \cdot (1 - r_{\min})
   $$

   对每条支路 $(i,j) \in E$ 和每个时段 $t$：

   **if** $\Delta^{k+1}_{t,i,j} > \epsilon$ **then**

   - 识别当前解的方向：
     $$
     \hat{\mathbf{d}}^{k+1} = \frac{(P^{k+1}_{t,i,j}, Q^{k+1}_{t,i,j})^T}{\|(P^{k+1}_{t,i,j}, Q^{k+1}_{t,i,j})\|_2}
     $$

   - 在 $\hat{\mathbf{d}}^{k+1}$ 附近生成新方向 $\mathbf{d}^{new}$（例如，微扰 $\hat{\mathbf{d}}^{k+1}$ 或在其邻近角度）

   - 添加新的方向割约束到 $\mathcal{L}$：
     $$
     \frac{(\mathbf{d}^{new})^T \mathbf{x}_{ij}}{\|\mathbf{d}^{new}\|_2} \geq r_\ell \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}}
     $$

   **end if**

   ---

   **步骤 3.5：判断层内停滞**

   **if** $|\mathcal{H}| \geq W_{\text{stag}}$ **then**

   - 提取最近 $W_{\text{stag}}$ 次迭代的间隙改进：
     $$
     \text{improvements} = \{\mathcal{H}[k-i] - \mathcal{H}[k-i+1] \mid i = 1, 2, \ldots, W_{\text{stag}}\}
     $$

   - **if** $\max(\text{improvements}) < \epsilon_{\text{stag}}$ **then**

     输出："层 $\ell$ 内搜索停滞"

     **步骤 3.5.1：层次提升判断**

     **if** $\ell < L_{\max}$ **then**

     - 输出："提升到层 $\ell + 1$"

     - $\ell \leftarrow \ell + 1$

     - 计算新层的收缩因子：
       $$
       r_\ell = r_{\min} + \frac{\ell-1}{L_{\max}-1} \cdot (1 - r_{\min})
       $$

     - 生成第 $\ell$ 层全局方向集 $\mathcal{D}_\ell$（更密集的方向）

     - 对所有支路 $(i,j) \in E$、所有时段 $t$ 和所有方向 $\mathbf{d} \in \mathcal{D}_\ell$：
       - 添加约束 $\frac{\mathbf{d}^T \mathbf{x}_{ij}}{\|\mathbf{d}\|_2} \geq r_\ell \cdot \sqrt{l_{t,i,j} \cdot v_{t,i}}$ 到 $\mathcal{L}$

     - 清空历史记录：$\mathcal{H} \leftarrow \emptyset$（开始新层的停滞判断）

     **else**

     - 输出："已达到最大层数 $L_{\max}$，无法进一步提升"

     - **break**（退出主循环）

     **end if**

     **end if**

   **end if**

   ---

   **步骤 3.6：更新迭代计数**

   $k \leftarrow k + 1$

**end while**

---

**第三阶段：终止处理**

4. **终止条件判断**：

   **if** $k \geq K_{\max}$ **then**

   输出："达到最大迭代次数 $K_{\max}$，算法终止"

   输出："最终最大间隙 $\Delta^{k}_{\max}$"

   **end if**

5. **返回结果**：

   返回最后一次迭代的解 $(P^{k}_{t,i,j}, Q^{k}_{t,i,j}, l^{k}_{t,i,j}, v^{k}_{t,i})$ 及终止状态

---

**算法说明**：

1. **层内迭代 vs 层间提升**：
   - **层内**（步骤 3.4）：固定 $r_\ell$，通过自适应选择方向 $\mathbf{d}$ 来逐步逼近锥面
   - **层间**（步骤 3.5.1）：只有当层内搜索停滞时，才增大 $r$ 并添加新的全局方向集

2. **停滞判断**（步骤 3.5）：
   - 使用滑动窗口监测最近 $W_{\text{stag}}$ 次迭代的改进
   - 如果所有改进都小于阈值 $\epsilon_{\text{stag}}$，则认为停滞

3. **收缩因子策略**：
   - 在同一层内，$r$ 保持不变
   - 层间递增：$r_1 < r_2 < \cdots < r_{L_{\max}} = 1$

4. **方向选择策略**：
   - **全局方向集** $\mathcal{D}_\ell$：每层初始化时添加，覆盖全空间
   - **自适应方向**：每次迭代根据当前解自适应选择，针对性强

---

---

### 7.5 算法执行示例

为了更好地理解算法的执行逻辑，下面给出一个典型的执行流程示例。

**场景**：配电网 SOCP 松弛解存在较小的锥松弛间隙，需要锥还原。

**参数设置**：
- 初始收缩因子：$r_{\min} = 0.5$
- 最大层数：$L_{\max} = 3$
- 收敛容差：$\epsilon = 10^{-6}$
- 停滞容差：$\epsilon_{\text{stag}} = 10^{-7}$
- 停滞窗口：$W_{\text{stag}} = 3$

---

**执行过程**：

**初始化阶段**：
- 层数：$\ell = 1$
- 收缩因子：$r_1 = 0.5$
- 添加全局方向集 $\mathcal{D}_1 = \{(1,1)^T, (-1,1)^T, (-1,-1)^T, (1,-1)^T\}$ 的割约束

---

**层 1（$\ell = 1, r_1 = 0.5$）：粗略约束阶段**

| 迭代 | 操作 | 最大间隙 $\Delta_{\max}$ | 改进量 | 说明 |
|------|------|-------------------------|--------|------|
| $k=1$ | 求解 + 识别方向 + 添加新割（$r=0.5$） | $5.0 \times 10^{-4}$ | - | 初始求解 |
| $k=2$ | 求解 + 识别方向 + 添加新割（$r=0.5$） | $3.2 \times 10^{-4}$ | $1.8 \times 10^{-4}$ | 显著改进 |
| $k=3$ | 求解 + 识别方向 + 添加新割（$r=0.5$） | $2.1 \times 10^{-4}$ | $1.1 \times 10^{-4}$ | 显著改进 |
| $k=4$ | 求解 + 识别方向 + 添加新割（$r=0.5$） | $1.5 \times 10^{-4}$ | $0.6 \times 10^{-4}$ | 改进减缓 |
| $k=5$ | 求解 + 识别方向 + 添加新割（$r=0.5$） | $1.45 \times 10^{-4}$ | $0.05 \times 10^{-4}$ | 改进很小 |
| $k=6$ | 求解 + 识别方向 + 添加新割（$r=0.5$） | $1.44 \times 10^{-4}$ | $0.01 \times 10^{-4}$ | 改进很小 |
| $k=7$ | 求解 + 识别方向 + 添加新割（$r=0.5$） | $1.43 \times 10^{-4}$ | $0.01 \times 10^{-4}$ | 改进很小 |

**停滞判断**：
- 最近3次迭代的改进：$\{0.05 \times 10^{-4}, 0.01 \times 10^{-4}, 0.01 \times 10^{-4}\}$
- 最大改进 $0.05 \times 10^{-4} < 10^{-7}$
- **判定：层 1 停滞**

**层次提升**：
- $\ell \leftarrow 2$
- 计算新收缩因子：$r_2 = 0.5 + \frac{1}{2} \cdot (1 - 0.5) = 0.75$
- 添加全局方向集 $\mathcal{D}_2$ 的割约束（更密集的8个方向）

---

**层 2（$\ell = 2, r_2 = 0.75$）：中等约束阶段**

| 迭代 | 操作 | 最大间隙 $\Delta_{\max}$ | 改进量 | 说明 |
|------|------|-------------------------|--------|------|
| $k=8$ | 添加全局方向集 + 求解 | $8.5 \times 10^{-5}$ | $5.8 \times 10^{-5}$ | 全局方向集带来显著改进 |
| $k=9$ | 求解 + 识别方向 + 添加新割（$r=0.75$） | $5.2 \times 10^{-5}$ | $3.3 \times 10^{-5}$ | 显著改进 |
| $k=10$ | 求解 + 识别方向 + 添加新割（$r=0.75$） | $3.1 \times 10^{-5}$ | $2.1 \times 10^{-5}$ | 显著改进 |
| $k=11$ | 求解 + 识别方向 + 添加新割（$r=0.75$） | $1.8 \times 10^{-5}$ | $1.3 \times 10^{-5}$ | 改进减缓 |
| $k=12$ | 求解 + 识别方向 + 添加新割（$r=0.75$） | $1.2 \times 10^{-5}$ | $0.6 \times 10^{-5}$ | 改进减缓 |
| $k=13$ | 求解 + 识别方向 + 添加新割（$r=0.75$） | $1.15 \times 10^{-5}$ | $0.05 \times 10^{-5}$ | 改进很小 |
| $k=14$ | 求解 + 识别方向 + 添加新割（$r=0.75$） | $1.14 \times 10^{-5}$ | $0.01 \times 10^{-5}$ | 改进很小 |
| $k=15$ | 求解 + 识别方向 + 添加新割（$r=0.75$） | $1.13 \times 10^{-5}$ | $0.01 \times 10^{-5}$ | 改进很小 |

**停滞判断**：
- 最近3次迭代的改进：$\{0.05 \times 10^{-5}, 0.01 \times 10^{-5}, 0.01 \times 10^{-5}\}$
- 最大改进 $0.05 \times 10^{-5} < 10^{-7}$
- **判定：层 2 停滞**

**层次提升**：
- $\ell \leftarrow 3$
- 计算新收缩因子：$r_3 = 0.5 + \frac{2}{2} \cdot (1 - 0.5) = 1.0$
- 添加全局方向集 $\mathcal{D}_3$ 的割约束（最密集的16个方向）

---

**层 3（$\ell = 3, r_3 = 1.0$）：精确约束阶段**

| 迭代 | 操作 | 最大间隙 $\Delta_{\max}$ | 改进量 | 说明 |
|------|------|-------------------------|--------|------|
| $k=16$ | 添加全局方向集 + 求解 | $2.8 \times 10^{-6}$ | $8.5 \times 10^{-6}$ | 全局方向集带来显著改进 |
| $k=17$ | 求解 + 识别方向 + 添加新割（$r=1.0$） | $1.3 \times 10^{-6}$ | $1.5 \times 10^{-6}$ | 显著改进 |
| $k=18$ | 求解 + 识别方向 + 添加新割（$r=1.0$） | $6.2 \times 10^{-7}$ | $6.8 \times 10^{-7}$ | 逼近收敛阈值 |

**收敛判断**：
- 最大间隙 $\Delta_{\max} = 6.2 \times 10^{-7} < 10^{-6}$
- **判定：算法收敛！**

---

**算法终止**：
- 总迭代次数：$k = 18$
- 最终间隙：$\Delta_{\max} = 6.2 \times 10^{-7}$
- 使用层数：3 层
- 返回解：$(P^{18}_{t,i,j}, Q^{18}_{t,i,j}, l^{18}_{t,i,j}, v^{18}_{t,i,j})$

---

**关键观察**：

1. **层内充分搜索**：
   - 层 1：7 次迭代（$k=1$ 到 $k=7$）
   - 层 2：8 次迭代（$k=8$ 到 $k=15$）
   - 层 3：3 次迭代（$k=16$ 到 $k=18$）
   - 在每层内进行了多次自适应方向搜索

2. **渐进逼近**：
   - 层 1 从 $5.0 \times 10^{-4}$ 降至 $1.43 \times 10^{-4}$（降低 71%）
   - 层 2 从 $1.43 \times 10^{-4}$ 降至 $1.13 \times 10^{-5}$（降低 92%）
   - 层 3 从 $1.13 \times 10^{-5}$ 降至 $6.2 \times 10^{-7}$（降低 95%）

3. **层次提升时机**：
   - 只有当层内搜索停滞时才增大 $r$
   - 避免了过早提升导致的计算浪费

4. **全局方向集的作用**：
   - 每次层次提升时添加全局方向集，带来显著改进
   - 层 2 提升：改进 $5.8 \times 10^{-5}$
   - 层 3 提升：改进 $8.5 \times 10^{-6}$

---

### 7.6 配电网潮流优化中的应用

在实际配电网潮流优化问题中：

**情况 1：松弛自然紧致**（最常见）

对于辐射状配电网络，在满足以下条件时，二阶锥松弛通常自然紧致：
- 网络拓扑为树状结构（无环）
- 负荷节点无功率注入（纯负荷）
- 目标函数为最小化网损

此时 SOCP 求解结果无需还原，直接满足锥面约束。

**情况 2：需要锥还原**

在以下情况下可能需要锥还原：
- 存在环状网络结构
- 存在分布式发电机（DG）
- 目标函数非典型（如最小化电压偏差）
- 存在储能装置或可调负荷

**还原策略选择**：

| 方法 | 适用场景 | 优点 | 缺点 |
|------|---------|------|------|
| **方向割方法** | 间隙较小，需精确还原 | 保持凸性，全局最优 | 迭代次数可能较多 |
| **差凸规划** | 间隙较大，快速收敛 | 收敛快，实现简单 | 可能收敛到局部解 |
| **直接投影** | 间隙很小，快速修正 | 最快，无需迭代 | 可能违反其他约束 |

### 7.7 数值实现建议

本节提供基于 Python 和 Gurobi 的实现建议。

#### 7.7.1 参数设置

**容差设置**：
```python
# 锥松弛间隙检测容差
epsilon = 1e-6

# 停滞判断容差
epsilon_stag = 1e-7

# 迭代收敛容差（Gurobi求解器）
epsilon_conv = 1e-8
```

**算法参数推荐**：
```python
# 方向割方法参数
r_min = 0.5          # 初始收缩因子
L_max = 3            # 最大层数
K_max = 50           # 最大迭代次数
W_stag = 3           # 停滞判断窗口大小

# 第一层方向（二维情况，4个方向）
directions_layer1 = [
    (1, 1),   # 45°
    (-1, 1),  # 135°
    (-1, -1), # 225°
    (1, -1)   # 315°
]

# 第二层方向（更密集，8-12个方向）
directions_layer2 = [
    (1, 0), (0, 1), (-1, 0), (0, -1),      # 坐标轴
    (2, 1), (1, 2), (-2, 1), (-1, 2),      # 第1、2象限细分
    (-2, -1), (-1, -2), (2, -1), (1, -2)   # 第3、4象限细分
]
```

#### 7.7.2 核心函数实现

**间隙计算函数**：
```python
import numpy as np

def compute_cone_gap(P, Q, l, v):
    """
    计算锥松弛间隙

    参数:
        P, Q: 支路功率 (T, num_branches)
        l: 支路电流平方 (T, num_branches)
        v: 节点电压平方（扩展到支路）(T, num_branches)

    返回:
        gaps: 每条支路每个时段的间隙
        max_gap: 最大间隙
    """
    gaps = l * v - (P**2 + Q**2)
    max_gap = np.max(gaps)
    return gaps, max_gap


def identify_direction(P, Q):
    """识别当前解的方向"""
    norm = np.sqrt(P**2 + Q**2)
    if norm < 1e-10:
        return np.array([1.0, 0.0])
    return np.array([P / norm, Q / norm])


def check_stagnation(history, window_size, threshold):
    """检查层内搜索是否停滞"""
    if len(history) < window_size:
        return False

    recent = history[-window_size:]
    improvements = [recent[i] - recent[i+1]
                   for i in range(len(recent)-1)]

    return all(imp < threshold for imp in improvements)
```

**方向割约束添加函数**：
```python
import gurobipy as gp

def add_directional_cut(model, P_vars, Q_vars, l_vars, v_vars,
                       direction, r, branch_idx, time_idx):
    """
    添加方向割约束到Gurobi模型

    约束形式: d^T·x / ||d|| >= r * sqrt(l*v)

    使用旋转二阶锥约束实现
    """
    d_P, d_Q = direction
    norm_d = np.sqrt(d_P**2 + d_Q**2)

    # 方向投影
    projection = (d_P * P_vars[time_idx, branch_idx] +
                 d_Q * Q_vars[time_idx, branch_idx]) / norm_d

    # 创建辅助变量 t >= sqrt(l*v)
    t = model.addVar(lb=0, name=f"aux_t_{time_idx}_{branch_idx}")

    # 旋转二阶锥约束：l*v <= t^2
    model.addQConstr(
        l_vars[time_idx, branch_idx] * v_vars[time_idx, branch_idx] <= t * t,
        name=f"soc_aux_{time_idx}_{branch_idx}"
    )

    # 方向割约束：projection >= r * t
    constr = model.addConstr(
        projection >= r * t,
        name=f"dir_cut_t{time_idx}_b{branch_idx}"
    )

    return constr
```

**主算法实现**：
```python
from collections import deque

def hierarchical_cone_recovery(model, solution_init, params):
    """
    层次化锥还原主算法

    参数:
        model: Gurobi SOCP 模型
        solution_init: 初始解字典 {'P': ..., 'Q': ..., 'l': ..., 'v': ...}
        params: 算法参数

    返回:
        solution: 还原后的解
        info: 算法信息
    """
    # 提取参数
    epsilon = params['epsilon']
    epsilon_stag = params['epsilon_stag']
    W_stag = params['W_stag']
    L_max = params['L_max']
    K_max = params['K_max']
    r_min = params['r_min']

    # 初始化
    layer = 1
    k = 0
    history = deque(maxlen=W_stag + 1)

    # 计算初始间隙
    gaps, max_gap = compute_cone_gap(
        solution_init['P'], solution_init['Q'],
        solution_init['l'], solution_init['v']
    )
    history.append(max_gap)
    print(f"初始间隙: {max_gap:.2e}")

    # 第一层：添加全局方向集
    r_current = r_min
    print(f"\n层 {layer}: r = {r_current:.3f}")

    T, num_branches = solution_init['P'].shape
    for direction in directions_layer1:
        for t in range(T):
            for b in range(num_branches):
                if gaps[t, b] > epsilon:
                    add_directional_cut(
                        model, model._P, model._Q,
                        model._l, model._v,
                        direction, r_current, b, t
                    )

    # 主迭代循环
    while k < K_max:
        k += 1
        print(f"\n迭代 {k}")

        # 求解SOCP
        model.optimize()
        if model.status != gp.GRB.OPTIMAL:
            print(f"求解失败: {model.status}")
            break

        # 提取解
        solution = extract_solution(model)

        # 计算间隙
        gaps, max_gap = compute_cone_gap(
            solution['P'], solution['Q'],
            solution['l'], solution['v']
        )
        history.append(max_gap)
        print(f"间隙: {max_gap:.2e}")

        # 检查收敛
        if max_gap <= epsilon:
            print(f"\n收敛！")
            return solution, {
                'converged': True,
                'iterations': k,
                'layers': layer,
                'final_gap': max_gap
            }

        # 层内自适应方向搜索
        r_current = r_min + (layer-1)/(L_max-1) * (1-r_min)

        for t in range(T):
            for b in range(num_branches):
                if gaps[t, b] > epsilon:
                    d_hat = identify_direction(
                        solution['P'][t,b],
                        solution['Q'][t,b]
                    )
                    add_directional_cut(
                        model, model._P, model._Q,
                        model._l, model._v,
                        d_hat, r_current, b, t
                    )

        # 检查停滞
        if check_stagnation(list(history), W_stag, epsilon_stag):
            print(f"层 {layer} 停滞")

            if layer < L_max:
                layer += 1
                r_current = r_min + (layer-1)/(L_max-1) * (1-r_min)
                print(f"\n提升到层 {layer}: r = {r_current:.3f}")

                # 添加新层全局方向集
                global_dirs = (directions_layer2 if layer == 2
                              else directions_layer1)
                for direction in global_dirs:
                    for t in range(T):
                        for b in range(num_branches):
                            if gaps[t, b] > epsilon:
                                add_directional_cut(
                                    model, model._P, model._Q,
                                    model._l, model._v,
                                    direction, r_current, b, t
                                )

                history.clear()
                history.append(max_gap)
            else:
                print(f"达到最大层数")
                break

    return solution, {
        'converged': False,
        'iterations': k,
        'layers': layer,
        'final_gap': max_gap
    }
```

#### 7.7.3 完整求解流程

**算法选择逻辑**：
```python
def solve_socp_with_cone_recovery(model, params):
    """
    求解SOCP并在必要时进行锥还原
    """
    # 1. 求解初始SOCP松弛
    print("求解初始SOCP...")
    model.optimize()

    if model.status != gp.GRB.OPTIMAL:
        raise RuntimeError(f"SOCP求解失败: {model.status}")

    # 2. 提取解
    solution = extract_solution(model)

    # 3. 检测间隙
    _, max_gap = compute_cone_gap(
        solution['P'], solution['Q'],
        solution['l'], solution['v']
    )

    print(f"初始间隙: {max_gap:.2e}")

    # 4. 根据间隙选择策略
    if max_gap < 1e-6:
        print("松弛紧致，无需还原")
        return solution, {'method': 'none', 'gap': max_gap}

    elif max_gap < 1e-4:
        print("\n使用方向割方法...")
        solution, info = hierarchical_cone_recovery(
            model, solution, params
        )
        info['method'] = 'directional_cut'
        return solution, info

    else:
        print("间隙较大，建议使用差凸规划")
        return solution, {'method': 'warning', 'gap': max_gap}


# 使用示例
params = {
    'epsilon': 1e-6,
    'epsilon_stag': 1e-7,
    'W_stag': 3,
    'L_max': 3,
    'K_max': 50,
    'r_min': 0.5
}

solution, info = solve_socp_with_cone_recovery(model, params)

print("\n=== 结果 ===")
print(f"方法: {info['method']}")
print(f"收敛: {info.get('converged', 'N/A')}")
print(f"迭代: {info.get('iterations', 'N/A')}")
print(f"层数: {info.get('layers', 'N/A')}")
print(f"最终间隙: {info.get('final_gap', info.get('gap')):.2e}")
```

### 7.8 与相角恢复的关系

锥还原方法确保了变量满足锥面约束，这是相角恢复的**必要条件**：

**完整求解流程**：

```
1. 求解 SOCP 松弛问题
   ↓
2. 检测锥松弛间隙
   ↓
3. [若需要] 锥还原（本章方法）
   ↓
4. 相角恢复（第 5 节方法）
   ↓
5. 计算完整潮流结果
```

**锥面约束的物理意义**：

在潮流方程中，锥面约束 $l_{t,i,j} \cdot v_{t,i} = P_{t,i,j}^2 + Q_{t,i,j}^2$ 等价于：

$$
|I_{ij}|^2 \cdot |V_i|^2 = |S_{ij}|^2
$$

这正是支路功率与电压、电流的精确关系。松弛到锥内部意味着违反了这一物理约束，因此需要还原。

---

### 7.9 算法伪代码

为了便于理解和实现，下面给出第 7.4 节算法的简洁伪代码形式。

---

**算法 1**：层次化方向割锥还原算法

---

**输入**：
- SOCP 松弛解：$(P^*, Q^*, l^*, v^*)$
- 收敛容差：$\epsilon = 10^{-6}$
- 停滞容差：$\epsilon_{\text{stag}} = 10^{-7}$
- 停滞窗口：$W_{\text{stag}} = 3$
- 最大层数：$L_{\max} = 3$
- 最大迭代次数：$K_{\max} = 50$
- 初始收缩因子：$r_{\min} = 0.5$

**输出**：
- 还原解：$(\bar{P}, \bar{Q}, \bar{l}, \bar{v})$
- 收敛状态

---

**算法主体**：

```
1:  procedure HIERARCHICAL_CONE_RECOVERY(P*, Q*, l*, v*, params)
2:
3:      ▷ 阶段一：初始化
4:      k ← 0, ℓ ← 1
5:      L ← ∅                                    ▷ 方向割约束集合
6:      H ← ∅                                    ▷ 间隙历史记录
7:      r₁ ← r_min
8:      D₁ ← {(1,1)ᵀ, (-1,1)ᵀ, (-1,-1)ᵀ, (1,-1)ᵀ}   ▷ 第一层全局方向集
9:
10:     ▷ 添加初始方向割约束
11:     for each (i,j) ∈ E, t ∈ [0,T-1], d ∈ D₁ do
12:         L ← L ∪ {dᵀ·x_ij / ||d|| ≥ r₁·√(l_{t,i,j}·v_{t,i})}
13:     end for
14:
15:     ▷ 阶段二：主迭代循环
16:     while k < K_max do
17:         k ← k + 1
18:
19:         ▷ 步骤 1：求解带约束的 SOCP 问题
20:         (Pᵏ, Qᵏ, lᵏ, vᵏ) ← SOLVE_SOCP(L)
21:
22:         ▷ 步骤 2：计算锥松弛间隙
23:         for each (i,j) ∈ E, t ∈ [0,T-1] do
24:             Δᵏ_{t,i,j} ← lᵏ_{t,i,j}·vᵏ_{t,i} - [(Pᵏ_{t,i,j})² + (Qᵏ_{t,i,j})²]
25:         end for
26:         Δᵏ_max ← max_{t,i,j} Δᵏ_{t,i,j}
27:         H ← H ∪ {Δᵏ_max}
28:
29:         ▷ 步骤 3：检查全局收敛
30:         if Δᵏ_max ≤ ε then
31:             return (Pᵏ, Qᵏ, lᵏ, vᵏ), CONVERGED
32:         end if
33:
34:         ▷ 步骤 4：层内自适应方向搜索（固定 rℓ）
35:         r_ℓ ← r_min + (ℓ-1)/(L_max-1)·(1 - r_min)
36:         for each (i,j) ∈ E, t ∈ [0,T-1] do
37:             if Δᵏ_{t,i,j} > ε then
38:                 d̂ ← (Pᵏ_{t,i,j}, Qᵏ_{t,i,j})ᵀ / ||(Pᵏ_{t,i,j}, Qᵏ_{t,i,j})||
39:                 L ← L ∪ {d̂ᵀ·x_ij / ||d̂|| ≥ r_ℓ·√(l_{t,i,j}·v_{t,i})}
40:             end if
41:         end for
42:
43:         ▷ 步骤 5：判断层内停滞
44:         if |H| ≥ W_stag then
45:             recent ← LAST_N_ELEMENTS(H, W_stag)
46:             improvements ← {recent[i] - recent[i+1] | i = 1,...,W_stag-1}
47:
48:             if max(improvements) < ε_stag then
49:                 print "层 ℓ 搜索停滞"
50:
51:                 ▷ 步骤 6：层次提升
52:                 if ℓ < L_max then
53:                     ℓ ← ℓ + 1
54:                     r_ℓ ← r_min + (ℓ-1)/(L_max-1)·(1 - r_min)
55:                     D_ℓ ← GENERATE_GLOBAL_DIRECTIONS(ℓ)
56:                     print "提升到层 ℓ，r = r_ℓ"
57:
58:                     ▷ 添加新层的全局方向集
59:                     for each (i,j) ∈ E, t ∈ [0,T-1], d ∈ D_ℓ do
60:                         L ← L ∪ {dᵀ·x_ij / ||d|| ≥ r_ℓ·√(l_{t,i,j}·v_{t,i})}
61:                     end for
62:
63:                     H ← ∅                    ▷ 清空历史，开始新层监测
64:                     H ← H ∪ {Δᵏ_max}
65:                 else
66:                     print "达到最大层数 L_max"
67:                     break
68:                 end if
69:             end if
70:         end if
71:     end while
72:
73:     ▷ 阶段三：终止处理
74:     if k ≥ K_max then
75:         print "达到最大迭代次数"
76:     end if
77:     return (Pᵏ, Qᵏ, lᵏ, vᵏ), NOT_FULLY_CONVERGED
78:
79: end procedure
```

---

**算法 2**：生成全局方向集

---

```
1:  function GENERATE_GLOBAL_DIRECTIONS(ℓ)
2:
3:      if ℓ = 1 then
4:          ▷ 第一层：4 个象限等分线方向
5:          D_ℓ ← {(1,1)ᵀ, (-1,1)ᵀ, (-1,-1)ᵀ, (1,-1)ᵀ}
6:
7:      else if ℓ = 2 then
8:          ▷ 第二层：12 个更密集的方向
9:          D_ℓ ← {(1,0)ᵀ, (0,1)ᵀ, (-1,0)ᵀ, (0,-1)ᵀ,
10:                (2,1)ᵀ, (1,2)ᵀ, (-2,1)ᵀ, (-1,2)ᵀ,
11:                (-2,-1)ᵀ, (-1,-2)ᵀ, (2,-1)ᵀ, (1,-2)ᵀ}
12:
13:     else if ℓ ≥ 3 then
14:         ▷ 第三层及以上：16 个或更密集的方向
15:         D_ℓ ← ANGULAR_SAMPLING(16)     ▷ 在 [0°, 360°) 均匀采样
16:     end if
17:
18:     return D_ℓ
19:
20: end function
```

---

### 算法复杂度分析

**时间复杂度（单次迭代）**：

| 操作 | 复杂度 | 说明 |
|------|--------|------|
| SOCP 求解（第 20 行） | $O(n^3)$ | $n$ 为决策变量数 |
| 间隙计算（第 23-27 行） | $O(T \cdot \|E\|)$ | 遍历所有支路和时段 |
| 方向割添加（第 36-41 行） | $O(T \cdot \|E\|)$ | 最多添加 $T \cdot \|E\|$ 个约束 |
| **总计** | **$O(n^3)$** | 由 SOCP 求解主导 |

**空间复杂度**：

| 数据结构 | 复杂度 | 说明 |
|---------|--------|------|
| 决策变量 | $O(T \cdot \|E\|)$ | 支路变量和节点变量 |
| 方向割约束集合 $\mathcal{L}$ | $O(K_{\max} \cdot T \cdot \|E\|)$ | 最坏情况（每次迭代添加 $T \cdot \|E\|$ 个） |
| **总计** | **$O(K_{\max} \cdot T \cdot \|E\|)$** | - |

**迭代次数**：
- 典型情况：10-30 次迭代
- 最坏情况：$K_{\max} = 50$ 次迭代
- 层数：通常 2-3 层

---

### 算法关键步骤说明

**1. 初始化（第 3-13 行）**
   - 设置初始层数 $\ell = 1$，收缩因子 $r_1 = r_{\min}$
   - 生成第一层全局方向集 $\mathcal{D}_1$（4 个象限等分线）
   - 为所有支路和时段添加初始方向割约束

**2. 层内自适应方向搜索（第 34-41 行）**
   - **核心特点**：在固定的 $r_\ell$ 下进行迭代
   - 根据当前解识别方向 $\hat{\mathbf{d}}$
   - 在识别方向上添加新的方向割约束
   - 重复此过程直至收敛或停滞

**3. 停滞判断（第 43-50 行）**
   - 使用滑动窗口监测最近 $W_{\text{stag}}$ 次迭代
   - 计算每次迭代的间隙改进量
   - 若最大改进小于阈值 $\epsilon_{\text{stag}}$，判定为停滞

**4. 层次提升（第 51-68 行）**
   - **触发条件**：当前层搜索停滞
   - **操作**：
     - 增大收缩因子：$r_{\ell+1} > r_\ell$
     - 生成更密集的全局方向集 $\mathcal{D}_{\ell+1}$
     - 添加新层的所有方向割约束
     - 清空历史记录，开始新层的停滞监测

**5. 清空历史（第 63-64 行）**
   - 防止上一层的历史数据影响新层的停滞判断
   - 确保每层独立监测收敛行为

---

### 算法与详细流程的对应关系

| 伪代码行号 | 详细算法（7.4 节） | 说明 |
|-----------|------------------|------|
| 3-13 | 步骤 1-2 | 初始化和第一层方向割 |
| 16-32 | 步骤 3.1-3.3 | 主循环：求解、计算间隙、检查收敛 |
| 34-41 | 步骤 3.4 | **层内自适应方向搜索**（固定 $r_\ell$） |
| 43-70 | 步骤 3.5-3.5.1 | 停滞判断和层次提升 |
| 73-77 | 步骤 4-5 | 终止处理和返回结果 |

---

## 8. 模型特点（修订）

1. **凸优化问题**：通过二阶锥松弛，将非凸的交流潮流问题转化为凸优化问题，可以高效求解并获得全局最优解

2. **精确性保证**：
   - 在辐射状配电网络中，二阶锥松弛通常自然紧致
   - 当松弛不紧时，通过锥还原方法确保解满足物理约束
   - 结合相角恢复技术，可获得原问题的高质量解

3. **多时段优化**：考虑多个时间断面的联合优化，支持时变的负荷、发电机出力限制等

4. **适用场景**：适用于辐射状配电网络的最优潮流计算，可用于配电网规划、运行优化等场景

5. **模型依据**：基于 Branch Flow Model (BFM) 和二阶锥松弛 (SOCP relaxation) 技术，并提供完整的锥还原求解方法

6. **鲁棒性**：
   - 方向割方法保持问题凸性，确保全局最优性
   - 差凸规划方法提供快速收敛的备选方案
   - 层次化策略平衡了求解效率和精度

---

**注**：该模型实现基于 Gurobi 优化器求解，使用了高精度参数设置（`MIPGap=1e-8`, `OptimalityTol=1e-8`, `FeasibilityTol=1e-8`）以确保求解精度。锥还原方法可在必要时自动激活，确保解的物理可行性。
