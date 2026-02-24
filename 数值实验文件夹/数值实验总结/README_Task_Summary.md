# 任务完成总结

## 已完成的三个任务

### 任务 1：将二阶锥约束改为等式约束，创建 ACOPF_Model ✅

**实现内容：**
- 创建了 `compare_models.py` 脚本，包含两个模型：
  - **SOCP-AC-OPF**：使用二阶锥松弛（不等式 ≤）
  - **AC-OPF**：使用等式约束（=）

**关键差异：**
```python
# SOCP 模型（约束 5-SOCP）
model.addQConstr(lhs <= rhs, name=f"SOCP_t{t}_{i}_{j}")

# AC 模型（约束 5-AC）
model.addQConstr(P^2 + Q^2 == l * v, name=f"AC_t{t}_{i}_{j}")
```

**数学表达：**
- SOCP: $(2P)^2 + (2Q)^2 + (l-v)^2 \leq (l+v)^2$ （凸优化）
- AC: $P^2 + Q^2 = l \cdot v$ （非凸优化）

**求解结果：**
- 两个模型都成功求解到最优
- SOCP 使用 Barrier 算法（10次迭代）
- AC 使用 MIP 算法（1个节点，180次迭代）

---

### 任务 2：优化输出，将每个时间断面的结果保存到 CSV 文件 ✅

**实现功能：**
- 创建 `save_results_to_csv()` 函数
- 为每个时间断面 $t$ 生成独立的 CSV 文件
- 删除了冗余的 print 输出

**生成的文件：**

#### SOCP 模型输出
- `output/SOCP_branch_t{t}.csv` - 支路潮流结果
  - 列：from_bus, to_bus, P_MW, Q_MVar, l_pu, I_magnitude_pu, r_ij, x_ij
- `output/SOCP_bus_t{t}.csv` - 节点结果
  - 列：bus, v_pu_squared, V_magnitude_pu, p_inject_pu, q_inject_pu, P_inject_MW, Q_inject_MVar
- `output/SOCP_gen_t{t}.csv` - 发电机结果
  - 列：gen_id, Pg_pu, Qg_pu, Pg_MW, Qg_MVar

#### AC 模型输出
- `output/AC_branch_t{t}.csv` - 支路潮流结果
- `output/AC_bus_t{t}.csv` - 节点结果
- `output/AC_gen_t{t}.csv` - 发电机结果

**文件格式示例：**
```csv
bus,v_pu_squared,V_magnitude_pu,p_inject_pu,q_inject_pu,P_inject_MW,Q_inject_MVar
1,1.0,1.0,0.3918,0.2435,3.9177,2.4351
2,0.9941,0.9970,-0.01,-0.006,-0.1,-0.06
...
```

---

### 任务 3：对比两个模型的松弛间隙 ✅

**实现功能：**
- 创建 `compare_models()` 函数
- 全面对比目标函数值、决策变量和松弛间隙

**对比结果（IEEE 33节点系统）：**

#### 目标函数对比
| 指标 | SOCP 模型 | AC 模型 | 差异 |
|------|----------|---------|------|
| 目标函数 $\sum l$ | 0.7902633383 | 0.7902625606 | 7.78×10⁻⁷ |
| 相对松弛间隙 | - | - | **0.000098%** |

#### 松弛间隙统计
| 统计量 | 数值 |
|--------|------|
| 最大松弛量 | 3.17×10⁻⁷ |
| 平均松弛量 | 2.08×10⁻⁸ |
| 总松弛量 | 6.66×10⁻⁷ |

**结论：二阶锥松弛是紧的 (tight relaxation)** ✅

所有支路的松弛量 $\text{slack} = l \cdot v - (P^2 + Q^2)$ 都小于 $10^{-6}$，表明 SOCP 模型的解与 AC 模型的精确解几乎相同。

#### 典型支路的松弛间隙
| 支路 | SOCP: $P^2+Q^2$ | SOCP: $l \cdot v$ | AC: $P^2+Q^2$ | AC: $l \cdot v$ | 松弛量 |
|------|----------------|------------------|--------------|---------------|-------|
| 1→2 | 0.2127810856 | 0.2127811256 | 0.2127810564 | 0.2127810564 | 3.99×10⁻⁸ |
| 20→21 | 0.0003888015 | 0.0003891184 | 0.0003888010 | 0.0003888010 | 3.17×10⁻⁷ |

---

## 创建的文件

### 1. 核心脚本
- **`compare_models.py`** - 主要对比脚本
  - 包含 `solve_model()` 函数：求解 SOCP 或 AC 模型
  - 包含 `save_results_to_csv()` 函数：保存结果到 CSV
  - 包含 `compare_models()` 函数：对比两个模型

### 2. 文档
- **`Model_Comparison.md`** - 详细的数学模型对比文档
  - 包含两个模型的完整数学描述
  - 包含决策变量、参数和约束的详细说明
  - 包含数值对比结果和结论

- **`SOC_ACOPF_Model.md`** - SOCP 模型的数学描述（已更新）

### 3. Notebook
- **`Model_Comparison.ipynb`** - 交互式对比 notebook
  - 展示如何使用 `compare_models.py`
  - 包含可视化代码
  - 包含结果分析和统计

### 4. 输出文件（在 `output/` 目录）
- SOCP 和 AC 模型的所有结果 CSV 文件

---

## 使用方法

### 方法 1：直接运行 Python 脚本
```bash
conda activate CETLAB3_11_4
python compare_models.py
```

### 方法 2：在 Jupyter Notebook 中运行
```python
%run compare_models.py
```

### 方法 3：使用交互式 Notebook
打开 `Model_Comparison.ipynb` 并逐个运行单元格

---

## 主要发现

### 1. 松弛紧性验证 ✅
- SOCP 松弛在辐射状配电网络中是**紧的**
- 相对松弛间隙仅为 0.000098%（远小于 1%）
- 最大松弛量 3.17×10⁻⁷（远小于数值精度 10⁻⁶）

### 2. 求解性能对比
| 特性 | SOCP 模型 | AC 模型 |
|------|----------|---------|
| 问题性质 | 凸优化 | 非凸优化 |
| 求解算法 | Barrier | MIP |
| 迭代次数 | 10 | 180 |
| 求解时间 | 0.01s | 0.01s |
| 全局最优 | ✅保证 | ⚠️可能局部最优 |

### 3. 实用性建议
- ✅ **推荐使用 SOCP 模型**：对于辐射状配电网
- 原因：
  - 保证全局最优
  - 求解效率高
  - 解的精度与 AC 模型相同
- ⚠️ **使用 AC 模型**：仅在需要精确等式约束或环网情况下

---

## 技术细节

### 模型特点
1. **多时段优化**：支持 T 个时间断面
2. **Branch Flow Model**：基于 BFM 潮流模型
3. **二阶锥松弛**：使用 SOCP 将非凸问题转为凸问题
4. **高精度求解**：MIPGap=1e-8, OptimalityTol=1e-8, FeasibilityTol=1e-8

### 约束说明
- 根节点电压约束 (1)
- 节点净注入功率约束 (2a-2b)
- 潮流平衡约束 (3a-3b) - Branch Flow Model
- 电压降落约束 (4) - DistFlow Model
- **关键约束 (5)**：
  - SOCP: 二阶锥不等式
  - AC: 等式约束
- 发电机/负荷约束 (6-7)
- 电压上下限约束 (8)

### 求解器设置
```python
model.setParam('MIPGap',         1e-8)
model.setParam('OptimalityTol',  1e-8)
model.setParam('FeasibilityTol', 1e-8)
model.setParam('NonConvex',      2)  # AC 模型需要
```

---

## 参考文献

1. Baran, M. E., & Wu, F. F. (1989). Network reconfiguration in distribution systems for loss reduction and load balancing. *IEEE Transactions on Power Delivery*, 4(2), 1401-1407.

2. Farivar, M., & Low, S. H. (2013). Branch flow model: Relaxations and convexification. *IEEE Transactions on Power Systems*, 28(3), 2554-2564.

3. Low, S. H. (2014). Convex relaxation of optimal power flow—Part I: Formulations and equivalence. *IEEE Transactions on Control of Network Systems*, 1(1), 15-27.

---

## 结论

✅ **三个任务全部完成**

1. ✅ 成功实现 ACOPF_Model（等式约束）
2. ✅ 优化输出，所有结果保存到 CSV 文件
3. ✅ 全面对比两个模型的松弛间隙

**核心发现：在 IEEE 33 节点辐射状配电网络中，二阶锥松弛是紧的，SOCP 模型可以作为 AC-OPF 的精确凸松弛。**

---

**生成时间**：2026-02-06
**求解器**：Gurobi 13.0.0
**测试系统**：IEEE 33-bus Radial Distribution Network
