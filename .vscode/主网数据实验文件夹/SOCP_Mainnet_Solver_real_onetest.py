"""
环网下数据中心承载力评估模型 - 单/少量场景测试版（onetest）
============================================================
基于 power_system_model_onetest.md 的环网承载力模型实现

主要特点：
1. 单场景确定性模型，无场景上标ω
2. 使用二阶锥松弛（SOCP）求解环网潮流
3. 目标：最大化全局承载力下界 hc = min_t Σ_i P^IDC_{i,t}
4. 环网相角参考约束（式7）：消除旋转不变性
5. 场景输入控制：SCENARIO_INDEX（单组）或 N_SCENARIOS（前N组）
6. 输出目录：test_output/output_onetest/
7. 式(5)松弛紧程度测试：计算并报告每条支路每时段的松弛间隙 δ

约束编号对应 power_system_model_onetest.md：
  式(1)：目标函数 - 最大化 hc
  式(2)：承载力下界约束
  式(3)(4)：节点功率平衡
  式(5)：二阶锥松弛（显式约束）
  式(6)：辅助变量对称性（显式约束）
  式(7)：相角参考约束（环网关键）
  式(8)(9)：对角线元素（显式约束）
  式(10-12)：运行约束（变量上下界实现）
  式(13)：参考节点电压（变量固定 lb=ub=1 实现）
  式(17-29)：数据中心任务调度
  注：式(14)不进入优化（SOCP松弛替代）；式(15)(16)后处理计算；
      式(22)(25)(27)通过变量 lb=0 实现
"""

import os
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # 非交互后端，线程安全
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from filelock import FileLock
from case33_mainnet import case33_mainnet
from tqdm import tqdm
from generate_scenarios import generate_scenario_combinations, read_scenario_data

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

current_dir = os.path.dirname(os.path.abspath(__file__))


# ===========================================================
# SECTION 1: 网络参数处理
# ===========================================================

def load_network_data(ppc):
    """
    从 case33_mainnet 提取网络拓扑和参数

    Returns
    -------
    T : int
    N : list
    E : list of tuple  [(from, to, r, x, Gij, Bij, ...), ...]
    G_gen : dict
    DC_data : dict
    loops : list
    baseMVA : float
    """
    baseMVA = ppc["baseMVA"]
    T = int(ppc["T"])

    num_buses = ppc["bus_data"]["num_buses"]
    N = list(range(num_buses))

    Vbase = ppc["bus_data"]["base_kV"] * 1e3
    Sbase = baseMVA * 1e6
    Zbase = Vbase**2 / Sbase

    E = []
    for br in ppc["branch"]:
        i, j = int(br[0]-1), int(br[1]-1)
        r_pu = br[2] / Zbase
        x_pu = br[3] / Zbase
        z_sq = r_pu**2 + x_pu**2
        Gij = r_pu / z_sq
        Bij = -x_pu / z_sq
        E.append((i, j, r_pu, x_pu, Gij, Bij, 0.0, 0.0, 0.0, 0.0))

    G_gen = {}
    for gen_row in ppc["gen"]:
        bus = int(gen_row[0] - 1)
        G_gen[bus] = {
            'Pmin': gen_row[9] / baseMVA,
            'Pmax': gen_row[8] / baseMVA,
            'Qmin': gen_row[4] / baseMVA,
            'Qmax': gen_row[3] / baseMVA,
        }

    DC_data = ppc.get("DC", None)
    if DC_data is not None:
        DC_data['bus_0idx'] = [int(b-1) for b in DC_data['bus']]
        DC_data['K'] = list(range(len(DC_data['bus_0idx'])))
        DC_data['D'] = int(np.max(DC_data['tau_max']))

    loops = ppc.get("loops", [])

    return T, N, E, G_gen, DC_data, loops, baseMVA


# ===========================================================
# SECTION 2: 单场景SOCP模型构建与求解
# ===========================================================

def build_and_solve_single(scenario_data, output_dir, scenario_id):
    """
    构建并求解单场景SOCP环网承载力模型（流水线：数据输入→求解→返回结果）

    Parameters
    ----------
    scenario_data : dict
        单个场景数据 {'name', 'dc0_lambda', 'dc1_lambda', 'loadlist18', 'loadlist25'}
    output_dir : str
        输出根目录
    scenario_id : int
        场景编号（用于文件命名）

    Returns
    -------
    tuple : (scenario_id, hc_kW, scale, status, solve_time) 或失败时返回 (scenario_id, None, None, status, solve_time)
    """
    t_start = time.time()

    ppc = case33_mainnet(scenario_data)
    T, N, E, G_gen, DC_data, loops, baseMVA = load_network_data(ppc)

    K = DC_data['K']
    IDC_nodes = DC_data['bus_0idx']
    arcs = DC_data.get('arcs', [])
    load_ts = ppc["bus_data"]["load_timeseries"]
    Qd_ratio = ppc["bus_data"]["Qd_ratio"]
    lambda_shape = ppc["DC"]['lambda_local']  # (K, T)

    # ---- 创建Gurobi模型 ----
    model = gp.Model(f"SOCP_s{scenario_id}")
    model.setParam("OutputFlag", 0)
    model.setParam("OptimalityTol", 1e-6)
    model.setParam("FeasibilityTol", 1e-6)
    model.setParam("NonConvex", 2)
    model.setParam("TimeLimit", 600)
    model.setParam("Threads", max(1, os.cpu_count() // 2))

    # ================================================================
    # 变量声明
    # 约束处理原则：
    #   式(10-12)：通过 lb/ub 实现，不添加独立约束
    #   式(13)：U_ref 固定为 lb=ub=1.0
    #   式(22)(25)(27)：变量 lb=0 实现
    #   式(14)：不进优化模型（SOCP松弛替代）
    #   式(15)(16)：后处理计算，不加入优化
    # ================================================================

    # 全局变量
    scale = model.addVar(lb=0.0, ub=GRB.INFINITY, name="scale")
    hc    = model.addVar(lb=0.0, ub=GRB.INFINITY, name="hc")

    # 环网潮流变量
    T_var = {}
    Z_var = {}
    U_var = {}
    P_G   = {}
    Q_G   = {}
    P_IDC = {}

    for t in range(T):
        # 支路辅助变量（正向+反向）
        for (i, j, *_) in E:
            for (a, b) in [(i, j), (j, i)]:
                T_var[a, b, t] = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"T_{a}_{b}_t{t}")
                Z_var[a, b, t] = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"Z_{a}_{b}_t{t}")

        for i in N:
            # 式(13)：参考节点电压固定（lb=ub=1.0），其余节点用式(10)上下界
            if i == 0:
                U_var[i, t] = model.addVar(lb=1.0, ub=1.0, name=f"U_{i}_t{t}")
            else:
                U_var[i, t] = model.addVar(lb=0.9**2, ub=1.1**2, name=f"U_{i}_t{t}")  # 式(10)

            # 式(8)(9)：对角线辅助变量
            T_var[i, i, t] = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"T_{i}_{i}_t{t}")
            Z_var[i, i, t] = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"Z_{i}_{i}_t{t}")

            # 式(11)(12)：发电机出力上下界
            if i in G_gen:
                P_G[i, t] = model.addVar(lb=G_gen[i]['Pmin'], ub=G_gen[i]['Pmax'], name=f"PG_{i}_t{t}")
                Q_G[i, t] = model.addVar(lb=G_gen[i]['Qmin'], ub=G_gen[i]['Qmax'], name=f"QG_{i}_t{t}")

            if i in IDC_nodes:
                P_IDC[i, t] = model.addVar(lb=0.0, ub=GRB.INFINITY, name=f"PIDC_{i}_t{t}")

    # 数据中心任务调度变量
    lambda_var = {}
    q_dc = {}
    r_dc = {}  # lb=0 实现式(22)
    p_dc = {}  # lb=0 实现式(25)
    y_mig = {} # lb=0 实现式(25)

    for k in K:
        tau_k = int(DC_data['tau_max'][k])
        for t in range(T):
            lambda_var[k, t] = model.addVar(lb=0.0, ub=GRB.INFINITY, name=f"lam_k{k}_t{t}")  # lb=0 实现式(27)
            for d in range(tau_k + 1):
                q_dc[k, t, d]  = model.addVar(lb=0.0, ub=GRB.INFINITY, name=f"q_k{k}_t{t}_d{d}")
                r_dc[k, t, d]  = model.addVar(lb=0.0, ub=GRB.INFINITY, name=f"r_k{k}_t{t}_d{d}")
                p_dc[k, t, d]  = model.addVar(lb=0.0, ub=GRB.INFINITY, name=f"p_k{k}_t{t}_d{d}")

    if arcs:
        for (k_from, k_to) in arcs:
            tau_from = int(DC_data['tau_max'][k_from])
            for t in range(T):
                for d in range(tau_from + 1):
                    y_mig[k_from, k_to, t, d] = model.addVar(lb=0.0, ub=GRB.INFINITY,
                                                               name=f"ymig_{k_from}_{k_to}_t{t}_d{d}")

    model.update()

    # ================================================================
    # 添加约束
    # ================================================================

    for t in range(T):
        # ---- 式(8)(9)：对角线元素 ----
        for i in N:
            model.addConstr(T_var[i, i, t] == U_var[i, t],  name=f"Tii_{i}_t{t}")
            model.addConstr(Z_var[i, i, t] == 0.0,          name=f"Zii_{i}_t{t}")

        # ---- 式(6)：辅助变量对称性 ----
        for (i, j, *_) in E:
            model.addConstr(T_var[i, j, t] ==  T_var[j, i, t], name=f"Tsym_{i}_{j}_t{t}")
            model.addConstr(Z_var[i, j, t] == -Z_var[j, i, t], name=f"Zasym_{i}_{j}_t{t}")

        # ---- 式(3)(4)：节点功率平衡 ----
        for i in N:
            P_L_i = 0.0 if i == 0 else -load_ts[i-1][t] / 1000.0 / baseMVA
            Q_L_i = 0.0 if i == 0 else -load_ts[i-1][t] * Qd_ratio[i-1] / 1000.0 / baseMVA

            P_G_i   = P_G.get((i, t), 0)
            Q_G_i   = Q_G.get((i, t), 0)
            P_IDC_i = P_IDC.get((i, t), 0) if i in IDC_nodes else 0

            sum_TG_ZB = gp.LinExpr()
            sum_ZG_TB = gp.LinExpr()
            for (ii, jj, _, _, Gij, Bij, *_) in E:
                if ii == i:
                    sum_TG_ZB += T_var[i, jj, t] * Gij + Z_var[i, jj, t] * Bij
                    sum_ZG_TB += Z_var[i, jj, t] * Gij - T_var[i, jj, t] * Bij
                elif jj == i:
                    sum_TG_ZB += T_var[ii, i, t] * Gij + Z_var[ii, i, t] * Bij
                    sum_ZG_TB += Z_var[ii, i, t] * Gij - T_var[ii, i, t] * Bij

            model.addConstr(P_G_i - P_L_i - P_IDC_i - sum_TG_ZB == 0, name=f"PBal_{i}_t{t}")
            model.addConstr(Q_G_i - Q_L_i - sum_ZG_TB == 0,            name=f"QBal_{i}_t{t}")

        # ---- 式(5)：二阶锥松弛约束 ----
        for (i, j, *_) in E:
            model.addConstr(
                4*T_var[i, j, t]*T_var[i, j, t] + 4*Z_var[i, j, t]*Z_var[i, j, t] +
                (U_var[i, t] - U_var[j, t])*(U_var[i, t] - U_var[j, t])
                <= (U_var[i, t] + U_var[j, t])*(U_var[i, t] + U_var[j, t]),
                name=f"SOCP_{i}_{j}_t{t}"
            )

        # ---- 式(2)：承载力下界约束 ----
        model.addConstr(
            hc <= gp.quicksum(P_IDC[i, t] for i in IDC_nodes),
            name=f"hc_lb_t{t}"
        )

    # ---- 式(7)：相角参考约束（环路约束）----
    for t in range(T):
        for loop_idx, loop_branches in enumerate(loops):
            angle_sum = gp.LinExpr()
            for (i, j) in loop_branches:
                # 平坦电压初值：T0=1, Z0=0 → 线性化后等于 Z_var[i,j,t]
                angle_sum += Z_var[i, j, t]
            model.addConstr(angle_sum == 0, name=f"AngleRef_loop{loop_idx}_t{t}")

    # ---- 数据中心任务调度约束（式17-29）----
    mu_max = DC_data['mu_max']
    alpha  = DC_data['alpha']
    beta   = DC_data['beta']
    Bmax   = DC_data.get('Bmax', {})
    s      = DC_data.get('s', 1.0)

    for k in K:
        tau_k = int(DC_data['tau_max'][k])
        bus_k = DC_data['bus_0idx'][k]
        alpha_pu = lambda t_: alpha[k, t_] / 1000.0 / baseMVA
        beta_pu  = lambda t_: beta[k, t_]  / 1000.0 / baseMVA

        for t in range(T):
            # 式(17)：λ_{t,u} = λ^shape_{t,u} · scale
            model.addConstr(lambda_var[k, t] == lambda_shape[k, t] * scale,
                            name=f"LamScale_k{k}_t{t}")

            # 式(18)：d=0 队列到达
            inflow_0 = gp.quicksum(y_mig[kv, k, t, 0]
                                   for (kv, ku) in arcs if ku == k) if arcs else 0
            model.addConstr(q_dc[k, t, 0] == lambda_var[k, t] + inflow_0,
                            name=f"Q0_k{k}_t{t}")

            # 式(19)(20)：d≥1 队列状态更新
            for d in range(1, tau_k + 1):
                inflow_d = gp.quicksum(y_mig[kv, k, t, d]
                                       for (kv, ku) in arcs if ku == k) if arcs else 0
                if t == 0:
                    # 式(20)：初始条件
                    model.addConstr(q_dc[k, 0, d] == 0, name=f"QInit_k{k}_d{d}")
                else:
                    # 式(19)：状态转移
                    model.addConstr(q_dc[k, t, d] == r_dc[k, t-1, d-1] + inflow_d,
                                    name=f"Q_k{k}_t{t}_d{d}")

            # 式(21)：剩余量更新
            for d in range(tau_k + 1):
                outflow_d = gp.quicksum(y_mig[k, kv, t, d]
                                        for (ku, kv) in arcs if ku == k) if arcs else 0
                model.addConstr(
                    r_dc[k, t, d] == q_dc[k, t, d] - p_dc[k, t, d] - outflow_d,
                    name=f"Res_k{k}_t{t}_d{d}"
                )

            # 式(23)：最高级队列清空
            model.addConstr(r_dc[k, t, tau_k] == 0, name=f"Deadline_k{k}_t{t}")

            # 式(28)：处理能力约束（多变量线性不等式，不能用变量界代替）
            model.addConstr(
                gp.quicksum(p_dc[k, t, d] for d in range(tau_k + 1)) <= mu_max[k],
                name=f"Cap_k{k}_t{t}"
            )

            # 式(29)：功耗映射（转换为标幺值）
            model.addConstr(
                P_IDC[bus_k, t] ==
                alpha_pu(t) * gp.quicksum(p_dc[k, t, d] for d in range(tau_k + 1)) + beta_pu(t),
                name=f"PwrMap_k{k}_t{t}"
            )

        # 式(26)：周期性闭环约束
        for d in range(tau_k + 1):
            model.addConstr(q_dc[k, T-1, d] == q_dc[k, 0, d], name=f"Periodic_k{k}_d{d}")

    # 式(24)：带宽约束
    if arcs:
        for (k_from, k_to) in arcs:
            tau_from = int(DC_data['tau_max'][k_from])
            for t in range(T):
                model.addConstr(
                    s * gp.quicksum(y_mig[k_from, k_to, t, d] for d in range(tau_from + 1))
                    <= Bmax.get((k_from, k_to), 0),
                    name=f"BW_{k_from}_{k_to}_t{t}"
                )

    # ================================================================
    # 目标函数（式1）
    # ================================================================
    model.setObjective(-hc, GRB.MINIMIZE)
    model.optimize()

    solve_time = time.time() - t_start
    status = model.Status

    if model.SolCount == 0:
        if status in (GRB.INFEASIBLE, GRB.INF_OR_UNBD):
            os.makedirs(output_dir, exist_ok=True)
            model.computeIIS()
            model.write(os.path.join(output_dir, f"infeasible_s{scenario_id}.ilp"))
        return scenario_id, None, None, status, solve_time

    hc_val    = hc.X
    scale_val = scale.X

    # ---- 收集变量最优值（.X 值，不再保留 Gurobi 变量对象）----
    result_arrays = {}

    U_arr = np.array([[U_var[i, t].X for t in range(T)] for i in N])           # (N, T)
    result_arrays['U'] = U_arr

    PG_buses = sorted({i for (i, _) in P_G})
    PG_arr = np.array([[P_G.get((i, t), None) and P_G[i, t].X or 0.0
                        for t in range(T)] for i in PG_buses])
    result_arrays['PG_buses'] = np.array(PG_buses)
    result_arrays['PG'] = PG_arr

    QG_arr = np.array([[Q_G.get((i, t), None) and Q_G[i, t].X or 0.0
                        for t in range(T)] for i in PG_buses])
    result_arrays['QG'] = QG_arr

    PIDC_arr = np.array([[P_IDC[i, t].X for t in range(T)] for i in IDC_nodes])  # (K, T)
    result_arrays['PIDC'] = PIDC_arr

    lam_arr = np.array([[lambda_var[k, t].X for t in range(T)] for k in K])      # (K, T)
    result_arrays['lambda'] = lam_arr

    # 队列变量：shape (K, T, D+1)
    max_tau = DC_data['D']
    q_arr = np.zeros((len(K), T, max_tau + 1))
    r_arr = np.zeros((len(K), T, max_tau + 1))
    p_arr = np.zeros((len(K), T, max_tau + 1))
    for k in K:
        tau_k = int(DC_data['tau_max'][k])
        for t in range(T):
            for d in range(tau_k + 1):
                q_arr[k, t, d] = q_dc[k, t, d].X
                r_arr[k, t, d] = r_dc[k, t, d].X
                p_arr[k, t, d] = p_dc[k, t, d].X
    result_arrays['q'] = q_arr
    result_arrays['r'] = r_arr
    result_arrays['p'] = p_arr

    # ================================================================
    # 式(5) SOCP 松弛紧程度测试
    # δ_{ij,t} = (rhs - lhs) / rhs，越接近 0 表示松弛越紧
    # ================================================================
    SOCP_GAP_THRESHOLD = 1e-4
    socp_gap_arr = np.zeros((len(E), T))   # (n_branches, T)
    branch_labels = []

    for br_idx, (i, j, *_) in enumerate(E):
        branch_labels.append((i, j))
        for t in range(T):
            Ti = T_var[i, j, t].X
            Zi = Z_var[i, j, t].X
            Ui = U_var[i, t].X
            Uj = U_var[j, t].X
            lhs = 4*Ti**2 + 4*Zi**2 + (Ui - Uj)**2
            rhs = (Ui + Uj)**2
            socp_gap_arr[br_idx, t] = (rhs - lhs) / rhs if rhs > 1e-12 else 0.0

    result_arrays['socp_gap'] = socp_gap_arr
    result_arrays['branch_labels'] = np.array(branch_labels)  # (n_branches, 2)

    # 打印松紧度报告
    gap_max  = socp_gap_arr.max()
    gap_mean = socp_gap_arr.mean()
    n_over   = int((socp_gap_arr > SOCP_GAP_THRESHOLD).any(axis=1).sum())
    n_br     = len(E)
    if gap_max < SOCP_GAP_THRESHOLD:
        conclusion = "松弛严格紧，SOCP 解对应物理可行的交流潮流解"
    elif gap_max < 1e-2:
        conclusion = "松弛近似紧，需进一步核查超阈值支路"
    else:
        conclusion = "松弛偏松，存在较大间隙，需重点排查"

    print(f"\n{'='*50}")
    print(f"  式(5) SOCP 松弛紧程度报告  [场景 {scenario_id}: {scenario_data['name']}]")
    print(f"{'='*50}")
    print(f"  最大松弛间隙 (全支路全时段): {gap_max:.3e}")
    print(f"  平均松弛间隙               : {gap_mean:.3e}")
    print(f"  超过阈值({SOCP_GAP_THRESHOLD:.0e})的支路数量   : {n_over} / {n_br}")
    print(f"  结论: {conclusion}")
    if n_over > 0:
        print(f"  超阈值支路详情（最差时段）:")
        for br_idx, (bi, bj) in enumerate(branch_labels):
            br_gap = socp_gap_arr[br_idx]
            if br_gap.max() > SOCP_GAP_THRESHOLD:
                worst_t = int(br_gap.argmax())
                print(f"    支路{br_idx:>3d} ({bi+1:>2d}-{bj+1:>2d}): max_gap={br_gap.max():.3e}  最差时段 t={worst_t}")
    print(f"{'='*50}\n")

    result_meta = {
        'scenario_id':   scenario_id,
        'scenario_name': scenario_data['name'],
        'hc': hc_val,
        'scale': scale_val,
        'T': T,
        'N_buses': len(N),
        'K': K,
        'IDC_nodes': IDC_nodes,
        'baseMVA': baseMVA,
        'status': status,
        'solve_time': solve_time,
    }

    model.dispose()

    return scenario_id, hc_val, scale_val, status, solve_time, result_arrays, result_meta


# ===========================================================
# SECTION 3: 结果保存模块
# ===========================================================

def save_results(scenario_id, result_arrays, result_meta, output_dir):
    """
    保存单场景结果到磁盘

    - 变量数组 → output_dir/scenario_{id}/results.npz
    - 标量汇总 → output_dir/summary.csv（filelock 保护）
    """
    scenario_str = result_meta['scenario_name']
    sc_dir = os.path.join(output_dir, f"scenario_{scenario_str}")
    os.makedirs(sc_dir, exist_ok=True)

    # 保存 npz（压缩，快速）
    np.savez_compressed(
        os.path.join(sc_dir, "results.npz"),
        **result_arrays,
        meta_keys=list(result_meta.keys()),
        meta_vals=[str(v) for v in result_meta.values()]
    )

    # 追加写入汇总 CSV（进程间文件锁）
    summary_path = os.path.join(output_dir, "summary.csv")
    lock_path = summary_path + ".lock"
    row = {
        'scenario_id':   scenario_id,
        'scenario_name': result_meta['scenario_name'],
        'hc_kW':         result_meta['hc'] * 1000.0 * result_meta['baseMVA'],
        'scale':         result_meta['scale'],
        'status':        result_meta['status'],
        'solve_time_s':  result_meta['solve_time'],
    }
    with FileLock(lock_path):
        write_header = not os.path.exists(summary_path)
        df_row = pd.DataFrame([row])
        df_row.to_csv(summary_path, mode='a', header=write_header, index=False)


# ===========================================================
# SECTION 4: 结果读取模块
# ===========================================================

def load_results(scenario_id, output_dir, scenario_str=None):
    """
    从磁盘加载单场景结果

    Parameters
    ----------
    scenario_id  : int
    output_dir   : str
    scenario_str : str or None
        场景字符串（文件夹名后缀）。若为 None 则回退到 scenario_{scenario_id}。

    Returns
    -------
    result_arrays : dict  (numpy arrays)
    result_meta   : dict  (标量信息)
    """
    folder = f"scenario_{scenario_str}" if scenario_str is not None else f"scenario_{scenario_id}"
    sc_dir = os.path.join(output_dir, folder)
    npz = np.load(os.path.join(sc_dir, "results.npz"), allow_pickle=True)

    result_arrays = {k: npz[k] for k in npz.files
                     if k not in ('meta_keys', 'meta_vals')}

    meta_keys = list(npz['meta_keys'])
    meta_vals = list(npz['meta_vals'])
    result_meta = {}
    for k, v in zip(meta_keys, meta_vals):
        try:
            result_meta[k] = int(v)
        except ValueError:
            try:
                result_meta[k] = float(v)
            except ValueError:
                result_meta[k] = v

    return result_arrays, result_meta


# ===========================================================
# SECTION 5: 可视化模块
# ===========================================================

def visualize_results(scenario_id, result_arrays, result_meta, output_dir):
    """
    对单场景结果生成标准图表并保存（使用 Agg 后端，线程安全）
    """
    scenario_str = result_meta['scenario_name']
    fig_dir = os.path.join(output_dir, f"scenario_{scenario_str}", "figures")
    os.makedirs(fig_dir, exist_ok=True)

    T        = result_meta['T']
    baseMVA  = result_meta['baseMVA']
    hc_val   = result_meta['hc'] * 1000.0 * baseMVA
    IDC_nodes = result_meta['IDC_nodes']
    K        = result_meta['K']
    t_ax     = np.arange(T)

    U_arr    = result_arrays['U']       # (N, T)
    PIDC_arr = result_arrays['PIDC']    # (K, T)
    lam_arr  = result_arrays['lambda']  # (K, T)
    p_arr    = result_arrays['p']       # (K, T, D+1)

    # ---- 图1：节点电压幅值时序 ----
    fig, ax = plt.subplots(figsize=(12, 4))
    plot_nodes = [0, 10, 17, 23, 32] if U_arr.shape[0] > 32 else list(range(min(5, U_arr.shape[0])))
    for i in plot_nodes:
        V_vals = np.sqrt(np.maximum(U_arr[i], 0))
        ax.plot(t_ax, V_vals, marker='o', ms=2, lw=1.2, label=f"节点{i+1}")
    ax.axhline(0.9, color='red', ls='--', lw=1.2, label='下限0.9')
    ax.axhline(1.1, color='red', ls='--', lw=1.2, label='上限1.1')
    ax.set(xlabel='时段 t', ylabel='电压幅值 (p.u.)',
           title=f'节点电压时序 [场景{scenario_id}: {result_meta["scenario_name"]}]')
    ax.legend(ncol=4, fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "voltage.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---- 图2：IDC承载力时序 ----
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    # 子图1：总功率 vs hc
    ax = axes[0]
    total_kW = PIDC_arr.sum(axis=0) * 1000.0 * baseMVA
    ax.plot(t_ax, total_kW, lw=2, label='总IDC功率')
    ax.axhline(hc_val, color='red', ls='--', lw=1.5, label=f'hc = {hc_val:.1f} kW')
    ax.set(ylabel='功率 (kW)', title='IDC承载力时序')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.ticklabel_format(useOffset=False, style='plain', axis='y')
    # 子图2：各DC分解
    ax = axes[1]
    for idx, i in enumerate(IDC_nodes):
        power_kW = PIDC_arr[idx] * 1000.0 * baseMVA
        ax.plot(t_ax, power_kW, marker='o', ms=3, lw=1.5, label=f"DC{K[idx]}(节点{i+1})")
    ax.set(xlabel='时段 t', ylabel='功率 (kW)', title='各IDC功率分解')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.ticklabel_format(useOffset=False, style='plain', axis='y')
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "capacity.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---- 图3：任务到达量λ时序 ----
    fig, ax = plt.subplots(figsize=(12, 4))
    for idx, k in enumerate(K):
        ax.plot(t_ax, lam_arr[idx], marker='s', ms=3, lw=1.5, label=f"DC{k}")
    ax.set(xlabel='时段 t', ylabel='任务量 (job/period)',
           title='最优任务到达量 λ')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "lambda.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---- 图4：任务处理量堆叠图（各DC）----
    for idx, k in enumerate(K):
        fig, ax = plt.subplots(figsize=(12, 4))
        D_plus1 = p_arr.shape[2]
        bottom = np.zeros(T)
        for d in range(D_plus1):
            vals = p_arr[idx, :, d]
            if vals.max() > 1e-9:
                ax.bar(t_ax, vals, bottom=bottom, label=f"d={d}", alpha=0.8)
                bottom += vals
        ax.set(xlabel='时段 t', ylabel='处理量 (job)',
               title=f'DC{k} 任务处理量（按延迟级堆叠）')
        ax.legend(fontsize=7, ncol=4)
        ax.grid(alpha=0.3, axis='y')
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, f"task_dc{k}.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)

    # ---- 图5：式(5) SOCP 松弛间隙热力图 ----
    if 'socp_gap' in result_arrays:
        socp_gap = result_arrays['socp_gap']   # (n_branches, T)
        n_br = socp_gap.shape[0]
        fig, ax = plt.subplots(figsize=(max(8, T * 0.4), max(5, n_br * 0.25)))
        im = ax.imshow(socp_gap, aspect='auto', cmap='YlOrRd',
                       vmin=0, vmax=max(socp_gap.max(), 1e-8))
        plt.colorbar(im, ax=ax, label='松弛间隙 δ')
        ax.set(xlabel='时段 t', ylabel='支路编号',
               title=f'式(5) SOCP 松弛间隙 δ [场景{scenario_id}: {result_meta["scenario_name"]}]')
        ax.set_yticks(range(n_br))
        if 'branch_labels' in result_arrays:
            bl = result_arrays['branch_labels']
            ax.set_yticklabels([f"{int(bl[b,0])+1}-{int(bl[b,1])+1}" for b in range(n_br)],
                               fontsize=max(4, 8 - n_br // 10))
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "socp_gap.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)


# ===========================================================
# SECTION 6: 场景遍历主流程（并行求解 + 异步可视化）
# ===========================================================

def _worker(args):
    """进程池工作函数（顶层函数，可被 pickle）"""
    scenario_data, output_dir, scenario_id = args
    return build_and_solve_single(scenario_data, output_dir, scenario_id)


def run_all_scenarios(all_scenarios, output_dir, max_workers=None, viz_workers=4):
    """
    场景遍历求解主流程

    Parameters
    ----------
    all_scenarios : list of dict
    output_dir : str
    max_workers : int or None
        进程池大小，默认 cpu_count // 2（每个进程自用剩余核）
    viz_workers : int
        可视化线程池大小
    """
    os.makedirs(output_dir, exist_ok=True)
    N_total = len(all_scenarios)
    if max_workers is None:
        max_workers = max(1, (os.cpu_count() or 2) // 2)

    print(f"\n{'='*70}")
    print(f"场景遍历求解")
    print(f"  总场景数: {N_total}")
    print(f"  求解进程数: {max_workers}")
    print(f"  可视化线程数: {viz_workers}")
    print(f"{'='*70}\n")

    args_list = [(sd, output_dir, idx) for idx, sd in enumerate(all_scenarios)]

    # 线程池用于异步可视化（I/O密集，不阻塞求解）
    viz_pool = ThreadPoolExecutor(max_workers=viz_workers)
    viz_futures = []

    with tqdm(total=N_total, desc="场景求解进度", unit="场景",
              dynamic_ncols=True, leave=True) as progress:
        with ProcessPoolExecutor(max_workers=max_workers) as proc_pool:
            future_map = {proc_pool.submit(_worker, args): args[2] for args in args_list}

            for future in as_completed(future_map):
                sid = future_map[future]
                try:
                    ret = future.result()
                except Exception as exc:
                    tqdm.write(f"[ERR] 场景{sid} 求解异常: {exc}")
                    progress.update(1)
                    continue

                if ret[1] is None:
                    tqdm.write(f"[FAIL] 场景{sid} 无可行解 (status={ret[3]}, t={ret[4]:.1f}s)")
                    progress.update(1)
                    continue

                scenario_id, hc_val, scale_val, _, solve_time, result_arrays, result_meta = ret
                hc_kW = hc_val * 1000.0 * result_meta['baseMVA']
                tqdm.write(f"[OK] 场景{scenario_id:>3d} | hc={hc_kW:.1f} kW | scale={scale_val:.4f} | t={solve_time:.1f}s")

                # 同步保存结果（轻量，主进程执行）
                save_results(scenario_id, result_arrays, result_meta, output_dir)

                # 异步提交可视化任务（不阻塞）
                vf = viz_pool.submit(visualize_results, scenario_id, result_arrays, result_meta, output_dir)
                viz_futures.append(vf)

                progress.update(1)

    # 等待所有可视化任务完成（带进度条）
    n_viz = len(viz_futures)
    tqdm.write(f"\n等待 {n_viz} 个可视化任务完成...")
    with tqdm(total=n_viz, desc="可视化进度", unit="场景",
              dynamic_ncols=True, leave=True) as viz_progress:
        for vf in as_completed(viz_futures):
            try:
                vf.result()
            except Exception as exc:
                tqdm.write(f"[VIZ ERR] {exc}")
            viz_progress.update(1)

    viz_pool.shutdown(wait=False)
    tqdm.write("所有场景求解与可视化完成。\n")


# ===========================================================
# SECTION 7: 跨场景汇总分析
# ===========================================================

def analyze_all_scenarios(output_dir):
    """
    所有场景求解完毕后，读取结果并进行跨场景汇总分析

    Parameters
    ----------
    output_dir : str
    all_scenarios : list of dict  (用于获取场景名称等元信息，当前由 summary.csv 提供）
    """
    summary_path = os.path.join(output_dir, "summary.csv")
    if not os.path.exists(summary_path):
        print("[ANALYSIS] 未找到 summary.csv，跳过汇总分析")
        return

    analysis_dir = os.path.join(output_dir, "analysis")
    os.makedirs(analysis_dir, exist_ok=True)

    df_sum = pd.read_csv(summary_path)
    solved = df_sum[df_sum['hc_kW'].notna()].copy()

    print(f"\n{'='*70}")
    print(f"跨场景汇总分析")
    print(f"  有效场景数: {len(solved)} / {len(df_sum)}")
    print(f"  hc_kW 均值: {solved['hc_kW'].mean():.2f}")
    print(f"  hc_kW 最小: {solved['hc_kW'].min():.2f}")
    print(f"  hc_kW 最大: {solved['hc_kW'].max():.2f}")
    print(f"{'='*70}\n")

    # ---- 图1：hc 分布直方图 ----
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.hist(solved['hc_kW'], bins=20, edgecolor='black', alpha=0.8)
    ax.axvline(solved['hc_kW'].mean(), color='red', ls='--', lw=1.5,
               label=f"均值={solved['hc_kW'].mean():.1f} kW")
    ax.set(xlabel='承载力 hc (kW)', ylabel='场景数', title='全场景承载力分布')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(analysis_dir, "hc_distribution.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---- 图2：hc 时序（按场景排序）----
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(solved['scenario_id'], solved['hc_kW'], alpha=0.8)
    ax.set(xlabel='场景编号', ylabel='hc (kW)', title='各场景承载力下界 hc')
    ax.grid(alpha=0.3, axis='y')
    ax.ticklabel_format(useOffset=False, style='plain', axis='y')
    fig.tight_layout()
    fig.savefig(os.path.join(analysis_dir, "hc_per_scenario.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---- 图3：scale 与 hc 散点图 ----
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(solved['scale'], solved['hc_kW'], alpha=0.7, s=20)
    ax.set(xlabel='scale', ylabel='hc (kW)', title='scale vs hc（线性关系验证）')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(analysis_dir, "scale_vs_hc.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---- 跨场景电压越限频次（需读取各场景 npz）----
    valid_ids = sorted(solved['scenario_id'].tolist())
    sid_to_name = dict(zip(solved['scenario_id'].tolist(), solved['scenario_name'].tolist()))
    N_buses = None
    voltage_min_all = []   # 各场景各节点最低电压

    for sid in valid_ids:
        try:
            arrays, _ = load_results(sid, output_dir, scenario_str=sid_to_name.get(sid))
        except Exception:
            continue
        if N_buses is None:
            N_buses = arrays['U'].shape[0]
        V_min_i = np.sqrt(np.maximum(arrays['U'], 0)).min(axis=1)  # (N,)
        voltage_min_all.append(V_min_i)

    if voltage_min_all:
        V_min_mat = np.stack(voltage_min_all, axis=0)   # (S, N)
        violation_rate = (V_min_mat < 0.9).mean(axis=0)  # (N,)

        fig, ax = plt.subplots(figsize=(12, 4))
        ax.bar(range(1, N_buses+1), violation_rate * 100, alpha=0.8)
        ax.set(xlabel='节点编号', ylabel='电压越限率 (%)',
               title='各节点低电压越限频次（全场景）')
        ax.grid(alpha=0.3, axis='y')
        fig.tight_layout()
        fig.savefig(os.path.join(analysis_dir, "voltage_violation.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)

    # 保存汇总报告
    solved.describe().to_csv(os.path.join(analysis_dir, "report.csv"))
    print(f"汇总分析完成，结果保存至 {analysis_dir}\n")


# ===========================================================
# MAIN
# ===========================================================

if __name__ == "__main__":

    # ================================================================
    # 场景输入控制参数（见 power_system_model_onetest.md § 场景输入控制）
    #
    #   SCENARIO_INDEX = 0      → 只跑第 0 组（单场景串行，最适合调试）
    #   SCENARIO_INDEX = 7      → 只跑第 7 组
    #   SCENARIO_INDEX = None   → 跑前 N_SCENARIOS 组（支持并行）
    # ================================================================
    SCENARIO_INDEX = 1     # int 或 None
    N_SCENARIOS    = 3     # 当 SCENARIO_INDEX 为 None 时生效

    # ================================================================
    # 输出目录（固定为脚本所在目录下 test_output/output_onetest/）
    # ================================================================
    output_dir = os.path.join(current_dir, "test_output", "output_onetest")
    os.makedirs(output_dir, exist_ok=True)

    # ---- 生成全量场景（笛卡尔积）----
    result_scenario = generate_scenario_combinations(
        '海发光伏电站',
        '九岗光伏电站',
        '智良站10kV备用线528',
        '云航站10kVF11线',
        5
    )
    print(f"总共生成 {len(result_scenario)} 种场景组合")

    csv_file_path = os.path.join(current_dir, 'test_input/光伏与数据中心典型日曲线.csv')

    all_scenarios = []
    for combo in tqdm(result_scenario, desc="处理场景组合", unit="个"):
        dc0, dc1, load18, load25, scenario_str = read_scenario_data(
            combo[2], combo[3], combo[0], combo[1], csv_file_path
        )
        all_scenarios.append({
            'name': scenario_str,
            'dc0_lambda': np.array(dc0),
            'dc1_lambda': np.array(dc1),
            'loadlist18': np.array(load18) * 0.1,
            'loadlist25': np.array(load25) * 0.1,
        })

    # ---- 按 SCENARIO_INDEX / N_SCENARIOS 截取待测场景 ----
    if SCENARIO_INDEX is not None:
        selected = [all_scenarios[SCENARIO_INDEX]]
        print(f"\n单场景模式：运行第 {SCENARIO_INDEX} 组场景 [{selected[0]['name']}]")
    else:
        selected = all_scenarios[:N_SCENARIOS]
        print(f"\n批量模式：运行前 {len(selected)} 组场景")

    print(f"输出目录: {output_dir}\n")

    # ---- 求解 ----
    if len(selected) == 1:
        # 单场景串行：直接在主进程运行，便于调试
        sd = selected[0]
        sid = SCENARIO_INDEX if SCENARIO_INDEX is not None else 0
        ret = build_and_solve_single(sd, output_dir, sid)
        if ret[1] is None:
            print(f"[FAIL] 场景{sid} 无可行解 (status={ret[3]}, t={ret[4]:.1f}s)")
        else:
            scenario_id, hc_val, scale_val, _, solve_time, result_arrays, result_meta = ret
            hc_kW = hc_val * 1000.0 * result_meta['baseMVA']
            print(f"[OK] 场景{scenario_id} | hc={hc_kW:.1f} kW | scale={scale_val:.4f} | t={solve_time:.1f}s")
            save_results(scenario_id, result_arrays, result_meta, output_dir)
            visualize_results(scenario_id, result_arrays, result_meta, output_dir)
            print(f"结果已保存至 {output_dir}/scenario_{result_meta['scenario_name']}/")
    else:
        # 多场景并行（ProcessPoolExecutor + 异步可视化）
        run_all_scenarios(selected, output_dir, max_workers=None, viz_workers=4)
        analyze_all_scenarios(output_dir)

    print(f"\n{'='*70}")
    print("任务完成！")
    print(f"  结果目录: {output_dir}")
    print(f"{'='*70}\n")
