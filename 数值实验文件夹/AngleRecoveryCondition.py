import networkx as nx
import numpy as np
import cmath

# Step 1: 构建树T（有向图）
def build_tree(ppc):
    G = nx.DiGraph()
    for idx, branch in enumerate(ppc["branch"]):
        fbus, tbus = int(branch[0]) - 1, int(branch[1]) - 1
        G.add_edge(fbus, tbus)
    balance_node = [i for i, row in enumerate(ppc["bus"]) if row[1] == 3][0]
    T = nx.dfs_tree(G, source=balance_node)
    tree_edges = list(T.edges())
    return T, tree_edges, G

# Step 2: 构建节点支路关联矩阵A
def build_adjacency_matrix(ppc, tree_edges, num_nodes):
    A = np.zeros((num_nodes, len(ppc["branch"])))
    tree_branch_indices = []
    for idx, branch in enumerate(ppc["branch"]):
        fbus, tbus = int(branch[0]) - 1, int(branch[1]) - 1
        if (fbus, tbus) in tree_edges or (tbus, fbus) in tree_edges:
            tree_branch_indices.append(idx)
            if (fbus, tbus) in tree_edges:
                A[fbus, idx] = 1
                A[tbus, idx] = -1
            else:
                A[fbus, idx] = -1
                A[tbus, idx] = 1
        else:
            A[fbus, idx] = 1
            A[tbus, idx] = -1
    non_tree_branch_indices = [i for i in range(len(ppc["branch"])) if i not in tree_branch_indices]
    A = A[:, tree_branch_indices + non_tree_branch_indices]
    A = A[1:, :]
    return A, tree_branch_indices

# Step 3: 分别计算A_t和A_l
def get_tree_and_loops_matrices(A, tree_edges,ppc):
    tree_branch_indices = []
    for idx, branch in enumerate(ppc["branch"]):
        fbus, tbus = int(branch[0]) - 1, int(branch[1]) - 1
        if (fbus, tbus) in tree_edges or (tbus, fbus) in tree_edges:
            tree_branch_indices.append(idx)
    A_t = A[:, :len(tree_edges)]
    A_l = A[:, len(tree_edges):]
    return A_t, A_l

# Step 4: 计算基本回路关联矩阵B_t
def compute_basic_loop_matrix(A_t, A_l):
    A_t_inv = np.linalg.inv(A_t)
    B_t = -np.dot(A_t_inv, A_l)
    return B_t.T

# Step 5: 计算基本回路关联矩阵B_f
def compute_full_loop_matrix(B_t, A_l):
    num_loops = A_l.shape[1]
    E_l = np.eye(num_loops)
    B_f = np.hstack([B_t, E_l])
    return B_f

# Step 6: 计算 beta 向量并排序
def compute_beta(E, P, Q, v, t_idx, tree_branch_indices):
    beta = {}
    for (start_bus, end_bus, r_ij, x_ij) in E:
        valP = P[(t_idx, start_bus, end_bus, r_ij, x_ij)].X
        valQ = Q[(t_idx, start_bus, end_bus, r_ij, x_ij)].X
        S_ij = complex(valP, valQ)
        z_ij_conj = complex(r_ij, -x_ij)
        zS = z_ij_conj * S_ij
        v_i = complex(v[t_idx, start_bus].X, 0)
        beta_ij = cmath.phase(v_i - zS) * 180 / np.pi
        beta_ij = (beta_ij + 180) % 360 - 180
        beta[(start_bus, end_bus)] = beta_ij

    beta_t = np.zeros(len(tree_branch_indices))
    beta_l = np.zeros(len(E) - len(tree_branch_indices))

    for i, idx in enumerate(tree_branch_indices):
        start_bus, end_bus, r_ij, x_ij = E[idx]
        beta_t[i] = beta[(start_bus, end_bus)]

    non_tree_branch_indices = [i for i in range(len(E)) if i not in tree_branch_indices]
    for i, idx in enumerate(non_tree_branch_indices):
        start_bus, end_bus, r_ij, x_ij = E[idx]
        beta_l[i] = beta[(start_bus, end_bus)]

    beta_combined = np.concatenate([beta_t, beta_l])

    return beta_combined,beta_t,beta_l

# 主程序：计算 B_f @ beta_combined
def compute_Bf_beta(ppc, E, P, Q, v, t_idx):
    T, tree_edges, G = build_tree(ppc)
    num_nodes = len(ppc["bus"])
    num_branches = len(ppc["branch"])

    # 构建节点支路关联矩阵A
    A, tree_branch_indices = build_adjacency_matrix(ppc, tree_edges, num_nodes)

    # 获取树支和余支的矩阵A_t, A_l
    A_t, A_l = get_tree_and_loops_matrices(A, tree_edges,ppc)

    # 计算基本回路关联矩阵B_t
    B_t = compute_basic_loop_matrix(A_t, A_l)

    # 计算最终的基本回路关联矩阵B_f
    B_f = compute_full_loop_matrix(B_t, A_l)

    # 计算并排列 beta 向量
    beta_combined,beta_t,beta_l = np.array(compute_beta(E, P, Q, v, t_idx, tree_branch_indices))

    # 返回 B_f @ beta_combined 的结果
    result = np.dot(B_f, beta_combined)
    if len(result) == 0:
        print( f"ppc可能是一个辐射网络，因此没有连支，结果为空")
    return result,beta_t,A_t


