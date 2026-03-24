"""
场景组合生成器
用于生成所有可能的场景组合
"""
from itertools import product
import pandas as pd
import numpy as np


def generate_scenario_combinations(str1, str2, str3, str4, N):
    """
    生成所有可能的场景组合

    参数:
        str1 (str): 第一个设备名称
        str2 (str): 第二个设备名称
        str3 (str): 第三个设备名称
        str4 (str): 第四个设备名称
        N (int): 场景数目

    返回:
        list: 包含所有场景组合的列表，每个组合是一个包含4个元素的列表

    示例:
        >>> generate_scenario_combinations('海发光伏电站', '九岗光伏电站',
                                          '智良站10kV备用线528', '云航站10kVF11线', 5)
        [['海发光伏电站_场景1', '九岗光伏电站_场景1', '智良站10kV备用线528_场景1', '云航站10kVF11线_场景1'],
         ['海发光伏电站_场景1', '九岗光伏电站_场景1', '智良站10kV备用线528_场景1', '云航站10kVF11线_场景2'],
         ...
        ]
    """
    # 设备名称列表
    devices = [str1, str2, str3, str4]

    # 场景范围 (1到N)
    scenarios = range(1, N + 1)

    # 生成所有可能的场景组合
    all_combinations = []

    # 使用itertools.product生成笛卡尔积
    for combo in product(scenarios, repeat=4):
        # combo是一个元组，如 (1, 1, 1, 1), (1, 1, 1, 2), ...
        combination = [f"{device}_场景{scene}" for device, scene in zip(devices, combo)]
        all_combinations.append(combination)

    return all_combinations


def read_scenario_data(dc0_name, dc1_name, loadlist18_name, loadlist25_name, csv_path):
    """
    从CSV文件中读取指定列的场景数据

    参数:
        dc0_name (str): 第一个数据中心列名（返回为numpy array）
        dc1_name (str): 第二个数据中心列名（返回为numpy array）
        loadlist18_name (str): 第一个负荷列名（返回为list）
        loadlist25_name (str): 第二个负荷列名（返回为list）
        csv_path (str): CSV文件路径

    返回:
        tuple: (dc0_data, dc1_data, loadlist18_data, loadlist25_data, scenario_str)
               - dc0_data: numpy.ndarray
               - dc1_data: numpy.ndarray
               - loadlist18_data: list
               - loadlist25_data: list
               - scenario_str: str (格式: "场景" + 四个参数名称各自的最后一个字符)

    示例:
        >>> dc0, dc1, load18, load25, scenario = read_scenario_data(
                '智良站10kV备用线528_场景1',
                '云航站10kVF11线_场景1',
                '海发光伏电站_场景1',
                '九岗光伏电站_场景1',
                'code/SCOP_HC/data/光伏与数据中心典型日曲线.csv'
            )
        >>> print(scenario)  # 输出: 场景1111
    """
    # 读取CSV文件，跳过第二行（概率行）
    # 先读取所有数据
    df_full = pd.read_csv(csv_path)

    # 删除第一行数据（即概率行，在pandas中是索引0）
    df = df_full.drop(index=0).reset_index(drop=True)

    # 提取对应列的数据
    dc0_data = np.array(df[dc0_name].values, dtype=float)
    dc1_data = np.array(df[dc1_name].values, dtype=float)
    loadlist18_data = df[loadlist18_name].tolist()
    loadlist25_data = df[loadlist25_name].tolist()

    # 生成场景字符串：场景 + 四个参数名称各自的最后一个字符
    scenario_str = "场景" + dc0_name[-1] + dc1_name[-1] + loadlist18_name[-1] + loadlist25_name[-1]

    return dc0_data, dc1_data, loadlist18_data, loadlist25_data, scenario_str


# # 示例使用
# if __name__ == "__main__":
#     # 示例1: 测试场景组合生成
#     print("=" * 60)
#     print("示例1: 生成场景组合")
#     print("=" * 60)
#     result = generate_scenario_combinations(
#         '海发光伏电站',
#         '九岗光伏电站',
#         '智良站10kV备用线528',
#         '云航站10kVF11线',
#         5
#     )

#     print(f"总共生成 {len(result)} 种场景组合\n")
#     print("前5个组合:")
#     for i, combo in enumerate(result[:5], 1):
#         print(f"{i}. {combo}")

#     print("\n...")
#     print("\n后5个组合:")
#     for i, combo in enumerate(result[-5:], len(result)-4):
#         print(f"{i}. {combo}")

#     # 示例2: 读取CSV数据
#     print("\n\n" + "=" * 60)
#     print("示例2: 读取场景数据")
#     print("=" * 60)

#     csv_file_path = 'code/SCOP_HC/data/光伏与数据中心典型日曲线.csv'

#     dc0, dc1, load18, load25 = read_scenario_data(
#         '智良站10kV备用线528_场景1',
#         '云航站10kVF11线_场景1',
#         '海发光伏电站_场景1',
#         '九岗光伏电站_场景1',
#         csv_file_path
#     )

#     print(f"\ndc0 (智良站10kV备用线528_场景1) 类型: {type(dc0)}")
#     print(f"dc0 形状: {dc0.shape}")
#     print(f"dc0 前5个值: {dc0[:5]}")

#     print(f"\ndc1 (云航站10kVF11线_场景1) 类型: {type(dc1)}")
#     print(f"dc1 形状: {dc1.shape}")
#     print(f"dc1 前5个值: {dc1[:5]}")

#     print(f"\nload18 (海发光伏电站_场景1) 类型: {type(load18)}")
#     print(f"load18 长度: {len(load18)}")
#     print(f"load18 前5个值: {load18[:5]}")

#     print(f"\nload25 (九岗光伏电站_场景1) 类型: {type(load25)}")
#     print(f"load25 长度: {len(load25)}")
#     print(f"load25 前5个值: {load25[:5]}")
