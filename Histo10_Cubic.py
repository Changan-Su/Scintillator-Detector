# -*- coding: utf-8 -*-
"""
Histo10_Cubic：立方体晶体六面 SiPM 光子计数批处理与热力图脚本
================================================================

针对 Geant4 输出的「每晶体 6 个面、每面 4×4 SiPM」的 PhotonFaceBlockEvent CSV，
批量合并多线程文件、写出合并 CSV 与聚合 CSV，并生成 6 面（2×3 子图）光子计数热力图。

使用方法
--------
  python Histo10_Cubic.py <ResultsDir> [选项]
  python Histo10_Cubic.py --results <ResultsDir> [选项]

参数
----
  results               Results 根目录（位置参数或 --results）
  --results PATH        Results 根目录（与位置参数二选一）
  --output PATH         输出根目录；不指定则为 当前目录/Output/SiPM6_Output_<时间戳>
  --cmap NAME           Matplotlib 色图（默认 viridis）
  --vmin F / --vmax F   热力图颜色范围（不指定则自动）

输入
----
  每个配置子目录下需有：
    AnaEx01_nt_PhotonFaceBlockEvent_t*.csv  （多线程）或
    AnaEx01_nt_PhotonFaceBlockEvent.csv     （单线程/合并后）
  列：EventID, CrystalID, iy, iz, Face, j, k, SiPMBlockID, PhotonCount

输出（每个配置子目录对应 Output 下同名子目录）
----------------------------------------------
  merged_event.csv       各线程事件行合并（原始 9 列）
  merged_face_jk.csv     按 (Face, j, k) 聚合后的计数表，用于绘图
  reconstructed_position.csv  多种重建算法结果（每算法一行）
  sipm_6faces_heatmap.png  6 面 4×4 热力图（Face 0..5：+X,-X,+Y,-Y,+Z,-Z）

示例
----
  python Histo10_Cubic.py Results
  python Histo10_Cubic.py Results --output ./Output/SiPM6_Out --cmap plasma
"""

import argparse
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 六面编号与几何方向对应（与 DetectorConstruction 中 face 定义一致）
FACE_NAMES = {0: "+X", 1: "-X", 2: "+Y", 3: "-Y", 4: "+Z", 5: "-Z"}
DEFAULT_HALF_LENGTHS_CM = {"x": 1.25, "y": 1.25, "z": 1.25}

# 每个面的热力图坐标语义：
# imshow 中横轴对应 k，纵轴对应 j；这里补充其在世界坐标中的方向含义。
FACE_AXIS_LABELS = {
    0: ("k ( +Z )", "j ( +Y )"),  # +X 面：offset=(0,u,v)
    1: ("k ( +Z )", "j ( +Y )"),  # -X 面：offset=(0,u,v)
    2: ("k ( +Z )", "j ( +X )"),  # +Y 面：offset=(u,0,v)
    3: ("k ( +Z )", "j ( +X )"),  # -Y 面：offset=(u,0,v)
    4: ("k ( +Y )", "j ( +X )"),  # +Z 面：offset=(u,v,0)
    5: ("k ( +Y )", "j ( +X )"),  # -Z 面：offset=(u,v,0)
}

# 每面 SiPM 网格边长（与 Geant4 中 SiPm_np 一致）
SIPM_GRID_N = 4

# ---------------------------------------------------------------------------
# 算法 half_side_ratio：用「沿某轴左右半空间」选中的 SiPM 光子数之和（非整面总和）
#
# 约定：
# - 矩阵与 build_face_matrix 一致：matrix[j, k]，j、k ∈ 0..SIPM_GRID_N-1
# - Face 0..5 含义见 FACE_NAMES / FACE_AXIS_LABELS
#
# X 轴示例（当前默认实现）：
# - +X 半空间：+X 面整面 4x4 + 四个侧面（±Y、±Z）上靠 +X 的 2x4 半面
# - -X 半空间：-X 面整面 4x4 + 四个侧面（±Y、±Z）上靠 -X 的 2x4 半面
#
# Y、Z 轴按同一思路对称推广：主面整面 4x4 + 四个相关侧面的 2x4 半面。
# 若你的物理定义不同，只需改下列表中的 (face, j, k)，不必改 reconstruct_half_side_ratio 的比值形式。
# ---------------------------------------------------------------------------

# X：+X / -X 半空间
HALF_SIDE_X_PLUS_CELLS: list[tuple[int, int, int]] = [
    # Face 0 (+X)：整面 4x4
    (0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3),
    (0, 1, 0), (0, 1, 1), (0, 1, 2), (0, 1, 3),
    (0, 2, 0), (0, 2, 1), (0, 2, 2), (0, 2, 3),
    (0, 3, 0), (0, 3, 1), (0, 3, 2), (0, 3, 3),
    # Face 2 (+Y)、3 (-Y)、4 (+Z)、5 (-Z)：靠 +X 的 2x4 半面（j=2,3）
    (2, 2, 0), (2, 2, 1), (2, 2, 2), (2, 2, 3),
    (2, 3, 0), (2, 3, 1), (2, 3, 2), (2, 3, 3),
    (3, 2, 0), (3, 2, 1), (3, 2, 2), (3, 2, 3),
    (3, 3, 0), (3, 3, 1), (3, 3, 2), (3, 3, 3),
    (4, 2, 0), (4, 2, 1), (4, 2, 2), (4, 2, 3),
    (4, 3, 0), (4, 3, 1), (4, 3, 2), (4, 3, 3),
    (5, 2, 0), (5, 2, 1), (5, 2, 2), (5, 2, 3),
    (5, 3, 0), (5, 3, 1), (5, 3, 2), (5, 3, 3),
]
HALF_SIDE_X_MINUS_CELLS: list[tuple[int, int, int]] = [
    # Face 1 (-X)：整面 4x4
    (1, 0, 0), (1, 0, 1), (1, 0, 2), (1, 0, 3),
    (1, 1, 0), (1, 1, 1), (1, 1, 2), (1, 1, 3),
    (1, 2, 0), (1, 2, 1), (1, 2, 2), (1, 2, 3),
    (1, 3, 0), (1, 3, 1), (1, 3, 2), (1, 3, 3),
    # Face 2 (+Y)、3 (-Y)、4 (+Z)、5 (-Z)：靠 -X 的 2x4 半面（j=0,1）
    (2, 0, 0), (2, 0, 1), (2, 0, 2), (2, 0, 3),
    (2, 1, 0), (2, 1, 1), (2, 1, 2), (2, 1, 3),
    (3, 0, 0), (3, 0, 1), (3, 0, 2), (3, 0, 3),
    (3, 1, 0), (3, 1, 1), (3, 1, 2), (3, 1, 3),
    (4, 0, 0), (4, 0, 1), (4, 0, 2), (4, 0, 3),
    (4, 1, 0), (4, 1, 1), (4, 1, 2), (4, 1, 3),
    (5, 0, 0), (5, 0, 1), (5, 0, 2), (5, 0, 3),
    (5, 1, 0), (5, 1, 1), (5, 1, 2), (5, 1, 3),
]

# Y：+Y / -Y 半空间（侧面为 ±X、±Z）
HALF_SIDE_Y_PLUS_CELLS: list[tuple[int, int, int]] = [
    # Face 2 (+Y)：整面 4x4
    (2, 0, 0), (2, 0, 1), (2, 0, 2), (2, 0, 3),
    (2, 1, 0), (2, 1, 1), (2, 1, 2), (2, 1, 3),
    (2, 2, 0), (2, 2, 1), (2, 2, 2), (2, 2, 3),
    (2, 3, 0), (2, 3, 1), (2, 3, 2), (2, 3, 3),
    # Face 0 (+X)、1 (-X)：靠 +Y 的 2x4 半面（j=2,3）
    (0, 2, 0), (0, 2, 1), (0, 2, 2), (0, 2, 3),
    (0, 3, 0), (0, 3, 1), (0, 3, 2), (0, 3, 3),
    (1, 2, 0), (1, 2, 1), (1, 2, 2), (1, 2, 3),
    (1, 3, 0), (1, 3, 1), (1, 3, 2), (1, 3, 3),
    # Face 4 (+Z)、5 (-Z)：靠 +Y 的 2x4 半面（k=2,3）
    (4, 0, 2), (4, 0, 3), (4, 1, 2), (4, 1, 3),
    (4, 2, 2), (4, 2, 3), (4, 3, 2), (4, 3, 3),
    (5, 0, 2), (5, 0, 3), (5, 1, 2), (5, 1, 3),
    (5, 2, 2), (5, 2, 3), (5, 3, 2), (5, 3, 3),
]
HALF_SIDE_Y_MINUS_CELLS: list[tuple[int, int, int]] = [
    # Face 3 (-Y)：整面 4x4
    (3, 0, 0), (3, 0, 1), (3, 0, 2), (3, 0, 3),
    (3, 1, 0), (3, 1, 1), (3, 1, 2), (3, 1, 3),
    (3, 2, 0), (3, 2, 1), (3, 2, 2), (3, 2, 3),
    (3, 3, 0), (3, 3, 1), (3, 3, 2), (3, 3, 3),
    # Face 0 (+X)、1 (-X)：靠 -Y 的 2x4 半面（j=0,1）
    (0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3),
    (0, 1, 0), (0, 1, 1), (0, 1, 2), (0, 1, 3),
    (1, 0, 0), (1, 0, 1), (1, 0, 2), (1, 0, 3),
    (1, 1, 0), (1, 1, 1), (1, 1, 2), (1, 1, 3),
    # Face 4 (+Z)、5 (-Z)：靠 -Y 的 2x4 半面（k=0,1）
    (4, 0, 0), (4, 0, 1), (4, 1, 0), (4, 1, 1),
    (4, 2, 0), (4, 2, 1), (4, 3, 0), (4, 3, 1),
    (5, 0, 0), (5, 0, 1), (5, 1, 0), (5, 1, 1),
    (5, 2, 0), (5, 2, 1), (5, 3, 0), (5, 3, 1),
]

# Z：+Z / -Z 半空间（侧面为 ±X、±Y）
HALF_SIDE_Z_PLUS_CELLS: list[tuple[int, int, int]] = [
    # Face 4 (+Z)：整面 4x4
    (4, 0, 0), (4, 0, 1), (4, 0, 2), (4, 0, 3),
    (4, 1, 0), (4, 1, 1), (4, 1, 2), (4, 1, 3),
    (4, 2, 0), (4, 2, 1), (4, 2, 2), (4, 2, 3),
    (4, 3, 0), (4, 3, 1), (4, 3, 2), (4, 3, 3),
    # Face 0 (+X)、1 (-X)：靠 +Z 的 2x4 半面（k=2,3）
    (0, 0, 2), (0, 0, 3), (0, 1, 2), (0, 1, 3),
    (0, 2, 2), (0, 2, 3), (0, 3, 2), (0, 3, 3),
    (1, 0, 2), (1, 0, 3), (1, 1, 2), (1, 1, 3),
    (1, 2, 2), (1, 2, 3), (1, 3, 2), (1, 3, 3),
    # Face 2 (+Y)、3 (-Y)：靠 +Z 的 2x4 半面（k=2,3）
    (2, 0, 2), (2, 0, 3), (2, 1, 2), (2, 1, 3),
    (2, 2, 2), (2, 2, 3), (2, 3, 2), (2, 3, 3),
    (3, 0, 2), (3, 0, 3), (3, 1, 2), (3, 1, 3),
    (3, 2, 2), (3, 2, 3), (3, 3, 2), (3, 3, 3),
]
HALF_SIDE_Z_MINUS_CELLS: list[tuple[int, int, int]] = [
    # Face 5 (-Z)：整面 4x4
    (5, 0, 0), (5, 0, 1), (5, 0, 2), (5, 0, 3),
    (5, 1, 0), (5, 1, 1), (5, 1, 2), (5, 1, 3),
    (5, 2, 0), (5, 2, 1), (5, 2, 2), (5, 2, 3),
    (5, 3, 0), (5, 3, 1), (5, 3, 2), (5, 3, 3),
    # Face 0 (+X)、1 (-X)：靠 -Z 的 2x4 半面（k=0,1）
    (0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
    (0, 2, 0), (0, 2, 1), (0, 3, 0), (0, 3, 1),
    (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1),
    (1, 2, 0), (1, 2, 1), (1, 3, 0), (1, 3, 1),
    # Face 2 (+Y)、3 (-Y)：靠 -Z 的 2x4 半面（k=0,1）
    (2, 0, 0), (2, 0, 1), (2, 1, 0), (2, 1, 1),
    (2, 2, 0), (2, 2, 1), (2, 3, 0), (2, 3, 1),
    (3, 0, 0), (3, 0, 1), (3, 1, 0), (3, 1, 1),
    (3, 2, 0), (3, 2, 1), (3, 3, 0), (3, 3, 1),
]


def sum_sipm_cells(
    matrices: dict[int, np.ndarray],
    cells: list[tuple[int, int, int]],
) -> float:
    """
    对给定 (Face, j, k) 列表求光子计数之和。
    缺失的面或越界下标按 0 计（与「该格无数据」一致）。
    """
    total = 0.0
    for face, j, k in cells:
        mat = matrices.get(face)
        if mat is None:
            continue
        if 0 <= j < SIPM_GRID_N and 0 <= k < SIPM_GRID_N:
            total += float(mat[j, k])
    return total


def reconstruct_half_side_ratio(
    matrices: dict[int, np.ndarray],
    half_lengths_cm: dict[str, float],
) -> tuple[float, float, float]:
    """
    算法：half_side_ratio

    用「沿 ±X / ±Y / ±Z 划分的半空间」内、按 HALF_SIDE_*_CELLS 选中的 SiPM 总光子数，
    做与 linear_scaled 同形的比值映射到 cm：

      axis_rec = half * (N(+) - N(-)) / (N(+) + N(-))

    其中 N(+)、N(-) 来自 sum_sipm_cells，不是整面 face_totals。

    具体包含哪些格子由模块级 HALF_SIDE_X_PLUS_CELLS 等列表定义，可自行修改。
    """
    n_xp = sum_sipm_cells(matrices, HALF_SIDE_X_PLUS_CELLS)
    n_xm = sum_sipm_cells(matrices, HALF_SIDE_X_MINUS_CELLS)
    n_yp = sum_sipm_cells(matrices, HALF_SIDE_Y_PLUS_CELLS)
    n_ym = sum_sipm_cells(matrices, HALF_SIDE_Y_MINUS_CELLS)
    n_zp = sum_sipm_cells(matrices, HALF_SIDE_Z_PLUS_CELLS)
    n_zm = sum_sipm_cells(matrices, HALF_SIDE_Z_MINUS_CELLS)

    hx, hy, hz = half_lengths_cm["x"], half_lengths_cm["y"], half_lengths_cm["z"]
    denom_x = n_xp + n_xm
    denom_y = n_yp + n_ym
    denom_z = n_zp + n_zm

    angle_x = np.arccos(1-2 * n_xp/denom_x)
    angle_y = np.arccos(1-2 *n_yp/denom_y)
    angle_z = np.arccos(1-2 *n_zp/denom_z)
    x_rec = -hx / np.tan(angle_x)
    y_rec = -hy / np.tan(angle_y)
    z_rec = -hz / np.tan(angle_z)
    return x_rec, y_rec, z_rec


# Geant4 PhotonFaceBlockEvent CSV 的列名（9 列）
EVENT_COLUMNS = ["EventID", "CrystalID", "iy", "iz", "Face", "j", "k", "SiPMBlockID", "PhotonCount"]


def parse_event_csv(csv_path: Path) -> pd.DataFrame:
    """
    读取单份 PhotonFaceBlockEvent CSV（无表头、# 注释行），返回规范 9 列 DataFrame。
    列：EventID, CrystalID, iy, iz, Face, j, k, SiPMBlockID, PhotonCount；
    Face/j/k/PhotonCount 转为数值并丢弃含无效值的行。
    """
    df = pd.read_csv(csv_path, comment="#", header=None, sep=",")
    if df.empty or df.shape[1] < 9:
        raise ValueError(f"Invalid or empty event CSV: {csv_path}")
    df = df.iloc[:, :9].copy()#按整数位置（行号、列号）选取数据，和列名无关。
    df.columns = EVENT_COLUMNS#列名赋值
    for col in EVENT_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")#将指定列转换为数值类型，遇到错误时替换为NaN
    df = df.dropna(subset=["Face", "j", "k", "PhotonCount"]).copy()#删除包含NaN的行
    df["Face"] = df["Face"].astype(int)#将指定列转换为整数类型
    df["j"] = df["j"].astype(int)#将指定列转换为整数类型
    df["k"] = df["k"].astype(int)#将指定列转换为整数类型
    return df#返回处理后的DataFrame 


def aggregate_for_heatmap(df_event: pd.DataFrame) -> pd.DataFrame:
    """
    按 (Face, j, k) 对事件表聚合光子计数，得到每格总计数。
    返回列：Face, j, k, Count。用于绘制 6 个 4×4 热力图。
    """
    agg = df_event.groupby(["Face", "j", "k"], as_index=False)["PhotonCount"].sum()
    agg = agg.rename(columns={"PhotonCount": "Count"})
    return agg


def merge_thread_event_csvs(config_dir: Path) -> pd.DataFrame:
    """
    合并同一配置目录下所有线程 CSV：AnaEx01_nt_PhotonFaceBlockEvent_t*.csv。
    若没有 _t* 文件则尝试无后缀的 AnaEx01_nt_PhotonFaceBlockEvent.csv。
    返回合并后的事件表（所有行纵向拼接）。
    """
    thread_files = sorted(config_dir.glob("AnaEx01_nt_PhotonFaceBlockEvent_t*.csv"))
    frames = []
    for tf in thread_files:
        try:
            frames.append(parse_event_csv(tf))
        except ValueError:
            continue

    if frames:
        return pd.concat(frames, ignore_index=True)

    fallback = config_dir / "AnaEx01_nt_PhotonFaceBlockEvent.csv"
    if fallback.exists():
        return parse_event_csv(fallback)
    raise FileNotFoundError(f"No usable PhotonFaceBlockEvent CSV in {config_dir}")


def build_face_matrix(agg: pd.DataFrame, face: int) -> np.ndarray:
    """
    从聚合表 agg（列 Face, j, k, Count）中抽取指定 face，填成 4×4 矩阵。
    j、k 为 0..3，矩阵[i,j] 对应 j 行 k 列（与热力图 imshow 一致）。
    """
    matrix = np.zeros((4, 4), dtype=float)
    face_df = agg[agg["Face"] == face]#根据Face列的值等于face的行，返回一个DataFrame
    for _, row in face_df.iterrows():#遍历face_df的每一行
        j, k = int(row["j"]), int(row["k"])
        if 0 <= j < 4 and 0 <= k < 4:
            matrix[j, k] = float(row["Count"])
    return matrix


def compute_reconstructed_position(face_totals: dict[int, float]) -> tuple[float, float, float]:
    """
    算法 0：最原始的“对面差分比值”重建。

    物理直觉：
    - 如果源更靠近 +X 面，那么 +X 面总光子数通常会比 -X 面更大；
    - 因此可用 N(+X)-N(-X) 表示“偏向哪一侧”；
    - 再除以 N(+X)+N(-X) 做归一化，去掉总光子数规模的影响。

    用六面总光子数重建归一化坐标：
      x = (N(+X)-N(-X)) / (N(+X)+N(-X))
      y = (N(+Y)-N(-Y)) / (N(+Y)+N(-Y))
      z = (N(+Z)-N(-Z)) / (N(+Z)+N(-Z))

    注意：
    - 这个算法更像“方向指标”，不是严格的 cm 坐标；
    - 它只使用六个面的总和，没有利用每个面内部 4x4 分布。
    """
    n_px, n_nx = face_totals.get(0, 0.0), face_totals.get(1, 0.0)
    n_py, n_ny = face_totals.get(2, 0.0), face_totals.get(3, 0.0)
    n_pz, n_nz = face_totals.get(4, 0.0), face_totals.get(5, 0.0)

    x_denom = n_px + n_nx
    y_denom = n_py + n_ny
    z_denom = n_pz + n_nz

    x_rec = (n_px - n_nx) / x_denom if x_denom > 0 else 0.0
    y_rec = (n_py - n_ny) / y_denom if y_denom > 0 else 0.0
    z_rec = (n_pz - n_nz) / z_denom if z_denom > 0 else 0.0
    return x_rec, y_rec, z_rec

def read_half_lengths_cm(metadata_path: Optional[Path]) -> dict[str, float]:
    """
    从 metadata.csv 读取晶体尺寸并转换为半长（cm）。
    若 metadata 不存在或读取失败，则回退到默认 25 mm 立方体。

    为什么这里要算“半长”：
    - 后面的重建算法通常都把输出约束在 [-half_length, +half_length]；
    - 例如 25 mm = 2.5 cm，全长的一半就是 1.25 cm；
    - 这样算法返回值可以直接解释为“晶体内的物理坐标（cm）”。
    """
    if metadata_path is None or not metadata_path.exists():
        return DEFAULT_HALF_LENGTHS_CM.copy()

    try:
        meta_df = pd.read_csv(metadata_path)
        meta = dict(zip(meta_df["Key"].astype(str), meta_df["Value"].astype(str)))
        size_xz_mm = float(meta.get("detector.crystalSize_mm", 25.0))
        size_y_mm = float(meta.get("detector.crystalSizeY_mm", size_xz_mm))
        return {
            "x": size_xz_mm / 20.0,
            "y": size_y_mm / 20.0,
            "z": size_xz_mm / 20.0,
        }
    except Exception:
        return DEFAULT_HALF_LENGTHS_CM.copy()


def reconstruct_linear_scaled(
    face_totals: dict[int, float],
    half_lengths_cm: dict[str, float],
) -> tuple[float, float, float]:
    """
    算法 1：linear_scaled

    每一轴独立：用一对对面的总光子数 N 做差分，再映射到 [-half, +half] cm。
    单轴公式（以 x 为例，+X/-X）：
      x_rec = half_x * (N(+X) - N(-X)) / (N(+X) + N(-X))
    分母为 0 表示该轴无信息，返回 0。
    """
    # +X / -X
    n_px = max(face_totals.get(0, 0.0), 0.0)
    n_nx = max(face_totals.get(1, 0.0), 0.0)
    hx = half_lengths_cm["x"]
    denom_x = n_px + n_nx
    x_rec = hx * (n_px - n_nx) / denom_x if denom_x > 0 else 0.0

    # +Y / -Y
    n_py = max(face_totals.get(2, 0.0), 0.0)
    n_ny = max(face_totals.get(3, 0.0), 0.0)
    hy = half_lengths_cm["y"]
    denom_y = n_py + n_ny
    y_rec = hy * (n_py - n_ny) / denom_y if denom_y > 0 else 0.0

    # +Z / -Z
    n_pz = max(face_totals.get(4, 0.0), 0.0)
    n_nz = max(face_totals.get(5, 0.0), 0.0)
    hz = half_lengths_cm["z"]
    denom_z = n_pz + n_nz
    z_rec = hz * (n_pz - n_nz) / denom_z if denom_z > 0 else 0.0

    return x_rec, y_rec, z_rec


def reconstruct_sqrt_scaled(
    face_totals: dict[int, float],
    half_lengths_cm: dict[str, float],
) -> tuple[float, float, float]:
    """
    算法 2：sqrt_scaled

    与 linear_scaled 相同结构，但对每边计数先取 sqrt，再差分：
      x_rec = half_x * (sqrt(N(+X)) - sqrt(N(-X))) / (sqrt(N(+X)) + sqrt(N(-X)))
    这样大计数被压扁，极端亮面不会完全主导比值。
    """
    n_px = max(face_totals.get(0, 0.0), 0.0)
    n_nx = max(face_totals.get(1, 0.0), 0.0)
    hx = half_lengths_cm["x"]
    spx, snx = np.sqrt(n_px), np.sqrt(n_nx)
    denom_x = spx + snx
    x_rec = hx * (spx - snx) / denom_x if denom_x > 0 else 0.0

    n_py = max(face_totals.get(2, 0.0), 0.0)
    n_ny = max(face_totals.get(3, 0.0), 0.0)
    hy = half_lengths_cm["y"]
    spy, sny = np.sqrt(n_py), np.sqrt(n_ny)
    denom_y = spy + sny
    y_rec = hy * (spy - sny) / denom_y if denom_y > 0 else 0.0

    n_pz = max(face_totals.get(4, 0.0), 0.0)
    n_nz = max(face_totals.get(5, 0.0), 0.0)
    hz = half_lengths_cm["z"]
    spz, snz = np.sqrt(n_pz), np.sqrt(n_nz)
    denom_z = spz + snz
    z_rec = hz * (spz - snz) / denom_z if denom_z > 0 else 0.0

    return x_rec, y_rec, z_rec


def reconstruct_log_scaled(
    face_totals: dict[int, float],
    half_lengths_cm: dict[str, float],
) -> tuple[float, float, float]:
    """
    算法 3：log_scaled

    与 linear_scaled 相同结构，但对每边计数用 log(1+N)：
      x_rec = half_x * (log(1+N(+X)) - log(1+N(-X))) / (log(1+N(+X)) + log(1+N(-X)))
    对超大计数压缩更强，结果通常比 sqrt 更靠中心。
    """
    n_px = max(face_totals.get(0, 0.0), 0.0)
    n_nx = max(face_totals.get(1, 0.0), 0.0)
    hx = half_lengths_cm["x"]
    lpx, lnx = np.log1p(n_px), np.log1p(n_nx)
    denom_x = lpx + lnx
    x_rec = hx * (lpx - lnx) / denom_x if denom_x > 0 else 0.0

    n_py = max(face_totals.get(2, 0.0), 0.0)
    n_ny = max(face_totals.get(3, 0.0), 0.0)
    hy = half_lengths_cm["y"]
    lpy, lny = np.log1p(n_py), np.log1p(n_ny)
    denom_y = lpy + lny
    y_rec = hy * (lpy - lny) / denom_y if denom_y > 0 else 0.0

    n_pz = max(face_totals.get(4, 0.0), 0.0)
    n_nz = max(face_totals.get(5, 0.0), 0.0)
    hz = half_lengths_cm["z"]
    lpz, lnz = np.log1p(n_pz), np.log1p(n_nz)
    denom_z = lpz + lnz
    z_rec = hz * (lpz - lnz) / denom_z if denom_z > 0 else 0.0

    return x_rec, y_rec, z_rec


def bin_centers_cm(half_length_cm: float) -> np.ndarray:
    """
    4 个等宽 bin 的中心位置，范围 [-half_length, +half_length]，单位 cm。

    因为每个面只有 4x4 个 SiPM block，所以面内坐标不是连续的，
    这里只能把每一格近似看成一个“离散采样点中心”。
    """
    return np.linspace(-0.75 * half_length_cm, 0.75 * half_length_cm, 4)


def face_centroid(
    matrix: np.ndarray,
    half_j_cm: float,
    half_k_cm: float,
) -> tuple[float, float] | None:
    """
    对单个 4×4 面矩阵计算质心，返回 (j_axis_cm, k_axis_cm)。

    思想和图像质心/重心完全一样：
    - 哪些格子计数高，就认为哪里“更亮”；
    - 用计数当权重，求出该面内的加权平均位置。

    这一步开始真正使用每个面内部的 4x4 分布信息，
    而不只是把整个面的计数加总。
    """
    total = float(np.sum(matrix))
    if total <= 0:
        return None

    j_centers = bin_centers_cm(half_j_cm)
    k_centers = bin_centers_cm(half_k_cm)
    # 先把 4x4 矩阵沿列/行压成一维分布，再做加权平均位置。
    j_coord = float(np.dot(np.sum(matrix, axis=1), j_centers) / total)
    k_coord = float(np.dot(np.sum(matrix, axis=0), k_centers) / total)
    return j_coord, k_coord


def reconstruct_centroid_mean(
    matrices: dict[int, np.ndarray],
    face_totals: dict[int, float],
    half_lengths_cm: dict[str, float],
) -> tuple[float, float, float]:
    """
    利用每个面的 4×4 分布质心重建位置。
    各轴用所有提供该轴信息的面做简单平均。

    对应关系：
    - +X/-X 面上看到的是 y-z 平面，因此它们给 y、z 提供信息；
    - +Y/-Y 面上看到的是 x-z 平面，因此它们给 x、z 提供信息；
    - +Z/-Z 面上看到的是 x-y 平面，因此它们给 x、y 提供信息。

    例如 x 轴的估计，不来自 +X/-X 面，而来自：
    - +Y 面上的 x 质心
    - -Y 面上的 x 质心
    - +Z 面上的 x 质心
    - -Z 面上的 x 质心

    然后把这些值简单平均。
    """
    axis_terms = {"x": [], "y": [], "z": []}

    for face, matrix in matrices.items():
        if face_totals.get(face, 0.0) <= 0:
            continue
        if face in (0, 1):
            centroid = face_centroid(matrix, half_lengths_cm["y"], half_lengths_cm["z"])
            if centroid is not None:
                axis_terms["y"].append(centroid[0])
                axis_terms["z"].append(centroid[1])
        elif face in (2, 3):
            centroid = face_centroid(matrix, half_lengths_cm["x"], half_lengths_cm["z"])
            if centroid is not None:
                axis_terms["x"].append(centroid[0])
                axis_terms["z"].append(centroid[1])
        elif face in (4, 5):
            centroid = face_centroid(matrix, half_lengths_cm["x"], half_lengths_cm["y"])
            if centroid is not None:
                axis_terms["x"].append(centroid[0])
                axis_terms["y"].append(centroid[1])

    def mean_or_zero(values: list[float]) -> float:
        return float(np.mean(values)) if values else 0.0

    return (
        mean_or_zero(axis_terms["x"]),
        mean_or_zero(axis_terms["y"]),
        mean_or_zero(axis_terms["z"]),
    )


def reconstruct_centroid_weighted(
    matrices: dict[int, np.ndarray],
    face_totals: dict[int, float],
    half_lengths_cm: dict[str, float],
) -> tuple[float, float, float]:
    """
    利用每个面的 4×4 分布质心重建位置。
    各轴按对应面的总光子数加权平均。

    和 centroid_mean 的区别只有一点：
    - centroid_mean：所有面一票同权
    - centroid_weighted：更亮的面权重更大

    这样做的直觉是：信号更强的面，质心通常更可信。
    """
    axis_terms = {"x": [], "y": [], "z": []}

    for face, matrix in matrices.items():
        weight = float(face_totals.get(face, 0.0))
        if weight <= 0:
            continue
        if face in (0, 1):
            centroid = face_centroid(matrix, half_lengths_cm["y"], half_lengths_cm["z"])
            if centroid is not None:
                axis_terms["y"].append((weight, centroid[0]))
                axis_terms["z"].append((weight, centroid[1]))
        elif face in (2, 3):
            centroid = face_centroid(matrix, half_lengths_cm["x"], half_lengths_cm["z"])
            if centroid is not None:
                axis_terms["x"].append((weight, centroid[0]))
                axis_terms["z"].append((weight, centroid[1]))
        elif face in (4, 5):
            centroid = face_centroid(matrix, half_lengths_cm["x"], half_lengths_cm["y"])
            if centroid is not None:
                axis_terms["x"].append((weight, centroid[0]))
                axis_terms["y"].append((weight, centroid[1]))

    def weighted_mean_or_zero(items: list[tuple[float, float]]) -> float:
        if not items:
            return 0.0
        weights = np.array([w for w, _ in items], dtype=float)
        values = np.array([v for _, v in items], dtype=float)
        return float(np.sum(weights * values) / np.sum(weights))

    return (
        weighted_mean_or_zero(axis_terms["x"]),
        weighted_mean_or_zero(axis_terms["y"]),
        weighted_mean_or_zero(axis_terms["z"]),
    )


def reconstruct_hybrid_sqrt_centroid(
    face_totals: dict[int, float],
    matrices: dict[int, np.ndarray],
    half_lengths_cm: dict[str, float],
) -> tuple[float, float, float]:
    """
    混合算法：sqrt 对面差分 + centroid weighted 平均，各占 50%。

    设计动机：
    - sqrt_scaled 更擅长反映“更靠近哪一边”；
    - centroid_weighted 更擅长利用面内分布提供横向位置细节；
    - 所以先用最简单的 50/50 平均把两者混合。

    这不是严格最优权重，只是一个容易理解、容易实验的起点。
    """
    x_sqrt, y_sqrt, z_sqrt = reconstruct_sqrt_scaled(face_totals, half_lengths_cm)
    x_cent, y_cent, z_cent = reconstruct_centroid_weighted(matrices, face_totals, half_lengths_cm)
    return (
        0.5 * (x_sqrt + x_cent),
        0.5 * (y_sqrt + y_cent),
        0.5 * (z_sqrt + z_cent),
    )


def compute_reconstruction_rows(
    face_totals: dict[int, float],
    matrices: dict[int, np.ndarray],
    half_lengths_cm: dict[str, float],
) -> list[dict[str, float | str]]:
    """
    统一收集所有可比较的重建算法输出。
    返回值可直接写入 reconstructed_position.csv。

    这里可以理解成“算法注册表”：
    - 每个 tuple 的第一个元素是算法名字；
    - 第二个元素是该算法算出来的 (x_rec, y_rec, z_rec)。

    以后如果你还想加新算法，最省事的方式就是：
    1. 先新写一个 reconstruct_xxx(...) 函数；
    2. 再把它塞进下面这个 algorithms 列表；
    3. 后续 CSV 和分析脚本都会自动识别到它。
    """
    rows = []
    algorithms = [
        ("linear_scaled", reconstruct_linear_scaled(face_totals, half_lengths_cm)),
        ("sqrt_scaled", reconstruct_sqrt_scaled(face_totals, half_lengths_cm)),
        ("log_scaled", reconstruct_log_scaled(face_totals, half_lengths_cm)),
        ("half_side_ratio", reconstruct_half_side_ratio(matrices, half_lengths_cm)),
        ("centroid_mean", reconstruct_centroid_mean(matrices, face_totals, half_lengths_cm)),
        ("centroid_weighted", reconstruct_centroid_weighted(matrices, face_totals, half_lengths_cm)),
        ("hybrid_sqrt_centroid", reconstruct_hybrid_sqrt_centroid(face_totals, matrices, half_lengths_cm)),
    ]
    for algorithm, (x_rec, y_rec, z_rec) in algorithms:
        rows.append(
            {
                "algorithm": algorithm,
                "x_rec": x_rec,
                "y_rec": y_rec,
                "z_rec": z_rec,
            }
        )
    return rows



def compute_face_totals(agg: pd.DataFrame) -> dict[int, float]:
    """
    从聚合表（Face, j, k, Count）计算六面总光子数。
    返回字典键为 Face 0..5，缺失面自动补 0。
    """
    totals_series = agg.groupby("Face")["Count"].sum()
    return {face: float(totals_series.get(face, 0.0)) for face in range(6)}


def plot_6_faces(
    agg: pd.DataFrame,
    output_path: Path,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    title: str,
    reference_text: str,
) -> None:
    """
    根据聚合表 agg 绘制 6 面热力图并保存。
    布局 2×3，每子图为 4×4 矩阵，共享 colorbar；格子内标注整数计数，
    根据背景亮度自动选白/黑文字。vmin/vmax 为 None 时按数据范围自动。
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    matrices = [build_face_matrix(agg, face) for face in range(6)]
    face_totals = compute_face_totals(agg)
    x_rec, y_rec, z_rec = compute_reconstructed_position(face_totals)
    stacked = np.stack(matrices)
    draw_vmin = float(np.min(stacked)) if vmin is None else vmin
    draw_vmax = float(np.max(stacked)) if vmax is None else vmax

    last_im = None
    for face, ax in enumerate(axes.flat):
        mat = matrices[face]
        last_im = ax.imshow(mat, cmap=cmap, vmin=draw_vmin, vmax=draw_vmax, origin="lower")
        ax.set_title(
            f"Face {face} ({FACE_NAMES.get(face, 'Unknown')}, normal) | Total={int(face_totals[face])}"
        )
        x_label, y_label = FACE_AXIS_LABELS.get(face, ("k", "j"))
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_xticks(range(4))
        ax.set_yticks(range(4))
        for j in range(4):
            for k in range(4):
                val = int(mat[j, k])
                color = "white" if val > (draw_vmin + draw_vmax) / 2 else "black"
                ax.text(k, j, f"{val}", ha="center", va="center", fontsize=8, color=color)

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes.ravel().tolist(), shrink=0.9)
        cbar.set_label("Photon Count")

    fig.text(
        0.5,
        0.01,
        reference_text,
        ha="center",
        va="bottom",
        fontsize=10,
    )
    fig.suptitle(title)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def process_one_config(
    config_dir: Path,
    output_dir: Path,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
) -> None:
    """
    处理单个配置目录：合并线程 CSV → 写 merged_event.csv、merged_face_jk.csv → 画热力图。
    config_dir：Results 下某一配置子目录；output_dir：Output 下同名子目录。
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_src = config_dir / "metadata.csv"
    metadata_dst = output_dir / "metadata.csv"
    metadata_copied = False
    if metadata_src.exists() and metadata_src.is_file():
        shutil.copy2(metadata_src, metadata_dst)
        metadata_copied = True
    merged_event = merge_thread_event_csvs(config_dir)
    merged_event_path = output_dir / "merged_event.csv"
    merged_event.to_csv(merged_event_path, index=False)

    agg = aggregate_for_heatmap(merged_event)
    merged_agg_path = output_dir / "merged_face_jk.csv"
    agg.to_csv(merged_agg_path, index=False)

    # 先准备所有算法都会用到的公共输入：
    # 1. half_lengths_cm：晶体半长（决定最终坐标量纲是 cm）
    # 2. face_totals：六个面的总光子数（给差分类算法使用）
    # 3. matrices：六个面的 4x4 分布矩阵（给质心类算法使用）
    half_lengths_cm = read_half_lengths_cm(metadata_src if metadata_copied else None)
    face_totals = compute_face_totals(agg)
    matrices = {face: build_face_matrix(agg, face) for face in range(6)}
    reconstructed_path = output_dir / "reconstructed_position.csv"
    reconstructed_df = pd.DataFrame(
        compute_reconstruction_rows(face_totals, matrices, half_lengths_cm)
    )
    reconstructed_df.to_csv(reconstructed_path, index=False)
    reference_row = reconstructed_df[reconstructed_df["algorithm"] == "hybrid_sqrt_centroid"]
    if reference_row.empty:
        reference_row = reconstructed_df.iloc[[0]]
    reference_row = reference_row.iloc[0]
    reference_text = (
        "Reference reconstruction (cm): "
        f"{reference_row['algorithm']} "
        f"(x,y,z)=({reference_row['x_rec']:.4f}, {reference_row['y_rec']:.4f}, {reference_row['z_rec']:.4f}) "
        "| full list in reconstructed_position.csv"
    )

    heatmap_path = output_dir / "sipm_6faces_heatmap.png"
    plot_6_faces(
        agg=agg,
        output_path=heatmap_path,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        title=f"6-Face SiPM Photon Heatmap - {config_dir.name}",
        reference_text=reference_text,
    )
    print(f"[OK] {config_dir.name}")
    print(f"     merged_event.csv:   {merged_event_path}")
    print(f"     merged_face_jk.csv: {merged_agg_path}")
    print(f"     reconstructed.csv:  {reconstructed_path}")
    if metadata_copied:
        print(f"     metadata.csv:       {metadata_dst}")
    print(f"     heatmap:            {heatmap_path}")


def main() -> None:
    """解析命令行，遍历 Results 下所有配置子目录并逐个调用 process_one_config。"""
    parser = argparse.ArgumentParser(
        description="Histo10_Cubic: 批量合并六面 SiPM 事件 CSV 并生成 6 面热力图",
    )
    parser.add_argument("results", nargs="?", help="Results 根目录路径")
    parser.add_argument("--results", dest="results_alt", help="Results 根目录路径（备用）")
    parser.add_argument("--output", default=None, help="输出根目录；默认 Output/SiPM6_Output_<时间戳>")
    parser.add_argument("--cmap", default="viridis", help="热力图的 Matplotlib 色图名")
    parser.add_argument("--vmin", type=float, default=None, help="热力图颜色最小值（不指定则自动）")
    parser.add_argument("--vmax", type=float, default=None, help="热力图颜色最大值（不指定则自动）")
    args = parser.parse_args()

    results_path = args.results or args.results_alt
    if not results_path:
        raise ValueError("请提供 Results 路径（位置参数或 --results）")

    results_dir = Path(results_path).resolve()
    if not results_dir.exists() or not results_dir.is_dir():
        raise NotADirectoryError(f"Results 目录不存在: {results_dir}")

    if args.output:
        output_root = Path(args.output).resolve()
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_root = Path.cwd() / "Output" / f"SiPM6_Output_{ts}"
    output_root.mkdir(parents=True, exist_ok=True)
    print(f"Output root: {output_root}")

    config_dirs = sorted([p for p in results_dir.iterdir() if p.is_dir()], key=lambda p: p.name)
    if not config_dirs:
        raise ValueError(f"在 {results_dir} 中未找到配置子目录")

    ok = 0
    for cfg in config_dirs:
        try:
            process_one_config(
                config_dir=cfg,
                output_dir=output_root / cfg.name,
                cmap=args.cmap,
                vmin=args.vmin,
                vmax=args.vmax,
            )
            ok += 1
        except Exception as exc:
            print(f"[SKIP] {cfg.name}: {exc}")

    print(f"[DONE] 已处理 {ok}/{len(config_dirs)} 个配置目录。")


if __name__ == "__main__":
    main()
