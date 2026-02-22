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
  sipm_6faces_heatmap.png  6 面 4×4 热力图（Face 0..5：+X,-X,+Y,-Y,+Z,-Z）

示例
----
  python Histo10_Cubic.py Results
  python Histo10_Cubic.py Results --output ./Output/SiPM6_Out --cmap plasma
"""

import argparse
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 六面编号与几何方向对应（与 DetectorConstruction 中 face 定义一致）
FACE_NAMES = {0: "+X", 1: "-X", 2: "+Y", 3: "-Y", 4: "+Z", 5: "-Z"}

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
    df = df.iloc[:, :9].copy()
    df.columns = EVENT_COLUMNS
    for col in EVENT_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["Face", "j", "k", "PhotonCount"]).copy()
    df["Face"] = df["Face"].astype(int)
    df["j"] = df["j"].astype(int)
    df["k"] = df["k"].astype(int)
    return df


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
    face_df = agg[agg["Face"] == face]
    for _, row in face_df.iterrows():
        j, k = int(row["j"]), int(row["k"])
        if 0 <= j < 4 and 0 <= k < 4:
            matrix[j, k] = float(row["Count"])
    return matrix


def plot_6_faces(
    agg: pd.DataFrame,
    output_path: Path,
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    title: str,
) -> None:
    """
    根据聚合表 agg 绘制 6 面热力图并保存。
    布局 2×3，每子图为 4×4 矩阵，共享 colorbar；格子内标注整数计数，
    根据背景亮度自动选白/黑文字。vmin/vmax 为 None 时按数据范围自动。
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    matrices = [build_face_matrix(agg, face) for face in range(6)]
    stacked = np.stack(matrices)
    draw_vmin = float(np.min(stacked)) if vmin is None else vmin
    draw_vmax = float(np.max(stacked)) if vmax is None else vmax

    last_im = None
    for face, ax in enumerate(axes.flat):
        mat = matrices[face]
        last_im = ax.imshow(mat, cmap=cmap, vmin=draw_vmin, vmax=draw_vmax, origin="lower")
        ax.set_title(f"Face {face} ({FACE_NAMES.get(face, 'Unknown')})")
        ax.set_xlabel("k")
        ax.set_ylabel("j")
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
    merged_event = merge_thread_event_csvs(config_dir)
    merged_event_path = output_dir / "merged_event.csv"
    merged_event.to_csv(merged_event_path, index=False)

    agg = aggregate_for_heatmap(merged_event)
    merged_agg_path = output_dir / "merged_face_jk.csv"
    agg.to_csv(merged_agg_path, index=False)

    heatmap_path = output_dir / "sipm_6faces_heatmap.png"
    plot_6_faces(
        agg=agg,
        output_path=heatmap_path,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        title=f"6-Face SiPM Photon Heatmap - {config_dir.name}",
    )
    print(f"[OK] {config_dir.name}")
    print(f"     merged_event.csv:   {merged_event_path}")
    print(f"     merged_face_jk.csv: {merged_agg_path}")
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
