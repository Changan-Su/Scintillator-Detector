#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_position_accuracy.py
============================
输入：一个包含多个单次测量子目录的总文件夹（如 Output/SiPM6_Output_20260314_220905）。
      每个子目录需含：
        metadata.csv             — 包含 source.fp_source_x/y/z_mm 等 Key,Value 行
        reconstructed_position.csv — 含 algorithm, x_rec, y_rec, z_rec 列

单位说明：
  metadata 中 source.fp_source_*_mm 为 mm，读入时统一转换为 cm；
  reconstructed_position.csv 中 rec 位置为 cm；bias/dx/dy/dz/d 单位均为 cm。

输出（总文件夹下）：
  accuracy_summary.csv       每行 (folder, algorithm, true_xyz(cm), rec_xyz(cm), dx/dy/dz/d(cm))
                             末尾按算法分组的 BIAS / STD（单位 cm）
  scatter_<algorithm>.png    每个算法一张三视图散点图（XY / XZ / YZ 面）

使用方式：
  python analyze_position_accuracy.py <总文件夹路径>
  python analyze_position_accuracy.py Output/SiPM6_Output_20260314_220905
  python analyze_position_accuracy.py Output/SiPM6_Output_20260314_220905 --output my_summary.csv
  python analyze_position_accuracy.py Output/... --no-plot   # 跳过散点图生成
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


METADATA_FILE = "metadata.csv"
RECON_FILE = "reconstructed_position.csv"
META_KEY_X = "source.fp_source_x_mm"
META_KEY_Y = "source.fp_source_y_mm"
META_KEY_Z = "source.fp_source_z_mm"


def read_metadata(meta_path: Path) -> Dict[str, str]:
    df = pd.read_csv(meta_path)
    if not {"Key", "Value"}.issubset(df.columns):
        raise ValueError(f"metadata.csv must have Key,Value columns: {meta_path}")
    return dict(zip(df["Key"].astype(str), df["Value"].astype(str)))


def read_recon(recon_path: Path) -> pd.DataFrame:
    df = pd.read_csv(recon_path)
    required = {"algorithm", "x_rec", "y_rec", "z_rec"}
    if not required.issubset(df.columns):
        raise ValueError(f"reconstructed_position.csv missing columns {required - set(df.columns)}: {recon_path}")
    for col in ("x_rec", "y_rec", "z_rec"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def process_subdir(subdir: Path) -> Optional[List[dict]]:
    meta_path = subdir / METADATA_FILE
    recon_path = subdir / RECON_FILE

    if not meta_path.exists():
        print(f"  [SKIP] no {METADATA_FILE}: {subdir.name}")
        return None
    if not recon_path.exists():
        print(f"  [SKIP] no {RECON_FILE}: {subdir.name}")
        return None

    try:
        meta = read_metadata(meta_path)
    except Exception as e:
        print(f"  [SKIP] bad metadata: {subdir.name} — {e}")
        return None

    for key in (META_KEY_X, META_KEY_Y, META_KEY_Z):
        if key not in meta:
            print(f"  [SKIP] missing key '{key}' in metadata: {subdir.name}")
            return None

    try:
        # metadata stores positions in mm; convert to cm to match rec positions
        true_x = float(meta[META_KEY_X]) / 10.0
        true_y = float(meta[META_KEY_Y]) / 10.0
        true_z = float(meta[META_KEY_Z]) / 10.0
    except ValueError as e:
        print(f"  [SKIP] non-numeric true position: {subdir.name} — {e}")
        return None

    try:
        recon_df = read_recon(recon_path)
    except Exception as e:
        print(f"  [SKIP] bad reconstructed_position: {subdir.name} — {e}")
        return None

    rows = []
    for _, row in recon_df.iterrows():
        algorithm = str(row["algorithm"])
        rx, ry, rz = float(row["x_rec"]), float(row["y_rec"]), float(row["z_rec"])
        dx = abs(rx - true_x)
        dy = abs(ry - true_y)
        dz = abs(rz - true_z)
        d = math.sqrt(dx**2 + dy**2 + dz**2)
        rows.append(
            {
                "folder": subdir.name,
                "algorithm": algorithm,
                "true_x": true_x,
                "true_y": true_y,
                "true_z": true_z,
                "rec_x": rx,
                "rec_y": ry,
                "rec_z": rz,
                "dx": dx,
                "dy": dy,
                "dz": dz,
                "d": d,
            }
        )
    return rows


def compute_stats(rows: List[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    stat_rows = []
    for algo, group in df.groupby("algorithm", sort=False):
        bias_x = group["dx"].mean()
        bias_y = group["dy"].mean()
        bias_z = group["dz"].mean()
        bias_d = group["d"].mean()
        std_x = group["dx"].std(ddof=1) if len(group) > 1 else 0.0
        std_y = group["dy"].std(ddof=1) if len(group) > 1 else 0.0
        std_z = group["dz"].std(ddof=1) if len(group) > 1 else 0.0
        std_d = group["d"].std(ddof=1) if len(group) > 1 else 0.0
        stat_rows.append(
            {
                "folder": f"[SUMMARY algorithm={algo}]",
                "algorithm": algo,
                "true_x": "",
                "true_y": "",
                "true_z": "",
                "rec_x": "",
                "rec_y": "",
                "rec_z": "",
                "dx": f"BIAS={bias_x:.6f}",
                "dy": f"BIAS={bias_y:.6f}",
                "dz": f"BIAS={bias_z:.6f}",
                "d": f"BIAS={bias_d:.6f}",
            }
        )
        stat_rows.append(
            {
                "folder": f"[SUMMARY algorithm={algo}]",
                "algorithm": algo,
                "true_x": "",
                "true_y": "",
                "true_z": "",
                "rec_x": "",
                "rec_y": "",
                "rec_z": "",
                "dx": f"STD={std_x:.6f}",
                "dy": f"STD={std_y:.6f}",
                "dz": f"STD={std_z:.6f}",
                "d": f"STD={std_d:.6f}",
            }
        )
    return pd.DataFrame(stat_rows)


def plot_scatter_views(data_df: pd.DataFrame, output_dir: Path) -> None:
    """
    For each algorithm, produce a grid figure:
      - 3 columns: XY plane, XZ plane, YZ plane
      - N rows: one row per depth slice (sorted descending), where depth is the
        axis perpendicular to that column's view plane
          XY column → depth = true Z
          XZ column → depth = true Y
          YZ column → depth = true X

    Each subplot shows only the points belonging to that depth slice:
      filled circles = reconstructed positions, crosses = true positions.
    """
    # (view name, h-axis, v-axis, rec_h col, rec_v col,
    #  true_h col, true_v col, depth true-col, depth axis label)
    views = [
        ("XY", "x", "y", "rec_x", "rec_y", "true_x", "true_y", "true_z", "Z"),
        ("XZ", "x", "z", "rec_x", "rec_z", "true_x", "true_z", "true_y", "Y"),
        ("YZ", "y", "z", "rec_y", "rec_z", "true_y", "true_z", "true_x", "X"),
    ]

    for algo, group in data_df.groupby("algorithm", sort=False):
        # Depth levels per row (one row per view), sorted descending (left = largest depth)
        depth_levels = [
            sorted(group[depth_col].dropna().unique(), reverse=True)
            for *_, depth_col, _ in views
        ]
        n_rows = len(views)
        n_cols = max(len(d) for d in depth_levels)

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(4 * n_cols, 4 * n_rows),
            constrained_layout=True,
            squeeze=False,
        )
        fig.suptitle(f"Reconstruction scatter — algorithm: {algo}", fontsize=13)

        for row_idx, (view_name, hlabel, vlabel, rec_h, rec_v, true_h, true_v, depth_col, depth_label) in enumerate(views):
            depths = depth_levels[row_idx]

            for col_idx in range(n_cols):
                ax = axes[row_idx, col_idx]

                if col_idx >= len(depths):
                    ax.set_visible(False)
                    continue

                d = depths[col_idx]
                sl = group[group[depth_col] == d]

                ax.axhline(0, color="gray", linewidth=0.6, linestyle="--")
                ax.axvline(0, color="gray", linewidth=0.6, linestyle="--")
                ax.scatter(
                    sl[rec_h], sl[rec_v],
                    s=18, alpha=0.75, color="steelblue", label="Reconstructed", zorder=3,
                )
                ax.scatter(
                    sl[true_h], sl[true_v],
                    s=50, marker="+", linewidths=1.5, color="crimson", label="True", zorder=4,
                )

                ax.set_title(f"{view_name}  |  {depth_label} = {d:.4g} cm", fontsize=10)
                ax.set_xlabel(f"{hlabel.upper()} (cm)", fontsize=9)
                ax.set_ylabel(f"{vlabel.upper()} (cm)", fontsize=9)
                ax.set_aspect("equal", adjustable="datalim")
                ax.grid(True, linewidth=0.4, alpha=0.5)
                if col_idx == 0:
                    ax.legend(fontsize=8, loc="upper right")

        out_path = output_dir / f"scatter_{algo}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"       scatter plot: {out_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize reconstruction accuracy across multiple measurement sub-directories."
    )
    parser.add_argument("input_dir", help="Top-level folder containing measurement sub-directories.")
    parser.add_argument("--output", default=None, help="Output CSV filename (default: accuracy_summary.csv inside input_dir).")
    parser.add_argument("--no-plot", action="store_true", help="Skip scatter plot generation.")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        print(f"[ERROR] Not a directory: {input_dir}")
        return 1

    output_path = Path(args.output).resolve() if args.output else input_dir / "accuracy_summary.csv"

    subdirs = sorted([p for p in input_dir.iterdir() if p.is_dir()], key=lambda p: p.name)
    if not subdirs:
        print(f"[ERROR] No sub-directories found in {input_dir}")
        return 1

    print(f"Input:  {input_dir}")
    print(f"Output: {output_path}")
    print(f"Found {len(subdirs)} sub-directories.")

    all_rows: List[dict] = []
    ok = 0
    for subdir in subdirs:
        rows = process_subdir(subdir)
        if rows:
            all_rows.extend(rows)
            ok += 1

    if not all_rows:
        print("[ERROR] No valid data found in any sub-directory.")
        return 1

    data_df = pd.DataFrame(all_rows)
    stat_df = compute_stats(all_rows)
    combined = pd.concat([data_df, stat_df], ignore_index=True)
    combined.to_csv(output_path, index=False, float_format="%.6f")

    print(f"\n[DONE] Processed {ok}/{len(subdirs)} sub-directories.")
    print(f"       CSV:    {output_path}")

    if not args.no_plot:
        plot_scatter_views(data_df, input_dir)

    # Print summary to terminal
    for algo, group in data_df.groupby("algorithm", sort=False):
        n = len(group)
        print(f"\n  Algorithm: {algo}  (n={n})")
        print(f"    bias_x={group['dx'].mean():.4f}  std_x={group['dx'].std(ddof=1) if n>1 else 0:.4f}")
        print(f"    bias_y={group['dy'].mean():.4f}  std_y={group['dy'].std(ddof=1) if n>1 else 0:.4f}")
        print(f"    bias_z={group['dz'].mean():.4f}  std_z={group['dz'].std(ddof=1) if n>1 else 0:.4f}")
        print(f"    bias_d={group['d'].mean():.4f}   std_d={group['d'].std(ddof=1) if n>1 else 0:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
