# -*- coding: utf-8 -*-
"""
predict3.py —— 批量推理 + 汇总 + 散点图

相对 predict2.py 的增量
=======================
  1. `--input-dir` 接受**总文件夹**（如 Output/SiPM6_Output_20260329_032534），
     自动遍历所有含 merged_event.csv 的子目录；
     也兼容单个 config 目录（如果本身含 merged_event.csv）。
  2. 每个子目录仍然走 predict2 的对称镜像推理，把 `nn_mlp_sym` 追加到
     该子目录的 reconstructed_position.csv。
  3. 在总文件夹下输出：
       accuracy_summary_nn.csv    （和 flow2/analyze_position_accuracy.py 同款格式）
       scatter_nn_mlp_sym.png     （XY / XZ / YZ 三视图 × 深度切片）

用法：
    uv run python workflow/flow3_NN/predict3.py \
        --input-dir "Output/SiPM6_Output_20260329_032534" \
        --artifacts workflow/flow3_NN/artifacts_p5000_norm \
        --normalize-counts
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import torch

from build_dataset import events_to_tensor, read_true_position_cm
from model import PositionMLP
from predict2 import (
    canonicalize_to_positive_octant,
    detect_octant,
    uncanonicalize_output,
)

ALGO_NAME = "nn_mlp_sym"


# ---------------------------------------------------------------------------
# 单个 config 目录推理
# ---------------------------------------------------------------------------
def predict_one(
    cfg_dir: Path,
    model: torch.nn.Module,
    x_mean: np.ndarray,
    x_std: np.ndarray,
    y_scale: float,
    device: torch.device,
    args,
) -> Optional[dict]:
    event_csv = cfg_dir / "merged_event.csv"
    if not event_csv.exists():
        return None

    X = events_to_tensor(
        event_csv,
        min_photons=1,
        photons_per_sample=args.photons_per_sample,
        normalize_counts=args.normalize_counts,
    )
    if X is None or len(X) == 0:
        print(f"  [SKIP] no events: {cfg_dir.name}")
        return None

    signs = detect_octant(X)
    X_canon = canonicalize_to_positive_octant(X, signs)
    X_n = (X_canon - x_mean) / x_std

    preds = []
    with torch.no_grad():
        for i in range(0, len(X_n), args.batch_size):
            xb = torch.from_numpy(X_n[i : i + args.batch_size].astype(np.float32)).to(device)
            preds.append(model(xb).cpu().numpy())
    y_pred_canon = np.concatenate(preds, axis=0) * y_scale
    y_pred = uncanonicalize_output(y_pred_canon, signs)

    pos_mean = y_pred.mean(axis=0)

    # 事件级 CSV
    ev_df = pd.DataFrame({
        "event_idx": np.arange(len(y_pred)),
        "x_rec_cm": y_pred[:, 0],
        "y_rec_cm": y_pred[:, 1],
        "z_rec_cm": y_pred[:, 2],
        "total_photons": X.sum(axis=1),
        "sign_x": signs[:, 0],
        "sign_y": signs[:, 1],
        "sign_z": signs[:, 2],
    })
    ev_df.to_csv(cfg_dir / "nn_predicted_events_sym.csv", index=False)

    # 追加到该 config 的 reconstructed_position.csv
    if args.append_to_recon:
        rec_csv = cfg_dir / "reconstructed_position.csv"
        row = {"algorithm": ALGO_NAME,
               "x_rec": pos_mean[0], "y_rec": pos_mean[1], "z_rec": pos_mean[2]}
        if rec_csv.exists():
            df = pd.read_csv(rec_csv)
            df = df[df["algorithm"] != ALGO_NAME]
            df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
        else:
            df = pd.DataFrame([row])
        df.to_csv(rec_csv, index=False)

    # 真值
    true_pos = read_true_position_cm(cfg_dir / "metadata.csv")
    if true_pos is None:
        print(f"  [WARN] no metadata truth: {cfg_dir.name}")
        return None

    dx = abs(pos_mean[0] - true_pos[0])
    dy = abs(pos_mean[1] - true_pos[1])
    dz = abs(pos_mean[2] - true_pos[2])
    d = math.sqrt(dx * dx + dy * dy + dz * dz)

    print(f"  {cfg_dir.name}: true=({true_pos[0]:+.3f},{true_pos[1]:+.3f},{true_pos[2]:+.3f}) "
          f"pred=({pos_mean[0]:+.3f},{pos_mean[1]:+.3f},{pos_mean[2]:+.3f}) |err|={d:.4f}")

    return {
        "folder": cfg_dir.name,
        "algorithm": ALGO_NAME,
        "true_x": float(true_pos[0]),
        "true_y": float(true_pos[1]),
        "true_z": float(true_pos[2]),
        "rec_x": float(pos_mean[0]),
        "rec_y": float(pos_mean[1]),
        "rec_z": float(pos_mean[2]),
        "dx": dx, "dy": dy, "dz": dz, "d": d,
    }


# ---------------------------------------------------------------------------
# 汇总统计（和 flow2/analyze_position_accuracy.py 同款）
# ---------------------------------------------------------------------------
def compute_stats(rows: List[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    stat_rows = []
    for algo, group in df.groupby("algorithm", sort=False):
        bx, by, bz, bd = (group[c].mean() for c in ("dx", "dy", "dz", "d"))
        if len(group) > 1:
            sx, sy, sz, sd = (group[c].std(ddof=1) for c in ("dx", "dy", "dz", "d"))
        else:
            sx = sy = sz = sd = 0.0

        def _blank_row():
            return {"folder": f"[SUMMARY algorithm={algo}]", "algorithm": algo,
                    "true_x": "", "true_y": "", "true_z": "",
                    "rec_x": "", "rec_y": "", "rec_z": ""}

        bias = _blank_row()
        bias.update({"dx": f"BIAS={bx:.6f}", "dy": f"BIAS={by:.6f}",
                     "dz": f"BIAS={bz:.6f}", "d": f"BIAS={bd:.6f}"})
        stat_rows.append(bias)

        std = _blank_row()
        std.update({"dx": f"STD={sx:.6f}", "dy": f"STD={sy:.6f}",
                    "dz": f"STD={sz:.6f}", "d": f"STD={sd:.6f}"})
        stat_rows.append(std)
    return pd.DataFrame(stat_rows)


# ---------------------------------------------------------------------------
# 散点图（和 flow2/analyze_position_accuracy.py 同款）
# ---------------------------------------------------------------------------
def plot_scatter_views(data_df: pd.DataFrame, output_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    views = [
        ("XY", "x", "y", "rec_x", "rec_y", "true_x", "true_y", "true_z", "Z"),
        ("XZ", "x", "z", "rec_x", "rec_z", "true_x", "true_z", "true_y", "Y"),
        ("YZ", "y", "z", "rec_y", "rec_z", "true_y", "true_z", "true_x", "X"),
    ]

    for algo, group in data_df.groupby("algorithm", sort=False):
        depth_levels = [
            sorted(group[v[7]].dropna().unique(), reverse=True)
            for v in views
        ]
        n_rows = len(views)
        n_cols = max(len(d) for d in depth_levels)

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(4 * n_cols, 4 * n_rows),
            constrained_layout=True, squeeze=False,
        )
        fig.suptitle(f"Reconstruction scatter — algorithm: {algo}", fontsize=13)

        for row_idx, (view_name, hl, vl, rec_h, rec_v, true_h, true_v, depth_col, depth_label) in enumerate(views):
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
                ax.scatter(sl[rec_h], sl[rec_v], s=18, alpha=0.75,
                           color="steelblue", label="Reconstructed", zorder=3)
                ax.scatter(sl[true_h], sl[true_v], s=50, marker="+",
                           linewidths=1.5, color="crimson", label="True", zorder=4)

                ax.set_title(f"{view_name}  |  {depth_label} = {d:.4g} cm", fontsize=10)
                ax.set_xlabel(f"{hl.upper()} (cm)", fontsize=9)
                ax.set_ylabel(f"{vl.upper()} (cm)", fontsize=9)
                ax.set_aspect("equal", adjustable="datalim")
                ax.grid(True, linewidth=0.4, alpha=0.5)
                if col_idx == 0:
                    ax.legend(fontsize=8, loc="upper right")

        out_path = output_dir / f"scatter_{algo}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"scatter plot: {out_path}")


def plot_scatter_overlay(data_df: pd.DataFrame, output_dir: Path) -> None:
    """所有深度叠加到同一张图 —— 每个算法输出 1 张图，3 个子图（XY / XZ / YZ）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    views = [
        ("XY", "X", "Y", "rec_x", "rec_y", "true_x", "true_y"),
        ("XZ", "X", "Z", "rec_x", "rec_z", "true_x", "true_z"),
        ("YZ", "Y", "Z", "rec_y", "rec_z", "true_y", "true_z"),
    ]

    for algo, group in data_df.groupby("algorithm", sort=False):
        fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
        fig.suptitle(f"Reconstruction scatter (all depths overlaid) — algorithm: {algo}",
                     fontsize=13)

        for ax, (view_name, hl, vl, rec_h, rec_v, true_h, true_v) in zip(axes, views):
            ax.axhline(0, color="gray", linewidth=0.6, linestyle="--")
            ax.axvline(0, color="gray", linewidth=0.6, linestyle="--")
            # 每个 (true_h, true_v) 位置用直线连到对应 rec，看误差方向
            for _, r in group.iterrows():
                ax.plot([r[true_h], r[rec_h]], [r[true_v], r[rec_v]],
                        color="lightgray", linewidth=0.5, zorder=1)
            ax.scatter(group[rec_h], group[rec_v], s=18, alpha=0.7,
                       color="steelblue", label="Reconstructed", zorder=3)
            ax.scatter(group[true_h], group[true_v], s=60, marker="+",
                       linewidths=1.5, color="crimson", label="True", zorder=4)

            ax.set_title(f"{view_name} (all depths)", fontsize=11)
            ax.set_xlabel(f"{hl} (cm)", fontsize=10)
            ax.set_ylabel(f"{vl} (cm)", fontsize=10)
            ax.set_aspect("equal", adjustable="datalim")
            ax.grid(True, linewidth=0.4, alpha=0.5)
            ax.legend(fontsize=8, loc="upper right")

        out_path = output_dir / f"scatter_{algo}_overlay.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"overlay scatter plot: {out_path}")


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def iter_config_dirs(root: Path):
    """如果 root 本身含 merged_event.csv 则只处理 root；否则遍历子目录。"""
    if (root / "merged_event.csv").exists():
        yield root
        return
    for child in sorted(root.iterdir(), key=lambda p: p.name):
        if child.is_dir() and (child / "merged_event.csv").exists():
            yield child


def run(args) -> None:
    artifacts: Path = args.artifacts.resolve()
    input_dir: Path = args.input_dir.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"not a directory: {input_dir}")

    # 加载一次模型
    norm = np.load(artifacts / "norm.npz")
    x_mean = norm["x_mean"]; x_std = norm["x_std"]; y_scale = float(norm["y_scale"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PositionMLP(in_dim=96, hidden=args.hidden, out_dim=3, dropout=0.0).to(device)
    model.load_state_dict(torch.load(artifacts / "best.pt", map_location=device))
    model.eval()
    print(f"device = {device}  artifacts = {artifacts}")

    cfg_dirs = list(iter_config_dirs(input_dir))
    print(f"found {len(cfg_dirs)} config dirs under {input_dir}\n")
    if not cfg_dirs:
        raise SystemExit("no config with merged_event.csv found")

    rows: List[dict] = []
    ok = 0
    for cfg_dir in cfg_dirs:
        r = predict_one(cfg_dir, model, x_mean, x_std, y_scale, device, args)
        if r is not None:
            rows.append(r)
            ok += 1

    if not rows:
        raise SystemExit("no valid config produced a prediction")

    # 汇总 CSV
    data_df = pd.DataFrame(rows)
    stat_df = compute_stats(rows)
    combined = pd.concat([data_df, stat_df], ignore_index=True)
    out_csv = Path(args.output).resolve() if args.output else input_dir / "accuracy_summary_nn.csv"
    combined.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"\n[DONE] processed {ok}/{len(cfg_dirs)} configs")
    print(f"       CSV: {out_csv}")

    # 终端汇总
    for algo, group in data_df.groupby("algorithm", sort=False):
        n = len(group)
        print(f"\n  Algorithm: {algo}  (n={n})")
        print(f"    bias_x={group['dx'].mean():.4f}  std_x={group['dx'].std(ddof=1) if n>1 else 0:.4f}")
        print(f"    bias_y={group['dy'].mean():.4f}  std_y={group['dy'].std(ddof=1) if n>1 else 0:.4f}")
        print(f"    bias_z={group['dz'].mean():.4f}  std_z={group['dz'].std(ddof=1) if n>1 else 0:.4f}")
        print(f"    bias_d={group['d'].mean():.4f}   std_d={group['d'].std(ddof=1) if n>1 else 0:.4f}")

    # 散点图：逐深度切片 + 全深度叠加
    if not args.no_plot:
        plot_scatter_views(data_df, input_dir)
        plot_scatter_overlay(data_df, input_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="NN 批量推理 + 精度汇总 + 散点图（全卦限对称镜像版）"
    )
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="总文件夹或单个 config 目录")
    parser.add_argument("--artifacts", type=Path,
                        default=Path("workflow/flow3_NN/artifacts_p5000_norm"))
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--photons-per-sample", type=int, default=5000)
    parser.add_argument("--normalize-counts", action="store_true",
                        help="必须和训练时一致（artifacts 名含 _norm 时要加）")
    parser.add_argument("--append-to-recon", action="store_true", default=True,
                        help="把 nn_mlp_sym 追加到每个 config 的 reconstructed_position.csv（默认 True）")
    parser.add_argument("--no-append-recon", dest="append_to_recon", action="store_false")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--output", default=None,
                        help="汇总 CSV 路径（默认 <input-dir>/accuracy_summary_nn.csv）")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
