# -*- coding: utf-8 -*-
"""
predict_double_batch.py —— flow4 双点批量推理 + 汇总 + 散点图

输入：一个总目录（如 Results_flow4_testdata），其下含若干 DP_*/ 子目录。
对每个 config 走和 predict_double.py 完全一致的推理流程，并把：
  - 每 config 的 mean prediction 与真值对齐后的误差
  - 汇总到 accuracy_summary_double_nn.csv（末尾追加 BIAS / STD 两行）
  - 6 子图散点（A/B 三轴 true vs pred）→ scatter_double_nn.png
  - separation true vs pred 散点 → separation_double_nn.png

设计动机：和 flow3_NN/predict3.py 的"批量 + 汇总 + 散点"对应起来，
方便把测试集 / 训练集自检结果一次性铺成表 + 图。

用法：
    uv run python workflow/flow4_NN_double/predict_double_batch.py \
        --input-dir Results_flow4_testdata \
        --photons-per-sample 2000 --normalize-counts
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import torch

from build_dataset_double import (
    aggregate_to_samples,
    load_config_photons,
    read_true_positions_mm,
)
from model_double import DoublePointMLP
from predict_double import best_match_to_truth


# ---------------------------------------------------------------------------
# 单 config 推理（精简版：只返回汇总行所需字段，不打日志）
# ---------------------------------------------------------------------------
def predict_one(
    cfg_dir: Path,
    model: torch.nn.Module,
    x_mean: np.ndarray,
    x_std: np.ndarray,
    y_scale: float,
    device: torch.device,
    args,
) -> dict | None:
    """返回该 config 的汇总 dict；缺数据返回 None。"""
    meta = read_true_positions_mm(cfg_dir / "metadata.csv")
    if meta is None:
        print(f"  [SKIP] no metadata: {cfg_dir.name}")
        return None
    truth_6, frac_a = meta

    loaded = load_config_photons(cfg_dir)
    if loaded is None:
        print(f"  [SKIP] no photon CSV: {cfg_dir.name}")
        return None
    flat_ch, counts, _ = loaded

    rng = np.random.default_rng(args.seed)
    X = aggregate_to_samples(
        flat_ch, counts,
        photons_per_sample=args.photons_per_sample,
        rng=rng,
        normalize_counts=args.normalize_counts,
    )
    if X is None or len(X) == 0:
        print(f"  [SKIP] not enough photons: {cfg_dir.name}")
        return None

    X_n = (X - x_mean) / x_std
    preds = []
    with torch.no_grad():
        for i in range(0, len(X_n), args.batch_size):
            xb = torch.from_numpy(X_n[i : i + args.batch_size].astype(np.float32)).to(device)
            preds.append(model(xb).cpu().numpy())
    y_pred = np.concatenate(preds, axis=0) * y_scale  # (N, 6) mm

    # 与真值最优匹配（逐样本）
    y_pred_aln = best_match_to_truth(y_pred, truth_6)

    # 可选写出每 config CSV
    if args.per_config_csv:
        out_df = pd.DataFrame({
            "sample_id": np.arange(len(y_pred_aln)),
            "ax_pred": y_pred_aln[:, 0], "ay_pred": y_pred_aln[:, 1], "az_pred": y_pred_aln[:, 2],
            "bx_pred": y_pred_aln[:, 3], "by_pred": y_pred_aln[:, 4], "bz_pred": y_pred_aln[:, 5],
        })
        out_df.to_csv(cfg_dir / "nn_double_predicted.csv", index=False, float_format="%.6f")

    # 取 sample 平均，抑制 statistical fluctuation
    mean_pred = y_pred_aln.mean(axis=0)

    dx_a = float(mean_pred[0] - truth_6[0])
    dy_a = float(mean_pred[1] - truth_6[1])
    dz_a = float(mean_pred[2] - truth_6[2])
    dist_a = float(np.sqrt(dx_a * dx_a + dy_a * dy_a + dz_a * dz_a))

    dx_b = float(mean_pred[3] - truth_6[3])
    dy_b = float(mean_pred[4] - truth_6[4])
    dz_b = float(mean_pred[5] - truth_6[5])
    dist_b = float(np.sqrt(dx_b * dx_b + dy_b * dy_b + dz_b * dz_b))

    sep_true = float(np.linalg.norm(truth_6[:3] - truth_6[3:]))
    sep_pred = float(np.linalg.norm(mean_pred[:3] - mean_pred[3:]))
    sep_err = sep_pred - sep_true

    return {
        "config_name": cfg_dir.name,
        "ax_true": float(truth_6[0]), "ay_true": float(truth_6[1]), "az_true": float(truth_6[2]),
        "bx_true": float(truth_6[3]), "by_true": float(truth_6[4]), "bz_true": float(truth_6[5]),
        "fraction_a": float(frac_a),
        "ax_pred": float(mean_pred[0]), "ay_pred": float(mean_pred[1]), "az_pred": float(mean_pred[2]),
        "bx_pred": float(mean_pred[3]), "by_pred": float(mean_pred[4]), "bz_pred": float(mean_pred[5]),
        "dx_a": dx_a, "dy_a": dy_a, "dz_a": dz_a, "dist_a": dist_a,
        "dx_b": dx_b, "dy_b": dy_b, "dz_b": dz_b, "dist_b": dist_b,
        "sep_true": sep_true, "sep_pred": sep_pred, "sep_err": sep_err,
        "n_samples": int(len(y_pred_aln)),
    }


# ---------------------------------------------------------------------------
# 汇总尾行（BIAS / STD），mimic flow3 predict3.py 的 footer 习惯
# ---------------------------------------------------------------------------
SUMMARY_NUM_COLS = (
    "dx_a", "dy_a", "dz_a", "dist_a",
    "dx_b", "dy_b", "dz_b", "dist_b",
    "sep_err",
)


def make_summary_rows(df: pd.DataFrame) -> List[dict]:
    """计算每个数值列的 BIAS（均值）和 STD（样本标准差），返回两行 dict。"""
    n = len(df)

    def _blank() -> dict:
        return {c: "" for c in df.columns}

    bias = _blank()
    bias["config_name"] = "[SUMMARY BIAS]"
    std = _blank()
    std["config_name"] = "[SUMMARY STD]"

    for c in SUMMARY_NUM_COLS:
        b = float(df[c].mean())
        s = float(df[c].std(ddof=1)) if n > 1 else 0.0
        bias[c] = f"BIAS={b:.6f}"
        std[c] = f"STD={s:.6f}"
    return [bias, std]


# ---------------------------------------------------------------------------
# 散点图
# ---------------------------------------------------------------------------
def plot_scatter_overlay(df: pd.DataFrame, out_path: Path) -> None:
    """3 子图（XY / XZ / YZ）三视图叠加，所有 config 画在同一张图上。

    每个 config 对应一对 (A, B) 双点：
      - true_A 用红色 `+`，pred_A 用红色实心圆
      - true_B 用 royalblue `+`，pred_B 用 royalblue 实心圆
      - 真值 → 预测之间用同色细线相连，便于追踪误差方向
    25 mm 立方体边界用灰色虚线方框示意（±12.5 mm）。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 三视图：(子图标题, 横轴标签, 纵轴标签, A 横列, A 纵列, B 横列, B 纵列)
    views = [
        ("XY view", "X", "Y",
         ("ax_true", "ay_true"), ("ax_pred", "ay_pred"),
         ("bx_true", "by_true"), ("bx_pred", "by_pred")),
        ("XZ view", "X", "Z",
         ("ax_true", "az_true"), ("ax_pred", "az_pred"),
         ("bx_true", "bz_true"), ("bx_pred", "bz_pred")),
        ("YZ view", "Y", "Z",
         ("ay_true", "az_true"), ("ay_pred", "az_pred"),
         ("by_true", "bz_true"), ("by_pred", "bz_pred")),
    ]

    half = 12.5  # 25 mm 立方体半边长
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    fig.suptitle("Double-point reconstruction (all configs overlaid)", fontsize=13)

    color_a = "red"
    color_b = "royalblue"

    for idx, (ax, (title, hl, vl, a_true_cols, a_pred_cols, b_true_cols, b_pred_cols)) in enumerate(zip(axes, views)):
        # 25 mm 立方体边界（灰色虚线方框）
        ax.plot([-half, half, half, -half, -half],
                [-half, -half, half, half, -half],
                color="gray", linestyle="--", linewidth=0.8, zorder=1)
        ax.axhline(0, color="gray", linewidth=0.4, linestyle=":", zorder=1)
        ax.axvline(0, color="gray", linewidth=0.4, linestyle=":", zorder=1)

        # 提取整列数据（所有 config）
        a_true_h = df[a_true_cols[0]].to_numpy(dtype=float)
        a_true_v = df[a_true_cols[1]].to_numpy(dtype=float)
        a_pred_h = df[a_pred_cols[0]].to_numpy(dtype=float)
        a_pred_v = df[a_pred_cols[1]].to_numpy(dtype=float)
        b_true_h = df[b_true_cols[0]].to_numpy(dtype=float)
        b_true_v = df[b_true_cols[1]].to_numpy(dtype=float)
        b_pred_h = df[b_pred_cols[0]].to_numpy(dtype=float)
        b_pred_v = df[b_pred_cols[1]].to_numpy(dtype=float)

        # true → pred 连线（同色细线，alpha 低，便于追踪 pair）
        for i in range(len(df)):
            ax.plot([a_true_h[i], a_pred_h[i]], [a_true_v[i], a_pred_v[i]],
                    color=color_a, alpha=0.3, lw=0.5, zorder=2)
            ax.plot([b_true_h[i], b_pred_h[i]], [b_true_v[i], b_pred_v[i]],
                    color=color_b, alpha=0.3, lw=0.5, zorder=2)

        # A 真值（+）与预测（实心圆）
        ax.scatter(a_true_h, a_true_v, s=50, marker="+", linewidths=1.5,
                   color=color_a, label="true A", zorder=4)
        ax.scatter(a_pred_h, a_pred_v, s=18, alpha=0.7,
                   color=color_a, label="pred A", zorder=3)
        # B 真值（+）与预测（实心圆）
        ax.scatter(b_true_h, b_true_v, s=50, marker="+", linewidths=1.5,
                   color=color_b, label="true B", zorder=4)
        ax.scatter(b_pred_h, b_pred_v, s=18, alpha=0.7,
                   color=color_b, label="pred B", zorder=3)

        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(f"{hl} (mm)", fontsize=10)
        ax.set_ylabel(f"{vl} (mm)", fontsize=10)
        ax.grid(True, linewidth=0.4, alpha=0.5)

        # 仅第一个子图画 legend，其它子图不重复
        if idx == 0:
            ax.legend(fontsize=8, loc="upper right")

    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_error_vs_separation(df: pd.DataFrame, out_path: Path) -> None:
    """单子图：横轴 sep_true，纵轴 A/B 两点的位置误差（mm），叠加两组散点。

    诊断"两点距离越近，pred 越塌陷成一坨"等失效模式。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sep_true = df["sep_true"].to_numpy(dtype=float)
    err_a = df["dist_a"].to_numpy(dtype=float)
    err_b = df["dist_b"].to_numpy(dtype=float)

    med_a = float(np.median(err_a))
    med_b = float(np.median(err_b))

    fig, ax = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)

    ax.scatter(sep_true, err_a, s=30, alpha=0.8, color="red",
               label=f"|err A| (median = {med_a:.2f} mm)", zorder=3)
    ax.scatter(sep_true, err_b, s=30, alpha=0.8, color="royalblue",
               label=f"|err B| (median = {med_b:.2f} mm)", zorder=3)

    # 中位数水平参考线（同色虚线，alpha 较低）
    ax.axhline(med_a, color="red", linestyle="--", linewidth=0.8, alpha=0.5, zorder=2)
    ax.axhline(med_b, color="royalblue", linestyle="--", linewidth=0.8, alpha=0.5, zorder=2)

    # 坐标范围：x 从 0 到稍大于 sep_true 最大值；y 从 0 到稍大于最大误差
    x_hi = float(sep_true.max()) if len(sep_true) else 1.0
    y_hi = float(max(err_a.max(), err_b.max())) if len(err_a) else 1.0
    ax.set_xlim(0.0, x_hi * 1.05 + 1e-6)
    ax.set_ylim(0.0, y_hi * 1.05 + 1e-6)

    ax.set_title("Per-point error vs true A-B separation", fontsize=12)
    ax.set_xlabel("True A-B separation (mm)", fontsize=10)
    ax.set_ylabel("Position error (mm)", fontsize=10)
    ax.grid(True, linewidth=0.4, alpha=0.5)
    ax.legend(fontsize=9, loc="upper right")

    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_separation(df: pd.DataFrame, out_path: Path) -> None:
    """单子图：sep_true vs sep_pred + 红色 y=x。诊断"两点是否被压成一坨"。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
    x = df["sep_true"].to_numpy(dtype=float)
    y = df["sep_pred"].to_numpy(dtype=float)
    lo = float(min(x.min(), y.min(), 0.0))
    hi = float(max(x.max(), y.max()))
    pad = 0.05 * (hi - lo + 1e-6)
    lo -= pad; hi += pad
    ax.plot([lo, hi], [lo, hi], color="crimson", linewidth=1.0, zorder=1, label="y = x")
    ax.scatter(x, y, s=30, alpha=0.8, color="seagreen", zorder=3)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("A-B separation: true vs pred (mm)", fontsize=12)
    ax.set_xlabel("sep_true (mm)", fontsize=10)
    ax.set_ylabel("sep_pred (mm)", fontsize=10)
    ax.grid(True, linewidth=0.4, alpha=0.5)
    ax.legend(fontsize=9, loc="upper left")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
COLUMN_ORDER = [
    "config_name",
    "ax_true", "ay_true", "az_true", "bx_true", "by_true", "bz_true", "fraction_a",
    "ax_pred", "ay_pred", "az_pred", "bx_pred", "by_pred", "bz_pred",
    "dx_a", "dy_a", "dz_a", "dist_a",
    "dx_b", "dy_b", "dz_b", "dist_b",
    "sep_true", "sep_pred", "sep_err", "n_samples",
]


def discover_dp_dirs(root: Path) -> List[Path]:
    return [
        p for p in sorted(root.iterdir(), key=lambda q: q.name)
        if p.is_dir() and p.name.startswith("DP_")
    ]


def run(args) -> None:
    input_dir: Path = args.input_dir.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"not a directory: {input_dir}")

    artifacts: Path = args.artifacts.resolve()
    norm = np.load(artifacts / "norm.npz")
    x_mean = norm["x_mean"]; x_std = norm["x_std"]; y_scale = float(norm["y_scale"])

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = DoublePointMLP(in_dim=96, hidden=args.hidden, out_dim=6, dropout=0.0).to(device)
    model.load_state_dict(torch.load(artifacts / "best.pt", map_location=device))
    model.eval()
    print(f"device = {device}  artifacts = {artifacts}")

    dp_dirs = discover_dp_dirs(input_dir)
    n_total = len(dp_dirs)
    print(f"在 {input_dir} 下发现 {n_total} 个 DP_* config\n")
    if n_total == 0:
        raise SystemExit("没有 DP_* 子目录可处理")

    rows: List[dict] = []
    for i, cfg_dir in enumerate(dp_dirs, 1):
        r = predict_one(cfg_dir, model, x_mean, x_std, y_scale, device, args)
        if r is None:
            continue
        rows.append(r)
        print(f"[{i}/{n_total}] {cfg_dir.name}  "
              f"d_a={r['dist_a']:.2f}mm d_b={r['dist_b']:.2f}mm sep_err={r['sep_err']:+.2f}mm")

    if not rows:
        raise SystemExit("无有效 config 产生预测")

    df = pd.DataFrame(rows)[COLUMN_ORDER]
    summary_rows = make_summary_rows(df)
    out_df = pd.concat([df, pd.DataFrame(summary_rows)], ignore_index=True)

    out_csv = input_dir / "accuracy_summary_double_nn.csv"
    out_df.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"\n[DONE] 处理 {len(rows)}/{n_total} 个 config")
    print(f"       CSV: {out_csv}")

    # 散点图
    scatter_path = input_dir / "scatter_double_nn.png"
    sep_path = input_dir / "separation_double_nn.png"
    err_sep_path = input_dir / "error_vs_separation_double_nn.png"
    plot_scatter_overlay(df, scatter_path)
    plot_separation(df, sep_path)
    plot_error_vs_separation(df, err_sep_path)
    print(f"       scatter: {scatter_path}")
    print(f"       separation: {sep_path}")
    print(f"       error-vs-separation: {err_sep_path}")

    # 终端汇总
    mean_da = float(df["dist_a"].mean())
    mean_db = float(df["dist_b"].mean())
    mean_se = float(df["sep_err"].abs().mean())
    print(f"\n  mean dist_a = {mean_da:.3f} mm")
    print(f"  mean dist_b = {mean_db:.3f} mm")
    print(f"  mean |sep_err| = {mean_se:.3f} mm")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="flow4 双点批量推理 + 精度汇总 + 散点图"
    )
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="总目录（含若干 DP_* 子目录），如 Results_flow4_testdata")
    parser.add_argument("--artifacts", type=Path,
                        default=Path("workflow/flow4_NN_double/artifacts"))
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--photons-per-sample", type=int, default=2000,
                        help="必须和训练时一致")
    parser.add_argument("--normalize-counts", action="store_true",
                        help="必须和训练时一致")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--per-config-csv", action="store_true", default=False,
                        help="每个 config 也写一份 nn_double_predicted.csv（默认关闭，避免散落几百个文件）")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
