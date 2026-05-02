# -*- coding: utf-8 -*-
"""
evaluate_double.py —— 双点模型测试集评估（permutation-aware 指标）

流程
====
  1. 读 dataset.npz + split_idx.npz + norm.npz + best.pt
  2. 测试集前向，反归一化回 mm
  3. 对每个样本，比较 identity-match 和 swap-match 总误差，取较小者
     作为"最佳匹配"，再算各项指标
  4. 输出 test_metrics.csv 和 2x3 hexbin 散点（true vs pred 6 个分量）
     另附一张点距 / 间距误差直方图

运行：
    uv run python workflow/flow4_NN_double/evaluate_double.py
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import torch

from model_double import DoublePointMLP


def best_match(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    """逐样本：返回经"最优 A↔B 匹配"对齐后的预测，shape 同 pred。"""
    pa, pb = pred[:, :3], pred[:, 3:]
    ya, yb = true[:, :3], true[:, 3:]
    err_id = np.linalg.norm(pa - ya, axis=1) + np.linalg.norm(pb - yb, axis=1)
    err_sw = np.linalg.norm(pa - yb, axis=1) + np.linalg.norm(pb - ya, axis=1)
    swap = err_sw < err_id
    out = pred.copy()
    out[swap, :3] = pb[swap]
    out[swap, 3:] = pa[swap]
    return out


def evaluate(args) -> None:
    artifacts: Path = args.artifacts.resolve()

    data = np.load(args.dataset, allow_pickle=True)
    X_all: np.ndarray = data["X"].astype(np.float32)
    y_all: np.ndarray = data["y"].astype(np.float32)   # mm

    split = np.load(artifacts / "split_idx.npz")
    idx_test: np.ndarray = split["idx_test"]

    norm = np.load(artifacts / "norm.npz")
    x_mean = norm["x_mean"]; x_std = norm["x_std"]; y_scale = float(norm["y_scale"])

    X_test = (X_all[idx_test] - x_mean) / x_std
    y_test = y_all[idx_test]                            # mm 真值

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = DoublePointMLP(in_dim=96, hidden=args.hidden, out_dim=6, dropout=0.0).to(device)
    state = torch.load(artifacts / "best.pt", map_location=device)
    model.load_state_dict(state)
    model.eval()

    preds = []
    with torch.no_grad():
        for i in range(0, len(X_test), args.batch_size):
            xb = torch.from_numpy(X_test[i : i + args.batch_size]).to(device)
            preds.append(model(xb).cpu().numpy())
    y_pred = np.concatenate(preds, axis=0) * y_scale     # 反归一化回 mm

    # 最优匹配（A↔B 对齐）
    y_pred_aln = best_match(y_pred, y_test)

    # 每点 3D 误差
    err_a = y_pred_aln[:, :3] - y_test[:, :3]
    err_b = y_pred_aln[:, 3:] - y_test[:, 3:]

    mae_a = np.mean(np.abs(err_a), axis=0)               # (3,)
    mae_b = np.mean(np.abs(err_b), axis=0)
    rmse_a = np.sqrt(np.mean(err_a ** 2, axis=0))
    rmse_b = np.sqrt(np.mean(err_b ** 2, axis=0))
    d_a = np.linalg.norm(err_a, axis=1)
    d_b = np.linalg.norm(err_b, axis=1)
    rmse3d_a = float(np.sqrt(np.mean(d_a ** 2)))
    rmse3d_b = float(np.sqrt(np.mean(d_b ** 2)))

    # 间距（separation）误差
    sep_true = np.linalg.norm(y_test[:, :3] - y_test[:, 3:], axis=1)
    sep_pred = np.linalg.norm(y_pred_aln[:, :3] - y_pred_aln[:, 3:], axis=1)
    sep_err = sep_pred - sep_true
    sep_mae = float(np.mean(np.abs(sep_err)))
    sep_rmse = float(np.sqrt(np.mean(sep_err ** 2)))

    print("=== Test Metrics (mm, A/B aligned by best match) ===")
    print(f"Point A  MAE  (x,y,z) = ({mae_a[0]:.3f}, {mae_a[1]:.3f}, {mae_a[2]:.3f})  3D RMSE = {rmse3d_a:.3f}")
    print(f"Point B  MAE  (x,y,z) = ({mae_b[0]:.3f}, {mae_b[1]:.3f}, {mae_b[2]:.3f})  3D RMSE = {rmse3d_b:.3f}")
    print(f"Separation: MAE = {sep_mae:.3f}  RMSE = {sep_rmse:.3f}")

    # 写 CSV
    with open(artifacts / "test_metrics.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric", "ax", "ay", "az", "bx", "by", "bz", "3d_A", "3d_B", "separation"])
        w.writerow(["MAE_mm",
                    f"{mae_a[0]:.6f}", f"{mae_a[1]:.6f}", f"{mae_a[2]:.6f}",
                    f"{mae_b[0]:.6f}", f"{mae_b[1]:.6f}", f"{mae_b[2]:.6f}",
                    "", "", f"{sep_mae:.6f}"])
        w.writerow(["RMSE_mm",
                    f"{rmse_a[0]:.6f}", f"{rmse_a[1]:.6f}", f"{rmse_a[2]:.6f}",
                    f"{rmse_b[0]:.6f}", f"{rmse_b[1]:.6f}", f"{rmse_b[2]:.6f}",
                    f"{rmse3d_a:.6f}", f"{rmse3d_b:.6f}", f"{sep_rmse:.6f}"])
    print(f"指标已写入 {artifacts / 'test_metrics.csv'}")

    # 2x3 hexbin scatter
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["ax", "ay", "az", "bx", "by", "bz"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 9), constrained_layout=True)
    limit = y_scale * 1.15
    for i, name in enumerate(labels):
        ax = axes.flat[i]
        hb = ax.hexbin(y_test[:, i], y_pred_aln[:, i], gridsize=60, mincnt=1, cmap="viridis")
        ax.plot([-limit, limit], [-limit, limit], "r--", linewidth=1)
        ax.set_xlim(-limit, limit); ax.set_ylim(-limit, limit)
        ax.set_xlabel(f"{name}_true (mm)")
        ax.set_ylabel(f"{name}_pred (mm)")
        ax.set_title(f"{name}")
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)
        fig.colorbar(hb, ax=ax, label="count")
    fig.suptitle(f"NN Double-Point — Test Set (N={len(y_test):,}, A/B aligned)", fontsize=13)
    out_png = artifacts / "test_scatter.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"图已保存：{out_png}")

    # 距离 / 间距误差直方图
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    axes[0].hist(d_a, bins=60, color="steelblue", edgecolor="white")
    axes[0].set_title(f"|err_A|  3D-RMSE={rmse3d_a:.3f} mm")
    axes[0].set_xlabel("mm"); axes[0].grid(True, alpha=0.3)
    axes[1].hist(d_b, bins=60, color="darkorange", edgecolor="white")
    axes[1].set_title(f"|err_B|  3D-RMSE={rmse3d_b:.3f} mm")
    axes[1].set_xlabel("mm"); axes[1].grid(True, alpha=0.3)
    axes[2].hist(sep_err, bins=60, color="seagreen", edgecolor="white")
    axes[2].axvline(0, color="red", linestyle="--")
    axes[2].set_title(f"separation err  MAE={sep_mae:.3f} mm")
    axes[2].set_xlabel("pred − true (mm)"); axes[2].grid(True, alpha=0.3)
    fig.savefig(artifacts / "test_error_hist.png", dpi=150)
    plt.close(fig)
    print(f"误差直方图：{artifacts / 'test_error_hist.png'}")


def main() -> None:
    parser = argparse.ArgumentParser(description="flow4 双点模型测试集评估")
    parser.add_argument("--dataset", type=Path,
                        default=Path("workflow/flow4_NN_double/artifacts/dataset.npz"))
    parser.add_argument("--artifacts", type=Path,
                        default=Path("workflow/flow4_NN_double/artifacts"))
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=4096)
    args = parser.parse_args()
    evaluate(args)


if __name__ == "__main__":
    main()
