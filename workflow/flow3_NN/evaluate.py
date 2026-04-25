# -*- coding: utf-8 -*-
"""
evaluate.py —— 在测试集上评估已训练模型

功能
====
  1. 读 dataset.npz + split_idx.npz + norm.npz + best.pt
  2. 在测试集上前向，反归一化回 cm
  3. 算 MAE / RMSE（逐轴 + 3D）
  4. 画 3 张 scatter：x_true vs x_pred, y, z
  5. 画一张 3D 误差直方图
  6. 把指标写进 test_metrics.csv

运行：
    uv run python workflow/flow3_NN/evaluate.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import torch

from model import PositionMLP


def evaluate(args) -> None:
    artifacts: Path = args.artifacts.resolve()

    # 1. 加载 dataset + split + norm
    data = np.load(args.dataset, allow_pickle=True)
    X_all: np.ndarray = data["X"].astype(np.float32)
    y_all: np.ndarray = data["y"].astype(np.float32)

    split = np.load(artifacts / "split_idx.npz")
    idx_test: np.ndarray = split["idx_test"]

    norm = np.load(artifacts / "norm.npz")
    x_mean = norm["x_mean"]; x_std = norm["x_std"]; y_scale = float(norm["y_scale"])

    X_test = (X_all[idx_test] - x_mean) / x_std   # 标准化
    y_test = y_all[idx_test]                       # 保持 cm，作为真值

    # 2. 加载模型权重
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PositionMLP(in_dim=96, hidden=args.hidden, out_dim=3, dropout=0.0).to(device)
    state = torch.load(artifacts / "best.pt", map_location=device)
    model.load_state_dict(state)
    model.eval()

    # 3. 推理（分批，防止 GPU 显存爆掉）
    preds = []
    with torch.no_grad():
        for i in range(0, len(X_test), args.batch_size):
            xb = torch.from_numpy(X_test[i : i + args.batch_size]).to(device)
            preds.append(model(xb).cpu().numpy())
    y_pred_n = np.concatenate(preds, axis=0)       # 归一化空间的预测
    y_pred = y_pred_n * y_scale                     # 反归一化回 cm

    # 4. 指标
    err = y_pred - y_test                           # (N, 3)
    mae_axis = np.mean(np.abs(err), axis=0)
    rmse_axis = np.sqrt(np.mean(err ** 2, axis=0))
    rmse_3d = float(np.sqrt(np.mean(np.sum(err ** 2, axis=1))))

    print("=== Test Metrics (cm) ===")
    print(f"MAE  (x, y, z) = ({mae_axis[0]:.4f}, {mae_axis[1]:.4f}, {mae_axis[2]:.4f})")
    print(f"RMSE (x, y, z) = ({rmse_axis[0]:.4f}, {rmse_axis[1]:.4f}, {rmse_axis[2]:.4f})")
    print(f"3D RMSE        = {rmse_3d:.4f}")

    # 5. 写 CSV
    import csv
    with open(artifacts / "test_metrics.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric", "x", "y", "z", "3d"])
        w.writerow(["MAE_cm",  f"{mae_axis[0]:.6f}",  f"{mae_axis[1]:.6f}",  f"{mae_axis[2]:.6f}",  ""])
        w.writerow(["RMSE_cm", f"{rmse_axis[0]:.6f}", f"{rmse_axis[1]:.6f}", f"{rmse_axis[2]:.6f}", f"{rmse_3d:.6f}"])
    print(f"指标已写入 {artifacts / 'test_metrics.csv'}")

    # 6. 画图
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 10), constrained_layout=True)
    axis_names = ["x", "y", "z"]
    limit = y_scale * 1.15

    for i, name in enumerate(axis_names):
        ax = axes.flat[i]
        # hexbin：预测点太多时，用六角热图比 scatter 快且清晰
        hb = ax.hexbin(y_test[:, i], y_pred[:, i], gridsize=60, mincnt=1, cmap="viridis")
        lims = (-limit, limit)
        ax.plot(lims, lims, "r--", linewidth=1, label="y = x (ideal)")
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel(f"{name}_true (cm)")
        ax.set_ylabel(f"{name}_pred (cm)")
        ax.set_title(f"{name} axis  MAE={mae_axis[i]:.3f} RMSE={rmse_axis[i]:.3f}")
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=8)
        fig.colorbar(hb, ax=ax, label="count")

    # 第 4 子图：3D 距离误差分布
    ax = axes.flat[3]
    dist = np.sqrt(np.sum(err ** 2, axis=1))
    ax.hist(dist, bins=60, color="steelblue", edgecolor="white")
    ax.axvline(rmse_3d, color="red", linestyle="--", label=f"RMSE_3D = {rmse_3d:.3f} cm")
    ax.set_xlabel("|pred − true|  (cm)")
    ax.set_ylabel("count")
    ax.set_title("3D 距离误差分布")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.suptitle(f"NN Position Reconstruction — Test Set (N={len(y_test):,})", fontsize=13)
    out_png = artifacts / "test_scatter.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"图已保存：{out_png}")


def main() -> None:
    parser = argparse.ArgumentParser(description="在测试集上评估 NN 模型")
    parser.add_argument("--dataset", type=Path, default=Path("workflow/flow3_NN/artifacts/dataset.npz"))
    parser.add_argument("--artifacts", type=Path, default=Path("workflow/flow3_NN/artifacts"))
    parser.add_argument("--hidden", type=int, default=256, help="必须和训练时一致")
    parser.add_argument("--batch-size", type=int, default=4096)
    args = parser.parse_args()
    evaluate(args)


if __name__ == "__main__":
    main()
