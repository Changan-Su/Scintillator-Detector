# -*- coding: utf-8 -*-
"""
train.py —— 训练主脚本

流程
====
  1. 读 dataset.npz → 得到 X (N,96), y (N,3)
  2. 标准化 X（每个通道减均值除以标准差）；y 除以 half_length 归一化到 [-1,1]
  3. 8:1:1 划分 train / val / test
  4. 封装成 PyTorch 的 Dataset + DataLoader（自动分 batch、shuffle）
  5. 构建模型、损失、优化器、学习率调度
  6. 循环 epoch：训练 → 验证 → 若 val 更好就保存 → 早停
  7. 产物：best.pt（最好模型权重）、norm.npz（标准化参数）、loss_curve.png、split_idx.npz

产物位置：workflow/flow3_NN/artifacts/

运行：
    uv run python workflow/flow3_NN/train.py
可选：
    --epochs 200        --batch-size 256
    --lr 1e-3           --hidden 256
    --seed 42           --device auto|cuda|cpu
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

# Python 知识点：相对 import 要求是"包"；同目录下脚本用下面这句即可
from model import PositionMLP, count_parameters


# ---------------------------------------------------------------------------
# 工具：设种子（保证可复现）
# ---------------------------------------------------------------------------
def set_seed(seed: int) -> None:
    """让 numpy / torch / Python 的随机数发生器都从同一个种子出发。"""
    import random

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# ---------------------------------------------------------------------------
# 数据加载 + 标准化 + 划分
# ---------------------------------------------------------------------------
def load_and_prepare(
    dataset_npz: Path,
    y_scale_cm: float,
    seed: int,
) -> dict:
    """
    读 .npz，返回一个包含 train/val/test 张量和归一化参数的大字典。

    参数
      dataset_npz : build_dataset.py 产出的 .npz
      y_scale_cm  : 晶体半长（cm），用来把标签归一化到 [-1,1]
    """
    data = np.load(dataset_npz, allow_pickle=True)
    X: np.ndarray = data["X"].astype(np.float32)   # (N, 96)
    y: np.ndarray = data["y"].astype(np.float32)   # (N, 3) in cm

    N = len(X)
    # 先洗牌下标，再按比例切。洗牌保证 train/val/test 分布一致。
    rng = np.random.default_rng(seed)
    idx = rng.permutation(N)
    n_test = int(N * 0.10)
    n_val = int(N * 0.10)
    n_train = N - n_val - n_test
    idx_train = idx[:n_train]
    idx_val = idx[n_train : n_train + n_val]
    idx_test = idx[n_train + n_val :]

    # 标准化参数"只能用训练集算"，否则泄露验证/测试集的信息到训练过程
    x_mean = X[idx_train].mean(axis=0)
    x_std = X[idx_train].std(axis=0) + 1e-6   # +epsilon 防除 0

    def norm_x(a: np.ndarray) -> np.ndarray:
        return (a - x_mean) / x_std

    def norm_y(a: np.ndarray) -> np.ndarray:
        return a / y_scale_cm

    X_train, y_train = norm_x(X[idx_train]), norm_y(y[idx_train])
    X_val,   y_val   = norm_x(X[idx_val]),   norm_y(y[idx_val])
    X_test,  y_test  = norm_x(X[idx_test]),  norm_y(y[idx_test])

    print(f"划分：train={n_train:,}  val={n_val:,}  test={n_test:,}")

    return {
        "X_train": X_train, "y_train": y_train,
        "X_val": X_val,     "y_val": y_val,
        "X_test": X_test,   "y_test": y_test,
        "idx_train": idx_train, "idx_val": idx_val, "idx_test": idx_test,
        "x_mean": x_mean, "x_std": x_std, "y_scale": y_scale_cm,
    }


def make_loader(X: np.ndarray, y: np.ndarray, batch_size: int, shuffle: bool) -> DataLoader:
    """把 numpy 数组包成 PyTorch 的 DataLoader。"""
    ds = TensorDataset(torch.from_numpy(X), torch.from_numpy(y))
    # num_workers=0：Windows 下多进程 DataLoader 有坑，先用 0 最稳
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, num_workers=0)


# ---------------------------------------------------------------------------
# 单轮训练 / 验证
# ---------------------------------------------------------------------------
def run_epoch(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.Module,
    optimizer: torch.optim.Optimizer | None,
    device: torch.device,
) -> float:
    """
    optimizer=None 表示"验证模式"（不更新权重）。
    返回该 epoch 的平均 loss（已按样本数加权）。
    """
    is_train = optimizer is not None
    model.train(is_train)           # 切换 Dropout/BN 的行为

    total_loss = 0.0
    total_n = 0
    # torch.no_grad：关掉自动求导，验证时省一半内存、稍快
    ctx = torch.enable_grad() if is_train else torch.no_grad()
    with ctx:
        for xb, yb in loader:
            xb = xb.to(device, non_blocking=True)
            yb = yb.to(device, non_blocking=True)

            pred = model(xb)
            loss = criterion(pred, yb)

            if is_train:
                optimizer.zero_grad()     # 清空上一轮梯度
                loss.backward()           # 自动微分反传
                optimizer.step()          # 按梯度更新权重

            bs = xb.size(0)
            total_loss += loss.item() * bs
            total_n += bs
    return total_loss / max(total_n, 1)


# ---------------------------------------------------------------------------
# 主训练函数
# ---------------------------------------------------------------------------
def train(args) -> None:
    set_seed(args.seed)
    artifacts_dir: Path = args.artifacts.resolve()
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    # 设备选择
    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    print(f"device = {device}")

    # 1. 数据
    pack = load_and_prepare(args.dataset, y_scale_cm=args.y_scale, seed=args.seed)
    train_loader = make_loader(pack["X_train"], pack["y_train"], args.batch_size, shuffle=True)
    val_loader   = make_loader(pack["X_val"],   pack["y_val"],   args.batch_size, shuffle=False)

    # 2. 模型 + 损失 + 优化器
    model = PositionMLP(in_dim=96, hidden=args.hidden, out_dim=3, dropout=args.dropout).to(device)
    print(f"可训练参数：{count_parameters(model):,}")

    # SmoothL1（又叫 Huber loss）：小误差时像 MSE，大误差时像 MAE
    # —— 比纯 MSE 更抗离群点，对物理数据特别合适
    criterion = nn.SmoothL1Loss()
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    # 学习率调度：验证 loss 不再下降时自动减半
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=8
    )

    # 3. 训练循环
    history = {"train_loss": [], "val_loss": [], "lr": []}
    best_val = float("inf")
    patience_left = args.early_stop_patience

    t0 = time.time()
    for epoch in range(1, args.epochs + 1):
        tl = run_epoch(model, train_loader, criterion, optimizer, device)
        vl = run_epoch(model, val_loader, criterion, None, device)
        cur_lr = optimizer.param_groups[0]["lr"]
        scheduler.step(vl)

        history["train_loss"].append(tl)
        history["val_loss"].append(vl)
        history["lr"].append(cur_lr)

        improved = vl < best_val - 1e-6
        flag = ""
        if improved:
            best_val = vl
            patience_left = args.early_stop_patience
            torch.save(model.state_dict(), artifacts_dir / "best.pt")
            flag = "  [saved]"
        else:
            patience_left -= 1

        if epoch == 1 or epoch % 5 == 0 or improved:
            print(f"epoch {epoch:3d} | train {tl:.5f} | val {vl:.5f} | lr {cur_lr:.2e}{flag}")

        if patience_left <= 0:
            print(f"early stop at epoch {epoch} (no improvement for {args.early_stop_patience} epochs)")
            break

    print(f"\n用时 {time.time() - t0:.1f}s  最佳 val loss = {best_val:.5f}")

    # 4. 保存归一化参数 + 划分下标
    np.savez(
        artifacts_dir / "norm.npz",
        x_mean=pack["x_mean"], x_std=pack["x_std"],
        y_scale=np.float32(pack["y_scale"]),
    )
    np.savez(
        artifacts_dir / "split_idx.npz",
        idx_train=pack["idx_train"], idx_val=pack["idx_val"], idx_test=pack["idx_test"],
    )
    print(f"已保存：best.pt / norm.npz / split_idx.npz 到 {artifacts_dir}")

    # 5. loss 曲线图
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(7, 4), constrained_layout=True)
        ax.plot(history["train_loss"], label="train")
        ax.plot(history["val_loss"], label="val")
        ax.set_xlabel("epoch")
        ax.set_ylabel("SmoothL1 loss (normalized)")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.savefig(artifacts_dir / "loss_curve.png", dpi=150)
        plt.close(fig)
        print(f"loss 曲线：{artifacts_dir / 'loss_curve.png'}")
    except Exception as exc:
        print(f"(loss 曲线保存失败：{exc})")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="NN 训练脚本")
    parser.add_argument(
        "--dataset", type=Path,
        default=Path("workflow/flow3_NN/artifacts/dataset.npz"),
    )
    parser.add_argument(
        "--artifacts", type=Path,
        default=Path("workflow/flow3_NN/artifacts"),
    )
    parser.add_argument("--y-scale", type=float, default=1.25, help="晶体半长（cm），把 y 归一化到 [-1,1]")
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--early-stop-patience", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto")
    args = parser.parse_args()
    train(args)


if __name__ == "__main__":
    main()
