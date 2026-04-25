# -*- coding: utf-8 -*-
"""
train2.py —— GPU 加速版训练脚本（大 batch + 数据常驻显存）

和 train.py 的区别
==================
  1. 整份数据一次性 .to(cuda)，之后再也不经过 CPU→GPU 拷贝
  2. 不用 DataLoader，手写张量切片，跳过 Python 调度开销
  3. 默认 batch_size=8192（train.py 是 512），让 GPU 真正忙起来
  4. 默认产物目录改成 artifacts_gpu/，和 train.py 产物隔离
     —— 评估/推理时加 --artifacts workflow/flow3_NN/artifacts_gpu 即可

什么时候用 train2.py 而不是 train.py？
  - 数据集整体 < GPU 显存（RTX 5060 8 GB，能吃 ~1.5 GB 的 X）
  - 想把 epoch 时间从几十秒压到 1–2 秒
  - 不需要复杂的 data augmentation

运行：
    uv run python workflow/flow3_NN/train2.py
可选：
    --batch-size 16384   --epochs 300
    --artifacts workflow/flow3_NN/artifacts_gpu
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from model import PositionMLP, count_parameters


def set_seed(seed: int) -> None:
    import random

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# ---------------------------------------------------------------------------
# 数据：一次性搬到 GPU
# ---------------------------------------------------------------------------
def load_and_prepare_gpu(
    dataset_npz: Path,
    y_scale_cm: float,
    seed: int,
    device: torch.device,
) -> dict:
    """
    读 .npz → 算归一化 → 切 train/val/test → 全部 .to(device)。
    之后训练循环里的张量都已经在 GPU 上，不再有 H2D 拷贝。
    """
    data = np.load(dataset_npz, allow_pickle=True)
    X: np.ndarray = data["X"].astype(np.float32)
    y: np.ndarray = data["y"].astype(np.float32)

    N = len(X)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(N)
    n_test = int(N * 0.10)
    n_val = int(N * 0.10)
    n_train = N - n_val - n_test
    idx_train = idx[:n_train]
    idx_val = idx[n_train : n_train + n_val]
    idx_test = idx[n_train + n_val :]

    x_mean = X[idx_train].mean(axis=0)
    x_std = X[idx_train].std(axis=0) + 1e-6

    Xn = (X - x_mean) / x_std
    yn = y / y_scale_cm

    # 估算占用：N×96×4B + N×3×4B ≈ N × 400B
    total_mb = Xn.nbytes / 1e6 + yn.nbytes / 1e6
    print(f"数据搬到 {device}（约 {total_mb:.0f} MB）...")

    # 整块搬家
    X_t = torch.from_numpy(Xn).to(device)
    y_t = torch.from_numpy(yn).to(device)

    # 把下标也转成 GPU 张量，后面切片更快
    idx_train_t = torch.from_numpy(idx_train).long().to(device)
    idx_val_t = torch.from_numpy(idx_val).long().to(device)
    idx_test_t = torch.from_numpy(idx_test).long().to(device)

    print(f"划分：train={n_train:,}  val={n_val:,}  test={n_test:,}")

    return {
        "X": X_t, "y": y_t,
        "idx_train": idx_train_t, "idx_val": idx_val_t, "idx_test": idx_test_t,
        "idx_train_np": idx_train, "idx_val_np": idx_val, "idx_test_np": idx_test,
        "x_mean": x_mean, "x_std": x_std, "y_scale": y_scale_cm,
    }


# ---------------------------------------------------------------------------
# 单轮：手写 batch 循环，不用 DataLoader
# ---------------------------------------------------------------------------
def run_epoch_gpu(
    model: nn.Module,
    X: torch.Tensor,
    y: torch.Tensor,
    idx: torch.Tensor,
    criterion: nn.Module,
    optimizer: torch.optim.Optimizer | None,
    batch_size: int,
    shuffle: bool,
) -> float:
    """
    在 GPU 张量上直接按下标切 batch。
    X/y 全量在 device 上；idx 是该 split 的样本下标；每个 epoch 只 shuffle idx。
    """
    is_train = optimizer is not None
    model.train(is_train)

    if shuffle:
        perm = torch.randperm(idx.numel(), device=idx.device)
        idx = idx[perm]

    total_loss = 0.0
    total_n = 0
    ctx = torch.enable_grad() if is_train else torch.no_grad()
    with ctx:
        for start in range(0, idx.numel(), batch_size):
            sel = idx[start : start + batch_size]
            xb = X.index_select(0, sel)
            yb = y.index_select(0, sel)

            pred = model(xb)
            loss = criterion(pred, yb)

            if is_train:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()

            bs = xb.size(0)
            total_loss += loss.item() * bs
            total_n += bs
    return total_loss / max(total_n, 1)


# ---------------------------------------------------------------------------
# 主训练
# ---------------------------------------------------------------------------
def train(args) -> None:
    set_seed(args.seed)
    artifacts_dir: Path = args.artifacts.resolve()
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    print(f"device = {device}")
    if device.type == "cuda":
        print(f"GPU = {torch.cuda.get_device_name(device)}  "
              f"capability={torch.cuda.get_device_capability(device)}")
        # cuDNN 自动寻找最快算法（固定 shape 时收益明显）
        torch.backends.cudnn.benchmark = True

    pack = load_and_prepare_gpu(args.dataset, args.y_scale, args.seed, device)

    model = PositionMLP(in_dim=96, hidden=args.hidden, out_dim=3, dropout=args.dropout).to(device)
    print(f"可训练参数：{count_parameters(model):,}")

    criterion = nn.SmoothL1Loss()
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=8
    )

    history = {"train_loss": [], "val_loss": [], "lr": []}
    best_val = float("inf")
    patience_left = args.early_stop_patience

    t0 = time.time()
    for epoch in range(1, args.epochs + 1):
        t_ep = time.time()
        tl = run_epoch_gpu(model, pack["X"], pack["y"], pack["idx_train"],
                           criterion, optimizer, args.batch_size, shuffle=True)
        vl = run_epoch_gpu(model, pack["X"], pack["y"], pack["idx_val"],
                           criterion, None, args.batch_size, shuffle=False)
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
            dt = time.time() - t_ep
            print(f"epoch {epoch:3d} | train {tl:.5f} | val {vl:.5f} | "
                  f"lr {cur_lr:.2e} | {dt:.2f}s{flag}")

        if patience_left <= 0:
            print(f"early stop at epoch {epoch} (no improvement for {args.early_stop_patience} epochs)")
            break

    print(f"\n总用时 {time.time() - t0:.1f}s  最佳 val loss = {best_val:.5f}")

    np.savez(
        artifacts_dir / "norm.npz",
        x_mean=pack["x_mean"], x_std=pack["x_std"],
        y_scale=np.float32(pack["y_scale"]),
    )
    np.savez(
        artifacts_dir / "split_idx.npz",
        idx_train=pack["idx_train_np"],
        idx_val=pack["idx_val_np"],
        idx_test=pack["idx_test_np"],
    )
    print(f"已保存：best.pt / norm.npz / split_idx.npz 到 {artifacts_dir}")

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


def main() -> None:
    parser = argparse.ArgumentParser(description="NN 训练脚本（GPU 大 batch 版）")
    parser.add_argument(
        "--dataset", type=Path,
        default=Path("workflow/flow3_NN/artifacts/dataset.npz"),
        help="和 train.py 共用同一份 dataset.npz",
    )
    parser.add_argument(
        "--artifacts", type=Path,
        default=Path("workflow/flow3_NN/artifacts_gpu"),
        help="产物目录（和 train.py 的 artifacts/ 隔离）",
    )
    parser.add_argument("--y-scale", type=float, default=1.25)
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--batch-size", type=int, default=8192,
                        help="大 batch 才能喂饱 GPU；显存不够就降到 4096")
    parser.add_argument("--epochs", type=int, default=300,
                        help="大 batch 下每 epoch 步数少，总 epoch 数适当调多")
    parser.add_argument("--lr", type=float, default=2e-3,
                        help="大 batch 常配大一点的 lr（经验法则：lr 随 batch 线性放大）")
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--early-stop-patience", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto")
    args = parser.parse_args()
    train(args)


if __name__ == "__main__":
    main()
