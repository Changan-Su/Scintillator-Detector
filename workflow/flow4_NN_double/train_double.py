# -*- coding: utf-8 -*-
"""
train_double.py —— 双点 MLP 训练（GPU 大 batch + permutation-invariant loss）

核心区别于 flow3 train2.py
==========================
  1. 标签 y 是 6 维 (ax, ay, az, bx, by, bz)，而非 3 维。
  2. **Loss 必须对 A↔B 排列不变**：
        loss_id   = SmoothL1(pa, ya) + SmoothL1(pb, yb)
        loss_swap = SmoothL1(pa, yb) + SmoothL1(pb, ya)
        loss      = mean( min(loss_id, loss_swap) )    # element-wise min in batch
     —— 模型若预测出 (B, A) 顺序，训练时不应被惩罚。
  3. y_scale 用 12.5 mm（= 晶体半边长），所有 6 维共用同一个 scale。
  4. 单位全程 mm（与 flow3 的 cm 不同，参考 flow4 README）。

CLI：
    uv run python workflow/flow4_NN_double/train_double.py --epochs 200
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from model_double import DoublePointMLP, count_parameters


def set_seed(seed: int) -> None:
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# ---------------------------------------------------------------------------
# permutation-invariant SmoothL1（per-sample reduction = sum over 3 dims）
# ---------------------------------------------------------------------------
def perm_invariant_loss(
    pred: torch.Tensor,            # (B, 6)
    target: torch.Tensor,          # (B, 6)
    beta: float = 1.0,
    return_match: bool = False,
):
    """min(identity-match, swap-match) SmoothL1，逐样本取较小者后再 mean。

    实现要点：
      - 用 reduction='none' 拿到逐元素 loss
      - 沿坐标维 sum → 得到逐样本 (B,) 的 loss_id 和 loss_swap
      - torch.minimum 元素级别取小者
      - 最后再 mean 成标量
    """
    pa, pb = pred[:, :3], pred[:, 3:]
    ya, yb = target[:, :3], target[:, 3:]

    # 每对 3D 距离的 SmoothL1，sum 而非 mean，确保和"两点损失"语义一致
    l_id_a = F.smooth_l1_loss(pa, ya, beta=beta, reduction="none").sum(dim=1)
    l_id_b = F.smooth_l1_loss(pb, yb, beta=beta, reduction="none").sum(dim=1)
    l_id = l_id_a + l_id_b

    l_sw_a = F.smooth_l1_loss(pa, yb, beta=beta, reduction="none").sum(dim=1)
    l_sw_b = F.smooth_l1_loss(pb, ya, beta=beta, reduction="none").sum(dim=1)
    l_sw = l_sw_a + l_sw_b

    per_sample = torch.minimum(l_id, l_sw)
    if return_match:
        return per_sample.mean(), (l_id <= l_sw)  # True = identity match
    return per_sample.mean()


# ---------------------------------------------------------------------------
# 数据加载（一次性 .to(device)，与 flow3 train2.py 同款）
# ---------------------------------------------------------------------------
def load_and_prepare_gpu(
    dataset_npz: Path,
    y_scale_mm: float,
    seed: int,
    device: torch.device,
) -> dict:
    data = np.load(dataset_npz, allow_pickle=True)
    X: np.ndarray = data["X"].astype(np.float32)
    y: np.ndarray = data["y"].astype(np.float32)   # mm
    assert X.shape[1] == 96, f"X 列数应为 96，实际 {X.shape[1]}"
    assert y.shape[1] == 6, f"y 列数应为 6，实际 {y.shape[1]}"

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
    yn = y / y_scale_mm   # 标签归一化到 ~[-1, 1]

    total_mb = Xn.nbytes / 1e6 + yn.nbytes / 1e6
    print(f"数据搬到 {device}（约 {total_mb:.0f} MB）...")
    X_t = torch.from_numpy(Xn).to(device)
    y_t = torch.from_numpy(yn).to(device)

    idx_train_t = torch.from_numpy(idx_train).long().to(device)
    idx_val_t = torch.from_numpy(idx_val).long().to(device)
    idx_test_t = torch.from_numpy(idx_test).long().to(device)
    print(f"划分：train={n_train:,}  val={n_val:,}  test={n_test:,}")

    return {
        "X": X_t, "y": y_t,
        "idx_train": idx_train_t, "idx_val": idx_val_t, "idx_test": idx_test_t,
        "idx_train_np": idx_train, "idx_val_np": idx_val, "idx_test_np": idx_test,
        "x_mean": x_mean, "x_std": x_std, "y_scale": y_scale_mm,
    }


# ---------------------------------------------------------------------------
# 单 epoch
# ---------------------------------------------------------------------------
def run_epoch(
    model: nn.Module,
    X: torch.Tensor,
    y: torch.Tensor,
    idx: torch.Tensor,
    optimizer: torch.optim.Optimizer | None,
    batch_size: int,
    shuffle: bool,
) -> float:
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
            loss = perm_invariant_loss(pred, yb)

            if is_train:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
            bs = xb.size(0)
            total_loss += loss.item() * bs
            total_n += bs
    return total_loss / max(total_n, 1)


# ---------------------------------------------------------------------------
# permutation-invariance sanity check
# ---------------------------------------------------------------------------
def sanity_check_perm_invariance(device: torch.device) -> None:
    """构造随机 (pred, target)，验证 swap target 的 A/B 后 loss 不变。"""
    torch.manual_seed(0)
    pred = torch.randn(64, 6, device=device)
    target = torch.randn(64, 6, device=device)
    target_swapped = torch.cat([target[:, 3:], target[:, :3]], dim=1)
    l1 = perm_invariant_loss(pred, target).item()
    l2 = perm_invariant_loss(pred, target_swapped).item()
    print(f"[sanity] perm-invariance: loss={l1:.6f}  swapped-target loss={l2:.6f}  "
          f"diff={abs(l1 - l2):.2e}  (应≈0)")


# ---------------------------------------------------------------------------
# 主训练
# ---------------------------------------------------------------------------
def train(args) -> None:
    set_seed(args.seed)
    artifacts: Path = args.artifacts.resolve()
    artifacts.mkdir(parents=True, exist_ok=True)

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    print(f"device = {device}")
    if device.type == "cuda":
        print(f"GPU = {torch.cuda.get_device_name(device)}")
        torch.backends.cudnn.benchmark = True

    sanity_check_perm_invariance(device)

    pack = load_and_prepare_gpu(args.dataset, args.y_scale, args.seed, device)

    model = DoublePointMLP(
        in_dim=96, hidden=args.hidden, out_dim=6, dropout=args.dropout,
    ).to(device)
    print(f"可训练参数：{count_parameters(model):,}")

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.lr, weight_decay=args.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=8,
    )

    history = {"train_loss": [], "val_loss": [], "lr": []}
    best_val = float("inf")
    patience_left = args.early_stop_patience

    t0 = time.time()
    for epoch in range(1, args.epochs + 1):
        t_ep = time.time()
        tl = run_epoch(model, pack["X"], pack["y"], pack["idx_train"],
                       optimizer, args.batch_size, shuffle=True)
        vl = run_epoch(model, pack["X"], pack["y"], pack["idx_val"],
                       None, args.batch_size, shuffle=False)
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
            torch.save(model.state_dict(), artifacts / "best.pt")
            flag = "  [saved]"
        else:
            patience_left -= 1

        if epoch == 1 or epoch % 5 == 0 or improved:
            dt = time.time() - t_ep
            print(f"epoch {epoch:3d} | train {tl:.5f} | val {vl:.5f} | "
                  f"lr {cur_lr:.2e} | {dt:.2f}s{flag}")

        if patience_left <= 0:
            print(f"early stop at epoch {epoch}")
            break

    print(f"\n总用时 {time.time() - t0:.1f}s   最佳 val loss = {best_val:.5f}")

    np.savez(
        artifacts / "norm.npz",
        x_mean=pack["x_mean"], x_std=pack["x_std"],
        y_scale=np.float32(pack["y_scale"]),
    )
    np.savez(
        artifacts / "split_idx.npz",
        idx_train=pack["idx_train_np"],
        idx_val=pack["idx_val_np"],
        idx_test=pack["idx_test_np"],
    )
    print(f"已保存：best.pt / norm.npz / split_idx.npz 到 {artifacts}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(7, 4), constrained_layout=True)
        ax.plot(history["train_loss"], label="train")
        ax.plot(history["val_loss"], label="val")
        ax.set_xlabel("epoch")
        ax.set_ylabel("perm-invariant SmoothL1 (normalized)")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.savefig(artifacts / "loss_curve.png", dpi=150)
        plt.close(fig)
        print(f"loss 曲线：{artifacts / 'loss_curve.png'}")
    except Exception as exc:
        print(f"(loss 曲线保存失败：{exc})")


def main() -> None:
    parser = argparse.ArgumentParser(description="flow4 双点 MLP 训练")
    parser.add_argument(
        "--dataset", type=Path,
        default=Path("workflow/flow4_NN_double/artifacts/dataset.npz"),
    )
    parser.add_argument(
        "--artifacts", type=Path,
        default=Path("workflow/flow4_NN_double/artifacts"),
    )
    parser.add_argument("--y-scale", type=float, default=12.5,
                        help="标签归一化常数（mm），默认 12.5 = 晶体半边长")
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--batch-size", type=int, default=8192)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--lr", type=float, default=2e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--early-stop-patience", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto")
    args = parser.parse_args()
    train(args)


if __name__ == "__main__":
    main()
