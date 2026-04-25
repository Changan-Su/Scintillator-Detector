# -*- coding: utf-8 -*-
"""
predict2.py —— 对称镜像版推理

动机
====
训练数据只覆盖了 +X/+Y/+Z 卦限（位置坐标都在 [0, +10] mm），
但晶体本身关于三个轴都是对称的（立方体几何 + 均匀 SiPM 分布）。

本脚本利用这一点：对任意卦限的事件，先把输入 96 维 SiPM 向量**镜像**
到 +++ 卦限的等价形式，用模型推理，再把输出坐标反镜像回去。

判卦限的依据：
  sign_x = sign(N(+X face) - N(-X face))       # +X 面更亮 → 源在 +X 半空间
  sign_y = sign(N(+Y face) - N(-Y face))
  sign_z = sign(N(+Z face) - N(-Z face))

SiPM 通道布局（来自 DetectorConstruction）：
  - Face 0=+X, 1=-X, 2=+Y, 3=-Y, 4=+Z, 5=-Z
  - Face 0/1：j 沿 +Y，k 沿 +Z
  - Face 2/3：j 沿 +X，k 沿 +Z
  - Face 4/5：j 沿 +X，k 沿 +Y
  - 扁平下标：flat = face*16 + j*4 + k

镜像规则（以 X 轴为例）：
  - 交换 face 0 ↔ face 1（两个对面互换）
  - face 2,3,4,5 上 j 方向翻转（j → 3-j），因为 j 轴沿 +X

用法
====
  uv run python workflow/flow3_NN/predict2.py \
      --config-dir "Output/.../S0p3_Xm6p25_Ym6p25_Zm12p5_..." \
      --artifacts workflow/flow3_NN/artifacts_p5000_norm \
      --photons-per-sample 5000 \
      --normalize-counts
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from build_dataset import events_to_tensor, read_true_position_cm
from model import PositionMLP


N_FACES = 6
N_J = 4
N_K = 4
N_CHANNELS = N_FACES * N_J * N_K   # 96


# ---------------------------------------------------------------------------
# 预计算三个轴的镜像通道重排（96 维索引置换）
# ---------------------------------------------------------------------------
def _build_mirror_perm(axis: str) -> np.ndarray:
    """
    返回长度 96 的 int 数组 perm，使得 X_mirrored[i] = X[perm[i]]。

    几何：
      axis='x': 源在 ±X 对称 → 交换 face 0↔1；face 2/3/4/5 的 j 翻转
      axis='y': 源在 ±Y 对称 → 交换 face 2↔3；face 0/1 的 j 翻转；face 4/5 的 k 翻转
      axis='z': 源在 ±Z 对称 → 交换 face 4↔5；face 0/1/2/3 的 k 翻转
    """
    perm = np.empty(N_CHANNELS, dtype=np.int64)
    for new_idx in range(N_CHANNELS):
        face = new_idx // 16
        j = (new_idx % 16) // 4
        k = new_idx % 4

        if axis == "x":
            if face == 0:
                old_face, old_j, old_k = 1, j, k
            elif face == 1:
                old_face, old_j, old_k = 0, j, k
            else:   # 2,3,4,5：j 翻转
                old_face, old_j, old_k = face, 3 - j, k
        elif axis == "y":
            if face == 2:
                old_face, old_j, old_k = 3, j, k
            elif face == 3:
                old_face, old_j, old_k = 2, j, k
            elif face in (0, 1):         # j 沿 +Y → 翻 j
                old_face, old_j, old_k = face, 3 - j, k
            else:                         # face 4,5：k 沿 +Y → 翻 k
                old_face, old_j, old_k = face, j, 3 - k
        elif axis == "z":
            if face == 4:
                old_face, old_j, old_k = 5, j, k
            elif face == 5:
                old_face, old_j, old_k = 4, j, k
            else:                         # face 0,1,2,3：k 沿 +Z → 翻 k
                old_face, old_j, old_k = face, j, 3 - k
        else:
            raise ValueError(f"bad axis: {axis}")

        perm[new_idx] = old_face * 16 + old_j * 4 + old_k
    return perm


PERM_X = _build_mirror_perm("x")
PERM_Y = _build_mirror_perm("y")
PERM_Z = _build_mirror_perm("z")


def _face_sum(X: np.ndarray, face: int) -> np.ndarray:
    """对每个样本，返回指定面 16 个 SiPM 的总计数，shape=(N,)"""
    return X[:, face * 16 : (face + 1) * 16].sum(axis=1)


def detect_octant(X: np.ndarray) -> np.ndarray:
    """
    对每个样本判源所在卦限。
    返回 shape=(N, 3)，元素 ∈ {+1, -1}，对应 (sign_x, sign_y, sign_z)。
    """
    sign_x = np.where(_face_sum(X, 0) >= _face_sum(X, 1), 1, -1)   # +X vs -X
    sign_y = np.where(_face_sum(X, 2) >= _face_sum(X, 3), 1, -1)
    sign_z = np.where(_face_sum(X, 4) >= _face_sum(X, 5), 1, -1)
    return np.stack([sign_x, sign_y, sign_z], axis=1).astype(np.int8)


def canonicalize_to_positive_octant(X: np.ndarray, signs: np.ndarray) -> np.ndarray:
    """
    把每个样本镜像到 +++ 卦限：
      - signs[i, a] = -1 表示第 a 轴原本在负半空间，需要对该样本做该轴镜像。
    返回同 shape=(N, 96) 的新数组。

    实现技巧：逐样本按掩码应用 3 个置换。虽然是循环，但 N 通常不大（几千以内）。
    想完全向量化可以按 2^3 = 8 种 sign 组合分组处理，这里用朴素写法便于读。
    """
    X_out = X.copy()
    mask_x = signs[:, 0] == -1
    mask_y = signs[:, 1] == -1
    mask_z = signs[:, 2] == -1

    if mask_x.any():
        X_out[mask_x] = X_out[mask_x][:, PERM_X]
    if mask_y.any():
        X_out[mask_y] = X_out[mask_y][:, PERM_Y]
    if mask_z.any():
        X_out[mask_z] = X_out[mask_z][:, PERM_Z]
    return X_out


def uncanonicalize_output(y_pred: np.ndarray, signs: np.ndarray) -> np.ndarray:
    """
    把 +++ 卦限里的预测坐标翻回原卦限：sign=-1 的轴 → 取负。
    y_pred : (N, 3)    signs : (N, 3)
    """
    return y_pred * signs.astype(y_pred.dtype)


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def predict(args) -> None:
    artifacts: Path = args.artifacts.resolve()
    cfg_dir: Path = args.config_dir.resolve()

    norm = np.load(artifacts / "norm.npz")
    x_mean = norm["x_mean"]; x_std = norm["x_std"]; y_scale = float(norm["y_scale"])

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PositionMLP(in_dim=96, hidden=args.hidden, out_dim=3, dropout=0.0).to(device)
    model.load_state_dict(torch.load(artifacts / "best.pt", map_location=device))
    model.eval()

    # 读事件并聚合（和训练口径一致）
    event_csv = cfg_dir / "merged_event.csv"
    X = events_to_tensor(
        event_csv,
        min_photons=1,
        photons_per_sample=args.photons_per_sample,
        normalize_counts=args.normalize_counts,
    )
    if X is None or len(X) == 0:
        raise SystemExit(f"没有可用事件：{event_csv}")

    # ---- 对称处理核心 ----
    signs = detect_octant(X)                                 # (N, 3) ∈ {-1,+1}
    X_canon = canonicalize_to_positive_octant(X, signs)      # 镜像到 +++ 卦限
    X_n = (X_canon - x_mean) / x_std                          # 用训练时的 norm

    preds = []
    with torch.no_grad():
        for i in range(0, len(X_n), args.batch_size):
            xb = torch.from_numpy(X_n[i : i + args.batch_size].astype(np.float32)).to(device)
            preds.append(model(xb).cpu().numpy())
    y_pred_canon = np.concatenate(preds, axis=0) * y_scale    # cm，在 +++ 卦限
    y_pred = uncanonicalize_output(y_pred_canon, signs)       # 翻回原卦限

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
    ev_out = cfg_dir / "nn_predicted_events_sym.csv"
    ev_df.to_csv(ev_out, index=False)
    print(f"每事件预测（对称版）已保存：{ev_out}")

    # 卦限投票统计（健康检查：同一 config 内所有事件大多应落在同一卦限）
    oct_counts = pd.Series([tuple(s) for s in signs]).value_counts()
    print("\n卦限分布（同一 config 内事件理论上应高度集中）：")
    for oct_key, cnt in oct_counts.head(4).items():
        print(f"  {oct_key}: {cnt}  ({100*cnt/len(signs):.1f}%)")

    pos_mean = y_pred.mean(axis=0)
    pos_median = np.median(y_pred, axis=0)
    print(f"\n事件平均 (x,y,z) = ({pos_mean[0]:+.4f}, {pos_mean[1]:+.4f}, {pos_mean[2]:+.4f}) cm")
    print(f"事件中位数       = ({pos_median[0]:+.4f}, {pos_median[1]:+.4f}, {pos_median[2]:+.4f}) cm")

    true_pos = read_true_position_cm(cfg_dir / "metadata.csv")
    if true_pos is not None:
        err = pos_mean - true_pos
        print(f"真值             = ({true_pos[0]:+.4f}, {true_pos[1]:+.4f}, {true_pos[2]:+.4f}) cm")
        print(f"误差             = ({err[0]:+.4f}, {err[1]:+.4f}, {err[2]:+.4f}) cm  "
              f"|err|={np.linalg.norm(err):.4f}")

    if args.append_to_recon:
        rec_csv = cfg_dir / "reconstructed_position.csv"
        row = {"algorithm": "nn_mlp_sym", "x_rec": pos_mean[0],
               "y_rec": pos_mean[1], "z_rec": pos_mean[2]}
        if rec_csv.exists():
            df = pd.read_csv(rec_csv)
            df = df[df["algorithm"] != "nn_mlp_sym"]
            df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
        else:
            df = pd.DataFrame([row])
        df.to_csv(rec_csv, index=False)
        print(f"已追加到 {rec_csv}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="NN 重建（对称镜像版）—— 覆盖所有 8 个卦限"
    )
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--artifacts", type=Path,
                        default=Path("workflow/flow3_NN/artifacts_p5000_norm"))
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--photons-per-sample", type=int, default=5000)
    parser.add_argument("--normalize-counts", action="store_true")
    parser.add_argument("--append-to-recon", action="store_true")
    args = parser.parse_args()
    predict(args)


if __name__ == "__main__":
    main()
