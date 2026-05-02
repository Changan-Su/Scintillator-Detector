# -*- coding: utf-8 -*-
"""
predict_double.py —— 单 config 双点推理

输入：一个 DP_*/ 目录（含 metadata.csv + AnaEx01_nt_PhotonFaceBlockEvent_t*.csv）
流程：
  1. 用与训练完全一致的方式聚合光子（photons-per-sample / normalize-counts）
  2. 加载 best.pt + norm.npz，前向得到 (N, 6) 预测（mm）
  3. 与真值最优匹配后保存 nn_double_predicted.csv
  4. 终端打印 "mean pred A/B" vs "ground truth A/B"

运行：
    uv run python workflow/flow4_NN_double/predict_double.py \
        --config-dir Results_flow4_trainingdata/DP_0001_S0p3_...
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from build_dataset_double import (
    aggregate_to_samples,
    load_config_photons,
    read_true_positions_mm,
)
from model_double import DoublePointMLP


def best_match_to_truth(pred_6: np.ndarray, truth_6: np.ndarray) -> np.ndarray:
    """逐样本：把 (pa, pb) 重排为与 truth (ya, yb) 总误差更小的顺序。"""
    pa, pb = pred_6[:, :3], pred_6[:, 3:]
    ya, yb = truth_6[:3], truth_6[3:]
    err_id = np.linalg.norm(pa - ya, axis=1) + np.linalg.norm(pb - yb, axis=1)
    err_sw = np.linalg.norm(pa - yb, axis=1) + np.linalg.norm(pb - ya, axis=1)
    swap = err_sw < err_id
    out = pred_6.copy()
    out[swap, :3] = pb[swap]
    out[swap, 3:] = pa[swap]
    return out


def predict(args) -> None:
    cfg_dir: Path = args.config_dir.resolve()
    if not cfg_dir.is_dir():
        raise SystemExit(f"not a dir: {cfg_dir}")

    artifacts: Path = args.artifacts.resolve()
    norm = np.load(artifacts / "norm.npz")
    x_mean = norm["x_mean"]; x_std = norm["x_std"]; y_scale = float(norm["y_scale"])

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = DoublePointMLP(in_dim=96, hidden=args.hidden, out_dim=6, dropout=0.0).to(device)
    model.load_state_dict(torch.load(artifacts / "best.pt", map_location=device))
    model.eval()
    print(f"device = {device}  artifacts = {artifacts}")

    # ---- 与训练一致的样本聚合 ----
    loaded = load_config_photons(cfg_dir)
    if loaded is None:
        raise SystemExit(f"no PhotonFaceBlockEvent CSVs in {cfg_dir}")
    flat_ch, counts, _ = loaded
    rng = np.random.default_rng(args.seed)
    X = aggregate_to_samples(
        flat_ch, counts,
        photons_per_sample=args.photons_per_sample,
        rng=rng,
        normalize_counts=args.normalize_counts,
    )
    if X is None or len(X) == 0:
        raise SystemExit(f"no samples after aggregation (N photons not enough?)")
    print(f"聚合得到 {len(X)} 个样本（每样本 {args.photons_per_sample} 个光子）")

    X_n = (X - x_mean) / x_std
    preds = []
    with torch.no_grad():
        for i in range(0, len(X_n), args.batch_size):
            xb = torch.from_numpy(X_n[i : i + args.batch_size].astype(np.float32)).to(device)
            preds.append(model(xb).cpu().numpy())
    y_pred = np.concatenate(preds, axis=0) * y_scale     # mm

    # ---- 真值（如果有）→ 最优匹配 ----
    truth = read_true_positions_mm(cfg_dir / "metadata.csv")
    if truth is not None:
        truth_6, frac_a = truth
        y_pred_aln = best_match_to_truth(y_pred, truth_6)
    else:
        truth_6, frac_a = None, None
        y_pred_aln = y_pred

    # ---- 写 CSV ----
    out_df = pd.DataFrame({
        "sample_id": np.arange(len(y_pred_aln)),
        "ax_pred": y_pred_aln[:, 0], "ay_pred": y_pred_aln[:, 1], "az_pred": y_pred_aln[:, 2],
        "bx_pred": y_pred_aln[:, 3], "by_pred": y_pred_aln[:, 4], "bz_pred": y_pred_aln[:, 5],
    })
    out_csv = cfg_dir / "nn_double_predicted.csv"
    out_df.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"已写：{out_csv}")

    # ---- 终端汇总 ----
    mean_pred = y_pred_aln.mean(axis=0)
    print("\n--- mean prediction (mm) ---")
    print(f"A_pred = ({mean_pred[0]:+.3f}, {mean_pred[1]:+.3f}, {mean_pred[2]:+.3f})")
    print(f"B_pred = ({mean_pred[3]:+.3f}, {mean_pred[4]:+.3f}, {mean_pred[5]:+.3f})")
    if truth_6 is not None:
        print(f"A_true = ({truth_6[0]:+.3f}, {truth_6[1]:+.3f}, {truth_6[2]:+.3f})")
        print(f"B_true = ({truth_6[3]:+.3f}, {truth_6[4]:+.3f}, {truth_6[5]:+.3f})")
        print(f"fraction_a = {frac_a:.3f}")
        err_a = np.linalg.norm(mean_pred[:3] - truth_6[:3])
        err_b = np.linalg.norm(mean_pred[3:] - truth_6[3:])
        print(f"|mean_A - A_true| = {err_a:.3f} mm   |mean_B - B_true| = {err_b:.3f} mm")


def main() -> None:
    parser = argparse.ArgumentParser(description="flow4 双点单 config 推理")
    parser.add_argument("--config-dir", type=Path, required=True,
                        help="DP_* 目录（含 metadata.csv + AnaEx01_nt_PhotonFaceBlockEvent_t*.csv）")
    parser.add_argument("--artifacts", type=Path,
                        default=Path("workflow/flow4_NN_double/artifacts"))
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--photons-per-sample", type=int, default=2000,
                        help="必须和训练时一致")
    parser.add_argument("--normalize-counts", action="store_true",
                        help="必须和训练时一致")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    predict(args)


if __name__ == "__main__":
    main()
