# -*- coding: utf-8 -*-
"""
predict.py —— 对单个配置目录做 NN 位置重建

用法：
    uv run python workflow/flow3_NN/predict.py \
        --config-dir "Output/SiPM6_Output_20260404_234359/<某 config>"

会做两件事：
  1. 事件级推理：读该 config 的 merged_event.csv，对每个事件预测 (x,y,z)，
     保存为 nn_predicted_events.csv
  2. 事件平均 → 一个 "per-config 重建位置"，打印出来并可追加到
     reconstructed_position.csv（作为新的 "nn_mlp" 算法行）
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from build_dataset import events_to_tensor, read_true_position_cm
from model import PositionMLP


def predict(args) -> None:
    artifacts: Path = args.artifacts.resolve()
    cfg_dir: Path = args.config_dir.resolve()

    # 1. 载入归一化参数 + 模型
    norm = np.load(artifacts / "norm.npz")
    x_mean = norm["x_mean"]; x_std = norm["x_std"]; y_scale = float(norm["y_scale"])

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = PositionMLP(in_dim=96, hidden=args.hidden, out_dim=3, dropout=0.0).to(device)
    model.load_state_dict(torch.load(artifacts / "best.pt", map_location=device))
    model.eval()

    # 2. 读这个 config 的事件
    #    - photons_per_sample：把该 config 的单光子行合并成"一个γ事件样本"，
    #      和训练时的口径对齐（训练用 5000 光子/样本就写 5000）
    #    - normalize_counts：转成光子分数，让模型对光子总数免疫
    #      这样哪怕真实事件是 3万 光子，推理依然准确
    event_csv = cfg_dir / "merged_event.csv"
    X = events_to_tensor(
        event_csv,
        min_photons=1,
        photons_per_sample=args.photons_per_sample,
        normalize_counts=args.normalize_counts,
    )
    if X is None or len(X) == 0:
        raise SystemExit(f"没有可用事件：{event_csv}")
    X_n = (X - x_mean) / x_std

    # 3. 分批推理
    preds = []
    with torch.no_grad():
        for i in range(0, len(X_n), args.batch_size):
            xb = torch.from_numpy(X_n[i : i + args.batch_size].astype(np.float32)).to(device)
            preds.append(model(xb).cpu().numpy())
    y_pred = np.concatenate(preds, axis=0) * y_scale   # cm

    # 4. 事件级结果 → CSV
    ev_df = pd.DataFrame({
        "event_idx": np.arange(len(y_pred)),
        "x_rec_cm": y_pred[:, 0],
        "y_rec_cm": y_pred[:, 1],
        "z_rec_cm": y_pred[:, 2],
        "total_photons": X.sum(axis=1),
    })
    ev_out = cfg_dir / "nn_predicted_events.csv"
    ev_df.to_csv(ev_out, index=False)
    print(f"每事件预测已保存：{ev_out}")

    # 5. 事件平均 → 单个配置位置
    pos_mean = y_pred.mean(axis=0)
    pos_median = np.median(y_pred, axis=0)
    print(f"\n事件平均 (x,y,z) = ({pos_mean[0]:.4f}, {pos_mean[1]:.4f}, {pos_mean[2]:.4f}) cm")
    print(f"事件中位数       = ({pos_median[0]:.4f}, {pos_median[1]:.4f}, {pos_median[2]:.4f}) cm")

    true_pos = read_true_position_cm(cfg_dir / "metadata.csv")
    if true_pos is not None:
        err = pos_mean - true_pos
        print(f"真值             = ({true_pos[0]:.4f}, {true_pos[1]:.4f}, {true_pos[2]:.4f}) cm")
        print(f"误差             = ({err[0]:+.4f}, {err[1]:+.4f}, {err[2]:+.4f}) cm  "
              f"|err|={np.linalg.norm(err):.4f}")

    # 6. 可选：追加到 reconstructed_position.csv
    if args.append_to_recon:
        rec_csv = cfg_dir / "reconstructed_position.csv"
        row = {"algorithm": "nn_mlp", "x_rec": pos_mean[0], "y_rec": pos_mean[1], "z_rec": pos_mean[2]}
        if rec_csv.exists():
            df = pd.read_csv(rec_csv)
            df = df[df["algorithm"] != "nn_mlp"]        # 去掉老的（如果有）
            df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
        else:
            df = pd.DataFrame([row])
        df.to_csv(rec_csv, index=False)
        print(f"已追加到 {rec_csv}")


def main() -> None:
    parser = argparse.ArgumentParser(description="对单个 config 做 NN 重建")
    parser.add_argument("--config-dir", type=Path, required=True, help="形如 Output/.../S0p3_X0_Y0_Z0_..._batch_0001")
    parser.add_argument("--artifacts", type=Path, default=Path("workflow/flow3_NN/artifacts"))
    parser.add_argument("--hidden", type=int, default=256, help="必须和训练一致")
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--append-to-recon", action="store_true",
                        help="把 nn_mlp 结果追加到该 config 的 reconstructed_position.csv")
    parser.add_argument("--photons-per-sample", type=int, default=5000,
                        help="聚合每 N 个光子行为一个样本（和训练时口径一致）")
    parser.add_argument("--normalize-counts", action="store_true",
                        help="转成光子分数（必须和训练时一致）")
    args = parser.parse_args()
    predict(args)


if __name__ == "__main__":
    main()
