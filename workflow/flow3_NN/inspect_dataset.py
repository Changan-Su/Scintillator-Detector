# -*- coding: utf-8 -*-
"""快速检查 dataset.npz 的数据分布，用来诊断训练结果。"""
import argparse
from pathlib import Path
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path,
                        default=Path("workflow/flow3_NN/artifacts/dataset.npz"))
    args = parser.parse_args()

    d = np.load(args.dataset, allow_pickle=True)
    X, y = d["X"], d["y"]

    print(f"样本数: {len(X):,}")
    print(f"输入维度: {X.shape[1]}")

    tot = X.sum(axis=1)
    print(f"\n每事件光子总数:")
    print(f"  min    = {tot.min():.1f}")
    print(f"  5%     = {np.percentile(tot, 5):.1f}")
    print(f"  median = {np.median(tot):.1f}")
    print(f"  mean   = {tot.mean():.1f}")
    print(f"  95%    = {np.percentile(tot, 95):.1f}")
    print(f"  max    = {tot.max():.1f}")

    nnz = (X > 0).sum(axis=1)
    print(f"\n每事件非零 SiPM 数（共 96 个）:")
    print(f"  median = {np.median(nnz):.1f}   max = {nnz.max()}")

    print(f"\ny 分布（cm）:")
    for i, name in enumerate("xyz"):
        col = y[:, i]
        print(f"  {name}: min={col.min():.3f}  max={col.max():.3f}  "
              f"mean={col.mean():.3f}  std={col.std():.3f}")

    # 唯一真值位置数
    uniq = np.unique(np.round(y, 3), axis=0)
    print(f"\n唯一真值位置数（round 到 0.001 cm）: {len(uniq)}")
    if len(uniq) < 30:
        print("前 20 个真值位置:")
        for pos in uniq[:20]:
            print(f"  ({pos[0]:+.3f}, {pos[1]:+.3f}, {pos[2]:+.3f})")


if __name__ == "__main__":
    main()
