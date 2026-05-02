# -*- coding: utf-8 -*-
"""
build_dataset_double.py —— 双点光源训练数据集构建脚本

功能
====
扫描一个或多个 flow4 双点光源结果根目录（例如 Results_flow4_trainingdata*/）
下的每个 DP_* 子目录：
  - metadata.csv                              读取 (A, B) 两点真值位置 + fraction_a
  - AnaEx01_nt_PhotonFaceBlockEvent_t*.csv    多线程 per-event 光子命中（无表头、# 注释）

把每个 config 内的所有光子命中行随机分组，每 N（= --photons-per-sample）个光子
合并成一个样本：X (96 维) + y (6 维：ax, ay, az, bx, by, bz, 单位 mm)。
不跨 config 混合 —— 每个聚合样本严格属于一个 config，标签明确。

为什么单位用 mm（而不是 cm）？
    flow4 的 metadata 字段后缀就是 _mm，README 里也明确"所有坐标都是 mm"，
    所以 flow4_NN_double 整条管线统一用 mm，省去无谓的 ×10 / ÷10 转换。
    （flow3_NN 用 cm 是历史遗留，不要混用。）

用法
====
    # 单个根目录
    uv run python workflow/flow4_NN_double/build_dataset_double.py \
        --input-roots Results_flow4_trainingdata \
        --photons-per-sample 2000 --normalize-counts

    # 多个根目录一起训练（覆盖不同 separation regime）
    uv run python workflow/flow4_NN_double/build_dataset_double.py \
        --input-roots Results_flow4_trainingdata \
                      Results_flow4_trainingdata_sep5mm \
                      Results_flow4_trainingdata_sep3mm \
        --photons-per-sample 2000 --normalize-counts
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# 与 Geant4 / flow3 一致的常量
# ---------------------------------------------------------------------------
N_FACES = 6
N_J = 4
N_K = 4
N_CHANNELS = N_FACES * N_J * N_K   # 96

EVENT_COLUMNS = [
    "EventID", "CrystalID", "iy", "iz", "Face", "j", "k", "SiPMBlockID", "PhotonCount",
]

POS_KEYS_A = ("source.fp_source_x_mm",   "source.fp_source_y_mm",   "source.fp_source_z_mm")
POS_KEYS_B = ("source.fp_source_b_x_mm", "source.fp_source_b_y_mm", "source.fp_source_b_z_mm")
KEY_FRACTION_A = "source.fraction_a"


# ---------------------------------------------------------------------------
# 1. 真值（A, B, fraction_a）解析
# ---------------------------------------------------------------------------
def read_true_positions_mm(metadata_csv: Path) -> tuple[np.ndarray, float] | None:
    """从 metadata.csv 读取 (A, B) 位置（mm）和 fraction_a。

    返回 (y6, frac_a) 其中 y6=[ax, ay, az, bx, by, bz] 都是 mm；任一字段缺失返回 None。
    """
    try:
        df = pd.read_csv(metadata_csv)
    except Exception:
        return None
    kv = dict(zip(df["Key"].astype(str), df["Value"].astype(str)))
    try:
        a = [float(kv[k]) for k in POS_KEYS_A]
        b = [float(kv[k]) for k in POS_KEYS_B]
        frac_a = float(kv.get(KEY_FRACTION_A, "1.0"))
    except (KeyError, ValueError):
        return None
    y6 = np.array(a + b, dtype=np.float32)  # mm
    return y6, frac_a


# ---------------------------------------------------------------------------
# 2. 多线程 per-event CSV 直接读取并聚合（不依赖 merged_event.csv）
# ---------------------------------------------------------------------------
def _parse_event_csv(csv_path: Path) -> pd.DataFrame | None:
    """读单份 PhotonFaceBlockEvent CSV（无表头、# 注释行）→ 9 列 DataFrame。"""
    try:
        df = pd.read_csv(csv_path, comment="#", header=None, sep=",")
    except Exception:
        return None
    if df.empty or df.shape[1] < 9:
        return None
    df = df.iloc[:, :9].copy()
    df.columns = EVENT_COLUMNS
    for col in EVENT_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["Face", "j", "k", "PhotonCount"]).copy()
    return df


def load_config_photons(cfg_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """合并 DP_* 目录下所有 _t*.csv，返回 (face, jk_flat_channel, count)。

    返回每行一个光子命中的扁平数组（int32 / float32）。无文件返回 None。
    """
    thread_files = sorted(cfg_dir.glob("AnaEx01_nt_PhotonFaceBlockEvent_t*.csv"))
    if not thread_files:
        fallback = cfg_dir / "AnaEx01_nt_PhotonFaceBlockEvent.csv"
        if fallback.exists():
            thread_files = [fallback]
    if not thread_files:
        return None

    frames = []
    for tf in thread_files:
        df = _parse_event_csv(tf)
        if df is not None and not df.empty:
            frames.append(df)
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)

    face = df["Face"].to_numpy(dtype=np.int32)
    j = df["j"].to_numpy(dtype=np.int32)
    k = df["k"].to_numpy(dtype=np.int32)
    counts = df["PhotonCount"].to_numpy(dtype=np.float32)
    valid = (face >= 0) & (face < N_FACES) & (j >= 0) & (j < N_J) & (k >= 0) & (k < N_K)
    if not valid.any():
        return None
    face, j, k, counts = face[valid], j[valid], k[valid], counts[valid]
    flat_ch = face * (N_J * N_K) + j * N_K + k
    return flat_ch.astype(np.int32), counts, face  # face 仅返回作占位，未必使用


def aggregate_to_samples(
    flat_ch: np.ndarray,
    counts: np.ndarray,
    photons_per_sample: int,
    rng: np.random.Generator,
    normalize_counts: bool,
) -> np.ndarray | None:
    """把扁平光子命中按 N 个一组合并成 (N_samples, 96) 的样本张量。"""
    n_rows = len(flat_ch)
    n_samples = n_rows // photons_per_sample
    if n_samples <= 0:
        return None

    perm = rng.permutation(n_rows)
    use = perm[: n_samples * photons_per_sample]
    sample_idx = np.repeat(np.arange(n_samples), photons_per_sample)

    X = np.zeros((n_samples, N_CHANNELS), dtype=np.float32)
    np.add.at(X, (sample_idx, flat_ch[use]), counts[use])

    if normalize_counts:
        totals = X.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1.0
        X = X / totals
    return X


# ---------------------------------------------------------------------------
# 3. 遍历多个根目录
# ---------------------------------------------------------------------------
def iter_dp_dirs(roots: Iterable[Path]) -> Iterable[Path]:
    for root in roots:
        if not root.is_dir():
            print(f"  [WARN] not a dir, skip: {root}")
            continue
        for child in sorted(root.iterdir(), key=lambda p: p.name):
            if not child.is_dir():
                continue
            # 兼容: 直接是 DP_*，或子目录里有 metadata.csv
            if (child / "metadata.csv").exists():
                yield child


def build(
    roots: List[Path],
    photons_per_sample: int,
    normalize_counts: bool,
    max_configs: int | None,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)

    X_parts: list[np.ndarray] = []
    y_parts: list[np.ndarray] = []
    cfg_ids: list[np.ndarray] = []
    frac_a_per_sample: list[np.ndarray] = []
    seps_per_sample: list[np.ndarray] = []
    cfg_names: list[str] = []

    n_ok = 0
    n_skip = 0
    for cfg_dir in iter_dp_dirs(roots):
        if max_configs is not None and n_ok >= max_configs:
            break

        meta = read_true_positions_mm(cfg_dir / "metadata.csv")
        if meta is None:
            n_skip += 1
            continue
        y6_mm, frac_a = meta

        loaded = load_config_photons(cfg_dir)
        if loaded is None:
            n_skip += 1
            continue
        flat_ch, counts, _ = loaded

        X_cfg = aggregate_to_samples(
            flat_ch, counts,
            photons_per_sample=photons_per_sample,
            rng=rng,
            normalize_counts=normalize_counts,
        )
        if X_cfg is None or len(X_cfg) == 0:
            n_skip += 1
            continue

        n_s = len(X_cfg)
        y_cfg = np.broadcast_to(y6_mm, (n_s, 6)).astype(np.float32).copy()
        sep = float(np.linalg.norm(y6_mm[:3] - y6_mm[3:]))
        cfg_id_idx = len(cfg_names)
        cfg_names.append(cfg_dir.name)

        X_parts.append(X_cfg)
        y_parts.append(y_cfg)
        cfg_ids.append(np.full(n_s, cfg_id_idx, dtype=np.int32))
        frac_a_per_sample.append(np.full(n_s, frac_a, dtype=np.float32))
        seps_per_sample.append(np.full(n_s, sep, dtype=np.float32))

        n_ok += 1
        if n_ok % 25 == 0:
            print(f"  [{n_ok}] {cfg_dir.name}  samples={n_s}  sep={sep:.2f}mm  fA={frac_a:.3f}")

    if not X_parts:
        raise RuntimeError("没有有效样本，请检查 --input-roots 是否包含 DP_* 子目录。")

    X = np.concatenate(X_parts, axis=0)
    y = np.concatenate(y_parts, axis=0)
    config_id = np.concatenate(cfg_ids, axis=0)
    fraction_a = np.concatenate(frac_a_per_sample, axis=0)
    separations = np.concatenate(seps_per_sample, axis=0)

    print(f"\n汇总：成功 {n_ok} 个 config，跳过 {n_skip} 个。")
    print(f"总样本 N = {len(X):,}   X.shape = {X.shape}   y.shape = {y.shape}")
    print(f"separation 分布: min={separations.min():.2f}  median={np.median(separations):.2f}  "
          f"max={separations.max():.2f}  (mm)")
    print(f"fraction_a 分布: min={fraction_a.min():.3f}  median={np.median(fraction_a):.3f}  "
          f"max={fraction_a.max():.3f}")

    return dict(
        X=X, y=y,
        config_id=config_id,
        fraction_a=fraction_a,
        separations=separations,
        config_names=np.array(cfg_names, dtype=object),
    )


# ---------------------------------------------------------------------------
# 4. CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="flow4 双点光源训练数据集构建")
    parser.add_argument(
        "--input-roots", nargs="+", required=True, type=Path,
        help="一个或多个 Results_flow4_trainingdata* 根目录（每个含若干 DP_* 子目录）",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("workflow/flow4_NN_double/artifacts/dataset.npz"),
        help="输出 .npz 路径",
    )
    parser.add_argument("--photons-per-sample", type=int, default=2000,
                        help="每 N 个光子命中合并成 1 个样本（推荐 1000-5000）")
    parser.add_argument("--normalize-counts", action="store_true",
                        help="每个样本除以总光子数（推荐打开，模型对总光子数免疫）")
    parser.add_argument("--max-configs", type=int, default=None,
                        help="只处理前 N 个 config（调试用，跨所有 root 累计）")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    roots = [p.resolve() for p in args.input_roots]
    print("扫描根目录：")
    for r in roots:
        print(f"  - {r}")

    pack = build(
        roots=roots,
        photons_per_sample=args.photons_per_sample,
        normalize_counts=args.normalize_counts,
        max_configs=args.max_configs,
        seed=args.seed,
    )

    out_path: Path = args.out.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, **pack)
    print(f"\n已保存：{out_path}")
    print(f"文件大小：{out_path.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
