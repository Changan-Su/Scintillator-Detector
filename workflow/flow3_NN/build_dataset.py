# -*- coding: utf-8 -*-
"""
build_dataset.py —— 从 Geant4 仿真的 Output 目录构建神经网络训练数据集

功能
====
扫描 Output/<SiPM6_Output_*>/<config>/ 下的：
  - metadata.csv       读取真值入射位置 (x, y, z, 单位 mm → cm)
  - merged_event.csv   事件级 SiPM 光子计数（每行一条 SiPM hit）

把每一个 Event 展开成一个 96 维向量（6 面 × 4 × 4 SiPM）作为一个训练样本，
对应的真值位置作为标签。所有样本拼成两个大数组后保存成 .npz。

为什么用"事件级"样本（不是整个 config 聚合）？
    - 同一个 config 里所有事件的真值位置相同（源固定），但每次事件的光子分布不一样
      —— 这正是 NN 能学到的"统计不确定性"的本质。
    - 事件级样本数瞬间从 O(配置数) 变成 O(配置数 × 每配置事件数)，
      从几千膨胀到几百万，NN 才有"素材"去学。

使用方法
========
    uv run python workflow/flow3_NN/build_dataset.py \
        --output-root Output/SiPM6_Output_20260404_234359 \
        --out workflow/flow3_NN/artifacts/dataset.npz

可选：
    --min-photons  每个事件总光子数下限，过滤掉噪声事件（默认 1）
    --max-configs  只处理前 N 个配置目录，方便先小规模跑通（默认全部）
"""

from __future__ import annotations  # 允许用新式类型注解（Python 3.10+ 本来就支持，这行向前兼容）

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# 常量：和 Geant4 侧 DetectorConstruction 里一致的 SiPM 网格参数
# ---------------------------------------------------------------------------
N_FACES = 6                # 立方体 6 个面
N_J = 4                    # 每面 j 方向 4 个 SiPM
N_K = 4                    # 每面 k 方向 4 个 SiPM
N_CHANNELS = N_FACES * N_J * N_K   # = 96，即每个事件的输入维度

# metadata.csv 里真值位置的 Key（你的 metadata 用 Key/Value 两列结构）
POS_KEY_X = "source.fp_source_x_mm"
POS_KEY_Y = "source.fp_source_y_mm"
POS_KEY_Z = "source.fp_source_z_mm"


# ---------------------------------------------------------------------------
# 1. 读真值位置
# ---------------------------------------------------------------------------
def read_true_position_cm(metadata_csv: Path) -> np.ndarray | None:
    """
    从 metadata.csv 解析源位置，返回 shape=(3,) 的 np.float32 数组，单位 cm。
    若关键 Key 缺失则返回 None（调用方决定跳过该 config）。

    metadata.csv 形如：
        Key,Value
        source.fp_source_x_mm,0
        source.fp_source_y_mm,0
        source.fp_source_z_mm,0
        ...
    """
    try:
        df = pd.read_csv(metadata_csv)
    except Exception:
        return None

    # 把两列 DataFrame 转成 {key: value} 字典
    kv = dict(zip(df["Key"].astype(str), df["Value"].astype(str)))

    try:
        x_mm = float(kv[POS_KEY_X])
        y_mm = float(kv[POS_KEY_Y])
        z_mm = float(kv[POS_KEY_Z])
    except (KeyError, ValueError):
        return None

    # 除以 10 从 mm 转 cm，与 Histo10_Cubic 里各重建算法约定一致
    return np.array([x_mm, y_mm, z_mm], dtype=np.float32) / 10.0


# ---------------------------------------------------------------------------
# 2. 把事件 CSV 转成 (N_events, 96) 张量
# ---------------------------------------------------------------------------
def events_to_tensor(
    event_csv: Path,
    min_photons: int = 1,
    max_events: int | None = None,
    rng: np.random.Generator | None = None,
    photons_per_sample: int | None = None,
    normalize_counts: bool = False,
) -> np.ndarray | None:
    """
    读 merged_event.csv，返回 shape=(N_samples, 96) 的 float32 数组。
    列：EventID, CrystalID, iy, iz, Face, j, k, SiPMBlockID, PhotonCount

    两种模式
    --------
    (A) photons_per_sample=None（默认，旧行为）
        每一行 CSV → 一个样本（适合 CSV 本身就是事件级聚合的情况）。

    (B) photons_per_sample=N（>0）
        把这个 config 的所有光子命中**随机分组**成每 N 个一批，
        每批 sum 成一个 96 维样本 —— 用来**把"一个γ事件被 split 成 N 个光子行"
        的数据重新拼回"一个γ事件 = 一个多光子样本"**。

    min_photons：样本总光子数 < 该阈值会被过滤（模式 B 下基本用不上）
    max_events ：最多保留多少样本，None=不限
    rng        ：np.random.Generator，外部传入以保证可复现
    """
    try:
        df = pd.read_csv(event_csv)
    except Exception:
        return None

    if df.empty:
        return None

    # 基本清洗：列名统一、数值列转 int
    required_cols = {"EventID", "Face", "j", "k", "PhotonCount"}
    if not required_cols.issubset(df.columns):
        return None

    # Face/j/k 是整数网格坐标，PhotonCount 可能是浮点 —— 统一成 int
    face = df["Face"].to_numpy(dtype=np.int32)
    j = df["j"].to_numpy(dtype=np.int32)
    k = df["k"].to_numpy(dtype=np.int32)
    counts = df["PhotonCount"].to_numpy(dtype=np.float32)
    event_id_raw = df["EventID"].to_numpy()

    # 合法性：face ∈ [0,5], j ∈ [0,3], k ∈ [0,3]
    valid = (face >= 0) & (face < N_FACES) & (j >= 0) & (j < N_J) & (k >= 0) & (k < N_K)
    if not valid.any():
        return None
    face, j, k, counts, event_id_raw = face[valid], j[valid], k[valid], counts[valid], event_id_raw[valid]

    # 通道下标：face*16 + j*4 + k  → 0..95
    flat_ch = face * (N_J * N_K) + j * N_K + k

    _rng = rng if rng is not None else np.random.default_rng()

    if photons_per_sample is not None and photons_per_sample > 0:
        # ---- 模式 B：把所有光子行随机分组，每 N 个合并成一个样本 ----
        n_rows = len(flat_ch)
        n_samples = n_rows // photons_per_sample
        if n_samples == 0:
            return None

        # 随机打乱光子顺序，再按 N 分块（确保每个样本是随机抽取的一批光子）
        perm = _rng.permutation(n_rows)
        use = perm[: n_samples * photons_per_sample]       # 丢弃余下不足 N 的尾巴
        sample_idx = np.repeat(np.arange(n_samples), photons_per_sample)

        X = np.zeros((n_samples, N_CHANNELS), dtype=np.float32)
        np.add.at(X, (sample_idx, flat_ch[use]), counts[use])
    else:
        # ---- 模式 A：EventID 当成样本 ID（CSV 本来就是事件级时用）----
        event_idx, _unique_ids = pd.factorize(event_id_raw)
        n_events = int(event_idx.max()) + 1 if len(event_idx) > 0 else 0
        if n_events == 0:
            return None

        X = np.zeros((n_events, N_CHANNELS), dtype=np.float32)
        np.add.at(X, (event_idx, flat_ch), counts)

    # 过滤总光子数过低的事件
    if min_photons > 0:
        keep = X.sum(axis=1) >= min_photons
        X = X[keep]

    # 每 config 样本数限流：超过阈值就随机采样
    if max_events is not None and len(X) > max_events:
        pick = _rng.choice(len(X), size=max_events, replace=False)
        X = X[pick]

    # 归一化成光子分数：每行除以总光子数 → 每行 sum=1，模型对光子总数免疫
    # 部署时真实 γ 事件不管有几千还是几万光子，只要分布相似，预测就稳
    if normalize_counts and len(X) > 0:
        totals = X.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1.0
        X = X / totals

    return X if len(X) > 0 else None


# ---------------------------------------------------------------------------
# 3. 遍历所有 config 子目录，组装 (X, y)
# ---------------------------------------------------------------------------
def iter_config_dirs(output_root: Path) -> Iterable[Path]:
    """返回 output_root 下所有直系子目录（跳过非目录文件）。"""
    for p in sorted(output_root.iterdir()):
        if p.is_dir():
            yield p


def build_dataset(
    output_root: Path,
    min_photons: int = 1,
    max_configs: int | None = None,
    max_events_per_config: int | None = None,
    photons_per_sample: int | None = None,
    normalize_counts: bool = False,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    主流程：扫描 → 每个 config 转成若干样本 → 拼接。

    返回：
      X : (N_total, 96)  float32
      y : (N_total, 3)   float32   cm
      names : list[str]  每个样本来自哪个 config（方便后续做分组评估）
    """
    X_parts: list[np.ndarray] = []
    y_parts: list[np.ndarray] = []
    names: list[str] = []
    rng = np.random.default_rng(seed)   # 固定种子，保证每次采样一样

    n_ok, n_skip = 0, 0
    for i, cfg_dir in enumerate(iter_config_dirs(output_root)):
        if max_configs is not None and n_ok >= max_configs:
            break

        meta_csv = cfg_dir / "metadata.csv"
        event_csv = cfg_dir / "merged_event.csv"
        if not (meta_csv.exists() and event_csv.exists()):
            n_skip += 1
            continue

        pos_cm = read_true_position_cm(meta_csv)
        if pos_cm is None:
            n_skip += 1
            continue

        X_cfg = events_to_tensor(
            event_csv,
            min_photons=min_photons,
            max_events=max_events_per_config,
            rng=rng,
            photons_per_sample=photons_per_sample,
            normalize_counts=normalize_counts,
        )
        if X_cfg is None:
            n_skip += 1
            continue

        # 同一个 config 里所有事件共享同一个真值位置
        # np.broadcast_to 不复制内存，只是 view；保险起见 .copy()
        y_cfg = np.broadcast_to(pos_cm, (len(X_cfg), 3)).astype(np.float32).copy()

        X_parts.append(X_cfg)
        y_parts.append(y_cfg)
        names.extend([cfg_dir.name] * len(X_cfg))

        n_ok += 1
        if n_ok % 50 == 0:
            print(f"  [{n_ok}] {cfg_dir.name}  events={len(X_cfg)}")

    if not X_parts:
        raise RuntimeError(f"没有找到有效样本（output_root={output_root}）。")

    X = np.concatenate(X_parts, axis=0)
    y = np.concatenate(y_parts, axis=0)

    print(f"\n汇总：成功 {n_ok} 个 config，跳过 {n_skip} 个。")
    print(f"总样本数 N = {len(X):,}   输入维度 = {X.shape[1]}   标签维度 = {y.shape[1]}")
    print(f"y 范围：x∈[{y[:,0].min():.3f}, {y[:,0].max():.3f}]  "
          f"y∈[{y[:,1].min():.3f}, {y[:,1].max():.3f}]  "
          f"z∈[{y[:,2].min():.3f}, {y[:,2].max():.3f}]  (cm)")
    return X, y, names


# ---------------------------------------------------------------------------
# 4. CLI 入口
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="从 Geant4 Output 目录构建 NN 训练数据集"
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="形如 Output/SiPM6_Output_20260404_234359 的根目录",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("workflow/flow3_NN/artifacts/dataset.npz"),
        help="输出 .npz 路径",
    )
    parser.add_argument("--min-photons", type=int, default=1, help="事件总光子数下限")
    parser.add_argument("--max-configs", type=int, default=None, help="只处理前 N 个 config（调试用）")
    parser.add_argument(
        "--max-events-per-config",
        type=int,
        default=5000,
        help="每个 config 最多保留的事件数（控制内存；默认 5000 ≈ 3 GB / 1800 configs）",
    )
    parser.add_argument("--seed", type=int, default=42, help="随机采样种子")
    parser.add_argument(
        "--photons-per-sample",
        type=int,
        default=None,
        help=("把每 N 个光子命中合并成一个样本，用于'单光子行 → 一个γ事件'的数据。"
              "None = 不合并（每行一个样本）。典型值 500 / 1000 / 2000。"),
    )
    parser.add_argument(
        "--normalize-counts",
        action="store_true",
        help=("每样本除以总光子数 → 转成光子分数（每行 sum=1）。"
              "打开后模型对真实事件的光子总数免疫，推荐部署用。"),
    )
    args = parser.parse_args()

    output_root: Path = args.output_root.resolve()
    if not output_root.is_dir():
        raise SystemExit(f"目录不存在：{output_root}")

    print(f"扫描：{output_root}")
    X, y, names = build_dataset(
        output_root,
        min_photons=args.min_photons,
        max_configs=args.max_configs,
        max_events_per_config=args.max_events_per_config,
        photons_per_sample=args.photons_per_sample,
        normalize_counts=args.normalize_counts,
        seed=args.seed,
    )

    out_path: Path = args.out.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # np.savez_compressed：保存多个数组到一个 .npz 文件，自动用 zip 压缩
    np.savez_compressed(
        out_path,
        X=X,
        y=y,
        config_names=np.array(names, dtype=object),
    )
    print(f"\n已保存：{out_path}")
    print(f"文件大小：{out_path.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
