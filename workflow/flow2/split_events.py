#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
split_events.py
===============
将 Results 中每个配置子目录的大批量仿真事件 CSV 按 EventID 分块，
拆分为若干等大批次，输出到独立子目录，供后续 Histo10_Cubic.py 分析。

典型用途：
  一次跑 10 万个事件，拆成 10 组各 1 万，
  然后对 10 组分别重建，观察统计误差分布。

使用方法：
  python split_events.py Results
  python split_events.py Results --batch-size 10000
  python split_events.py Results --batch-size 10000 --output Results_split
  python split_events.py Results --batch-size 5000 --dry-run

参数：
  results            Results 根目录路径（位置参数或 --results）
  --batch-size N     每批次包含的独立事件数（默认 10000）
  --output PATH      拆分后输出根目录（默认：当前目录下的 Results_split）
  --dry-run          只打印计划，不实际写文件

数据流：
  Results/<config>/AnaEx01_nt_PhotonFaceBlockEvent_t*.csv
    → 合并所有线程 CSV
    → 按 EventID 顺序切分为若干组，每组 batch-size 个独立事件
    → Results_split/<config>_batch_0001/AnaEx01_nt_PhotonFaceBlockEvent.csv
    → Results_split/<config>_batch_0001/metadata.csv   (从原始目录复制)
    → Results_split/<config>_batch_0002/...
    → ...
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

EVENT_COLUMNS = ["EventID", "CrystalID", "iy", "iz", "Face", "j", "k", "SiPMBlockID", "PhotonCount"]
OUTPUT_CSV_NAME = "AnaEx01_nt_PhotonFaceBlockEvent.csv"
METADATA_NAME = "metadata.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="按 EventID 将大批量仿真 CSV 拆分成多个等大批次子目录。"
    )
    parser.add_argument("results", nargs="?", help="Results 根目录路径")
    parser.add_argument("--results", dest="results_alt", help="Results 根目录路径（备用）")
    parser.add_argument(
        "--batch-size", type=int, default=1000,
        help="每批次包含的事件数（默认 1000）"
    )
    parser.add_argument(
        "--output", default="Results_split",
        help="输出根目录（默认 Results_split）"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="只打印计划，不实际写文件"
    )
    return parser.parse_args()


def load_merged_events(config_dir: Path) -> pd.DataFrame:
    """
    合并同一配置目录下所有线程 CSV，返回规范 9 列 DataFrame。
    若无多线程文件则尝试单一合并文件。
    """
    thread_files = sorted(config_dir.glob("AnaEx01_nt_PhotonFaceBlockEvent_t*.csv"))
    frames = []
    for tf in thread_files:
        try:
            df = pd.read_csv(tf, comment="#", header=None, sep=",")
            if df.empty or df.shape[1] < 9:
                continue
            df = df.iloc[:, :9].copy()
            df.columns = EVENT_COLUMNS
            for col in EVENT_COLUMNS:
                df[col] = pd.to_numeric(df[col], errors="coerce")
            df = df.dropna(subset=["EventID"]).copy()
            frames.append(df)
        except Exception as e:
            print(f"  [WARN] Failed to read {tf.name}: {e}")

    if frames:
        return pd.concat(frames, ignore_index=True)

    fallback = config_dir / "AnaEx01_nt_PhotonFaceBlockEvent.csv"
    if fallback.exists():
        df = pd.read_csv(fallback, comment="#", header=None, sep=",")
        if df.empty or df.shape[1] < 9:
            raise ValueError(f"Empty or invalid fallback CSV: {fallback}")
        df = df.iloc[:, :9].copy()
        df.columns = EVENT_COLUMNS
        for col in EVENT_COLUMNS:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.dropna(subset=["EventID"]).copy()
        return df

    raise FileNotFoundError(f"No usable PhotonFaceBlockEvent CSV in {config_dir}")


def split_config(
    config_dir: Path,
    output_root: Path,
    batch_size: int,
    dry_run: bool,
) -> int:
    """
    拆分单个配置目录的事件数据，返回生成的批次数。
    每批次包含 batch_size 个独立 EventID 对应的所有行。
    """
    try:
        df = load_merged_events(config_dir)
    except Exception as e:
        print(f"  [SKIP] {config_dir.name}: {e}")
        return 0

    event_ids = sorted(df["EventID"].dropna().astype(int).unique().tolist())
    n_events = len(event_ids)
    n_batches = (n_events + batch_size - 1) // batch_size

    print(f"  Config: {config_dir.name}")
    print(f"    Total unique events : {n_events}")
    print(f"    Batch size          : {batch_size}")
    print(f"    Number of batches   : {n_batches}")

    for i in range(n_batches):
        start = i * batch_size
        end = min(start + batch_size, n_events)
        batch_event_ids = set(event_ids[start:end])
        batch_name = f"{config_dir.name}_batch_{i + 1:04d}"

        if dry_run:
            print(
                f"    [DRY] {batch_name}/ "
                f"events {event_ids[start]}..{event_ids[end - 1]}, n={end - start}"
            )
            continue

        batch_dir = output_root / batch_name
        batch_dir.mkdir(parents=True, exist_ok=True)

        batch_df = df[df["EventID"].astype(int).isin(batch_event_ids)].copy()
        batch_df.to_csv(batch_dir / OUTPUT_CSV_NAME, index=False)

        meta_src = config_dir / METADATA_NAME
        if meta_src.exists():
            shutil.copy2(meta_src, batch_dir / METADATA_NAME)

        print(
            f"    batch_{i + 1:04d}: events {event_ids[start]}..{event_ids[end - 1]} "
            f"({end - start} events, {len(batch_df)} rows) → {batch_name}/"
        )

    return n_batches


def main() -> int:
    args = parse_args()

    results_path = args.results or args.results_alt
    if not results_path:
        print("[ERROR] 请提供 Results 路径（位置参数或 --results）")
        return 1

    results_dir = Path(results_path).resolve()
    if not results_dir.is_dir():
        print(f"[ERROR] Results 目录不存在: {results_dir}")
        return 1

    if args.batch_size <= 0:
        print("[ERROR] --batch-size 必须为正整数")
        return 1

    output_root = Path(args.output).resolve()
    if not args.dry_run:
        output_root.mkdir(parents=True, exist_ok=True)

    config_dirs = sorted(
        [p for p in results_dir.iterdir() if p.is_dir()],
        key=lambda p: p.name,
    )
    if not config_dirs:
        print(f"[ERROR] 在 {results_dir} 中未找到配置子目录")
        return 1

    print("=" * 60)
    print(f"Results    : {results_dir}")
    print(f"Output     : {output_root}")
    print(f"Batch size : {args.batch_size} events per batch")
    if args.dry_run:
        print("[DRY RUN] No files will be written.")
    print("=" * 60)

    total_batches = 0
    for cfg in config_dirs:
        total_batches += split_config(cfg, output_root, args.batch_size, args.dry_run)
        print()

    print("=" * 60)
    print(f"Done. Total batches: {total_batches}")
    if not args.dry_run:
        print(f"Output directory  : {output_root}")
        print(f"Next step         : python workflow/flow2/Histo10_Cubic.py {output_root.name}")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
