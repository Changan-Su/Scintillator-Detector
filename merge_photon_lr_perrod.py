
# -*- coding: utf-8 -*-
"""
Merge AnaEx01_nt_PhotonLRPerRod_t*.csv into one file.

Usage:
  python merge_photon_lr_perrod.py
  python merge_photon_lr_perrod.py --dir ./out --out merged.csv
  python merge_photon_lr_perrod.py --pattern "AnaEx01_nt_PhotonLRPerRod*.csv" --dedup

Notes:
- Designed for per-rod CSVs that contain columns:
    EventID, iz, iy, Left, Right
- Adds a "__thread__" column parsed from filename suffix "_tN".
- Tolerant to header/no-header. If no header, assumes 5 columns in the order above.
"""

import argparse
import re
from pathlib import Path
import sys
import pandas as pd
import numpy as np

EXPECTED_COLS = ["EventID", "iz", "iy", "Left", "Right"]

def parse_thread_id(name: str):
    m = re.search(r"_t(\d+)\.csv$", name)
    return int(m.group(1)) if m else None

def read_perrod_csv(path: Path) -> pd.DataFrame:
    """
    Robust reader:
      - Try header=0; if columns can be mapped (case-insensitive), use them.
      - Else, read header=None and assume first 5 columns correspond to EXPECTED_COLS.
    """
    # 1) try header=0
    try:
        df = pd.read_csv(path, comment="#", header=0)
        # map columns case-insensitively
        canon = {c.lower(): c for c in df.columns}
        rename = {}
        for want in EXPECTED_COLS:
            lw = want.lower()
            if lw in canon:
                rename[canon[lw]] = want
        if rename:
            df = df.rename(columns=rename)
        # If after rename we still don't have all, fall back
        if not set(EXPECTED_COLS).issubset(df.columns):
            raise ValueError("missing columns after rename")
        df = df[EXPECTED_COLS].copy()
        return df
    except Exception:
        pass

    # 2) header=None fallback
    df = pd.read_csv(path, comment="#", header=None, engine="python")
    if df.shape[1] < 5:
        raise ValueError(f"{path.name}: expected at least 5 columns but got {df.shape[1]}.")
    df = df.iloc[:, :5].copy()
    df.columns = EXPECTED_COLS
    return df

def coerce_types(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["EventID"] = pd.to_numeric(out["EventID"], errors="coerce").astype("Int64")
    out["iz"]      = pd.to_numeric(out["iz"], errors="coerce").astype("Int64")
    out["iy"]      = pd.to_numeric(out["iy"], errors="coerce").astype("Int64")
    out["Left"]    = pd.to_numeric(out["Left"], errors="coerce").astype(float)
    out["Right"]   = pd.to_numeric(out["Right"], errors="coerce").astype(float)
    out = out.dropna(subset=["EventID","iz","iy","Left","Right"]).copy()
    out["EventID"] = out["EventID"].astype(int)
    out["iz"]      = out["iz"].astype(int)
    out["iy"]      = out["iy"].astype(int)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=str, default=".", help="directory to search (default: current directory)")
    ap.add_argument("--pattern", type=str, default="AnaEx01_nt_PhotonLRPerRod*.csv", help="glob pattern of input files")
    ap.add_argument("--out", type=str, default="AnaEx01_nt_PhotonLRPerRod_merged.csv", help="output csv filename")
    ap.add_argument("--dedup", action="store_true", help="drop duplicate rows by (EventID, iz, iy) keeping first")
    ap.add_argument("--sort", action="store_true", help="sort output by (iz, iy, EventID)")
    args = ap.parse_args()

    root = Path(args.dir).resolve()
    files = sorted(root.glob(args.pattern))
    if not files:
        print(f"[ERROR] No files matched: {args.pattern} under {root}", file=sys.stderr)
        sys.exit(1)

    frames = []
    print(f"Found {len(files)} files:")
    for p in files:
        print("  -", p.name)
        try:
            df = read_perrod_csv(p)
        except Exception as e:
            print(f"[WARN] Skip {p.name}: {e}", file=sys.stderr)
            continue
        df = coerce_types(df)
        df["__thread__"] = parse_thread_id(p.name)
        frames.append(df)

    if not frames:
        print("[ERROR] No valid CSVs to merge.", file=sys.stderr)
        sys.exit(2)

    merged = pd.concat(frames, ignore_index=True)

    # optional de-duplication
    if args.dedup:
        before = len(merged)
        merged = merged.drop_duplicates(subset=["EventID","iz","iy"], keep="first")
        print(f"Dedup: {before} -> {len(merged)} rows (by EventID, iz, iy).")

    # optional sort
    if args.sort:
        merged = merged.sort_values(["iz","iy","EventID"], kind="mergesort").reset_index(drop=True)

    out_path = (root / args.out).resolve()
    merged.to_csv(out_path, index=False)
    print(f"\nMerged rows: {len(merged)}")
    print(f"Saved to: {out_path}")

    # Quick per-thread summary
    try:
        summ = merged.groupby("__thread__").size().reset_index(name="rows")
        print("\nRows per thread:")
        for _, r in summ.iterrows():
            print(f"  t{int(r['__thread__']) if pd.notnull(r['__thread__']) else 'NA'} : {int(r['rows'])}")
    except Exception:
        pass

if __name__ == "__main__":
    main()
