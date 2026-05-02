#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild double_point_ground_truth.csv from each DP_*/metadata.csv in a target dir.

Use when the per-config rows didn't get logged during the original sweep
(e.g. because Python looked for output dirs in the wrong root). The C++
metadata.csv has every truth field we need.

Usage:
  python rebuild_summary.py <target_dir>
  # default target_dir = Results_flow4_trainingdata
"""
from __future__ import annotations
import csv
import math
import sys
from pathlib import Path


def parse_metadata(meta_path: Path) -> dict:
    out = {}
    with meta_path.open("r", newline="") as fh:
        for i, row in enumerate(csv.reader(fh)):
            if i == 0 or len(row) < 2:
                continue
            out[row[0]] = row[1]
    return out


def main() -> int:
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("Results_flow4_trainingdata")
    target = target.resolve()
    if not target.exists():
        print(f"[ERROR] {target} not found")
        return 1

    summary_path = target / "double_point_ground_truth.csv"
    rows = []
    bad = []

    dp_dirs = sorted(p for p in target.iterdir() if p.is_dir() and p.name.startswith("DP_"))
    for d in dp_dirs:
        meta = d / "metadata.csv"
        if not meta.exists():
            bad.append(d.name + " (no metadata.csv)")
            continue
        try:
            m = parse_metadata(meta)
            ax, ay, az = float(m["source.fp_source_x_mm"]), float(m["source.fp_source_y_mm"]), float(m["source.fp_source_z_mm"])
            bx, by, bz = float(m["source.fp_source_b_x_mm"]), float(m["source.fp_source_b_y_mm"]), float(m["source.fp_source_b_z_mm"])
            fa = float(m["source.fraction_a"])
            sigma = m.get("detector.surfaceSigma", "")
            # config_name = directory minus trailing _<timestamp>; timestamp is 8+_+6+_+3 (e.g. _20260425_234217_887)
            name = d.name
            parts = name.rsplit("_", 3)  # split off the 3 trailing date/time/ms tokens
            config_name = parts[0] if len(parts) >= 4 else name

            beam_on_total = 35000  # from sweep CLI; not in metadata
            n_a_expected = int(round(fa * beam_on_total))
            n_b_expected = beam_on_total - n_a_expected
            sep = math.sqrt((ax-bx)**2 + (ay-by)**2 + (az-bz)**2)

            rows.append([
                config_name,
                f"{ax:g}", f"{ay:g}", f"{az:g}",
                f"{bx:g}", f"{by:g}", f"{bz:g}",
                f"{fa:.6f}", n_a_expected, n_b_expected, f"{sep:.4f}",
                sigma, beam_on_total, d.name,
            ])
        except Exception as e:
            bad.append(f"{d.name} ({e})")

    with summary_path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow([
            "config_name", "ax_mm", "ay_mm", "az_mm", "bx_mm", "by_mm", "bz_mm",
            "fraction_a", "expected_n_a", "expected_n_b", "separation_mm",
            "sigma", "beam_on_total", "results_subdir",
        ])
        w.writerows(rows)

    print(f"Wrote {summary_path} with {len(rows)} rows ({len(bad)} skipped).")
    for b in bad:
        print(f"  - {b}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
