#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch sweep for the **double-point** branch (single-process per config).

Each config = TWO source points A and B. The total beamOn (default 35000)
is split between them by a random fraction in [1/3, 2/3] (per-point ratio
1:2 ~ 2:1). The C++ side picks A vs B per primary by `/source/fraction_a`,
so a single `/run/beamOn` covers both points and the standard CSV output
needs no merging.

Two modes:
  --mode random   (default) sample N pairs of random positions.
  --mode explicit pass --point-a x y z --point-b x y z (single config).

Examples:
  # 5 random configs, default cube ±10 mm, default ratio 1:2..2:1
  python workflow/flow4_NN_double/run_batch_double_point.py --random-count 5

  # constrain the two random points to be within 8 mm of each other
  python workflow/flow4_NN_double/run_batch_double_point.py --random-count 100 \\
      --max-separation 8 --min-separation 2 --random-seed 42

  # explicit two-point config (positions in mm)
  python workflow/flow4_NN_double/run_batch_double_point.py --mode explicit \\
      --point-a 3 -2 5 --point-b -4 6 -1 --fraction-a 0.4

Units: all positions and distances are in **mm** (Geant4 internal length unit).
The /source/fp_source* messengers parse raw doubles → interpreted as mm.
"""

from __future__ import annotations

import argparse
import csv
import random
import shutil
import subprocess
import sys
import time
from decimal import Decimal
from pathlib import Path
from typing import List, Optional, Tuple


# ----------------------------- args -----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Two-point source batch (single beamOn per config; ratio 1:2..2:1)."
    )
    p.add_argument("--exe", default=r"build\Release\exampleB1.exe")
    p.add_argument("--geometry-mac", default="geometry.mac")
    p.add_argument("--source-mode", default="optical")
    p.add_argument("--source-distribution", default="Point")

    p.add_argument("--beam-on-total", type=int, default=35000,
                   help="Total beamOn for each config (split A+B in C++).")
    p.add_argument("--ratio-min", type=float, default=1.0 / 3.0,
                   help="Minimum fraction assigned to A (1/3 = ratio 1:2).")
    p.add_argument("--ratio-max", type=float, default=2.0 / 3.0,
                   help="Maximum fraction assigned to A (2/3 = ratio 2:1).")

    p.add_argument("--sigma", default="0.3", help="surfaceSigma fixed value.")

    # Mode selection.
    p.add_argument("--mode", choices=["random", "explicit"], default="random",
                   help="random: sample --random-count pairs; explicit: use --point-a/--point-b.")

    # --- random mode params ---
    p.add_argument("--random-count", type=int, default=10,
                   help="Number of double-point configs in random mode.")
    p.add_argument("--random-min", type=float, default=-10.0, help="Cube min (mm).")
    p.add_argument("--random-max", type=float, default=10.0, help="Cube max (mm).")
    p.add_argument("--random-precision", type=int, default=4)
    p.add_argument("--random-seed", type=int, default=None)
    p.add_argument("--min-separation", type=float, default=0.0,
                   help="Minimum 3D distance between A and B (mm). 0 = no constraint.")
    p.add_argument("--max-separation", type=float, default=0.0,
                   help="Maximum 3D distance between A and B (mm). 0 = no constraint.")
    p.add_argument("--max-sep-uniform", nargs=2, type=float, metavar=("MIN", "MAX"),
                   help="When set, sample per-config max_separation uniformly from [MIN, MAX] mm. "
                        "Overrides --max-separation. Useful for test sets spanning many separation "
                        "regimes in one sweep.")

    # --- explicit mode params ---
    p.add_argument("--point-a", nargs=3, type=float, metavar=("X", "Y", "Z"),
                   help="Explicit point A (mm). Required with --mode explicit.")
    p.add_argument("--point-b", nargs=3, type=float, metavar=("X", "Y", "Z"),
                   help="Explicit point B (mm). Required with --mode explicit.")
    p.add_argument("--fraction-a", type=float, default=None,
                   help="Explicit fraction A in [0,1]. If omitted in explicit mode, "
                        "sampled uniformly from [--ratio-min, --ratio-max].")

    p.add_argument("--prefix-base", default="DP",
                   help="Prefix prepended to every config name.")
    p.add_argument("--results-dir", default="Results",
                   help="Final output root (one subdir per config). HistoManager always "
                        "writes to ./Results/ first; if --results-dir is different, this "
                        "script moves each new subdir there after the run finishes.")
    p.add_argument("--summary-csv", default=None,
                   help="Optional ground-truth summary CSV path. "
                        "Default: <results-dir>/double_point_ground_truth.csv")

    p.add_argument("--dry-run", action="store_true",
                   help="Print plan only, do not run simulation.")
    return p.parse_args()


# ----------------------------- helpers -----------------------------

def fmt_decimal(v: Decimal) -> str:
    out = format(v.normalize(), "f")
    if "." in out:
        out = out.rstrip("0").rstrip(".")
    return out or "0"


def sanitize_token(text: str) -> str:
    return text.replace("-", "m").replace(".", "p")


Triplet = Tuple[Decimal, Decimal, Decimal]


def _quantize(precision: int) -> Decimal:
    return Decimal(10) ** -max(precision, 0)


def random_point(rng: random.Random, lo: float, hi: float, precision: int) -> Triplet:
    q = _quantize(precision)
    return tuple(Decimal(rng.uniform(lo, hi)).quantize(q) for _ in range(3))  # type: ignore[return-value]


def to_triplet(values: List[float], precision: int) -> Triplet:
    q = _quantize(precision)
    return tuple(Decimal(v).quantize(q) for v in values)  # type: ignore[return-value]


def distance(a: Triplet, b: Triplet) -> float:
    return float(sum((float(ai) - float(bi)) ** 2 for ai, bi in zip(a, b)) ** 0.5)


def sample_pair(rng: random.Random, lo: float, hi: float, precision: int,
                min_sep: float, max_sep: float, max_tries: int = 5000) -> Tuple[Triplet, Triplet]:
    """Sample two points in [lo, hi]^3 satisfying separation constraints.

    When max_sep > 0 we use a directed sampler (sample A then B in a spherical
    shell around A) — independent uniform sampling is hopelessly slow for
    small max_sep relative to the cube size (e.g. max_sep=3mm in 25mm cube
    has ~0.7% acceptance).
    """
    q = _quantize(precision)
    eff_min = max(min_sep, 0.0)
    if max_sep > 0 and eff_min > max_sep:
        raise ValueError(f"min_sep ({min_sep}) > max_sep ({max_sep}).")

    if max_sep <= 0:
        # No upper bound → independent uniform sampling is fine.
        for _ in range(max_tries):
            a = random_point(rng, lo, hi, precision)
            b = random_point(rng, lo, hi, precision)
            if eff_min <= 0 or distance(a, b) >= eff_min:
                return a, b
        raise RuntimeError(
            f"Failed to sample two points with min_sep={min_sep} after {max_tries} tries."
        )

    # Directed sampler: A uniform, B = A + r*direction, r∈[eff_min, max_sep].
    for _ in range(max_tries):
        a = random_point(rng, lo, hi, precision)
        # uniform direction on unit sphere via cube-rejection
        ux = uy = uz = 0.0
        for _ in range(50):
            ux = rng.uniform(-1, 1); uy = rng.uniform(-1, 1); uz = rng.uniform(-1, 1)
            mag = (ux * ux + uy * uy + uz * uz) ** 0.5
            if 1e-9 < mag <= 1.0:
                ux, uy, uz = ux / mag, uy / mag, uz / mag
                break
        # uniform-volume radius in shell [eff_min, max_sep]: r = (u·(R³ − r³) + r³)^(1/3)
        u = rng.random()
        r3 = u * (max_sep ** 3 - eff_min ** 3) + eff_min ** 3
        r = r3 ** (1.0 / 3.0)
        bx = float(a[0]) + r * ux
        by = float(a[1]) + r * uy
        bz = float(a[2]) + r * uz
        if not (lo <= bx <= hi and lo <= by <= hi and lo <= bz <= hi):
            continue
        b = (Decimal(bx).quantize(q), Decimal(by).quantize(q), Decimal(bz).quantize(q))
        d = distance(a, b)
        # Quantization can nudge d slightly; accept a small tolerance.
        if d > max_sep + float(q) * 2 or d < eff_min - float(q) * 2:
            continue
        return a, b
    raise RuntimeError(
        f"Failed to sample two points with min_sep={min_sep}, max_sep={max_sep} "
        f"after {max_tries} tries (directed sampler)."
    )


def write_macro(macro_path: Path, geometry_mac: str, sigma: str, prefix: str,
                source_mode: str, source_distribution: str,
                a: Triplet, b: Triplet, fraction_a: float, beam_on: int) -> None:
    macro_path.write_text(
        "\n".join([
            f"/control/execute {geometry_mac}",
            "/control/verbose 0",
            "/run/verbose 1",
            "/event/verbose 0",
            "/tracking/verbose 0",
            "/run/printProgress 100",
            "/vis/disable",
            f"/detector/surfaceSigma {sigma}",
            f"/results/prefix {prefix}",
            "/run/initialize",
            f"/source/mode {source_mode}",
            f"/source/distribution {source_distribution}",
            f"/source/fp_source {fmt_decimal(a[0])} {fmt_decimal(a[1])} {fmt_decimal(a[2])}",
            f"/source/fp_source_b {fmt_decimal(b[0])} {fmt_decimal(b[1])} {fmt_decimal(b[2])}",
            f"/source/fraction_a {fraction_a:.6f}",
            f"/run/beamOn {beam_on}",
            "",
        ]),
        encoding="utf-8",
    )


def latest_dir_with_prefix(results_dir: Path, prefix: str) -> Optional[Path]:
    if not results_dir.exists():
        return None
    cands = [p for p in results_dir.iterdir() if p.is_dir() and p.name.startswith(prefix)]
    if not cands:
        return None
    return max(cands, key=lambda p: p.stat().st_mtime)


def run_one(exe_path: Path, macro_path: Path) -> int:
    return subprocess.run([str(exe_path), str(macro_path)], check=False).returncode


# ----------------------------- main -----------------------------

def build_config_list(args: argparse.Namespace, rng: random.Random) -> List[Tuple[Triplet, Triplet, float]]:
    """Returns a list of (A, B, fraction_a) tuples ready to run."""
    out: List[Tuple[Triplet, Triplet, float]] = []

    if args.mode == "explicit":
        if args.point_a is None or args.point_b is None:
            raise SystemExit("[ERROR] --mode explicit requires --point-a and --point-b.")
        a = to_triplet(args.point_a, args.random_precision)
        b = to_triplet(args.point_b, args.random_precision)
        if args.fraction_a is not None:
            fa = float(args.fraction_a)
            if not 0.0 < fa < 1.0:
                raise SystemExit("[ERROR] --fraction-a must be strictly in (0, 1).")
        else:
            fa = rng.uniform(args.ratio_min, args.ratio_max)
        out.append((a, b, fa))
        return out

    # random mode
    if args.random_count <= 0:
        raise SystemExit("[ERROR] --random-count must be > 0 in random mode.")
    if args.random_max <= args.random_min:
        raise SystemExit("[ERROR] --random-max must be greater than --random-min.")
    if args.max_separation > 0 and args.min_separation > args.max_separation:
        raise SystemExit("[ERROR] --min-separation cannot exceed --max-separation.")

    use_varying = args.max_sep_uniform is not None
    if use_varying:
        ms_lo, ms_hi = args.max_sep_uniform
        if ms_lo <= 0 or ms_hi <= ms_lo:
            raise SystemExit(f"[ERROR] --max-sep-uniform requires 0 < MIN < MAX (got {ms_lo}, {ms_hi}).")
    for _ in range(args.random_count):
        max_sep = rng.uniform(ms_lo, ms_hi) if use_varying else args.max_separation
        a, b = sample_pair(rng, args.random_min, args.random_max, args.random_precision,
                           args.min_separation, max_sep)
        fa = rng.uniform(args.ratio_min, args.ratio_max)
        out.append((a, b, fa))
    return out


def main() -> int:
    args = parse_args()

    if not (0 < args.ratio_min < args.ratio_max < 1):
        print(f"[ERROR] ratio range invalid: [{args.ratio_min}, {args.ratio_max}], "
              f"must satisfy 0 < min < max < 1.")
        return 2

    root = Path.cwd()
    exe_path = (root / args.exe).resolve()
    # HistoManager hardcodes ./Results/ as the output root; we always look there
    # for the freshly-created config dir, then move it to results_dir if different.
    geant4_results_root = (root / "Results").resolve()
    geant4_results_root.mkdir(parents=True, exist_ok=True)
    results_dir = (root / args.results_dir).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    if not exe_path.exists() and not args.dry_run:
        print(f"[ERROR] Executable not found: {exe_path}")
        return 1

    rng = random.Random(args.random_seed)
    configs = build_config_list(args, rng)

    summary_path = Path(args.summary_csv) if args.summary_csv else (results_dir / "double_point_ground_truth.csv")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    write_summary_header = not summary_path.exists()
    summary_fh = summary_path.open("a", newline="")
    summary_w = csv.writer(summary_fh, lineterminator="\n")
    if write_summary_header:
        summary_w.writerow([
            "config_name", "ax_mm", "ay_mm", "az_mm", "bx_mm", "by_mm", "bz_mm",
            "fraction_a", "expected_n_a", "expected_n_b", "separation_mm",
            "sigma", "beam_on_total", "results_subdir",
        ])

    tmp_macro = root / "run_tmp_double_point.mac"

    print("=" * 64)
    print(f"Double-point sweep  mode={args.mode}  configs={len(configs)}")
    print(f"  total beamOn={args.beam_on_total}  ratio∈[{args.ratio_min:.3f}, {args.ratio_max:.3f}]")
    if args.mode == "random":
        sep = []
        if args.min_separation > 0:
            sep.append(f"≥{args.min_separation}")
        if args.max_separation > 0:
            sep.append(f"≤{args.max_separation}")
        sep_str = " AND ".join(sep) if sep else "none"
        print(f"  cube=[{args.random_min}, {args.random_max}] mm  "
              f"sep constraint: {sep_str}  seed={args.random_seed}")
    print(f"  output → {results_dir}")
    print("=" * 64)

    success = 0
    failures: List[str] = []

    try:
        for i, (a, b, fa) in enumerate(configs, start=1):
            sep = distance(a, b)
            n_a_expected = int(round(fa * args.beam_on_total))
            n_b_expected = args.beam_on_total - n_a_expected

            base_name = (
                f"{args.prefix_base}_{i:04d}_S{sanitize_token(args.sigma)}"
                f"_A{sanitize_token(fmt_decimal(a[0]))}_{sanitize_token(fmt_decimal(a[1]))}_{sanitize_token(fmt_decimal(a[2]))}"
                f"_B{sanitize_token(fmt_decimal(b[0]))}_{sanitize_token(fmt_decimal(b[1]))}_{sanitize_token(fmt_decimal(b[2]))}"
            )

            print(f"[{i}/{len(configs)}] {base_name}")
            print(f"   A=({fmt_decimal(a[0])}, {fmt_decimal(a[1])}, {fmt_decimal(a[2])})")
            print(f"   B=({fmt_decimal(b[0])}, {fmt_decimal(b[1])}, {fmt_decimal(b[2])})  d={sep:.3f} mm")
            print(f"   fraction_A={fa:.4f}  expected N_A≈{n_a_expected}  N_B≈{n_b_expected}")

            if args.dry_run:
                continue

            write_macro(tmp_macro, args.geometry_mac, args.sigma, base_name,
                        args.source_mode, args.source_distribution,
                        a, b, fa, args.beam_on_total)

            t0 = time.time()
            rc = run_one(exe_path, tmp_macro)
            time.sleep(0.05)
            if rc != 0:
                print(f"   ! simulation failed (rc={rc})")
                failures.append(base_name)
                continue

            # HistoManager wrote the dir under ./Results/. Find it there and move
            # to results_dir if the user asked for a different destination.
            src_dir = latest_dir_with_prefix(geant4_results_root, base_name)
            if src_dir is None:
                print(f"   ! output dir not found under {geant4_results_root}")
                failures.append(base_name + "(missing)")
                continue
            if results_dir == geant4_results_root:
                out_dir = src_dir
            else:
                dest = results_dir / src_dir.name
                if dest.exists():
                    dest = results_dir / (src_dir.name + "_dup")
                shutil.move(str(src_dir), str(dest))
                out_dir = dest

            summary_w.writerow([
                base_name,
                fmt_decimal(a[0]), fmt_decimal(a[1]), fmt_decimal(a[2]),
                fmt_decimal(b[0]), fmt_decimal(b[1]), fmt_decimal(b[2]),
                f"{fa:.6f}", n_a_expected, n_b_expected, f"{sep:.4f}",
                args.sigma, args.beam_on_total, out_dir.name,
            ])
            summary_fh.flush()

            print(f"   ok ({time.time() - t0:.1f}s) → {out_dir.relative_to(root)}")
            success += 1

    finally:
        summary_fh.close()
        if tmp_macro.exists():
            tmp_macro.unlink()

    print("=" * 64)
    print(f"Done. success={success}/{len(configs)}  failed={len(failures)}")
    for f in failures:
        print(f"  - {f}")
    print(f"Ground truth summary: {summary_path}")
    print("=" * 64)
    return 0 if not failures else 3


if __name__ == "__main__":
    sys.exit(main())
