#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch sweep script (Python): only source position + surfaceSigma.

Examples:
  python run_batch_sigma_position.py --dry-run
  python run_batch_sigma_position.py --sigma-start 0.3 --sigma-end 0.7 --sigma-step 0.1
  python run_batch_sigma_position.py --x-start 8 --x-end 12 --x-step 2 --beam-on 50000
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep only /detector/surfaceSigma and /source/fp_source(x,y,z)."
    )
    parser.add_argument("--exe", default=r"build\Release\exampleB1.exe", help="Path to exampleB1 executable.")
    parser.add_argument("--geometry-mac", default="geometry.mac", help="Geometry macro to execute first.")
    parser.add_argument("--source-mode", default="optical", help="Value for /source/mode.")
    parser.add_argument("--source-distribution", default="Point", help="Value for /source/distribution.")
    parser.add_argument("--beam-on", type=int, default=100000, help="Value for /run/beamOn.")
    parser.add_argument("--prefix-base", default="", help="Optional fixed prefix before auto sweep suffix.")
    parser.add_argument("--results-dir", default="Results", help="Results directory.")

    parser.add_argument("--sigma-start", default="0.1")
    parser.add_argument("--sigma-end", default="0.9")
    parser.add_argument("--sigma-step", default="0.1")

    parser.add_argument("--x-start", default="10.0")
    parser.add_argument("--x-end", default="10.0")
    parser.add_argument("--x-step", default="1.0")

    parser.add_argument("--y-start", default="-7.5")
    parser.add_argument("--y-end", default="-7.5")
    parser.add_argument("--y-step", default="1.0")

    parser.add_argument("--z-start", default="3.0")
    parser.add_argument("--z-end", default="3.0")
    parser.add_argument("--z-step", default="1.0")

    parser.add_argument("--dry-run", action="store_true", help="Print combinations only, do not run simulation.")
    return parser.parse_args()


def decimal_range(start: str, end: str, step: str) -> List[Decimal]:
    try:
        s = Decimal(start)
        e = Decimal(end)
        st = Decimal(step)
    except InvalidOperation as exc:
        raise ValueError(f"Invalid decimal range: start={start}, end={end}, step={step}") from exc

    if st == 0:
        raise ValueError("Step cannot be zero.")

    values: List[Decimal] = []
    v = s
    eps = Decimal("1e-12")
    if st > 0:
        while v <= e + eps:
            values.append(v)
            v += st
    else:
        while v >= e - eps:
            values.append(v)
            v += st
    return values


def fmt_decimal(v: Decimal) -> str:
    out = format(v.normalize(), "f")
    if "." in out:
        out = out.rstrip("0").rstrip(".")
    return out or "0"


def sanitize_token(text: str) -> str:
    return text.replace("-", "m").replace(".", "p")


def build_prefix(prefix_base: str, sigma: Decimal, x: Decimal, y: Decimal, z: Decimal) -> str:
    s = sanitize_token(fmt_decimal(sigma))
    xs = sanitize_token(fmt_decimal(x))
    ys = sanitize_token(fmt_decimal(y))
    zs = sanitize_token(fmt_decimal(z))
    suffix = f"S{s}_X{xs}_Y{ys}_Z{zs}"
    return f"{prefix_base}_{suffix}" if prefix_base else suffix


def latest_result_dir(results_dir: Path) -> str:
    dirs = [p for p in results_dir.iterdir() if p.is_dir()]
    if not dirs:
        return ""
    return max(dirs, key=lambda p: p.stat().st_mtime).name


def write_macro(
    macro_path: Path,
    geometry_mac: str,
    sigma: Decimal,
    prefix: str,
    source_mode: str,
    source_distribution: str,
    x: Decimal,
    y: Decimal,
    z: Decimal,
    beam_on: int,
) -> None:
    macro_path.write_text(
        "\n".join(
            [
                f"/control/execute {geometry_mac}",
                "/control/verbose 0",
                "/run/verbose 1",
                "/event/verbose 0",
                "/tracking/verbose 0",
                "/run/printProgress 100",
                "/vis/disable",
                f"/detector/surfaceSigma {fmt_decimal(sigma)}",
                f"/results/prefix {prefix}",
                "/run/initialize",
                f"/source/mode {source_mode}",
                f"/source/distribution {source_distribution}",
                f"/source/fp_source {fmt_decimal(x)} {fmt_decimal(y)} {fmt_decimal(z)}",
                f"/run/beamOn {beam_on}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def run_one(exe_path: Path, macro_path: Path) -> int:
    proc = subprocess.run([str(exe_path), str(macro_path)], check=False)
    return proc.returncode


def main() -> int:
    args = parse_args()

    root = Path.cwd()
    exe_path = (root / args.exe).resolve()
    results_dir = (root / args.results_dir).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    if not exe_path.exists():
        print(f"[ERROR] Executable not found: {exe_path}")
        return 1

    sigmas = decimal_range(args.sigma_start, args.sigma_end, args.sigma_step)
    xs = decimal_range(args.x_start, args.x_end, args.x_step)
    ys = decimal_range(args.y_start, args.y_end, args.y_step)
    zs = decimal_range(args.z_start, args.z_end, args.z_step)

    tmp_macro = root / "run_tmp_sigma_position.mac"
    run_count = 0

    print("=" * 60)
    print("Sweep start: surfaceSigma + source fp_source")
    print("=" * 60)

    try:
        for sigma in sigmas:
            for x in xs:
                for y in ys:
                    for z in zs:
                        run_count += 1
                        prefix = build_prefix(args.prefix_base, sigma, x, y, z)
                        print(
                            f"[{run_count}] sigma={fmt_decimal(sigma)} "
                            f"pos=({fmt_decimal(x)}, {fmt_decimal(y)}, {fmt_decimal(z)}) "
                            f"prefix={prefix}"
                        )

                        write_macro(
                            macro_path=tmp_macro,
                            geometry_mac=args.geometry_mac,
                            sigma=sigma,
                            prefix=prefix,
                            source_mode=args.source_mode,
                            source_distribution=args.source_distribution,
                            x=x,
                            y=y,
                            z=z,
                            beam_on=args.beam_on,
                        )

                        if args.dry_run:
                            print("    DRY_RUN=true, skip simulation.")
                            print()
                            continue

                        code = run_one(exe_path, tmp_macro)
                        if code != 0:
                            print(f"    ERROR: simulation failed (exit code {code}).")
                        else:
                            latest = latest_result_dir(results_dir)
                            if latest:
                                print(f"    Latest: Results/{latest}")
                        print()
    finally:
        if tmp_macro.exists():
            tmp_macro.unlink()

    print("=" * 60)
    print(f"Sweep done. Total runs: {run_count}")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())

