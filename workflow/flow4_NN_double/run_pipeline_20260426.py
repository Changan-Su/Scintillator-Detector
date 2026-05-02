#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master pipeline for the 2026-04-26 double-point experiment.

9 phases, sequential, ;-style continuation (no abort on phase error):
  1-6. Six sweeps × 1000 configs at max-sep 0.1 / 0.5 / 1 / 3 / 5 / 10 mm
  7.   One big sweep × 6000 configs at max-sep 10 mm (single regime)
  8.   Test sweep × 100 configs with per-config max-sep ∈ [0.01, 10] mm uniform
  9.   Build dataset A (6 mixed sweeps) → artifacts_3/dataset.npz
  10.  Train model A → artifacts_3/best.pt
  11.  Build dataset B (10 mm only 6k) → artifacts_4/dataset.npz
  12.  Train model B → artifacts_4/best.pt
  13.  Predict test set with model A → data/20260426/Output3/
  14.  Predict test set with model B → data/20260426/Output4/

Per-phase log files in repo root: pipeline_20260426_phase<NN>.txt
Master timing log: pipeline_20260426_master.log

Run from project root:
    uv run python workflow/flow4_NN_double/run_pipeline_20260426.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PYTHON = [sys.executable]  # uv-managed interpreter currently running
WF = "workflow/flow4_NN_double"
MASTER_LOG = ROOT / "pipeline_20260426_master.log"


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with MASTER_LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def run_phase(phase: int, name: str, cmd: list[str]) -> int:
    log_file = ROOT / f"pipeline_20260426_phase{phase:02d}.txt"
    log(f"=== Phase {phase}: {name} ===")
    log(f"    cmd: {' '.join(cmd)}")
    log(f"    log: {log_file.name}")
    t0 = time.time()
    with log_file.open("w", encoding="utf-8") as fh:
        rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t0
    status = "OK" if rc == 0 else f"FAIL(rc={rc})"
    log(f"    {status} in {dt/60:.1f} min")
    return rc


def sweep_cmd(count: int, max_sep: float, results_dir: str, seed: int, *, vary=None) -> list[str]:
    cmd = PYTHON + [
        f"{WF}/run_batch_double_point.py",
        "--random-count", str(count),
        "--random-min", "-12.5", "--random-max", "12.5",
        "--beam-on-total", "35000",
        "--results-dir", results_dir,
        "--random-seed", str(seed),
    ]
    if vary is not None:
        cmd += ["--max-sep-uniform", str(vary[0]), str(vary[1])]
    else:
        cmd += ["--max-separation", str(max_sep)]
    return cmd


def build_cmd(input_roots: list[str], out_npz: str) -> list[str]:
    return PYTHON + [
        f"{WF}/build_dataset_double.py",
        "--input-roots", *input_roots,
        "--photons-per-sample", "5000",
        "--normalize-counts",
        "--out", out_npz,
    ]


def train_cmd(dataset: str, artifacts: str) -> list[str]:
    return PYTHON + [
        f"{WF}/train_double.py",
        "--dataset", dataset,
        "--artifacts", artifacts,
    ]


def predict_cmd(input_dir: str, artifacts: str) -> list[str]:
    return PYTHON + [
        f"{WF}/predict_double_batch.py",
        "--input-dir", input_dir,
        "--artifacts", artifacts,
        "--photons-per-sample", "5000",
        "--normalize-counts",
    ]


def copy_outputs(src_dir: Path, dst_dir: Path) -> None:
    """Copy predict_double_batch outputs from input-dir to data/20260426/OutputN/."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    files = [
        "accuracy_summary_double_nn.csv",
        "scatter_double_nn.png",
        "error_vs_separation_double_nn.png",
        "separation_double_nn.png",
    ]
    for fn in files:
        src = src_dir / fn
        if src.exists():
            shutil.copy2(src, dst_dir / fn)
            log(f"    copied {fn} → {dst_dir.name}/")
        else:
            log(f"    [WARN] missing {src}")


def main() -> int:
    MASTER_LOG.write_text("", encoding="utf-8")  # truncate
    pipeline_t0 = time.time()
    log(f"Pipeline 2026-04-26 START   root={ROOT}")

    # ----- Phase 1-6: six 1000-config sweeps at fixed separations -----
    seps_dirs_seeds = [
        (0.1,  "Results_flow4_train1k_sep0p1mm", 2042),
        (0.5,  "Results_flow4_train1k_sep0p5mm", 2142),
        (1.0,  "Results_flow4_train1k_sep1mm",   2242),
        (3.0,  "Results_flow4_train1k_sep3mm",   2342),
        (5.0,  "Results_flow4_train1k_sep5mm",   2442),
        (10.0, "Results_flow4_train1k_sep10mm",  2542),
    ]
    for i, (sep, results, seed) in enumerate(seps_dirs_seeds, start=1):
        run_phase(i, f"sweep 1000 × max-sep {sep} mm → {results}",
                  sweep_cmd(1000, sep, results, seed))

    # ----- Phase 7: 6000-config 10mm-only sweep -----
    big_dir = "Results_flow4_train6k_sep10mm"
    run_phase(7, f"sweep 6000 × max-sep 10 mm → {big_dir}",
              sweep_cmd(6000, 10.0, big_dir, 2742))

    # ----- Phase 8: 100-config test set with varying max-sep -----
    test_dir = "Results_flow4_test100_uniform"
    run_phase(8, f"sweep 100 × max-sep uniform [0.01, 10] mm → {test_dir}",
              sweep_cmd(100, 0.0, test_dir, 2842, vary=(0.01, 10.0)))

    # ----- Phase 9: build dataset A from the 6 mixed sweeps -----
    art3 = ROOT / WF / "artifacts_3"
    art3.mkdir(parents=True, exist_ok=True)
    train_roots_A = [d for _, d, _ in seps_dirs_seeds]
    run_phase(9, "build dataset A (6 mixed seps) → artifacts_3",
              build_cmd(train_roots_A, f"{WF}/artifacts_3/dataset.npz"))

    # ----- Phase 10: train model A -----
    run_phase(10, "train model A → artifacts_3",
              train_cmd(f"{WF}/artifacts_3/dataset.npz", f"{WF}/artifacts_3"))

    # ----- Phase 11: build dataset B from the 6000-config 10mm sweep -----
    art4 = ROOT / WF / "artifacts_4"
    art4.mkdir(parents=True, exist_ok=True)
    run_phase(11, "build dataset B (10 mm 6k) → artifacts_4",
              build_cmd([big_dir], f"{WF}/artifacts_4/dataset.npz"))

    # ----- Phase 12: train model B -----
    run_phase(12, "train model B → artifacts_4",
              train_cmd(f"{WF}/artifacts_4/dataset.npz", f"{WF}/artifacts_4"))

    # ----- Phase 13: predict test set with model A → Output3 -----
    run_phase(13, "predict test set with model A",
              predict_cmd(test_dir, f"{WF}/artifacts_3"))
    copy_outputs(ROOT / test_dir, ROOT / "data" / "20260426" / "Output3")

    # ----- Phase 14: predict test set with model B → Output4 -----
    run_phase(14, "predict test set with model B",
              predict_cmd(test_dir, f"{WF}/artifacts_4"))
    copy_outputs(ROOT / test_dir, ROOT / "data" / "20260426" / "Output4")

    total = (time.time() - pipeline_t0) / 60
    log(f"Pipeline 2026-04-26 DONE in {total:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
