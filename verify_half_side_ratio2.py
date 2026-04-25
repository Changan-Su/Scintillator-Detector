#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_half_side_ratio2.py
==========================
Pure-math verification of the proposed reconstruct_half_side_ratio2 algorithm.

Uses exact solid-angle computation on a unit cube [0,1]^3 (6 faces, 4×4 cells
per face = 96 cells) to simulate noiseless SiPM photon counts.

Compares three reconstruction approaches:
  1. linear_scaled     — face-total ratio × half-length (baseline)
  2. half_side_ratio   — fixed X/Y/Z half-space split + arccos/arctan
  3. half_side_ratio2  — adaptive axis (weighted centroid → split → arccos/arctan)

Usage:
  python verify_half_side_ratio2.py
"""

import numpy as np

# ── Solid-angle kernel (from cube_azimuth_interactive.py) ──────────────

def rect_solid_angle(x1, x2, y1, y2, d):
    """Exact solid angle of an axis-aligned rectangle at normal distance |d|."""
    d = abs(d)
    if d < 1e-15:
        return 0.0

    def f(x, y):
        return np.arctan2(x * y, d * np.sqrt(x * x + y * y + d * d))

    return abs(f(x2, y2) - f(x1, y2) - f(x2, y1) + f(x1, y1))


# ── Build cell data for 6-face cube ───────────────────────────────────

def build_cells(sx, sy, sz, n=4):
    """
    Return list of (x, y, z, omega) for all 6·n·n cells.
    Cube is [0,1]^3, source at (sx, sy, sz).
    """
    edges = np.linspace(0.0, 1.0, n + 1)
    centers = (edges[:-1] + edges[1:]) / 2.0
    cells = []

    for j in range(n):
        for k in range(n):
            u1, u2 = edges[k], edges[k + 1]
            v1, v2 = edges[j], edges[j + 1]
            cu, cv = centers[k], centers[j]

            # x=0 face (-X): u→y, v→z, distance=sx
            om = rect_solid_angle(u1 - sy, u2 - sy, v1 - sz, v2 - sz, sx)
            cells.append((0.0, cu, cv, om))
            # x=1 face (+X): u→y, v→z, distance=1-sx
            om = rect_solid_angle(u1 - sy, u2 - sy, v1 - sz, v2 - sz, 1 - sx)
            cells.append((1.0, cu, cv, om))
            # y=0 face (-Y): u→x, v→z, distance=sy
            om = rect_solid_angle(u1 - sx, u2 - sx, v1 - sz, v2 - sz, sy)
            cells.append((cu, 0.0, cv, om))
            # y=1 face (+Y): u→x, v→z, distance=1-sy
            om = rect_solid_angle(u1 - sx, u2 - sx, v1 - sz, v2 - sz, 1 - sy)
            cells.append((cu, 1.0, cv, om))
            # z=0 face (-Z): u→x, v→y, distance=sz
            om = rect_solid_angle(u1 - sx, u2 - sx, v1 - sy, v2 - sy, sz)
            cells.append((cu, cv, 0.0, om))
            # z=1 face (+Z): u→x, v→y, distance=1-sz
            om = rect_solid_angle(u1 - sx, u2 - sx, v1 - sy, v2 - sy, 1 - sz)
            cells.append((cu, cv, 1.0, om))

    return cells


# ── Reconstruction algorithms ─────────────────────────────────────────

CENTER = 0.5
HALF = 0.5


def _arctan_pos(n_plus, n_minus, half_len):
    """arccos/arctan ratio → displacement from center."""
    denom = n_plus + n_minus
    if denom <= 0:
        return 0.0
    ratio = np.clip(n_plus / denom, 1e-12, 1 - 1e-12)
    angle = np.arccos(1.0 - 2.0 * ratio)
    t = np.tan(angle)
    if abs(t) < 1e-12:
        return np.sign(n_plus - n_minus) * half_len * 0.99
    return -half_len / t


def reconstruct_linear_scaled(cells):
    """Face-total ratio × half-length, each axis independent."""
    face_sums = {"x0": 0.0, "x1": 0.0, "y0": 0.0, "y1": 0.0, "z0": 0.0, "z1": 0.0}
    for x, y, z, om in cells:
        if x < 0.01:
            face_sums["x0"] += om
        elif x > 0.99:
            face_sums["x1"] += om
        if y < 0.01:
            face_sums["y0"] += om
        elif y > 0.99:
            face_sums["y1"] += om
        if z < 0.01:
            face_sums["z0"] += om
        elif z > 0.99:
            face_sums["z1"] += om

    def axis_rec(np_, nm_):
        d = np_ + nm_
        return HALF * (np_ - nm_) / d if d > 0 else 0.0

    return (
        CENTER + axis_rec(face_sums["x1"], face_sums["x0"]),
        CENTER + axis_rec(face_sums["y1"], face_sums["y0"]),
        CENTER + axis_rec(face_sums["z1"], face_sums["z0"]),
    )


def reconstruct_half_side_ratio(cells):
    """Fixed-axis half-space split + arccos/arctan, each axis independent."""
    n_xp, n_xm = 0.0, 0.0
    n_yp, n_ym = 0.0, 0.0
    n_zp, n_zm = 0.0, 0.0

    for x, y, z, om in cells:
        if x >= CENTER:
            n_xp += om
        else:
            n_xm += om
        if y >= CENTER:
            n_yp += om
        else:
            n_ym += om
        if z >= CENTER:
            n_zp += om
        else:
            n_zm += om

    return (
        CENTER + _arctan_pos(n_xp, n_xm, HALF),
        CENTER + _arctan_pos(n_yp, n_ym, HALF),
        CENTER + _arctan_pos(n_zp, n_zm, HALF),
    )


def reconstruct_half_side_ratio2(cells):
    """Adaptive axis: weighted centroid → split → arccos/arctan."""
    total = sum(om for _, _, _, om in cells)
    if total <= 0:
        return CENTER, CENTER, CENTER

    # Step 1: photon-weighted centroid → axis direction
    wx = sum(x * om for x, y, z, om in cells) / total
    wy = sum(y * om for x, y, z, om in cells) / total
    wz = sum(z * om for x, y, z, om in cells) / total

    d = np.array([wx - CENTER, wy - CENTER, wz - CENTER])
    d_len = np.linalg.norm(d)
    if d_len < 1e-15:
        return CENTER, CENTER, CENTER
    u = d / d_len

    # Step 2: project all cells onto axis, split ±
    n_plus, n_minus = 0.0, 0.0
    for x, y, z, om in cells:
        proj = (x - CENTER) * u[0] + (y - CENTER) * u[1] + (z - CENTER) * u[2]
        if proj >= 0:
            n_plus += om
        else:
            n_minus += om

    # Step 3: effective half-length along u (ray–box intersection)
    comps = []
    if abs(u[0]) > 1e-12:
        comps.append(HALF / abs(u[0]))
    if abs(u[1]) > 1e-12:
        comps.append(HALF / abs(u[1]))
    if abs(u[2]) > 1e-12:
        comps.append(HALF / abs(u[2]))
    half_len = min(comps) if comps else HALF

    r = _arctan_pos(n_plus, n_minus, half_len)

    # Step 4: map back to xyz
    rec = np.array([CENTER, CENTER, CENTER]) + r * u
    return float(rec[0]), float(rec[1]), float(rec[2])


# ── Helpers ────────────────────────────────────────────────────────────

def err3d(rec, true):
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(rec, true)))


def classify(sx, sy, sz):
    off = sum(1 for c in (sx, sy, sz) if abs(c - 0.5) > 0.01)
    if off <= 1:
        return "on-axis"
    elif off == 2:
        return "off-2D"
    else:
        return "off-3D"


# ── Main ───────────────────────────────────────────────────────────────

def main():
    n = 4

    # Build test positions
    vals = [0.25, 0.35, 0.5, 0.65, 0.75]
    positions = set()

    # On-axis
    for v in vals:
        positions.add((v, 0.5, 0.5))
        positions.add((0.5, v, 0.5))
        positions.add((0.5, 0.5, v))

    # Off-axis 2D & 3D
    off_vals = [0.3, 0.5, 0.7]
    for a in off_vals:
        for b in off_vals:
            for c in off_vals:
                if not (a == 0.5 and b == 0.5 and c == 0.5):
                    positions.add((a, b, c))

    # More extreme off-axis
    positions.add((0.75, 0.75, 0.75))
    positions.add((0.25, 0.25, 0.25))
    positions.add((0.75, 0.25, 0.75))
    positions.add((0.6, 0.7, 0.8))
    positions.add((0.35, 0.65, 0.75))

    positions = sorted(positions)

    # Compute
    results = []
    for sx, sy, sz in positions:
        cells = build_cells(sx, sy, sz, n)
        r_lin = reconstruct_linear_scaled(cells)
        r_fix = reconstruct_half_side_ratio(cells)
        r_adp = reconstruct_half_side_ratio2(cells)

        results.append({
            "true": (sx, sy, sz),
            "cat": classify(sx, sy, sz),
            "lin": r_lin,
            "fix": r_fix,
            "adp": r_adp,
            "e_lin": err3d(r_lin, (sx, sy, sz)),
            "e_fix": err3d(r_fix, (sx, sy, sz)),
            "e_adp": err3d(r_adp, (sx, sy, sz)),
        })

    # Print detailed table per category
    for cat_name in ["on-axis", "off-2D", "off-3D"]:
        group = [r for r in results if r["cat"] == cat_name]
        if not group:
            continue

        print(f"\n{'=' * 30} {cat_name} ({len(group)} points) {'=' * 30}")
        print(
            f"  {'True position':>22s}   "
            f"{'err(linear)':>11s}  {'err(fix-axis)':>13s}  {'err(adaptive)':>13s}   note"
        )
        print("  " + "-" * 90)

        for r in group:
            sx, sy, sz = r["true"]
            flag = ""
            if r["e_adp"] < r["e_fix"] * 0.95:
                flag = "adaptive wins"
            elif r["e_fix"] < r["e_adp"] * 0.95:
                flag = "fixed wins"
            print(
                f"  ({sx:.2f}, {sy:.2f}, {sz:.2f})   "
                f"{r['e_lin']:11.6f}  {r['e_fix']:13.6f}  {r['e_adp']:13.6f}   {flag}"
            )

        e_lin = np.mean([r["e_lin"] for r in group])
        e_fix = np.mean([r["e_fix"] for r in group])
        e_adp = np.mean([r["e_adp"] for r in group])
        print(f"  {'MEAN':>22s}   {e_lin:11.6f}  {e_fix:13.6f}  {e_adp:13.6f}")

    # Overall summary
    print(f"\n{'=' * 30} Overall Summary ({len(results)} points) {'=' * 30}")
    for name, key in [("linear_scaled", "e_lin"), ("half_side_ratio", "e_fix"), ("half_side_ratio2", "e_adp")]:
        errs = [r[key] for r in results]
        print(f"  {name:20s}  mean={np.mean(errs):.6f}  max={np.max(errs):.6f}  std={np.std(errs):.6f}")

    # On-axis vs off-axis breakdown
    print(f"\n{'=' * 30} On-axis vs Off-axis {'=' * 30}")
    for cat_name in ["on-axis", "off-2D", "off-3D"]:
        group = [r for r in results if r["cat"] == cat_name]
        if not group:
            continue
        e_fix = np.mean([r["e_fix"] for r in group])
        e_adp = np.mean([r["e_adp"] for r in group])
        ratio = e_adp / e_fix if e_fix > 0 else float("inf")
        print(f"  {cat_name:8s}  fix_mean={e_fix:.6f}  adp_mean={e_adp:.6f}  adp/fix={ratio:.3f}")


if __name__ == "__main__":
    main()
