#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
flow4_NN_double 训练数据可视化脚本

读取一个 flow4 双点扫描产物文件夹（包含 100+ 个 DP_*/metadata.csv，或者只有
一份 double_point_ground_truth.csv），把所有 (A, B) 真值位置画到 3D 立方体里。

用法:
    uv run python workflow/flow4_NN_double/visualize_points.py <folder>

例:
    uv run python workflow/flow4_NN_double/visualize_points.py Results_flow4_trainingdata
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib

# 默认 Agg 后端，--no-show 时不弹窗；交互模式下下面 main() 会切回默认后端。
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401, E402  仅用于注册 3D 投影


# ---------------------------------------------------------------------------
# 数据解析
# ---------------------------------------------------------------------------

# metadata.csv 里需要的 key 名
_KEY_AX = "source.fp_source_x_mm"
_KEY_AY = "source.fp_source_y_mm"
_KEY_AZ = "source.fp_source_z_mm"
_KEY_BX = "source.fp_source_b_x_mm"
_KEY_BY = "source.fp_source_b_y_mm"
_KEY_BZ = "source.fp_source_b_z_mm"
_KEY_FRAC = "source.fraction_a"
_KEY_CRYSTAL = "detector.crystalSize_mm"


def _parse_metadata_csv(path: Path) -> Optional[dict]:
    """解析单份 Key,Value 形式的 metadata.csv，返回 dict[str,str]。"""
    try:
        with path.open("r", newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader, None)
            if header != ["Key", "Value"]:
                # 还是尽量兼容，按两列读
                pass
            kv = {}
            for row in reader:
                if len(row) < 2:
                    continue
                kv[row[0].strip()] = row[1].strip()
        return kv
    except OSError as exc:
        print(f"[warn] 无法读取 {path}: {exc}", file=sys.stderr)
        return None


def _extract_pair(kv: dict, source_label: str) -> Optional[Tuple[float, float, float, float, float, float, float]]:
    """从 metadata kv 中抽 (ax, ay, az, bx, by, bz, fraction_a)。失败返回 None。"""
    needed = (_KEY_AX, _KEY_AY, _KEY_AZ, _KEY_BX, _KEY_BY, _KEY_BZ)
    if not all(k in kv for k in needed):
        missing = [k for k in needed if k not in kv]
        print(f"[warn] 跳过 {source_label}: 缺少字段 {missing}", file=sys.stderr)
        return None
    try:
        ax = float(kv[_KEY_AX])
        ay = float(kv[_KEY_AY])
        az = float(kv[_KEY_AZ])
        bx = float(kv[_KEY_BX])
        by = float(kv[_KEY_BY])
        bz = float(kv[_KEY_BZ])
        frac = float(kv.get(_KEY_FRAC, "1.0"))
    except ValueError as exc:
        print(f"[warn] 跳过 {source_label}: 数值解析失败 {exc}", file=sys.stderr)
        return None
    return (ax, ay, az, bx, by, bz, frac)


def load_pairs_from_metadata(folder: Path) -> Tuple[List[tuple], Optional[float]]:
    """优先：扫描 folder/DP_*/metadata.csv。返回 (pairs, crystal_size_mm or None)。"""
    pairs: List[tuple] = []
    crystal_size: Optional[float] = None
    sub_csvs = sorted(folder.glob("DP_*/metadata.csv"))
    if not sub_csvs:
        return pairs, crystal_size

    for csv_path in sub_csvs:
        kv = _parse_metadata_csv(csv_path)
        if kv is None:
            continue
        if crystal_size is None and _KEY_CRYSTAL in kv:
            try:
                crystal_size = float(kv[_KEY_CRYSTAL])
            except ValueError:
                pass
        pair = _extract_pair(kv, csv_path.parent.name)
        if pair is not None:
            pairs.append(pair)
    return pairs, crystal_size


def load_pairs_from_summary(folder: Path) -> List[tuple]:
    """回退：从 double_point_ground_truth.csv 读所有真值。"""
    summary = folder / "double_point_ground_truth.csv"
    if not summary.exists():
        return []
    pairs: List[tuple] = []
    with summary.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                pairs.append((
                    float(row["ax_mm"]),
                    float(row["ay_mm"]),
                    float(row["az_mm"]),
                    float(row["bx_mm"]),
                    float(row["by_mm"]),
                    float(row["bz_mm"]),
                    float(row.get("fraction_a", "1.0")),
                ))
            except (KeyError, ValueError) as exc:
                print(f"[warn] 跳过 summary 一行: {exc}", file=sys.stderr)
    return pairs


# ---------------------------------------------------------------------------
# 绘图
# ---------------------------------------------------------------------------

def _draw_cube_wireframe(ax, half: float, color: str = "gray", lw: float = 0.6) -> None:
    """画一个边长 = 2*half、以原点为中心的立方体线框。"""
    h = half
    corners = [
        (-h, -h, -h), (h, -h, -h), (h, h, -h), (-h, h, -h),
        (-h, -h,  h), (h, -h,  h), (h, h,  h), (-h, h,  h),
    ]
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),  # 底面
        (4, 5), (5, 6), (6, 7), (7, 4),  # 顶面
        (0, 4), (1, 5), (2, 6), (3, 7),  # 立柱
    ]
    for i, j in edges:
        xs = [corners[i][0], corners[j][0]]
        ys = [corners[i][1], corners[j][1]]
        zs = [corners[i][2], corners[j][2]]
        ax.plot(xs, ys, zs, color=color, lw=lw, alpha=0.7)


def _set_cube_limits(ax, half: float) -> None:
    pad = half * 0.05
    ax.set_xlim(-half - pad, half + pad)
    ax.set_ylim(-half - pad, half + pad)
    ax.set_zlim(-half - pad, half + pad)
    try:
        ax.set_box_aspect([1, 1, 1])
    except Exception:
        pass


def plot_3d(pairs: List[tuple], cube_size: float, title: str, alpha: float,
            draw_lines: bool, output_path: Path) -> None:
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    half = cube_size / 2.0
    _draw_cube_wireframe(ax, half)

    ax_pts = [(p[0], p[1], p[2]) for p in pairs]
    bx_pts = [(p[3], p[4], p[5]) for p in pairs]

    if ax_pts:
        xs, ys, zs = zip(*ax_pts)
        ax.scatter(xs, ys, zs, c="crimson", s=14, alpha=alpha, depthshade=True,
                   label=f"A points (n={len(ax_pts)})")
    if bx_pts:
        xs, ys, zs = zip(*bx_pts)
        ax.scatter(xs, ys, zs, c="royalblue", s=14, alpha=alpha, depthshade=True,
                   label=f"B points (n={len(bx_pts)})")
    if draw_lines:
        for p in pairs:
            ax.plot([p[0], p[3]], [p[1], p[4]], [p[2], p[5]],
                    color="gray", lw=0.4, alpha=alpha * 0.6)

    _set_cube_limits(ax, half)
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_zlabel("Z (mm)")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=9)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    print(f"[ok] 3D 图已保存: {output_path}")


def plot_projections(pairs: List[tuple], cube_size: float, title_base: str,
                     alpha: float, draw_lines: bool, base_output: Path) -> List[Path]:
    """生成 XY / XZ / YZ 三张 2D 投影图，文件名在 base_output 基础上加后缀。"""
    half = cube_size / 2.0
    pad = half * 0.05
    saved: List[Path] = []

    plane_specs = [
        ("xy", 0, 1, "X (mm)", "Y (mm)"),
        ("xz", 0, 2, "X (mm)", "Z (mm)"),
        ("yz", 1, 2, "Y (mm)", "Z (mm)"),
    ]

    for tag, ai, bi, xlabel, ylabel in plane_specs:
        fig, ax = plt.subplots(figsize=(7, 7))
        # 立方体边界框
        ax.add_patch(plt.Rectangle((-half, -half), cube_size, cube_size,
                                   fill=False, edgecolor="gray", lw=0.8))

        a_x = [p[ai] for p in pairs]
        a_y = [p[bi] for p in pairs]
        b_x = [p[3 + ai] for p in pairs]
        b_y = [p[3 + bi] for p in pairs]

        if draw_lines:
            for ax_v, ay_v, bx_v, by_v in zip(a_x, a_y, b_x, b_y):
                ax.plot([ax_v, bx_v], [ay_v, by_v],
                        color="gray", lw=0.4, alpha=alpha * 0.6)

        ax.scatter(a_x, a_y, c="crimson", s=16, alpha=alpha,
                   label=f"A points (n={len(pairs)})")
        ax.scatter(b_x, b_y, c="royalblue", s=16, alpha=alpha,
                   label=f"B points (n={len(pairs)})")

        ax.set_xlim(-half - pad, half + pad)
        ax.set_ylim(-half - pad, half + pad)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title_base} - {tag.upper()} projection")
        ax.legend(loc="upper left", fontsize=9)
        ax.grid(True, alpha=0.3)

        out = base_output.with_name(f"{base_output.stem}_{tag}{base_output.suffix}")
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        saved.append(out)
        print(f"[ok] {tag.upper()} 投影已保存: {out}")
        plt.close(fig)

    return saved


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="可视化 flow4_NN_double 扫描结果的 (A, B) 双点真值分布",
    )
    p.add_argument("folder", help="flow4 扫描产物文件夹（含 DP_*/metadata.csv 或 double_point_ground_truth.csv）")
    p.add_argument("--output", default=None,
                   help="输出 PNG 路径（默认：<folder>/double_point_visualization.png）")
    p.add_argument("--no-show", action="store_true", help="不弹出交互窗口（批量 / 无显示环境）")
    p.add_argument("--cube-size", type=float, default=None,
                   help="立方体边长（mm）。默认从 metadata 自动读取，缺失时退回 25")
    p.add_argument("--alpha", type=float, default=0.6, help="散点 / 连线透明度（默认 0.6）")
    p.add_argument("--no-lines", action="store_true", help="不画 A→B 连线")
    p.add_argument("--projections", action="store_true",
                   help="额外保存 XY / XZ / YZ 三张 2D 投影 PNG")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    folder = Path(args.folder)
    if not folder.exists() or not folder.is_dir():
        print(f"[err] 文件夹不存在或不是目录: {folder}", file=sys.stderr)
        return 2

    # 解析数据
    pairs, detected_cube = load_pairs_from_metadata(folder)
    used_path = "DP_*/metadata.csv"
    if not pairs:
        pairs = load_pairs_from_summary(folder)
        used_path = "double_point_ground_truth.csv"
        # summary 模式拿不到 crystal size
    if not pairs:
        print(f"[err] 在 {folder} 没找到任何双点真值数据 "
              f"（既无 DP_*/metadata.csv 也无 double_point_ground_truth.csv）",
              file=sys.stderr)
        return 1

    cube_size = args.cube_size if args.cube_size else (detected_cube or 25.0)

    # 输出路径
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_path = folder / "double_point_visualization.png"

    # no-show 时强制 Agg 后端
    if args.no_show:
        matplotlib.use("Agg", force=True)

    title = f"{folder.name}  |  N = {len(pairs)} pairs  |  cube = {cube_size:g} mm"

    plot_3d(
        pairs=pairs,
        cube_size=cube_size,
        title=title,
        alpha=args.alpha,
        draw_lines=not args.no_lines,
        output_path=out_path,
    )

    if args.projections:
        plot_projections(
            pairs=pairs,
            cube_size=cube_size,
            title_base=folder.name,
            alpha=args.alpha,
            draw_lines=not args.no_lines,
            base_output=out_path,
        )

    print(f"Plotted {len(pairs)} pairs from {folder} (source: {used_path}) -> {out_path}")

    if not args.no_show:
        plt.show()
    plt.close("all")
    return 0


if __name__ == "__main__":
    sys.exit(main())
