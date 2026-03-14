# -*- coding: utf-8 -*-
"""
python Histo9.py .\Results --doi-mode diff_over_sum


Histo9: 连续晶体稳健版 DOI 分析脚本。

与 Histo8 的主要区别：
1) 峰检测基于“未加权计数直方图”，避免 L+R 权重放大高光子事件。
2) 默认不强制目标峰数，不做“降阈值凑峰”。
3) 支持 total(L+R) 事件门限，减少低统计噪声点导致的假峰。
4) 图中同时展示归一化计数曲线与归一化加权曲线，便于对比。

使用方法
--------
  python Histo9.py <ResultsDir> [选项]
  python Histo9.py --results <ResultsDir> [选项]

参数
----
  results               Results 文件夹路径（位置参数或 --results）
  --results PATH        Results 文件夹路径（与位置参数二选一）
  --output PATH         输出文件夹；不指定则默认为当前目录下 Histo9_Output_<时间戳>

峰检测与直方图
--------------
  --n-bins N            直方图 bin 数（默认 300）
  --sigma-smooth F      高斯平滑 sigma（默认 1.2）
  --distance-ratio F    峰最小间距占比，相对 bin 数（默认 0.04）
  --prominence-ratio F  峰 prominence 比例，相对平滑峰高（默认 0.015）
  --min-prominence-abs F  峰 prominence 绝对下限（计数，默认 5.0）
  --min-width-bins N    峰最小宽度（bin，默认 2）
  --max-peaks N         保留峰数量上限；0 表示不限（默认 0）

事件与 DOI
----------
  --min-total-photons F 事件门限：仅保留 L+R >= 此值（默认 10.0）
  --min-events N        单 rod 最少事件数（过滤后，默认 30）
  --doi-mode MODE       DOI 定义：r_over_sum | r_over_l | diff_over_sum | log_r_over_l（默认 r_over_sum）
  --doi-min F           DOI 直方图区间下限（不指定则按 doi-mode 取默认）
  --doi-max F           DOI 直方图区间上限（不指定则按 doi-mode 取默认）

示例
----
  python Histo9.py Results
  python Histo9.py Results --output ./Histo9_Out --doi-mode r_over_sum --n-bins 400
  python Histo9.py --results Results --min-total-photons 20 --min-events 50
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths


def parse_geometry_mac(mac_path: Path) -> Dict[str, float]:
    config: Dict[str, float] = {}
    if not mac_path.exists():
        return config
    with open(mac_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            match = re.match(r"/detector/(\w+)\s+([\d.]+)", line)
            if not match:
                continue
            key, value = match.group(1), match.group(2)
            try:
                config[key] = float(value)
            except ValueError:
                pass
    return config


def merge_thread_csvs(folder: Path, pattern: str = "AnaEx01_nt_PhotonLRPerRod") -> Optional[Path]:
    merged_path = folder / f"{pattern}.csv"
    if merged_path.exists():
        try:
            with open(merged_path, "r", encoding="utf-8", errors="ignore") as f:
                has_data = any(line.strip() and not line.startswith("#") for line in f.readlines())
            if has_data:
                return merged_path
        except Exception:
            pass

    thread_files = sorted(folder.glob(f"{pattern}_t*.csv"))
    if not thread_files:
        return None

    frames: List[pd.DataFrame] = []
    print(f"  [Merge] Found {len(thread_files)} thread files in {folder.name}")
    for tf in thread_files:
        try:
            df = pd.read_csv(tf, comment="#", header=None, sep=",")
            if df.shape[1] >= 5:
                df = df.iloc[:, :5].copy()
                df.columns = ["EventID", "iz", "iy", "Left", "Right"]
                frames.append(df)
        except Exception as e:
            print(f"  [Warn] Failed to read {tf.name}: {e}")

    if not frames:
        return None

    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(merged_path, index=False)
    print(f"  [Merge] Created {merged_path.name} ({len(merged)} rows)")
    return merged_path


def read_photon_lr(csv_path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(csv_path, comment="#", header=0)
        if {"EventID", "iz", "iy", "Left", "Right"} <= set(df.columns):
            return df[["EventID", "iz", "iy", "Left", "Right"]]
    except Exception:
        pass

    try:
        df = pd.read_csv(csv_path, comment="#", header=None, sep=",")
        if df.shape[1] >= 5:
            df = df.iloc[:, :5].copy()
            df.columns = ["EventID", "iz", "iy", "Left", "Right"]
            return df
    except Exception:
        pass

    raise ValueError(f"Unable to parse CSV: {csv_path}")


def safe_norm(y: np.ndarray) -> np.ndarray:
    ymax = float(np.max(y)) if y.size > 0 else 0.0
    if ymax <= 0.0:
        return np.zeros_like(y, dtype=float)
    return y.astype(float) / ymax


def compute_doi(L: np.ndarray, R: np.ndarray, mode: str) -> np.ndarray:
    eps = 1e-6
    if mode == "r_over_sum":
        return R / (L + R + eps)
    if mode == "r_over_l":
        return R / (L + eps)
    if mode == "diff_over_sum":
        return (R - L) / (R + L + eps)
    if mode == "log_r_over_l":
        return np.log((R + eps) / (L + eps))
    raise ValueError(f"Unsupported DOI mode: {mode}")


def get_doi_label(mode: str) -> str:
    labels = {
        "r_over_sum": "R / (L + R)",
        "r_over_l": "R / L",
        "diff_over_sum": "(R - L) / (R + L)",
        "log_r_over_l": "log(R / L)",
    }
    return labels.get(mode, mode)


def get_default_doi_range(mode: str) -> tuple[float, float]:
    defaults = {
        "r_over_sum": (0.0, 1.0),
        "r_over_l": (0.0, 4.0),
        "diff_over_sum": (-1.0, 1.0),
        "log_r_over_l": (-3.0, 3.0),
    }
    return defaults.get(mode, (0.0, 1.0))


def analyze_one_rod(df: pd.DataFrame, iy: int, iz: int, args: argparse.Namespace) -> Dict:
    result = {
        "iy": iy,
        "iz": iz,
        "n_events": 0,
        "n_events_used": 0,
        "n_peaks": 0,
        "peak_positions": [],
        "avg_fwhm": np.nan,
        "min_pvr": np.nan,
        "min_pvr_db": np.nan,
        "hist_x": None,
        "hist_count_raw": None,
        "hist_count_smooth": None,
        "hist_weight_smooth": None,
        "hist_count_norm": None,
        "hist_weight_norm": None,
        "peaks": np.array([], dtype=int),
        "success": False,
        "error_msg": "",
    }

    sub = df[(df["iy"] == iy) & (df["iz"] == iz)].copy()
    if sub.empty:
        result["error_msg"] = "No data for this rod"
        return result

    L = pd.to_numeric(sub["Left"], errors="coerce").to_numpy(dtype=float)
    R = pd.to_numeric(sub["Right"], errors="coerce").to_numpy(dtype=float)
    total_photons = L + R
    doi = compute_doi(L, R, args.doi_mode)

    valid = np.isfinite(doi) & np.isfinite(total_photons) & (total_photons > 0)
    doi = doi[valid]
    total_photons = total_photons[valid]
    result["n_events"] = int(len(doi))
    if len(doi) < args.min_events:
        result["error_msg"] = "Insufficient events"
        return result

    strong = total_photons >= args.min_total_photons
    doi_u = doi[strong]
    total_u = total_photons[strong]
    result["n_events_used"] = int(len(doi_u))
    if len(doi_u) < args.min_events:
        result["error_msg"] = f"Insufficient events after total>={args.min_total_photons}"
        return result

    hist_count_raw, bins = np.histogram(doi_u, bins=args.n_bins, range=(args.doi_min, args.doi_max))
    hist_weight_raw, _ = np.histogram(doi_u, weights=total_u, bins=args.n_bins, range=(args.doi_min, args.doi_max))
    hist_x = (bins[:-1] + bins[1:]) / 2.0

    hist_count_smooth = gaussian_filter1d(hist_count_raw.astype(float), sigma=args.sigma_smooth)
    hist_weight_smooth = gaussian_filter1d(hist_weight_raw.astype(float), sigma=args.sigma_smooth)
    hist_count_norm = safe_norm(hist_count_smooth)
    hist_weight_norm = safe_norm(hist_weight_smooth)

    peak_distance = max(1, int(args.distance_ratio * len(hist_x)))
    prom_abs = max(args.min_prominence_abs, args.prominence_ratio * float(np.max(hist_count_smooth)))
    peaks, props = find_peaks(
        hist_count_smooth,
        distance=peak_distance,
        prominence=(prom_abs, None),
        width=(args.min_width_bins, None),
        wlen=int(0.15 * len(hist_x)) if len(hist_x) >= 10 else None,
    )

    if args.max_peaks > 0 and len(peaks) > args.max_peaks:
        order = np.argsort(props["prominences"])[::-1][: args.max_peaks]
        peaks = np.sort(peaks[order])

    result["hist_x"] = hist_x
    result["hist_count_raw"] = hist_count_raw
    result["hist_count_smooth"] = hist_count_smooth
    result["hist_weight_smooth"] = hist_weight_smooth
    result["hist_count_norm"] = hist_count_norm
    result["hist_weight_norm"] = hist_weight_norm
    result["peaks"] = peaks
    result["n_peaks"] = int(len(peaks))
    result["peak_positions"] = [float(hist_x[i]) for i in peaks]
    result["doi_label"] = get_doi_label(args.doi_mode)
    result["doi_min"] = args.doi_min
    result["doi_max"] = args.doi_max

    if len(peaks) > 0:
        widths_res = peak_widths(hist_count_smooth, peaks, rel_height=0.5)
        width_bins = widths_res[0]
        bin_size = 1.0 / float(args.n_bins)
        fwhm_values = width_bins * bin_size
        if len(fwhm_values) > 0:
            result["avg_fwhm"] = float(np.mean(fwhm_values))

    if len(peaks) >= 2:
        pvr_values: List[float] = []
        for i in range(len(peaks) - 1):
            p1, p2 = peaks[i], peaks[i + 1]
            valley = np.min(hist_count_smooth[p1 : p2 + 1])
            peak_min = min(hist_count_smooth[p1], hist_count_smooth[p2])
            if valley > 0:
                pvr_values.append(float(peak_min / valley))
        if pvr_values:
            result["min_pvr"] = float(np.min(pvr_values))
            result["min_pvr_db"] = float(10.0 * np.log10(result["min_pvr"]))

    result["success"] = True
    return result


def plot_single_rod(result: Dict, config: Dict, config_name: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6))
    if not result["success"] or result["hist_x"] is None:
        ax.text(
            0.5,
            0.5,
            f"Rod (iy={result['iy']}, iz={result['iz']})\n{result['error_msg']}",
            ha="center",
            va="center",
            fontsize=12,
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    else:
        x = result["hist_x"]
        cnt_norm = result["hist_count_norm"]
        w_norm = result["hist_weight_norm"]
        peaks = result["peaks"]

        ax.plot(x, cnt_norm, label="Count-based (smoothed, norm)", linewidth=2.0, color="tab:blue")
        ax.plot(x, w_norm, label="Weighted by L+R (smoothed, norm)", linewidth=1.5, alpha=0.8, color="tab:orange")
        if len(peaks) > 0:
            ax.scatter(x[peaks], cnt_norm[peaks], s=45, zorder=5, color="red", marker="x", label="Detected Peaks")
        ax.set_xlabel(f"DOI ({result['doi_label']})", fontsize=11)
        ax.set_ylabel("Normalized Intensity", fontsize=11)
        ax.set_xlim(result["doi_min"], result["doi_max"])
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=9)

    cfg = f"Nx={int(config.get('arrayNx', 0))} Ny={int(config.get('arrayNy', 0))} Nz={int(config.get('arrayNz', 0))}"
    stats = f"Peaks={result['n_peaks']} | Used={result['n_events_used']}/{result['n_events']}"
    if not np.isnan(result["avg_fwhm"]):
        stats += f" | AvgFWHM={result['avg_fwhm']:.4f}"
    if not np.isnan(result["min_pvr"]):
        stats += f" | MinPVR={result['min_pvr']:.2f} ({result['min_pvr_db']:.1f} dB)"
    ax.set_title(f"{config_name} | Rod(iy={result['iy']}, iz={result['iz']}) | {cfg}\n{stats}", fontsize=10, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_summary_grid(results: List[Dict], config: Dict, config_name: str, output_path: Path) -> None:
    ok_results = [r for r in results if r["success"] and r["hist_x"] is not None]
    if not ok_results:
        return
    fig, ax = plt.subplots(figsize=(16, 10))
    colors = plt.colormaps.get_cmap("tab20").resampled(len(ok_results))
    for idx, res in enumerate(ok_results):
        x = res["hist_x"]
        y = res["hist_count_norm"]
        iy, iz = res["iy"], res["iz"]
        ax.plot(x, y, linewidth=1.5, alpha=0.8, color=colors(idx), label=f"Rod(iy={iy},iz={iz}) P={res['n_peaks']}")
        peaks = res["peaks"]
        if len(peaks) > 0:
            ax.scatter(x[peaks], y[peaks], s=28, zorder=5, color=colors(idx), marker="x", alpha=0.9)

    cfg = f"Nx={int(config.get('arrayNx', 0))} Ny={int(config.get('arrayNy', 0))} Nz={int(config.get('arrayNz', 0))}"
    ax.set_title(f"{config_name} - All Rods (count-based normalized) | {cfg}", fontsize=14, fontweight="bold")
    doi_label = ok_results[0].get("doi_label", "DOI")
    doi_min = ok_results[0].get("doi_min", 0.0)
    doi_max = ok_results[0].get("doi_max", 1.0)
    ax.set_xlabel(f"DOI ({doi_label})", fontsize=12)
    ax.set_ylabel("Normalized Intensity", fontsize=12)
    ax.set_xlim(doi_min, doi_max)
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=9, framealpha=0.9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def process_config_folder(folder: Path, output_folder: Path, args: argparse.Namespace) -> List[Dict]:
    config_name = folder.name
    print(f"\n[Processing] {config_name}")
    mac_path = folder / "geometry.mac"
    config = parse_geometry_mac(mac_path)
    if not config:
        print("  [Warn] No geometry.mac found or failed to parse")
        return []
    print(f"  [Config] Nx={config.get('arrayNx', '?')} Ny={config.get('arrayNy', '?')} Nz={config.get('arrayNz', '?')}")

    csv_path = merge_thread_csvs(folder)
    if csv_path is None:
        print("  [Skip] No PhotonLRPerRod CSV found")
        return []

    try:
        df = read_photon_lr(csv_path)
        for c in ["EventID", "iz", "iy", "Left", "Right"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=["iz", "iy", "Left", "Right"]).copy()
        df["iz"] = df["iz"].astype(int)
        df["iy"] = df["iy"].astype(int)
    except Exception as e:
        print(f"  [Error] Failed to read CSV: {e}")
        return []

    ny = int(config.get("arrayNy", 1))
    nz = int(config.get("arrayNz", 1))
    config_output = output_folder / config_name
    config_output.mkdir(parents=True, exist_ok=True)

    all_results: List[Dict] = []
    for iy in range(ny):
        for iz in range(nz):
            print(f"  [Analyze] Rod iy={iy}, iz={iz}", end=" ... ")
            res = analyze_one_rod(df, iy, iz, args)
            all_results.append(res)
            if res["success"]:
                print(f"OK (peaks={res['n_peaks']}, used={res['n_events_used']})")
            else:
                print(f"FAILED ({res['error_msg']})")
            rod_output = config_output / f"rod_iy{iy}_iz{iz}.png"
            plot_single_rod(res, config, config_name, rod_output)

    summary_output = config_output / "summary.png"
    print("  [Summary] Generating overlay plot...")
    plot_summary_grid(all_results, config, config_name, summary_output)

    for res in all_results:
        res["config_name"] = config_name
        res["config"] = config
    return all_results


def main() -> None:
    ap = argparse.ArgumentParser(description="Histo9: 稳健 DOI 批量分析（连续晶体友好）")
    ap.add_argument("results", nargs="?", type=str, default=None, help="Results 文件夹路径")
    ap.add_argument("--results", dest="results_alt", type=str, default=None, help="Results 文件夹路径（备用）")
    ap.add_argument("--output", type=str, default=None, help="输出文件夹路径")
    ap.add_argument("--n-bins", type=int, default=300, help="直方图 bin 数")
    ap.add_argument("--sigma-smooth", type=float, default=1.2, help="高斯平滑 sigma")
    ap.add_argument("--distance-ratio", type=float, default=0.04, help="峰最小间距占比（相对 bin 数）")
    ap.add_argument("--prominence-ratio", type=float, default=0.015, help="峰 prominence 比例（相对平滑峰高）")
    ap.add_argument("--min-prominence-abs", type=float, default=5.0, help="峰 prominence 绝对下限（计数）")
    ap.add_argument("--min-width-bins", type=int, default=2, help="峰最小宽度（bin）")
    ap.add_argument("--min-total-photons", type=float, default=10.0, help="事件门限：仅保留 L+R >= 此值")
    ap.add_argument("--min-events", type=int, default=30, help="单 rod 最少事件数（过滤后）")
    ap.add_argument("--max-peaks", type=int, default=0, help="保留峰数量上限；0 表示不限")
    ap.add_argument(
        "--doi-mode",
        type=str,
        default="r_over_sum",
        choices=["r_over_sum", "r_over_l", "diff_over_sum", "log_r_over_l"],
        help="DOI 定义模式",
    )
    ap.add_argument("--doi-min", type=float, default=None, help="DOI 直方图区间下限")
    ap.add_argument("--doi-max", type=float, default=None, help="DOI 直方图区间上限")
    args = ap.parse_args()

    if args.doi_min is None or args.doi_max is None:
        dmin, dmax = get_default_doi_range(args.doi_mode)
        if args.doi_min is None:
            args.doi_min = dmin
        if args.doi_max is None:
            args.doi_max = dmax
    if args.doi_min >= args.doi_max:
        print(f"Error: invalid DOI range [{args.doi_min}, {args.doi_max}]")
        sys.exit(1)

    results_path = args.results or args.results_alt
    if results_path is None:
        print("Error: 请提供 Results 文件夹路径")
        sys.exit(1)
    results_dir = Path(results_path).resolve()
    if not results_dir.exists():
        print(f"Error: Results 文件夹不存在: {results_dir}")
        sys.exit(1)

    if args.output:
        output_dir = Path(args.output).resolve()
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path.cwd() / f"Histo9_Output_{ts}"
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")
    print(f"DOI mode: {args.doi_mode} ({get_doi_label(args.doi_mode)}), range=({args.doi_min}, {args.doi_max})")

    config_folders = sorted([f for f in results_dir.iterdir() if f.is_dir() and not f.name.startswith(".")], key=lambda x: x.name)
    if not config_folders:
        print(f"Error: 在 {results_dir} 中未找到配置子文件夹")
        sys.exit(1)
    print(f"Found {len(config_folders)} configuration folders")

    all_rows: List[Dict] = []
    for folder in config_folders:
        results = process_config_folder(folder, output_dir, args)
        for res in results:
            if not res["success"]:
                continue
            all_rows.append(
                {
                    "Config": res["config_name"],
                    "iy": res["iy"],
                    "iz": res["iz"],
                    "Nx": int(res["config"].get("arrayNx", 0)),
                    "Ny": int(res["config"].get("arrayNy", 0)),
                    "Nz": int(res["config"].get("arrayNz", 0)),
                    "Gap_mm": res["config"].get("crystalGap", np.nan),
                    "NumEvents": res["n_events"],
                    "NumEventsUsed": res["n_events_used"],
                    "NumPeaks": res["n_peaks"],
                    "DOIMode": args.doi_mode,
                    "DOIRangeMin": args.doi_min,
                    "DOIRangeMax": args.doi_max,
                    "PeakPositions": ",".join(f"{x:.4f}" for x in res["peak_positions"]),
                    "AvgFWHM": res["avg_fwhm"],
                    "MinPVR": res["min_pvr"],
                    "MinPVR_dB": res["min_pvr_db"],
                }
            )

    if all_rows:
        summary_csv = output_dir / "summary.csv"
        pd.DataFrame(all_rows).to_csv(summary_csv, index=False)
        print(f"\n[Done] Global summary saved to: {summary_csv}")
        print(f"Total rods analyzed: {len(all_rows)}")
    else:
        print("\n[Warn] No successful analyses to summarize")
    print(f"\n[Complete] All outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
