# -*- coding: utf-8 -*-
"""
批量处理 Results 文件夹中的所有配置，生成每个 rod 的 DOI 直方图。

功能：
1. 自动扫描 Results 文件夹下所有配置子文件夹
2. 解析 geometry.mac 获取晶体配置参数
3. 自动合并多线程 CSV 文件（如果需要）
4. 为每个 rod 生成带标注的直方图（峰数、FWHM、PVR）
5. 为每个配置生成汇总大图（所有 rod 的子图网格）
6. 生成全局汇总 CSV 表格

用法：
    python Histo8.py D:\\path\\to\\Results
    python Histo8.py --results D:\\path\\to\\Results --output ./Histo8_Output
    
输出结构：
    Histo8_Output_YYYYMMDD_HHMMSS/
    ├── 001_Ny1/
    │   ├── rod_iy0_iz0.png
    │   ├── rod_iy0_iz1.png
    │   └── summary.png
    ├── 002_Ny2/
    │   └── ...
    └── summary.csv
"""

import argparse
import re
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# ===================== 参数区（可按需微调） =====================
n_bins = 300
sigma_smooth = 1.0
distance_ratio = 0.035
prominence_ratio = 0.01
target_n = 11                   # 期望峰数；不限定可设 None
sigma_init_guess = 0.010
mu_window = 0.030
sigma_bounds = (0.003, 0.050)
# ===============================================================


def parse_geometry_mac(mac_path: Path) -> Dict[str, float]:
    """解析 geometry.mac 文件，提取晶体配置参数。"""
    config = {}
    if not mac_path.exists():
        return config
    
    with open(mac_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            
            # 匹配 /detector/arrayNx 11 格式
            match = re.match(r'/detector/(\w+)\s+([\d.]+)', line)
            if match:
                key = match.group(1)
                value = match.group(2)
                try:
                    config[key] = float(value)
                except ValueError:
                    pass
    
    return config


def merge_thread_csvs(folder: Path, pattern: str = "AnaEx01_nt_PhotonLRPerRod") -> Optional[Path]:
    """自动合并多线程 CSV 文件。"""
    merged_path = folder / f"{pattern}.csv"
    
    # 如果合并文件已存在且非空，直接返回
    if merged_path.exists():
        try:
            # 检查文件是否有数据行（不只是注释）
            with open(merged_path, 'r') as f:
                has_data = any(line.strip() and not line.startswith('#') 
                              for line in f.readlines())
            if has_data:
                return merged_path
            else:
                print(f"  [Merge] Existing merged file is empty, will regenerate")
        except Exception:
            pass
    
    # 查找所有 _t*.csv 文件
    thread_files = sorted(folder.glob(f"{pattern}_t*.csv"))
    if not thread_files:
        return None
    
    print(f"  [Merge] Found {len(thread_files)} thread files in {folder.name}")
    
    # 读取并合并（使用逗号分隔符）
    frames = []
    for tf in thread_files:
        try:
            df = pd.read_csv(tf, comment="#", header=None, sep=",")
            if df.shape[1] >= 5:
                df = df.iloc[:, :5].copy()
                df.columns = ["EventID","iz","iy","Left","Right"]
                frames.append(df)
        except Exception as e:
            print(f"  [Warn] Failed to read {tf.name}: {e}")
            continue
    
    if not frames:
        return None
    
    merged = pd.concat(frames, ignore_index=True)
    
    # 保存合并文件（不带注释头，直接保存数据）
    merged.to_csv(merged_path, index=False)
    print(f"  [Merge] Created {merged_path.name} ({len(merged)} rows)")
    
    return merged_path


def read_photon_lr(csv_path: Path) -> pd.DataFrame:
    """鲁棒读取器：优先尝试带表头；不行则按无表头强制列名。"""
    # 先尝试用逗号分隔符读取（Geant4 默认输出格式）
    try:
        df0 = pd.read_csv(csv_path, comment="#", header=None, sep=",")
        if df0.shape[1] < 5:
            raise ValueError(f"CSV 列数不足5列，实际 {df0.shape[1]}")
        df0 = df0.iloc[:, :5].copy()
        df0.columns = ["EventID","iz","iy","Left","Right"]
        # 转换数据类型
        for col in df0.columns:
            df0[col] = pd.to_numeric(df0[col], errors="coerce")
        df0 = df0.dropna()
        if len(df0) > 0:
            return df0
    except Exception as e:
        pass
    
    # 尝试带表头读取
    try:
        df0 = pd.read_csv(csv_path, comment="#", header=0)
        lower_map = {c.lower(): c for c in df0.columns}
        rename = {}
        want_list = ["eventid", "iz", "iy", "left", "right"]
        for want in want_list:
            if want in lower_map:
                src = lower_map[want]
                dst = "EventID" if want == "eventid" else want.capitalize() if want in ("left","right") else want
                rename[src] = dst
        if rename:
            df0 = df0.rename(columns=rename)
        if {"EventID","iz","iy","Left","Right"} <= set(df0.columns):
            return df0[["EventID","iz","iy","Left","Right"]]
    except Exception:
        pass

    # 回退：无表头读取，尝试自动分隔符
    try:
        df = pd.read_csv(csv_path, comment="#", header=None, sep=None, engine="python")
        if df.shape[1] < 5:
            raise ValueError(f"CSV 列数不足5列，实际 {df.shape[1]}")
        
        df = df.iloc[:, :5].copy()
        df.columns = ["EventID","iz","iy","Left","Right"]
        return df
    except Exception as e:
        raise ValueError(f"无法读取 CSV 文件 {csv_path.name}: {str(e)}")


def multi_gaussian(x, *params):
    """多峰高斯函数。"""
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        A, mu, sigma = params[i], params[i+1], params[i+2]
        y += A * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2))
    return y


def analyze_one_rod(df: pd.DataFrame, iy: int, iz: int) -> Dict:
    """
    分析单个 rod，返回统计信息字典。
    
    返回：
        {
            'iy': int,
            'iz': int,
            'n_events': int,
            'n_peaks': int,
            'fwhm_list': List[float],
            'avg_fwhm': float,
            'min_pvr': float,
            'min_pvr_db': float,
            'hist_x': np.ndarray,
            'hist_y': np.ndarray,
            'hist_y_smooth': np.ndarray,
            'fitted_y': Optional[np.ndarray],
            'peaks': np.ndarray,
            'success': bool,
            'error_msg': str
        }
    """
    result = {
        'iy': iy,
        'iz': iz,
        'n_events': 0,
        'n_peaks': 0,
        'fwhm_list': [],
        'avg_fwhm': np.nan,
        'min_pvr': np.nan,
        'min_pvr_db': np.nan,
        'hist_x': None,
        'hist_y': None,
        'hist_y_smooth': None,
        'fitted_y': None,
        'peaks': np.array([]),
        'success': False,
        'error_msg': ''
    }
    
    sub = df[(df["iy"] == iy) & (df["iz"] == iz)].copy()
    if sub.empty:
        result['error_msg'] = 'No data for this rod'
        return result
    
    L = pd.to_numeric(sub["Left"],  errors="coerce").to_numpy(dtype=float)
    R = pd.to_numeric(sub["Right"], errors="coerce").to_numpy(dtype=float)
    
    eps = 1e-6
    total = L + R
    asym  = R / (L + R + eps)
    
    mask = np.isfinite(asym) & np.isfinite(total) & (total > 0)
    asym  = asym[mask]
    total = total[mask]
    
    result['n_events'] = len(asym)
    
    if len(asym) < 10:
        result['error_msg'] = 'Insufficient data'
        return result
    
    hist_y, bins = np.histogram(asym, weights=total, bins=n_bins, range=(0.0, 1.0))
    hist_x = (bins[:-1] + bins[1:]) / 2.0
    hist_y_smooth = gaussian_filter1d(hist_y.astype(float), sigma=sigma_smooth)
    
    result['hist_x'] = hist_x
    result['hist_y'] = hist_y
    result['hist_y_smooth'] = hist_y_smooth
    
    # 峰检测
    peak_distance = max(1, int(distance_ratio * len(hist_x)))
    prom0 = prominence_ratio * float(np.max(hist_y_smooth)) if np.max(hist_y_smooth) > 0 else 0.0
    
    def detect_peaks(prom):
        return find_peaks(
            hist_y_smooth,
            distance=peak_distance,
            prominence=(prom, None),
            width=(2, None),
            wlen=int(0.15 * len(hist_x)) if len(hist_x) >= 10 else None
        )
    
    prom_candidates = [prom0,
                       0.008*np.max(hist_y_smooth),
                       0.006*np.max(hist_y_smooth),
                       0.004*np.max(hist_y_smooth),
                       0.003*np.max(hist_y_smooth),
                       0.002*np.max(hist_y_smooth),
                       0.001*np.max(hist_y_smooth)]
    
    peaks, props = np.array([], dtype=int), {}
    for prom in prom_candidates:
        peaks, props = detect_peaks(prom)
        if target_n is None:
            if len(peaks) >= 3:
                break
        else:
            if len(peaks) >= target_n:
                if len(peaks) > target_n:
                    order = np.argsort(props["prominences"])[::-1][:target_n]
                    peaks = peaks[order]
                break
    if len(peaks) == 0:
        peaks, props = detect_peaks(prom_candidates[-1])
    peaks = np.sort(peaks)
    
    result['peaks'] = peaks
    result['n_peaks'] = len(peaks)
    
    # 多峰高斯拟合
    init_params, lb, ub = [], [], []
    for pk in peaks:
        A_init = max(hist_y_smooth[pk], 1.0)
        mu_init = hist_x[pk]
        s_init  = sigma_init_guess
        init_params += [A_init, mu_init, s_init]
        lb += [0.0, mu_init - mu_window, sigma_bounds[0]]
        ub += [np.inf, mu_init + mu_window, sigma_bounds[1]]
    
    fitted_y, popt = None, None
    if len(init_params) >= 3:
        try:
            popt, _ = curve_fit(
                multi_gaussian, hist_x, hist_y_smooth,
                p0=init_params, bounds=(lb, ub), maxfev=20000
            )
            fitted_y = multi_gaussian(hist_x, *popt)
        except Exception as e1:
            # 放宽 sigma 上界重试
            try:
                ub2 = ub[:]
                for i in range(2, len(ub2), 3):
                    ub2[i] = max(0.08, ub2[i])
                popt, _ = curve_fit(
                    multi_gaussian, hist_x, hist_y_smooth,
                    p0=init_params, bounds=(lb, ub2), maxfev=40000
                )
                fitted_y = multi_gaussian(hist_x, *popt)
            except Exception as e2:
                result['error_msg'] = f'Fit failed: {str(e2)}'
                result['success'] = False
                return result
    
    if popt is not None:
        result['fitted_y'] = fitted_y
        
        # 计算 FWHM
        fwhm_list = []
        for i in range(0, len(popt), 3):
            sigma = popt[i+2]
            fwhm = 2.354820045 * sigma
            fwhm_list.append(fwhm)
        
        result['fwhm_list'] = fwhm_list
        result['avg_fwhm'] = float(np.mean(fwhm_list)) if fwhm_list else np.nan
        
        # 计算 PVR（简化版：相邻峰之间的 valley）
        if len(peaks) >= 2:
            pvr_list = []
            for k in range(len(peaks) - 1):
                pk1, pk2 = peaks[k], peaks[k+1]
                valley_region = fitted_y[pk1:pk2+1] if pk2 > pk1 else [fitted_y[pk1]]
                if len(valley_region) > 0:
                    valley = np.min(valley_region)
                    peak_min = min(fitted_y[pk1], fitted_y[pk2])
                    if valley > 0:
                        pvr = peak_min / valley
                        pvr_list.append(pvr)
            
            if pvr_list:
                result['min_pvr'] = float(np.min(pvr_list))
                result['min_pvr_db'] = 10.0 * np.log10(result['min_pvr'])
        
        result['success'] = True
    
    return result


def plot_single_rod(result: Dict, config: Dict, config_name: str, output_path: Path):
    """绘制单个 rod 的直方图。"""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    if not result['success'] or result['hist_x'] is None:
        ax.text(0.5, 0.5, f"Rod (iy={result['iy']}, iz={result['iz']})\n{result['error_msg']}", 
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    else:
        hist_x = result['hist_x']
        hist_y = result['hist_y']
        hist_y_smooth = result['hist_y_smooth']
        fitted_y = result['fitted_y']
        peaks = result['peaks']
        
        ax.plot(hist_x, hist_y, label="Raw Histogram", alpha=0.35, linewidth=1)
        ax.plot(hist_x, hist_y_smooth, label="Smoothed", linewidth=1.5)
        if fitted_y is not None:
            ax.plot(hist_x, fitted_y, label="Total Fit", linewidth=2.0, color='red')
        if len(peaks) > 0:
            ax.scatter(hist_x[peaks], hist_y_smooth[peaks], s=50, zorder=5, 
                      label="Detected Peaks", color='orange', marker='x')
        
        ax.set_xlabel("DOI (Asymmetry: R / (L + R))", fontsize=11)
        ax.set_ylabel("Weighted Counts (by L+R)", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=9)
    
    # 标题：包含配置信息和统计
    config_str = f"Nx={int(config.get('arrayNx', 0))} Ny={int(config.get('arrayNy', 0))} Nz={int(config.get('arrayNz', 0))}"
    stats_str = f"Peaks={result['n_peaks']}"
    if not np.isnan(result['avg_fwhm']):
        stats_str += f" | Avg FWHM={result['avg_fwhm']:.4f}"
    if not np.isnan(result['min_pvr']):
        stats_str += f" | Min PVR={result['min_pvr']:.2f} ({result['min_pvr_db']:.1f} dB)"
    
    title = f"{config_name} | Rod (iy={result['iy']}, iz={result['iz']}) | {config_str}\n{stats_str} | Events={result['n_events']}"
    ax.set_title(title, fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_summary_grid(results: List[Dict], config: Dict, config_name: str, output_path: Path):
    """绘制配置汇总大图（所有 rod 的直方图叠加显示）。"""
    if not results:
        return
    
    # 创建叠加图
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # 使用颜色映射为不同 rod 分配颜色
    colors = plt.colormaps.get_cmap('tab20').resampled(len(results))
    
    # 记录图例信息
    legend_handles = []
    legend_labels = []
    
    # 绘制所有 rod 的直方图
    for idx, res in enumerate(results):
        if not res['success'] or res['hist_x'] is None:
            continue
        
        iy, iz = res['iy'], res['iz']
        color = colors(idx)
        
        hist_x = res['hist_x']
        hist_y_smooth = res['hist_y_smooth']
        fitted_y = res['fitted_y']
        peaks = res['peaks']
        
        # 归一化处理（可选，使不同 rod 可比较）
        # hist_y_norm = hist_y_smooth / np.max(hist_y_smooth) if np.max(hist_y_smooth) > 0 else hist_y_smooth
        
        # 绘制平滑曲线
        label = f"Rod(iy={iy},iz={iz}) P={res['n_peaks']}"
        line, = ax.plot(hist_x, hist_y_smooth, linewidth=1.5, alpha=0.7, 
                       color=color, label=label)
        legend_handles.append(line)
        legend_labels.append(label)
        
        # 可选：绘制拟合曲线
        # if fitted_y is not None:
        #     ax.plot(hist_x, fitted_y, linewidth=1, alpha=0.5, 
        #            color=color, linestyle='--')
        
        # 可选：标记峰位置
        if len(peaks) > 0:
            ax.scatter(hist_x[peaks], hist_y_smooth[peaks], s=30, zorder=5, 
                      color=color, marker='x', alpha=0.8)
    
    # 设置坐标轴和标题
    ax.set_xlabel("DOI (Asymmetry: R / (L + R))", fontsize=12)
    ax.set_ylabel("Weighted Counts (by L+R)", fontsize=12)
    ax.grid(True, alpha=0.3)
    
    config_str = f"Nx={int(config.get('arrayNx', 0))} Ny={int(config.get('arrayNy', 0))} Nz={int(config.get('arrayNz', 0))}"
    title = f"{config_name} - All Rods Overlaid | {config_str}"
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # 添加图例（放在右侧）
    if legend_handles:
        ax.legend(legend_handles, legend_labels, 
                 loc='center left', bbox_to_anchor=(1, 0.5),
                 fontsize=9, framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def process_config_folder(folder: Path, output_folder: Path) -> List[Dict]:
    """处理单个配置文件夹。"""
    config_name = folder.name
    print(f"\n[Processing] {config_name}")
    
    # 解析 geometry.mac
    mac_path = folder / "geometry.mac"
    config = parse_geometry_mac(mac_path)
    if not config:
        print(f"  [Warn] No geometry.mac found or failed to parse")
        return []
    
    print(f"  [Config] Nx={config.get('arrayNx', '?')} Ny={config.get('arrayNy', '?')} Nz={config.get('arrayNz', '?')}")
    
    # 合并或查找 CSV
    csv_path = merge_thread_csvs(folder)
    if csv_path is None:
        print(f"  [Skip] No PhotonLRPerRod CSV found")
        return []
    
    # 读取数据
    try:
        df = read_photon_lr(csv_path)
        for c in ["EventID","iz","iy","Left","Right"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna(subset=["iz","iy","Left","Right"]).copy()
        df["iz"] = df["iz"].astype(int)
        df["iy"] = df["iy"].astype(int)
    except Exception as e:
        print(f"  [Error] Failed to read CSV: {e}")
        return []
    
    # 确定所有 rod
    ny = int(config.get('arrayNy', 1))
    nz = int(config.get('arrayNz', 1))
    
    # 创建输出子文件夹
    config_output = output_folder / config_name
    config_output.mkdir(parents=True, exist_ok=True)
    
    # 分析所有 rod
    all_results = []
    for iy in range(ny):
        for iz in range(nz):
            print(f"  [Analyze] Rod iy={iy}, iz={iz}", end=" ... ")
            result = analyze_one_rod(df, iy, iz)
            all_results.append(result)
            
            if result['success']:
                print(f"OK (peaks={result['n_peaks']})")
            else:
                print(f"FAILED ({result['error_msg']})")
            
            # 绘制单个 rod 图
            rod_output = config_output / f"rod_iy{iy}_iz{iz}.png"
            plot_single_rod(result, config, config_name, rod_output)
    
    # 绘制汇总大图
    summary_output = config_output / "summary.png"
    print(f"  [Summary] Generating grid plot...")
    plot_summary_grid(all_results, config, config_name, summary_output)
    
    # 添加配置信息到结果
    for res in all_results:
        res['config_name'] = config_name
        res['config'] = config
    
    return all_results


def main():
    ap = argparse.ArgumentParser(
        description="批量处理 Results 文件夹，生成所有 rod 的 DOI 直方图"
    )
    ap.add_argument("results", nargs='?', type=str, default=None,
                    help="Results 文件夹路径（可拖入）")
    ap.add_argument("--results", dest="results_alt", type=str, default=None,
                    help="Results 文件夹路径（备用参数）")
    ap.add_argument("--output", type=str, default=None,
                    help="输出文件夹路径（默认：Histo8_Output_时间戳）")
    args = ap.parse_args()
    
    # 确定 Results 路径
    results_path = args.results or args.results_alt
    if results_path is None:
        print("Error: 请提供 Results 文件夹路径")
        print("\n用法:")
        print("  python Histo8.py D:\\path\\to\\Results")
        print("  python Histo8.py --results D:\\path\\to\\Results --output ./output")
        sys.exit(1)
    
    results_dir = Path(results_path).resolve()
    if not results_dir.exists():
        print(f"Error: Results 文件夹不存在: {results_dir}")
        sys.exit(1)
    
    # 确定输出路径
    if args.output:
        output_dir = Path(args.output).resolve()
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path.cwd() / f"Histo8_Output_{timestamp}"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    # 查找所有配置子文件夹
    config_folders = [f for f in results_dir.iterdir() 
                     if f.is_dir() and not f.name.startswith('.')]
    config_folders = sorted(config_folders, key=lambda x: x.name)
    
    if not config_folders:
        print(f"Error: 在 {results_dir} 中未找到配置子文件夹")
        sys.exit(1)
    
    print(f"Found {len(config_folders)} configuration folders")
    
    # 处理每个配置
    all_summary_data = []
    for folder in config_folders:
        results = process_config_folder(folder, output_dir)
        
        # 收集汇总数据
        for res in results:
            if res['success']:
                row = {
                    'Config': res['config_name'],
                    'iy': res['iy'],
                    'iz': res['iz'],
                    'Nx': int(res['config'].get('arrayNx', 0)),
                    'Ny': int(res['config'].get('arrayNy', 0)),
                    'Nz': int(res['config'].get('arrayNz', 0)),
                    'Gap_mm': res['config'].get('crystalGap', np.nan),
                    'Size_mm': res['config'].get('crystalSize', np.nan),
                    'SizeY_mm': res['config'].get('crystalSizeY', np.nan),
                    'NumEvents': res['n_events'],
                    'NumPeaks': res['n_peaks'],
                    'AvgFWHM': res['avg_fwhm'],
                    'MinPVR': res['min_pvr'],
                    'MinPVR_dB': res['min_pvr_db']
                }
                all_summary_data.append(row)
    
    # 生成全局汇总 CSV
    if all_summary_data:
        summary_csv = output_dir / "summary.csv"
        summary_df = pd.DataFrame(all_summary_data)
        summary_df.to_csv(summary_csv, index=False)
        print(f"\n[Done] Global summary saved to: {summary_csv}")
        print(f"Total rods analyzed: {len(all_summary_data)}")
    else:
        print("\n[Warn] No successful analyses to summarize")
    
    print(f"\n[Complete] All outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
