# -*- coding: utf-8 -*-
"""
按 (iy, iz) 选择单条 rod，画 DOI (R/(L+R)) 直方图 + 平滑 + 多峰高斯拟合。
数据源：AnaEx01_nt_PhotonLRPerRod.csv（应包含列：EventID, iz, iy, Left, Right）
用法：
    python Histo6.py --iy 5 --iz 7
    python Histo6.py --iy 1 --iz 1 --csv AnaEx01_nt_PhotonLRPerRod.csv
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # ✅ 禁止弹窗（非交互环境）
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# ===================== 参数区（可按需微调） =====================
n_bins = 300
sigma_smooth = 1.0
distance_ratio = 0.035
prominence_ratio = 0.01
target_n = None                   # 期望峰数；不限定可设 None
sigma_init_guess = 0.010
mu_window = 0.030
sigma_bounds = (0.003, 0.050)
# ===============================================================

def read_photon_lr(csv_path: str) -> pd.DataFrame:
    """鲁棒读取器：优先尝试带表头；不行则按无表头强制列名。"""
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
            print("CSV columns detected (headered):", df0.columns.tolist())
            return df0[["EventID","iz","iy","Left","Right"]]
    except Exception:
        pass

    df = pd.read_csv(csv_path, comment="#", header=None, sep=None, engine="python")
    if df.shape[1] < 5:
        raise ValueError(f"CSV 列数不足5列，实际 {df.shape[1]}，请检查导出顺序/分隔符。")

    df = df.iloc[:, :5].copy()
    df.columns = ["EventID","iz","iy","Left","Right"]
    print("CSV columns forced:", df.columns.tolist())
    return df

def multi_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        A, mu, sigma = params[i], params[i+1], params[i+2]
        y += A * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2))
    return y

def analyze_one_rod(df: pd.DataFrame, iy: int, iz: int):
    sub = df[(df["iy"] == iy) & (df["iz"] == iz)].copy()
    if sub.empty:
        raise ValueError(f"选中的 rod(iy={iy}, iz={iz}) 在文件里没有记录。")

    L = pd.to_numeric(sub["Left"],  errors="coerce").to_numpy(dtype=float)
    R = pd.to_numeric(sub["Right"], errors="coerce").to_numpy(dtype=float)

    eps = 1e-6
    total = L + R
    asym  = R / (L + R + eps)

    mask = np.isfinite(asym) & np.isfinite(total) & (total > 0)
    asym  = asym[mask]
    total = total[mask]

    hist_y, bins = np.histogram(asym, weights=total, bins=n_bins, range=(0.0, 1.0))
    hist_x = (bins[:-1] + bins[1:]) / 2.0
    hist_y_smooth = gaussian_filter1d(hist_y.astype(float), sigma=sigma_smooth)

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
        except Exception:
            ub2 = ub[:]
            for i in range(2, len(ub2), 3):
                ub2[i] = max(0.08, ub2[i])
            popt, _ = curve_fit(
                multi_gaussian, hist_x, hist_y_smooth,
                p0=init_params, bounds=(lb, ub2), maxfev=40000
            )
        fitted_y = multi_gaussian(hist_x, *popt)

    # === 绘图部分 ===
    title = f"Rod (iy={iy}, iz={iz}) • Events={len(asym)}"
    plt.figure(figsize=(14, 5))
    plt.plot(hist_x, hist_y,        label="Raw Histogram", alpha=0.35)
    plt.plot(hist_x, hist_y_smooth, label="Smoothed",      linewidth=1.2)
    if fitted_y is not None:
        plt.plot(hist_x, fitted_y,  label="Total Fit",     linewidth=2.0)
    if len(peaks) > 0:
        plt.scatter(hist_x[peaks], hist_y_smooth[peaks], s=25, zorder=5, label="Detected Peaks")
    plt.xlabel("DOI (Asymmetry: R / (L + R))")
    plt.ylabel("Weighted Counts (by L+R)")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    # ✅ 自动保存图像文件，不弹窗
    output_file = f"DOI_iy{iy}_iz{iz}.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"[Saved] {output_file}")

    if popt is not None:
        print("\n== Gaussian Peak Fit Results ==")
        for i in range(0, len(popt), 3):
            A = popt[i]
            mu = popt[i+1]
            sigma = popt[i+2]
            fwhm = 2.354820045 * sigma
            print(f"Peak {i//3 + 1:2d}: μ = {mu:.6f}, A = {A:.1f}, σ = {sigma:.6f}, FWHM = {fwhm:.6f}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iy", type=int, required=True, help="rod 的 y 索引（0-based）")
    ap.add_argument("--iz", type=int, required=True, help="rod 的 z 索引（0-based）")
    ap.add_argument("--csv", type=str, default="AnaEx01_nt_PhotonLRPerRod_merged.csv",
                    help="PhotonLRPerRod CSV 路径")
    args = ap.parse_args()

    df = read_photon_lr(args.csv)
    for c in ["EventID","iz","iy","Left","Right"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["iz","iy","Left","Right"]).copy()
    df["iz"] = df["iz"].astype(int)
    df["iy"] = df["iy"].astype(int)
    analyze_one_rod(df, args.iy, args.iz)

if __name__ == "__main__":
    main()
