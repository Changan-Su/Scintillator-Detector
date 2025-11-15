# -*- coding: utf-8 -*-
"""
按 (iy, iz) 选择单条 rod，基于 DOI = R/(L+R) 的直方图（按 L+R 加权），
进行平滑与多峰高斯拟合，输出：
  - 相邻层之间的 Crosstalk（i→j、j→i、对称平均）
  - Peak-to-Valley Ratio (PVR) 及 PVR(dB)

数据源：AnaEx01_nt_PhotonLRPerRod.csv（至少包含列：EventID, iz, iy, Left, Right）

用法：
  python Histo6.py --iy 5 --iz 7
  python Histo6.py --iy 1 --iz 1 --csv AnaEx01_nt_PhotonLRPerRod.csv
  python Histo6.py --iy 1 --iz 1 --save-csv XT_PVR_iy1_iz1.csv
"""

import argparse
import numpy as np
import pandas as pd

# 不画图：不导入 pyplot
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import norm

# ===================== 可调参数 =====================
n_bins = 300                 # DOI 直方图的 bin 数
sigma_smooth = 1.0           # 直方图一维高斯平滑的 sigma（bin 单位）
distance_ratio = 0.035       # 峰间最小距离（相对 bins 的比例）
prominence_ratio = 0.01      # 初始 prominence 比例
target_n = 11                # 期望峰数；不限定可设 None
sigma_init_guess = 0.010     # 高斯初始 sigma
mu_window = 0.030            # 拟合时 μ 的搜索窗口（左右各 ±mu_window）
sigma_bounds = (0.003, 0.050)# sigma 的上下界
# ===================================================

def read_photon_lr(csv_path: str) -> pd.DataFrame:
    """鲁棒读取器：优先尝试带表头；失败则按无表头强制列名。"""
    # ① 尝试带表头
    try:
        df0 = pd.read_csv(csv_path, comment="#", header=0)
        lower_map = {c.lower(): c for c in df0.columns}
        rename = {}
        want = ["eventid", "iz", "iy", "left", "right"]
        for w in want:
            if w in lower_map:
                src = lower_map[w]
                if w == "eventid":
                    dst = "EventID"
                elif w in ("left", "right"):
                    dst = w.capitalize()
                else:
                    dst = w
                rename[src] = dst
        if rename:
            df0 = df0.rename(columns=rename)
        if {"EventID","iz","iy","Left","Right"} <= set(df0.columns):
            print("CSV columns detected (headered):", df0.columns.tolist())
            return df0[["EventID","iz","iy","Left","Right"]]
    except Exception:
        pass

    # ② 无表头读取（首行为数据）
    df = pd.read_csv(csv_path, comment="#", header=None, sep=None, engine="python")
    if df.shape[1] < 5:
        raise ValueError(f"CSV 列数不足5列，实际 {df.shape[1]}，请检查导出顺序/分隔符。")
    df = df.iloc[:, :5].copy()
    df.columns = ["EventID","iz","iy","Left","Right"]
    print("CSV columns forced:", df.columns.tolist())
    return df

def multi_gaussian(x, *params):
    """多个高斯之和：params = [A1, mu1, s1, A2, mu2, s2, ...]"""
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        A, mu, sigma = params[i], params[i+1], params[i+2]
        y += A * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2))
    return y

def analyze_one_rod(df: pd.DataFrame, iy: int, iz: int, csv_out: str = None):
    """对一根 rod(iy, iz) 执行：构建 DOI 直方图→平滑→找峰→多峰高斯拟合→计算 Crosstalk + PVR。"""
    sub = df[(df["iy"] == iy) & (df["iz"] == iz)].copy()
    if sub.empty:
        raise ValueError(f"选中的 rod(iy={iy}, iz={iz}) 在文件里没有记录。")

    # 1) 计算 DOI 并做加权直方图（权重 = L+R）
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

    # 2) 峰检测（自适应 prominence）
    peak_distance = max(1, int(distance_ratio * len(hist_x)))
    maxY = float(np.max(hist_y_smooth)) if hist_y_smooth.size else 0.0
    prom0 = prominence_ratio * maxY if maxY > 0 else 0.0

    def detect_peaks(prom):
        return find_peaks(
            hist_y_smooth,
            distance=peak_distance,
            prominence=(prom, None),
            width=(2, None),
            wlen=int(0.15 * len(hist_x)) if len(hist_x) >= 10 else None
        )

    prom_candidates = [prom0,
                       0.008*maxY, 0.006*maxY, 0.004*maxY,
                       0.003*maxY, 0.002*maxY, 0.001*maxY]

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

    # 3) 多峰高斯拟合
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
            # 放宽 sigma 上界重试
            ub2 = ub[:]
            for i in range(2, len(ub2), 3):
                ub2[i] = max(0.08, ub2[i])
            popt, _ = curve_fit(
                multi_gaussian, hist_x, hist_y_smooth,
                p0=init_params, bounds=(lb, ub2), maxfev=40000
            )
        fitted_y = multi_gaussian(hist_x, *popt)

    if popt is None:
        print("没有成功拟合到峰，无法计算 Crosstalk / PVR。")
        return

    # 4) 提取并按 μ 排序各个高斯分量
    comps = []
    for i in range(0, len(popt), 3):
        A, mu, s = popt[i], popt[i+1], popt[i+2]
        comps.append((A, mu, s))
    comps = sorted(comps, key=lambda t: t[1])

    # 5) 相邻峰对：找 valley（在两峰中心之间 fitted_y 最小点），并计算 Crosstalk + PVR
    out_rows = []
    for k in range(len(comps) - 1):
        A1, mu1, s1 = comps[k]
        A2, mu2, s2 = comps[k+1]

        # 在 [mu1, mu2] 内搜索最小拟合值点作为边界 x*
        seg = (hist_x >= mu1) & (hist_x <= mu2)
        if not np.any(seg):
            x_star = 0.5 * (mu1 + mu2)
            idx_v = np.argmin(np.abs(hist_x - x_star))
        else:
            idx_v_local = np.argmin(fitted_y[seg])
            idx_v = np.where(seg)[0][0] + idx_v_local
            x_star = hist_x[idx_v]

        # Crosstalk（用正态分布尾积分，独立于峰高）
        z1 = (x_star - mu1) / (s1 + 1e-12)
        z2 = (x_star - mu2) / (s2 + 1e-12)
        xt_i2j = (1.0 - norm.cdf(z1)) * 100.0  # i→j：右尾
        xt_j2i = (      norm.cdf(z2)) * 100.0  # j→i：左尾
        xt_sym = 0.5 * (xt_i2j + xt_j2i)

        # PVR：保守定义 = min(两峰峰高) / valley（在拟合曲线 total fit 上）
        valley = max(fitted_y[idx_v], 1e-12)
        peak_i = A1
        peak_j = A2
        pvr = min(peak_i, peak_j) / valley
        pvr_db = 10.0 * np.log10(pvr)

        out_rows.append({
            "Rod_iy": iy, "Rod_iz": iz,
            "Layer_i": k, "Layer_j": k+1,
            "Boundary": float(x_star),
            "XT_i2j_%": float(xt_i2j),
            "XT_j2i_%": float(xt_j2i),
            "XT_sym_%": float(xt_sym),
            "Valley": float(valley),
            "Peak_i": float(peak_i),
            "Peak_j": float(peak_j),
            "PVR": float(pvr),
            "PVR_dB": float(pvr_db),
        })

    # 6) 打印结果
    print(f"\n== Crosstalk & PVR for Rod (iy={iy}, iz={iz}) ==")
    for r in out_rows:
        print(f"层 {r['Layer_i']:02d} ↔ {r['Layer_j']:02d} | "
              f"XT(i→j)={r['XT_i2j_%']:6.2f}%  XT(j→i)={r['XT_j2i_%']:6.2f}%  "
              f"Sym={r['XT_sym_%']:6.2f}%  "
              f"PVR={r['PVR']:7.2f}  ({r['PVR_dB']:6.2f} dB)  "
              f"boundary x*={r['Boundary']:.5f}")

    # 7) 可选导出 CSV
    if csv_out:
        pd.DataFrame(out_rows).to_csv(csv_out, index=False)
        print(f"[Saved] Crosstalk & PVR → {csv_out}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iy", type=int, required=True, help="rod 的 y 索引（0-based）")
    ap.add_argument("--iz", type=int, required=True, help="rod 的 z 索引（0-based）")
    ap.add_argument("--csv", type=str, default="AnaEx01_nt_PhotonLRPerRod.csv",
                    help="PhotonLRPerRod CSV 路径")
    ap.add_argument("--save-csv", type=str, default=None,
                    help="将 Crosstalk & PVR 结果保存为 CSV 的路径（可选）")
    args = ap.parse_args()

    df = read_photon_lr(args.csv)
    # 统一为数值类型
    for c in ["EventID","iz","iy","Left","Right"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["iz","iy","Left","Right"]).copy()
    df["iz"] = df["iz"].astype(int)
    df["iy"] = df["iy"].astype(int)

    analyze_one_rod(df, args.iy, args.iz, csv_out=args.save_csv)

if __name__ == "__main__":
    main()
