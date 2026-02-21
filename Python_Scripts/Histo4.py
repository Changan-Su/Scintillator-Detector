import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# ==== 数据读取 ====
file_left  = './build/AnaEx01_nt_PhotonLeft.csv'
file_right = './build/AnaEx01_nt_PhotonRight.csv'
L = pd.read_csv(file_left, comment='#', header=0).iloc[:, 0].astype(float).to_numpy()
R = pd.read_csv(file_right, comment='#', header=0).iloc[:, 0].astype(float).to_numpy()

# ==== DOI 计算 ====
eps = 1e-6
total_photons = L + R
asymmetry = R / (L + R + eps)

# ==== 构建直方图 ====
hist_y, bins = np.histogram(asymmetry, weights=total_photons, bins=1000)
hist_x = (bins[:-1] + bins[1:]) / 2
hist_y_smooth = gaussian_filter1d(hist_y, sigma=2.0)

# ==== 自动寻找峰 ====
peaks, _ = find_peaks(hist_y_smooth, distance=30, height=np.max(hist_y_smooth) * 0.2)
n_peaks = len(peaks)

# ==== 多峰高斯模型 ====
def multi_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        A = params[i]
        mu = params[i + 1]
        sigma = params[i + 2]
        y += A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
    return y

# ==== 初始猜测参数 ====
init_params = []
for peak in peaks:
    A_init = hist_y_smooth[peak]
    mu_init = hist_x[peak]
    sigma_init = 0.01
    init_params += [A_init, mu_init, sigma_init]

# ==== 拟合 ====
popt, _ = curve_fit(multi_gaussian, hist_x, hist_y_smooth, p0=init_params)
fitted_y = multi_gaussian(hist_x, *popt)

# ==== 提取参数并计算 FWHM ====
fwhms = []
peak_positions = []
for i in range(n_peaks):
    A = popt[i * 3]
    mu = popt[i * 3 + 1]
    sigma = popt[i * 3 + 2]
    fwhm = 2.355 * sigma
    fwhms.append(fwhm)
    peak_positions.append(mu)

# ==== 峰谷比计算 ====
valleys, _ = find_peaks(-hist_y_smooth, distance=20)
pv_ratios = []
for peak_idx in peaks:
    left_valley = max([v for v in valleys if v < peak_idx], default=None)
    right_valley = min([v for v in valleys if v > peak_idx], default=None)
    if left_valley is not None and right_valley is not None:
        valley_min = min(hist_y_smooth[left_valley], hist_y_smooth[right_valley])
        if valley_min > 0:
            pv = hist_y_smooth[peak_idx] / valley_min
            pv_ratios.append((hist_x[peak_idx], pv))

# ==== 绘图 ====
plt.figure(figsize=(12, 6))
plt.plot(hist_x, hist_y, label="Raw Histogram", color='blue', alpha=0.4)
plt.plot(hist_x, hist_y_smooth, label="Smoothed", color='black')
plt.plot(hist_x, fitted_y, label="Total Fit", color='red', lw=2)
plt.xlabel("DOI (Asymmetry: R / (L + R))")
plt.ylabel("Total Photon Counts")
plt.title("Multi-Gaussian Fit with FWHM and P/V Ratios")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# ==== 打印 FWHM 和 P/V ====
print("\n== Gaussian Peak Fit Results ==")
for i, mu in enumerate(peak_positions):
    fwhm = fwhms[i]
    print(f"Peak {i+1:2d}: μ = {mu:.4f}, FWHM = {fwhm:.4f}")

print("\n== Peak-to-Valley Ratios ==")
for x, pv in pv_ratios:
    print(f"Peak at x = {x:.4f} → P/V = {pv:.2f}")
