import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import newton

# ==== 常数 ====
r  = 1.5
x2 = 30
eps = 1e-6
lambda_abs = 40.0  # 加入吸收长度


# ==== 读取数据 ====
file_left  = r'.\build\AnaEx01_nt_PhotonLeft.csv'
file_right = r'.\build\AnaEx01_nt_PhotonRight.csv'

L = pd.read_csv(file_left,  comment='#', header=0).iloc[:,0].astype(float).to_numpy()
R = pd.read_csv(file_right, comment='#', header=0).iloc[:,0].astype(float).to_numpy()

# ==== 先计算 total_photons，并对 L=R=0 的行做一次过滤 ====
total_photons = L + R
mask0 = ~((L == 0) & (R == 0))
L = L[mask0]
R = R[mask0]
total_photons = total_photons[mask0]

# ==== 定义 rho，并过滤掉 R=0 或 rho 太小的行 ====
rho_array = L / (R)
mask1 = (R > eps) & (rho_array > eps)
L = L[mask1]
R = R[mask1]
total_photons = total_photons[mask1]
rho_array = rho_array[mask1]

# ==== 待求根函数 F(x1; rho) ====
def F(x1, rho):
    # x1 可能被 newton 跑出界，先 clip 保证不过分发散
    x1 = np.clip(x1, 0.0, x2)
    exp_term = np.exp((x2 - 2*x1) / lambda_abs)
    u = 1 - x1/np.sqrt(r**2 + x1**2)
    v = 1 - (x2-x1)/np.sqrt(r**2 + (x2-x1)**2)
    return exp_term*(u/v) - rho

# ==== 用 Newton 逐点解 x1，失败的填 NaN ====
x1_array = np.full_like(rho_array, np.nan, dtype=float)
for i, rho in enumerate(rho_array):
    # —— 这里改成新的零级近似初值 ——  
    x0 = np.clip((x2 - np.log(rho)) / 2,
                 1e-6, x2 - 1e-6)
    try:
        sol = newton(lambda x: F(x, rho),
                     x0,
                     tol=1e-12,
                     maxiter=50)
        if 0.0 <= sol <= x2:
            x1_array[i] = sol
    except (RuntimeError, OverflowError):
        continue

# ==== 两种 DOI 计算方式 ====
asymmetry = ( R) / (L + R )
log_ratio = np.log((L ) / (R ))
doi2 = L / (L+R+eps) 

# ==== Asymmetry 直方图 ====
plt.figure(figsize=(8, 5))
plt.hist(asymmetry, weights=total_photons, bins=1000, color='blue', alpha=0.7)
plt.xlabel("DOI (Asymmetry: ( R)/(L + R))")
plt.ylabel("Total Photon Counts")
plt.title("Photon Counts vs DOI (Asymmetry)")
plt.grid(True)
plt.tight_layout()

# ==== Log-Ratio 直方图 ====
plt.figure(figsize=(8, 5))
plt.hist(log_ratio, weights=total_photons, bins=1000, color='green', alpha=0.7)
plt.xlabel("DOI (Log Ratio: log(L / R))")
plt.ylabel("Total Photon Counts")
plt.title("Photon Counts vs DOI (Log Ratio)")
plt.grid(True)
plt.tight_layout()

# #doi 2
# plt.figure(figsize=(8, 5))
# plt.hist(doi2, weights=total_photons, bins=1000, color='red', alpha=0.7)
# plt.xlabel("DOI : L/(L+R)")
# plt.ylabel("Total Photon Counts")
# plt.title("Photon Counts vs DOI (L/L+R)")
# plt.grid(True)
# plt.tight_layout()

# ==== 绘图：Photon Counts vs x1 ====
plt.figure(figsize=(8,5))
plt.hist(x1_array[~np.isnan(x1_array)],
         weights=total_photons[~np.isnan(x1_array)],
         bins=1000,
         color='purple',
         alpha=0.7)
plt.xlabel("x₁ (inverted from ρ = L/R)")
plt.ylabel("Total Photon Counts")
plt.title("Photon Counts vs x₁")
plt.grid(True)
plt.tight_layout()

from scipy.signal import find_peaks

# ==== 计算主峰位置并输出平均间隔 ====
hist_vals, bin_edges = np.histogram(
    x1_array[~np.isnan(x1_array)],
    bins=1000,
    weights=total_photons[~np.isnan(x1_array)]
)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
peaks, _ = find_peaks(hist_vals, prominence=5000)  # 可调参数
peak_positions = bin_centers[peaks]
peak_intervals = np.diff(peak_positions)

if len(peak_intervals) > 0:
    avg_interval = np.mean(peak_intervals)
    print(f"平均主峰间隔约为：{avg_interval:.3f} cm")
else:
    print("未能识别出多个峰，无法计算间隔")

plt.show()


