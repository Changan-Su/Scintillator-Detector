import pandas as pd
import matplotlib.pyplot as plt

filename = r'.\build\AnaEx01_nt_PhotonDepthNtuple.csv'

# 跳过前三行，手动设定列名
df = pd.read_csv(filename, skiprows=3, names=["DepthIndex"])

# 转换为整数类型
df["DepthIndex"] = pd.to_numeric(df["DepthIndex"], errors='coerce')  # 自动处理非数字
df = df.dropna().astype({"DepthIndex": int})  # 去除空值并转换为整数

# 绘图
plt.figure(figsize=(8, 5))
plt.hist(df["DepthIndex"], bins=range(df["DepthIndex"].min(), df["DepthIndex"].max() + 2), edgecolor='black', align='left')
plt.xlabel("Depth Index")
plt.ylabel("Photon Count")
plt.title("Histogram of Photon Depth Index")
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()
