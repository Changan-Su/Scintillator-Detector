# flow3_NN 上手指南

> 针对 **完全没学过神经网络、Python 也不熟** 的读者写。
> 我会把每个脚本在干什么、每段代码为什么这么写、常见问题怎么查，都讲一遍。

**目录**

1. [先决条件](#1-先决条件)
2. [脚本分工](#2-脚本分工)
3. [端到端跑一遍](#3-端到端跑一遍)
4. [每个脚本内部在做什么（逐步讲解）](#4-每个脚本内部在做什么逐步讲解)
5. [关键概念词典](#5-关键概念词典)
6. [排错 checklist](#6-排错-checklist)
7. [下一步可以怎么做](#7-下一步可以怎么做)

> **本文档最近一次大改：2026-04-18**
> 增加了 `train2.py`（GPU 大 batch 版）、`predict2.py`（对称镜像覆盖全卦限）、`predict3.py`（总文件夹批量 + accuracy_summary + 散点图），以及 `build_dataset.py` 的 `--photons-per-sample` / `--normalize-counts` 选项。设计动机见 `Document/Notes/Note-2026-04-18-NN-Position-Reconstruction.md` 末尾"事后复盘"一节。

---

## 1. 先决条件

- uv 环境已经建好，torch 能看到 GPU（见 `Note-2026-04-18-uv-Python-Manager.md`）
- 已经用 `Histo10_Cubic.py` 处理过一批 Output 目录，里面有 `merged_event.csv` 和 `metadata.csv`
- 在**项目根目录** (`D:\Geant4\Projects\Scintillator-Detector-Continuous`) 运行所有命令

验证 torch + GPU：

```powershell
uv run python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

期望：`True NVIDIA GeForce RTX 5060`

---

## 2. 脚本分工

| 脚本 | 作用 | 输入 | 产物 |
|---|---|---|---|
| `build_dataset.py`    | 仿真 CSV → `(X, y)`，支持**光子聚合**和**归一化** | `Output/SiPM6_Output_*/<config>/`      | `artifacts/dataset.npz` |
| `inspect_dataset.py`  | 快速体检 `dataset.npz`（样本数、光子分布、y 范围）| `dataset.npz`                          | 终端打印 |
| `train.py`            | 常规 DataLoader 训练                              | `dataset.npz`                          | `artifacts/best.pt` 等 |
| `train2.py`           | **数据常驻 GPU + 大 batch 版**（RTX 显卡约 10× 提速）| `dataset.npz`                       | `artifacts_gpu/best.pt` 等 |
| `evaluate.py`         | 测试集 MAE/RMSE + hexbin 散点图                   | 训练产物                               | `test_metrics.csv`、`test_scatter.png` |
| `predict.py`          | 对**单个** config 做推理（只在 +++ 卦限准）       | `best.pt` + `norm.npz` + config        | `nn_predicted_events.csv`（`nn_mlp`） |
| `predict2.py`         | 对称镜像版——覆盖**全部 8 个卦限**                 | 同上                                   | `nn_predicted_events_sym.csv`（`nn_mlp_sym`） |
| `predict3.py`         | **批量推理 + 精度汇总 + 三视图散点图**            | `best.pt` + 总文件夹                   | `accuracy_summary_nn.csv` + `scatter_nn_mlp_sym.png` |
| `model.py`            | 模型定义（被其它脚本 import）                     | -                                      | - |

```
build_dataset.py  →  dataset.npz
                         ↓
          train.py / train2.py  →  best.pt + norm.npz + split_idx.npz
                         ↓
    evaluate.py   predict.py   predict2.py   predict3.py(批量+汇总)
```

**train vs train2 怎么选**
- `train.py`：稳，和常见 PyTorch 教程一致，适合首次跑通和调试
- `train2.py`：数据一次性 `.to(cuda)` + 手写 batch 循环 + 默认 batch=8192。RTX 5060 上每 epoch ~2 s，总训练 < 1 分钟。默认产物目录 `artifacts_gpu/`，和 train.py 隔离

**predict / predict2 / predict3 怎么选**
- `predict.py`：单 config。训练只覆盖了 +X/+Y/+Z 卦限时，对其它卦限无效
- `predict2.py`：单 config。利用立方体对称性先把输入**镜像到 +++ 卦限**再推理，输出反镜像回来，覆盖全 8 卦限
- `predict3.py`：**批量版**。给总文件夹，遍历所有 config，每个走 predict2 的对称推理，并在总文件夹输出和 `flow2/analyze_position_accuracy.py` 完全同款的 `accuracy_summary_nn.csv` + `scatter_nn_mlp_sym.png`

---

## 3. 端到端跑一遍

### 3.0 完整应用流程（从仿真到可落地的位置重建模型）

你的最终目标 **不是** "跑完一次训练脚本"，而是得到一个 **"输入 SiPM 读数 → 输出 (x, y, z)"** 的可反复调用的模型。
这套 flow3_NN 在整个工作管线中的位置如下：

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  Stage A · 准备训练数据（一次性，耗时最长）                                 │
│  Geant4 仿真扫描       → Results/<config>/...t*.csv                         │
│  [run_batch 扫源位置]    （每个 config = 一个真值位置 + 几千事件）          │
│       ↓                                                                     │
│  Histo10_Cubic.py      → Output/SiPM6_Output_<时间戳>/<config>/             │
│  [合并线程、聚合]         merged_event.csv + metadata.csv                   │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ↓
┌─────────────────────────────────────────────────────────────────────────────┐
│  Stage B · 训练模型（flow3_NN 的四个脚本）                                  │
│  build_dataset.py      → artifacts/dataset.npz        ← 事件级 (X, y)       │
│       ↓                                                                     │
│  train.py              → artifacts/best.pt            ← 最好权重            │
│                           artifacts/norm.npz          ← 归一化参数          │
│                           artifacts/split_idx.npz     ← 测试集保留          │
│       ↓                                                                     │
│  evaluate.py           → artifacts/test_metrics.csv   ← 量化分数            │
│                           artifacts/test_scatter.png                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                          ┌───── 分数不满意？回去调参、扩数据 ──┐
                          │                                      │
                          ↓                                      ↑
┌─────────────────────────────────────────────────────────────────────────────┐
│  Stage C · 固化 & 应用（你的最终产出）                                      │
│  把 best.pt + norm.npz 打包成 "生产模型"                                    │
│       ↓                                                                     │
│  predict.py            → nn_predicted_events.csv      ← 每事件 (x,y,z)      │
│  [输入任意 merged_event.csv]   追加到 reconstructed_position.csv            │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 3.0.1 你需要的训练数据"覆盖度"

NN 对**没见过的位置**会外推得很烂。所以训练数据的关键要求：

- **覆盖整个晶体体积**：扫描网格要铺满你实际关心的 (x, y, z) 范围。如果你只在 z=0 平面扫了一层，NN 对 z≠0 的位置就不会定位。
- **网格密度合理**：立方体 25 mm 的典型扫描，每轴 5–11 个点（步长 2–5 mm）就够训 MLP。更密当然更好。
- **每个位置要有足够事件**：每个 config 至少几千事件，让 NN 看到同一位置下的"光子统计涨落"。这是 NN 超越解析公式的关键。
- **扫描其它变量时小心"泄漏"**：如果你同时改了 `surface sigma` 等物理参数且没把它放进输入，NN 会把不同 sigma 的事件"错认成"不同位置。建议先固定其它变量，只扫源位置。

### 3.0.2 训练 → 评估 → 迭代循环

第一次跑通后，你大概率会多次回到这里，直到分数满意：

| 如果 test RMSE_3D … | 该怎么办 |
|---|---|
| > 2 mm，远差于预期 | 数据覆盖不足 or 数据有 bug（先检查 `metadata.csv` 的位置单位、y 范围是否合理） |
| 1–2 mm，但 train/val 也差 | 模型太小或学习率不对，调 `--hidden 512 --lr 5e-4` |
| 1–2 mm，val 远差于 train | 过拟合，调 `--dropout 0.2 --weight-decay 1e-4` 或加数据 |
| < 1 mm，接近物理极限 | 进入"第 7 节下一步"阶段：对称性增强、CNN、不确定性 |

**每轮实验都要做到**：
- 固定随机种子（`--seed 42`）
- 复用同一个 `split_idx.npz` 评估（保证 test 集不变）
- 把 `test_metrics.csv` 另存到 `artifacts/exp_<日期>_<描述>/` 下备查

### 3.0.3 固化一个"生产模型"

当你得到满意的训练结果后：

1. **把当前 artifacts 整体拷贝**到一个版本化目录，例如：
   ```
   workflow/flow3_NN/artifacts_v1_20260418/
     ├── best.pt
     ├── norm.npz
     ├── split_idx.npz
     ├── test_metrics.csv
     ├── test_scatter.png
     └── EXPERIMENT.md   ← 手写一份记录：用哪份数据、超参数是什么、分数是多少
   ```
2. **EXPERIMENT.md 至少记三件事**：训练数据路径、超参数、test RMSE_3D。将来回头能复现。
3. **只要 `best.pt` + `norm.npz` 两个文件就够用**，`predict.py` 只依赖这两个。

### 3.0.4 在新数据上做真正的"位置重建"

这就是你的最终用途。场景分两类：

**场景 A · 新一批 Geant4 仿真数据**（验证泛化）

```powershell
# 对扫描里的一个新 config 做重建
uv run python workflow/flow3_NN/predict.py `
    --config-dir Output/SiPM6_Output_20260420_XXX/<某 config> `
    --append-to-recon
```
- 输出：每个事件的 (x, y, z) + 平均 → 和真值对比
- 作用：验证模型在新仿真条件（不同 sigma、不同入射角等）下是否还靠谱

**场景 B · 真实实验的 SiPM 读数**（真正的重建）

这是最终目标。你需要先把实验数据整理成和 `merged_event.csv` 一样的格式：
```
EventID, CrystalID, iy, iz, Face, j, k, SiPMBlockID, PhotonCount
```
然后 `predict.py` 走一样的路径。此时：
- **没有真值**，无法直接判断对错 → 靠 `test_scatter.png` 上的 test 集分数做背书
- **输入分布要和训练数据一致**：如果实验的光子计数量级、噪声水平和仿真差很大，要么重新训（加入真实数据），要么在 build 阶段加噪声让仿真更贴近实验（domain adaptation）

### 3.0.5 整体时间预算（参考）

| 阶段 | 单次耗时 | 频率 |
|---|---|---|
| Geant4 扫描（几百 config） | 数小时~一天 | 很少重做 |
| Histo10_Cubic 批处理 | 分钟级 | 每批数据一次 |
| build_dataset.py | 1–5 分钟 | 数据变了才跑 |
| train.py（RTX 5060） | 2–10 分钟 | 调参时反复跑 |
| evaluate.py | < 1 分钟 | 每次训完跑 |
| predict.py（单 config） | 秒级 | 可高频 |

**调超参数的时候只重跑 train+evaluate**，build_dataset.py 不用重跑。

---

下面是四步具体命令。

### Step 1 · 构建数据集

```powershell
uv run python workflow/flow3_NN/build_dataset.py `
    --output-root Output/SiPM6_Output_20260404_234359 `
    --photons-per-sample 5000 `
    --normalize-counts
```

> **⚠️ 2026-04-18 重要更新：两个必看的新选项**
>
> - `--photons-per-sample N`：我们这套仿真用了 `/source/mode optical`，每个 Geant4 event 只产生 1 个光子。所以 `merged_event.csv` 里一个 EventID 行仅 1 光子，**不能直接当训练样本**。该选项会把同一 config 内 N 个光子行随机组合成一个"γ 事件样本"（典型值 5000，和你单次 γ 事件的光子数量级匹配）。
> - `--normalize-counts`：把每个样本除以总光子数，转成"光子**分数**"（sum=1）。推理时真实事件光子数从几百到几万都有可能，加了这个以后模型对总光子数免疫，只看分布形状。
>
> **这两个参数是成对使用的**：训练时加了，推理时 predict/predict2/predict3 **必须同样加**，否则输入尺度不一致会让 Tanh 直接饱和，输出永远卡在 ±y_scale。

约定：加了 `--normalize-counts` 的 artifacts 目录后缀加 `_norm`，便于一眼看出来（如 `artifacts_p5000_norm/`）。

**第一次建议先小跑**（只用前 50 个 config，几秒钟出结果，验证路径和格式都对）：

```powershell
uv run python workflow/flow3_NN/build_dataset.py `
    --output-root Output/SiPM6_Output_20260404_234359 `
    --max-configs 50 `
    --photons-per-sample 5000 `
    --normalize-counts
```

构建完顺手体检一下：

```powershell
uv run python workflow/flow3_NN/inspect_dataset.py
```

看到类似输出就对了：

```
扫描：D:\...\SiPM6_Output_20260404_234359
  [50] S0p3_X0_Y0p9_Z0_..._batch_0050  events=1820
汇总：成功 50 个 config，跳过 0 个。
总样本数 N = 98,450   输入维度 = 96   标签维度 = 3
y 范围：x∈[-0.900, 0.900]  y∈[-0.900, 0.900]  z∈[-0.900, 0.900]  (cm)
已保存：workflow/flow3_NN/artifacts/dataset.npz
```

确认无误后把 `--max-configs 50` 去掉，跑全量。

### Step 2 · 训练

两个选项：

```powershell
# A. 常规版（稳，DataLoader）
uv run python workflow/flow3_NN/train.py

# B. GPU 大 batch 版（RTX 显卡推荐，约 10× 加速，产物到 artifacts_gpu/）
uv run python workflow/flow3_NN/train2.py `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm
```

train.py 默认 150 epochs、batch_size=512；train2.py 默认 300 epochs、batch_size=8192，数据一次性 `.to(cuda)` 避免 H2D 拷贝。RTX 5060 + 百万级样本：train.py 几分钟，train2.py 通常 < 1 分钟。

输出长这样：

```
device = cuda
划分：train=800,000  val=100,000  test=100,000
可训练参数：92,419
epoch   1 | train 0.08210 | val 0.05430 | lr 1.00e-03  ✓ saved
epoch   5 | train 0.02150 | val 0.02009 | lr 1.00e-03  ✓ saved
...
early stop at epoch 87 (no improvement for 20 epochs)
用时 124.3s  最佳 val loss = 0.00314
已保存：best.pt / norm.npz / split_idx.npz 到 .../artifacts
loss 曲线：.../artifacts/loss_curve.png
```

打开 `loss_curve.png`：**train 和 val 两条线都应该下降**，如果 val 在某处开始反弹 → 过拟合了，调 `--dropout 0.2` 或 `--weight-decay 1e-4`。

### Step 3 · 评估

```powershell
uv run python workflow/flow3_NN/evaluate.py
```

输出：

```
=== Test Metrics (cm) ===
MAE  (x, y, z) = (0.0412, 0.0398, 0.0421)
RMSE (x, y, z) = (0.0621, 0.0598, 0.0635)
3D RMSE        = 0.1074
```

打开 `test_scatter.png`：3 张 hexbin 图越贴近红色对角线越好。**RMSE_3D 大约 1 mm 量级就已经远好于 linear_scaled**。

### Step 4 · 推理（三种粒度）

**注意**：下面所有命令都必须带 `--normalize-counts`（如果训练时加了）。漏掉会让 Tanh 直接饱和，输出全是 ±y_scale。

**4.1 单 config，+++ 卦限**（训练数据覆盖到的那一象限）

```powershell
uv run python workflow/flow3_NN/predict.py `
    --config-dir Output/.../S0p3_X0p6_Y0_Z0_..._batch_0123 `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm `
    --normalize-counts --append-to-recon
```

**4.2 单 config，任意卦限**（对称镜像版）

```powershell
uv run python workflow/flow3_NN/predict2.py `
    --config-dir Output/.../S0p3_Xm6p25_Ym6p25_Zm6p25_... `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm `
    --normalize-counts --append-to-recon
```

predict2.py 会先从 `N(+X face) - N(-X face)` 等计数差判定输入落在哪个卦限，把输入通道重排"翻到"+++ 卦限等价形式，模型推理，再把坐标反翻回原卦限。原理见脚本顶部的 docstring。

**4.3 整个总文件夹（最常用）**

```powershell
uv run python workflow/flow3_NN/predict3.py `
    --input-dir Output/SiPM6_Output_20260329_032534 `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm `
    --normalize-counts
```

会在 `Output/SiPM6_Output_20260329_032534/` 下生成：
- `accuracy_summary_nn.csv`：每个 config 一行 `(true_x,true_y,true_z, rec_x,rec_y,rec_z, dx,dy,dz,d)`，末尾 BIAS/STD 汇总，格式和 `flow2/analyze_position_accuracy.py` 完全一致
- `scatter_nn_mlp_sym.png`：XY / XZ / YZ 三视图 × 各深度切片

各子 config 里同样会追加 `nn_mlp_sym` 行到 `reconstructed_position.csv`，方便和其它 7 个算法并排对比。

---

## 4. 每个脚本内部在做什么（逐步讲解）

### 4.1 `build_dataset.py`

核心问题：**怎么把"事件 CSV"转成 NN 能吃的 `(N_events, 96)` 矩阵？**

原始 CSV 每行是"某事件里某个 SiPM 收到了 n 个光子"：

```
EventID, Face, j, k, PhotonCount
  1773,    5,  0, 1,     1.0
  1773,    0,  1, 2,     3.0
  1774,    4,  2, 0,     2.0
  ...
```

目标：

```
X[0] = [Face0·j0·k0, Face0·j0·k1, ..., Face5·j3·k3]   ← 96 维，Event 1773 的光子分布
X[1] = [...]                                          ← Event 1774
...
```

实现关键：
- `face * 16 + j * 4 + k` 把 `(Face, j, k)` 压平成 0..95 的通道下标
- `pd.factorize(EventID)` 把稀疏的 EventID 紧凑化成 0..N-1
- `np.add.at(X, (event_idx, flat_ch), counts)` 是"在指定下标处累加"的**向量化操作**，比 `for row in df.iterrows()` 快几十倍

**为什么事件级，不是聚合级**：
- 同一个 config 里所有事件真值位置相同，但每次事件的 96 维向量都不一样（因为光子数是随机的）
- NN 需要看到"同一个 y 对应多种 X"才能学到"这种随机波动应该映射到同一个位置"，本质就是学"统计平均"
- 事件级样本数瞬间膨胀 1000 倍 → NN 有充分素材

### 4.2 `model.py` —— MLP

PyTorch 里定义模型的套路：

```python
class PositionMLP(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(     # 按顺序串一堆层
            nn.Linear(96, 256),       # 全连接：96 维 → 256 维
            nn.ReLU(),                # 激活函数（非线性）
            nn.Dropout(0.1),          # 训练时随机丢 10%
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(256, 3),        # 最后 256 → 3 (x, y, z)
            nn.Tanh(),                # 压到 (-1, 1)
        )

    def forward(self, x):             # 前向怎么算，PyTorch 会自动推反向
        return self.net(x)
```

**每一层在干什么**：

| 层 | 作用 | 直觉 |
|---|---|---|
| `Linear(a, b)` | `y = xW + b`，把 a 维变成 b 维 | 学一套线性变换 |
| `ReLU()` | `max(0, x)` | 没有它的话，堆多少层都等价于一层 |
| `Dropout(p)` | 训练时把 p 比例的神经元置 0 | 防"记忆"而非"理解"（正则化） |
| `Tanh()` | 压到 (-1, 1) | 输出层约束，配合归一化后的 y |

**参数量**：`96×256 + 256 + 256×256 + 256 + 256×3 + 3 = 91,395`。对晶体位置重建这种 3 维任务，这个规模刚好。

### 4.3 `train.py` —— 训练循环

这是最核心的一步。读懂这段，你就懂 80% 的 NN 训练代码。

#### 3 件必做的事

**a. 归一化（标准化）**

```python
x_mean = X_train.mean(axis=0)
x_std  = X_train.std(axis=0) + 1e-6
X_train = (X_train - x_mean) / x_std     # 每个通道变成均值 0、方差 1
y_train = y_train / y_scale              # y 从 cm 变成 [-1, 1]
```

为什么：
- 不同 SiPM 的光子数量级差很大（角落少、中心多）。不归一化的话，数值大的 SiPM 会"霸占"梯度
- y 归一到 [-1, 1] 配合 `Tanh`，网络不用硬学边界

**b. 训练/验证/测试三分**

```python
80% train  → 网络用它调参数
10% val    → 每 epoch 看一下，决定"训够了没"和"要不要降学习率"
10% test   → 训练完才能看，用来报告最终分数
```

**千万不能用 test 调超参数**，否则 test 分数就不客观了。

**c. 标准循环**

```python
for epoch in range(epochs):
    for xb, yb in train_loader:                 # 一个 mini-batch
        pred = model(xb)                        # 前向
        loss = criterion(pred, yb)              # 算误差
        optimizer.zero_grad()                   # 清零上轮梯度
        loss.backward()                         # 反向：自动算梯度
        optimizer.step()                        # 用梯度更新参数

    # 验证
    model.eval()
    with torch.no_grad():
        val_loss = mean_loss over val_loader

    if val_loss < best_val:
        torch.save(model.state_dict(), "best.pt")   # 只存"验证最好那次"
        best_val = val_loss
```

**为什么 mini-batch**：
- 一次把所有数据喂进去 → 显存炸
- 一次只喂一条 → 噪声太大，收敛慢
- 折中：每次喂 batch_size=512 条

**loss 函数选 SmoothL1 不选 MSE**：
- MSE 对离群点敏感：一个极端错误 → loss 平方放大 → 网络被带偏
- SmoothL1（= Huber loss）小误差时像 MSE，大误差时像 MAE，对物理数据更稳

**Adam 优化器**：
- 自带"按梯度历史自适应调学习率"的机制
- 不用手调 lr schedule 也能跑，新手首选
- `lr=1e-3`、`weight_decay=1e-5` 是合理起点

**学习率调度** (`ReduceLROnPlateau`)：
- 当 val loss 连续 8 个 epoch 不下降 → lr 砍半
- 模型卡住时能"再挪一点"

**早停** (early stopping)：
- 连续 20 个 epoch 没进步就停
- 省时间，也防过拟合

### 4.4 `evaluate.py` —— 测试集评估

核心反归一化：

```python
y_pred_n = model(X_test_n)    # 归一化空间 (-1, 1)
y_pred = y_pred_n * y_scale   # 乘回来就是 cm
err = y_pred - y_test         # 真值也用 cm
```

指标：
- **MAE** (Mean Absolute Error): `mean(|err|)`，物理单位 cm，直观
- **RMSE**: `sqrt(mean(err²))`，对大误差更敏感
- **3D RMSE**: `sqrt(mean(||err||²))`，总体"距离误差"

画图用 `hexbin`：
- scatter 点太多会堆成一坨
- hexbin 用六角格显示"密度"，可视化更清楚
- 红线 `y = x` 是完美预测的参考线

### 4.5 `predict.py` —— 对新 config 推理

复用 `build_dataset.events_to_tensor`，加载权重，前向。
支持把 `nn_mlp` 行追加到 `reconstructed_position.csv`，和你原来的 7 个算法并列。

---

## 5. 关键概念词典

| 词 | 意思 |
|---|---|
| **张量 (tensor)** | 多维数组，和 numpy 的 ndarray 几乎一样，但能在 GPU 上跑、能自动求导 |
| **epoch** | 把全部训练集过一遍 |
| **iteration / step** | 一个 mini-batch 的前向+反向+更新 |
| **batch_size** | 每个 iteration 喂多少条样本 |
| **loss (损失)** | 预测和真值的差距，越小越好 |
| **梯度 (gradient)** | loss 对每个参数的偏导，告诉我们"往哪个方向改参数能让 loss 变小" |
| **autograd** | PyTorch 的自动求导机制，你写前向，它帮你推反向 |
| **optimizer** | 拿梯度更新参数的算法（Adam/SGD 等） |
| **lr (learning rate)** | 每步往梯度方向挪多远 |
| **overfitting (过拟合)** | 把训练集背下来了，但验证集表现烂 |
| **dropout** | 训练时随机丢神经元，防过拟合 |
| **regularization (正则化)** | 一类限制模型"乱学"的技术，dropout 和 weight_decay 都算 |
| **标准化** | 把数据变成均值 0、方差 1 |
| **归一化** | 把数据压到固定区间（如 [-1, 1]） |
| **early stopping** | val loss 不再降就提前停 |
| **checkpoint** | 保存的模型权重文件（`.pt`） |
| **state_dict** | PyTorch 模型权重字典，`torch.save(model.state_dict(), ...)` |

---

## 6. 排错 checklist

### Loss 一直不下降

- ☐ 输入标准化了吗？`X.mean()` ≈ 0, `X.std()` ≈ 1?
- ☐ y 归一化了吗？`y ∈ [-1, 1]`?
- ☐ learning rate 是不是太小 (<1e-5) 或太大 (>1e-1)?
- ☐ 模型输出最后一层是不是多了个不该有的 sigmoid/relu 把范围夹死了？

### Loss 一开始就 NaN

- ☐ 学习率太大 → 调小 10 倍
- ☐ 输入有 inf / nan（某个 SiPM 计数是负数？）→ build 阶段就该过滤
- ☐ `X_std` 有 0 → 已经 `+1e-6`，但若某通道始终为 0 仍有问题，可以 log1p 预处理

### Train loss 降、val loss 涨（过拟合）

- ☐ 加大 dropout (0.2~0.3)
- ☐ 加大 weight_decay (1e-4)
- ☐ 减小模型（`--hidden 128`）
- ☐ 扩充数据（扫描更多源位置）
- ☐ 数据增强：利用立方体对称性做镜像翻转（同时翻 X 和 y）

### 测试集表现比训练差很多

- ☐ 训练时 val 也很好吗？如果 val 也差，就是过拟合
- ☐ 是不是推理时忘了用同一套 `x_mean`/`x_std`？脚本内都是复用 `norm.npz`，但你如果手写了别的推理逻辑就要注意
- ☐ 推理时忘了 `model.eval()` → dropout 还在 → 输出不稳

### 推理结果全是 ±y_scale（坐标饱和）

**症状**：predict/predict2/predict3 对**任何** config 的预测都是 (±1.25, ±1.25, ±1.25) 这种整齐的极值。

根因 99% 是**训练/推理预处理不一致**：
- ☐ 训练时用了 `--normalize-counts`（artifacts 目录名一般含 `_norm`），但 predict 命令**没加** `--normalize-counts` → 输入尺度差 1000× → Tanh 彻底饱和
- ☐ 训练时 `--photons-per-sample` 和推理时不一致 → 样本量级不同
- ☐ 推理的 config 是在**训练范围之外**（比如训练只到 10 mm，推理点在 12.5 mm）→ Tanh 外推到边界

**对比训练产物目录名来自检**：`artifacts_p5000_norm` 意味着 `--photons-per-sample 5000 --normalize-counts`，推理命令要对齐。

### predict 在负卦限全错

- ☐ 训练数据只覆盖了 +X/+Y/+Z 卦限（y ≥ 0），`predict.py` 没有镜像逻辑，对其它卦限必然失败 → 改用 `predict2.py` 或 `predict3.py`

### GPU 用不上 / 显存不足

- ☐ `torch.cuda.is_available()` 返回 False → 见 uv 笔记里的 CUDA 检查
- ☐ batch_size 太大 → 调小
- ☐ 模型太大 → `--hidden 128`

### Windows / PowerShell 特定

- ☐ 路径分隔符：PowerShell 里 `/` 也能用，但保险起见用 `\`
- ☐ 多行命令用 **反引号 `**（不是反斜杠）
- ☐ `num_workers=0`（脚本里已设）：Windows 多进程 DataLoader 不稳

---

## 7. 下一步可以怎么做

按"收益 / 工作量"从高到低：

1. **数据增强（对称性）** —— 立方体有 48 种对称操作，对应 48× 数据增强，几乎免费。实现要点：对输入做面 / j / k 重排的同时对 y 做同样的坐标变换。
2. **log1p 预处理** —— 光子数跨度大，`log1p(X)` 再标准化通常更稳。
3. **CNN 版** —— 把输入 reshape 成 `(6, 4, 4)`，用小 CNN。4×4 太小 CNN 收益有限，但可以试。
4. **不确定性估计** —— 改成预测 `(mu, sigma)`，用负对数似然作 loss。每个事件都会吐一个"我多有把握"的分数。
5. **前向模型做 warm start** —— 先用 linear_scaled 给个粗定位，再用 NN 做残差修正。
6. **端到端（而非事件级）** —— 训一个接受整个 config 聚合数据的模型，把事件平均的先验也学进去。

每做一次改动，都应该：
- 固定 seed、固定 split（用 `split_idx.npz`）
- 对比同一测试集上的 RMSE_3D
- 保留一份 `test_metrics.csv` 和 `test_scatter.png` 到 `artifacts/<实验名>/` 下

这样才能客观说"哪个改动真正带来了提升"。

---

## 附：最常见的命令

```powershell
# 1. 建数据集（推荐设置：光子聚合 + 归一化）
uv run python workflow/flow3_NN/build_dataset.py `
    --output-root Output/SiPM6_Output_20260404_234359 `
    --photons-per-sample 5000 --normalize-counts

# 2. 体检一下
uv run python workflow/flow3_NN/inspect_dataset.py

# 3. 训练（GPU 大 batch 版）
uv run python workflow/flow3_NN/train2.py `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm

# 4. 评估
uv run python workflow/flow3_NN/evaluate.py `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm

# 5. 单 config 推理（对称版，支持全部 8 卦限）
uv run python workflow/flow3_NN/predict2.py `
    --config-dir <某 config 路径> `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm `
    --normalize-counts --append-to-recon

# 6. 批量推理 + 汇总 + 散点图
uv run python workflow/flow3_NN/predict3.py `
    --input-dir Output/SiPM6_Output_20260329_032534 `
    --artifacts workflow/flow3_NN/artifacts_p5000_norm `
    --normalize-counts

# 7. 看模型结构 / 参数量
uv run python workflow/flow3_NN/model.py
```

**记住**：训练时传了 `--normalize-counts`，所有 predict* 命令都必须同样传。

祝玩得开心，第一次看到 `cuda_available: True` 和 loss 曲线哗啦啦掉下去，会很爽的。
