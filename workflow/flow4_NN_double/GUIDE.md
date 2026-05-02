# flow4_NN_double 上手指南

> 一句话定位：**双点光源**版本的 NN 位置重建管线。给定 96 维 SiPM 光子分布，
> 同时输出两个源 (A, B) 在晶体内的 3D 坐标 —— 共 6 维。
>
> 写给：**已经看过 `flow3_NN/GUIDE.md`、能跑通单点 NN，想理解"双点为什么不同"** 的读者。
> 如果你是第一次接触神经网络，请先读 `workflow/flow3_NN/GUIDE.md` —— 本文不重复"什么是 epoch / batch / dropout"这一类基础内容。

> **本文档最近一次大改：2026-04-25**
> 新增内容：5 个脚本的 build → train → eval → predict 端到端管线、**置换不变损失**详解、
> **多 sweep 文件夹合并训练**用法、与 flow3_NN 的逐项差异对照表。

---

## 目录

1. [先决条件](#1-先决条件)
2. [脚本分工](#2-脚本分工)
3. [端到端跑一遍](#3-端到端跑一遍)（含 [Stage 0：上游仿真](#30-stage-0--上游仿真)）
4. [每个脚本内部在做什么](#4-每个脚本内部在做什么逐步讲解)
5. [与 flow3_NN 的关键差异](#5-与-flow3_nn-的关键差异)
6. [关键概念词典](#6-关键概念词典)
7. [排错 checklist](#7-排错-checklist)
8. [下一步可以怎么做](#8-下一步可以怎么做)
9. [附：最常见的命令](#附最常见的命令)

---

## 1. 先决条件

- uv 环境 + torch + CUDA 可用（同 flow3，详见 `workflow/flow3_NN/GUIDE.md §1`）
- 已经至少跑完一次 flow4 上游仿真，根目录下有一个或多个 `Results_flow4_trainingdata*/`
  文件夹（每个里面是若干 `DP_*` config 子目录）

当前仓库里的 6 档分离距离训练数据：

```
Results_flow4_trainingdata           # 默认随机 (无 separation 限制)
Results_flow4_trainingdata_sep0p1mm  # ~0.1 mm 几乎重合
Results_flow4_trainingdata_sep1mm    # ~1 mm
Results_flow4_trainingdata_sep3mm    # ~3 mm
Results_flow4_trainingdata_sep5mm    # ~5 mm
```

验证 torch + GPU：

```powershell
uv run python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

期望：`True NVIDIA GeForce RTX 5060`

---

## 2. 脚本分工

| 脚本 | 作用 | 输入 | 产物 |
|---|---|---|---|
| `run_batch_double_point.py` | **Stage 0** · 批量双点仿真，单进程一次 beamOn 同时覆盖 (A, B) | exe + geometry.mac | `Results_flow4_trainingdata*/DP_*/...t*.csv` + `metadata.csv` + `double_point_ground_truth.csv` |
| `build_dataset_double.py` | **多 root 扫描** + 光子聚合 + 6 维标签 | 多个 `Results_flow4_trainingdata*/` | `artifacts/dataset.npz` |
| `model_double.py` | `DoublePointMLP` 96 → 256 → 256 → **6** + Tanh | - | - |
| `train_double.py` | GPU 大 batch + AdamW + **置换不变 SmoothL1 loss** | `dataset.npz` | `best.pt` / `norm.npz` / `split_idx.npz` / `loss_curve.png` |
| `evaluate_double.py` | 测试集 best-match 对齐后算 MAE/RMSE/separation 误差 | 训练产物 | `test_metrics.csv` / `test_scatter.png` / `test_error_hist.png` |
| `predict_double.py` | 单 config 推理 + 写 `nn_double_predicted.csv` | 训练产物 + 单 DP_* 目录 | `nn_double_predicted.csv` |

数据流：

```
Results_flow4_trainingdata/         ┐
Results_flow4_trainingdata_sep5mm/  ├──→  build_dataset_double.py  ──→  artifacts/dataset.npz
Results_flow4_trainingdata_sep3mm/  ┘                                                │
                                                                                     ↓
                                                            train_double.py  ──→  best.pt
                                                                                  norm.npz
                                                                                  split_idx.npz
                                                                                  loss_curve.png
                                                                                     │
                                                ┌────────────────────────────────────┴──────┐
                                                ↓                                           ↓
                                       evaluate_double.py                      predict_double.py
                                       test_metrics.csv                        nn_double_predicted.csv
                                       test_scatter.png                        （写到该 config 目录下）
                                       test_error_hist.png
```

> 与 flow3_NN 不同：flow4 **不需要** `split_events.py` 和 `Histo10_Cubic.py` —— `build_dataset_double.py`
> 直接读多线程 `_t*.csv` 原始光子命中行，跳过了"先 split → 再 merged_event.csv"两步。

---

## 3. 端到端跑一遍

### 3.0 Stage 0 · 上游仿真

完整双点仿真用法见 [`workflow/flow4_NN_double/README.md`](README.md)。这里只贴最常用的一行：

```powershell
# 100 个随机双点 config，A↔B 距离限制在 [2, 8] mm
uv run python workflow/flow4_NN_double/run_batch_double_point.py `
    --random-count 100 `
    --min-separation 2 --max-separation 8 `
    --random-seed 1234 `
    --results-dir Results_flow4_trainingdata_sep5mm
```

跑完后 `Results_flow4_trainingdata_sep5mm/` 下会有若干 `DP_*` 子目录，每个含
`metadata.csv` + 多线程 `AnaEx01_nt_PhotonFaceBlockEvent_t*.csv`。**这就是 NN 训练的原始输入**。

> NN 训练对**双点几何覆盖度**敏感。建议至少跑 3 档不同 separation：近重合（< 1 mm）、
> 中等（3–5 mm）、远分离（> 8 mm），这样模型才能学到"分得开 vs 分不开"的全谱。

### Step 1 · 构建数据集

**先小跑**（单 root + 前 20 个 config，几秒钟出结果，验证路径和格式都对）：

```powershell
uv run python workflow/flow4_NN_double/build_dataset_double.py `
    --input-roots Results_flow4_trainingdata `
    --max-configs 20 `
    --photons-per-sample 2000 --normalize-counts
```

确认无误后**多 root 全量**：

```powershell
uv run python workflow/flow4_NN_double/build_dataset_double.py `
    --input-roots Results_flow4_trainingdata `
                  Results_flow4_trainingdata_sep5mm `
                  Results_flow4_trainingdata_sep3mm `
                  Results_flow4_trainingdata_sep1mm `
                  Results_flow4_trainingdata_sep0p1mm `
    --photons-per-sample 2000 --normalize-counts
```

> ⚠️ **`--photons-per-sample` 和 `--normalize-counts` 必须和后续推理对齐**。
> 训练加了什么参数，predict_double.py 就要传什么参数。漏掉会让输入尺度差几个数量级，
> Tanh 直接饱和到 ±y_scale。和 flow3 的坑一模一样。

期望终端输出：

```
扫描根目录：
  - D:\...\Results_flow4_trainingdata
  - D:\...\Results_flow4_trainingdata_sep5mm
  - ...
  [25] DP_0025_..._batch  samples=17  sep=4.21mm  fA=0.512
  [50] DP_0050_..._batch  samples=17  sep=2.88mm  fA=0.486
  ...
汇总：成功 612 个 config，跳过 0 个。
总样本 N = 10,404   X.shape = (10404, 96)   y.shape = (10404, 6)
separation 分布: min=0.08  median=4.13  max=18.92  (mm)
fraction_a 分布: min=0.337  median=0.502  max=0.667
已保存：D:\...\workflow\flow4_NN_double\artifacts\dataset.npz
```

样本量级估算：每 config 35000 events × 概率合并后 ≈ 17 sample（35000 / 2000）。
**6 个 sweep × 100 config 大约 10000 样本**，比 flow3 单点（百万级）小 100×，所以
**调参时格外注意过拟合**（dropout / weight_decay 见 §7）。

### Step 2 · 训练

```powershell
uv run python workflow/flow4_NN_double/train_double.py --epochs 300
```

默认参数：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--y-scale` | `12.5` | 标签归一化常数（mm），= 晶体半边长 |
| `--hidden` | `256` | 隐层宽度 |
| `--dropout` | `0.1` | dropout 概率 |
| `--batch-size` | `8192` | GPU 大 batch（数据常驻 GPU，没有 H2D 拷贝开销） |
| `--epochs` | `300` | 最大轮数（一般早停） |
| `--lr` | `2e-3` | AdamW 起始学习率 |
| `--weight-decay` | `1e-5` | AdamW L2 正则 |
| `--early-stop-patience` | `30` | 连续 N 轮 val 不降就停 |
| `--seed` | `42` | 随机种子 |

期望输出（开头有一段 sanity check 验证 loss 排列不变性）：

```
device = cuda
GPU = NVIDIA GeForce RTX 5060
[sanity] perm-invariance: loss=0.831254  swapped-target loss=0.831254  diff=0.00e+00  (应≈0)
数据搬到 cuda（约 4 MB）...
划分：train=8,324  val=1,040  test=1,040
可训练参数：92,422
epoch   1 | train 0.45213 | val 0.34122 | lr 2.00e-03 | 0.03s  [saved]
epoch   5 | train 0.18044 | val 0.16782 | lr 2.00e-03 | 0.03s  [saved]
...
early stop at epoch 184
总用时 5.3s   最佳 val loss = 0.04211
```

打开 `loss_curve.png`：train 和 val 两条都应该平滑下降。如果 val 早早反弹 → 过拟合
（数据量小很容易碰到），调 `--dropout 0.2 --weight-decay 1e-4`。

### Step 3 · 评估

```powershell
uv run python workflow/flow4_NN_double/evaluate_double.py
```

输出：

```
=== Test Metrics (mm, A/B aligned by best match) ===
Point A  MAE  (x,y,z) = (0.512, 0.498, 0.521)  3D RMSE = 1.221
Point B  MAE  (x,y,z) = (0.523, 0.485, 0.510)  3D RMSE = 1.198
Separation: MAE = 0.612  RMSE = 0.873
指标已写入 .../artifacts/test_metrics.csv
图已保存：.../artifacts/test_scatter.png
误差直方图：.../artifacts/test_error_hist.png
```

可接受量级（参考，不是硬指标）：

| 3D-RMSE (A 或 B) | 评价 |
|---|---|
| > 5 mm | 模型几乎学不到东西，检查数据归一化、loss 是否在动 |
| 2–5 mm | 学到一些但训练数据不够，多铺几个 sweep |
| 1–2 mm | 不错，达到可用水平 |
| < 1 mm | 优秀 —— 但要警惕是否所有 separation regime 都准（小 sep 通常是难点） |

打开 `test_scatter.png`：6 张 hexbin 图（ax/ay/az/bx/by/bz）越贴近红色对角线越好。
**注意是 best-match 对齐后画的**，否则 A↔B 顺序歧义会让点散得一塌糊涂。

`test_error_hist.png`：`|err_A|`、`|err_B|` 距离直方图 + separation 误差直方图。
分布右尾（很大的误差）通常来自 fraction_a 极端（接近 0 或 1）的样本。

### Step 4 · 推理

单 config：

```powershell
uv run python workflow/flow4_NN_double/predict_double.py `
    --config-dir Results_flow4_trainingdata_sep5mm/DP_0001_S0p3_..._<ts> `
    --photons-per-sample 2000 --normalize-counts
```

终端输出：

```
device = cuda  artifacts = D:\...\artifacts
聚合得到 17 个样本（每样本 2000 个光子）
已写：D:\...\DP_0001_..._<ts>\nn_double_predicted.csv

--- mean prediction (mm) ---
A_pred = (+2.812, -3.901, +5.142)
B_pred = (-4.501, +5.821, -1.221)
A_true = (+3.000, -4.000, +5.000)
B_true = (-4.500, +6.000, -1.000)
fraction_a = 0.412
|mean_A - A_true| = 0.221 mm   |mean_B - B_true| = 0.247 mm
```

逐样本预测保存在该 config 目录下的 `nn_double_predicted.csv`：

```
sample_id, ax_pred, ay_pred, az_pred, bx_pred, by_pred, bz_pred
0,         2.811,   -3.902,  5.141,   -4.498,  5.823,   -1.222
1,         ...
```

`mean prediction` 是把 17 个样本预测取平均得到的（抑制 statistical fluctuation）。

### Step 4b · 批量推理

跑完单 config 验证流程对了之后，一般需要把整批测试数据一次性铺出来对比。
`predict_double_batch.py` 就是 flow3 `predict3.py` 的双点版本：扫一个总目录下所有 `DP_*/`，
每个 config 推理后把"mean pred vs truth"汇成一张表 + 散点图。

```powershell
uv run python workflow/flow4_NN_double/predict_double_batch.py `
    --input-dir Results_flow4_testdata `
    --photons-per-sample 2000 --normalize-counts
```

输出（都写在 `--input-dir` 里）：

- `accuracy_summary_double_nn.csv` —— 每行一个 config（true/pred A/B 各 3 维 + 误差 + separation），
  末尾两行是 `[SUMMARY BIAS]` 和 `[SUMMARY STD]`，方便快速看整体偏置和散布。表格视角看精度首选这里。
- `scatter_double_nn.png` —— 6 子图（A/B × x/y/z）的 true vs pred + 红色 y=x 对角线。
  视觉判断"哪些维度系统性偏移"最快。
- `separation_double_nn.png` —— 单图 sep_true vs sep_pred。重要诊断：
  模型有没有把两点压成一坨（点群偏在 y < x 一侧）？

加 `--per-config-csv` 可以同时在每个 DP_* 子目录写 `nn_double_predicted.csv`（默认关闭，
600 个文件太散）。其它 CLI 参数（`--photons-per-sample` / `--normalize-counts` / `--hidden`
/ `--seed`）必须和训练时完全一致，否则 Tanh 饱和。

> 对 600 套训练数据自检 (sanity check) 一下分布也很有用 —— 训练 RMSE 应该比测试 RMSE 略好；
> 如果反过来或者两边都很差，回头查 build/train 的 `--photons-per-sample` 和 `--normalize-counts`
> 是不是和这里对齐了。

---

## 4. 每个脚本内部在做什么（逐步讲解）

### 4.1 `build_dataset_double.py`

核心问题：**怎么把多个 sweep 文件夹下的双点仿真数据，合成一个 (X, y6) 训练集？**

**步骤 1：扫多个 root**

```python
def iter_dp_dirs(roots):
    for root in roots:
        for child in sorted(root.iterdir(), key=lambda p: p.name):
            if (child / "metadata.csv").exists():
                yield child
```

只接收 `--input-roots` 列表里**直接含 metadata.csv** 的子目录（即 DP_*）。
跨多 root 累计计数，方便用 `--max-configs` 调试时不会被某一档 sweep 的数量截断。

**步骤 2：从 metadata.csv 拿 6 维真值**

```python
POS_KEYS_A = ("source.fp_source_x_mm",   "source.fp_source_y_mm",   "source.fp_source_z_mm")
POS_KEYS_B = ("source.fp_source_b_x_mm", "source.fp_source_b_y_mm", "source.fp_source_b_z_mm")

a = [float(kv[k]) for k in POS_KEYS_A]
b = [float(kv[k]) for k in POS_KEYS_B]
y6 = np.array(a + b, dtype=np.float32)  # mm
```

注意字段后缀 `_mm` —— flow4 全程 mm，**不要**× 0.1 转 cm（那是 flow3 的事）。

**步骤 3：合并多线程 photon CSV**

每个 config 下有 `AnaEx01_nt_PhotonFaceBlockEvent_t0.csv ... _t7.csv`（线程数等于 OpenMP 设置）。
直接 concat 后 `Face × 16 + j × 4 + k` 压平成 0..95 的通道下标。
**和 flow3 完全一致**，唯一区别是这里不依赖 `merged_event.csv`，是直接吃原始线程文件。

**步骤 4：随机分组聚合（核心）**

```python
n_samples = n_rows // photons_per_sample
perm = rng.permutation(n_rows)
use = perm[: n_samples * photons_per_sample]
sample_idx = np.repeat(np.arange(n_samples), photons_per_sample)

X = np.zeros((n_samples, N_CHANNELS), dtype=np.float32)
np.add.at(X, (sample_idx, flat_ch[use]), counts[use])
```

- `rng.permutation` 把同 config 内的所有光子行打乱
- 截取整数倍的 N 个，按顺序切成 `n_samples` 块
- `np.add.at` 在 (sample_idx, channel) 二维下标处累加，向量化操作，比 for 循环快几十倍

**步骤 5：归一化（可选）**

```python
if normalize_counts:
    totals = X.sum(axis=1, keepdims=True)
    X = X / totals   # 每个样本变成 sum=1 的"分布"
```

打开后 NN 对总光子数完全免疫，只看 96 个通道的相对分布。**强烈推荐打开**。

**步骤 6：6 维标签广播 + 元数据保留**

```python
y_cfg = np.broadcast_to(y6_mm, (n_s, 6)).astype(np.float32).copy()
sep = float(np.linalg.norm(y6_mm[:3] - y6_mm[3:]))
```

同 config 内所有 sample 共享同一个 6 维标签。同时存了 `config_id` / `fraction_a` /
`separations` 三个辅助数组，留给将来做"按 separation 分桶分析"用。

> **不跨 config 混合**：每个 sample 严格属于一个 config，标签明确。如果跨 config 抽光子，
> 标签就没法定义了（半个 A 在某个 config，半个 A 在另一个 config）。

### 4.2 `model_double.py`

和 flow3 `model.py` 的代码 diff **只有一处**：输出层从 3 维改成 6 维。

```python
class DoublePointMLP(nn.Module):
    def __init__(self, in_dim=96, hidden=256, out_dim=6, dropout=0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, out_dim),    # ← flow3 是 nn.Linear(hidden, 3)
            nn.Tanh(),                     # ← Tanh 仍然把每维压到 (-1, 1)
        )
```

参数量：`96×256 + 256 + 256×256 + 256 + 256×6 + 6 = 92,422`。
比 flow3 的 91,395 多了 768（256×3 个新权重 + 3 个 bias）—— 几乎可以忽略。

**为什么模型本身不内嵌 permutation invariance？**

理论上可以做"对称化网络"：分别预测 (A, B) 后用 `min` 或 `sort` 强行规范顺序。
但这会破坏 MLP 的可微性，且推理阶段还是要做"和真值匹配"。所以这套方案的策略是：

- **模型保持简单 MLP**，自由输出 6 维
- **训练 loss 用 min(identity, swap)** 让模型自动学到"两点是个集合"
- **评估和推理时再做 best-match 对齐**

### 4.3 `train_double.py` —— 核心：置换不变损失

这是 flow4 唯一**真正的新概念**。其它部分（AdamW、ReduceLROnPlateau、早停、GPU 常驻数据）
和 flow3 `train2.py` 完全一致。

#### 4.3.1 问题：(A, B) 顺序怎么定？

仿真器 `metadata.csv` 里写 A 是某点、B 是另一点。但**模型怎么知道**应该把"A 槽位"输出为 A 还是 B？

物理上 A 和 B **完全对称**（都是相同的伽马源，只是位置不同）。如果训练时硬要求"模型输出
前 3 维必须等于真值 A"，那模型就得**记住每个具体 config 的 A↔B 标签**。这件事神经网络
学不会（因为输入只有光子分布，没有任何 A↔B 区分信息），梯度会乱跳。

#### 4.3.2 解决：min(identity, swap)

```python
def perm_invariant_loss(pred, target, beta=1.0):
    pa, pb = pred[:, :3], pred[:, 3:]
    ya, yb = target[:, :3], target[:, 3:]

    # identity match: pred A ↔ target A, pred B ↔ target B
    l_id_a = F.smooth_l1_loss(pa, ya, beta=beta, reduction="none").sum(dim=1)
    l_id_b = F.smooth_l1_loss(pb, yb, beta=beta, reduction="none").sum(dim=1)
    l_id = l_id_a + l_id_b

    # swap match: pred A ↔ target B, pred B ↔ target A
    l_sw_a = F.smooth_l1_loss(pa, yb, beta=beta, reduction="none").sum(dim=1)
    l_sw_b = F.smooth_l1_loss(pb, ya, beta=beta, reduction="none").sum(dim=1)
    l_sw = l_sw_a + l_sw_b

    per_sample = torch.minimum(l_id, l_sw)   # element-wise min
    return per_sample.mean()
```

公式版：

```
loss(pred, target) = mean_i [ min(
        SmoothL1(pa_i, ya_i) + SmoothL1(pb_i, yb_i),
        SmoothL1(pa_i, yb_i) + SmoothL1(pb_i, ya_i)
) ]
```

**关键技术点**：

- `reduction="none"` → 拿到逐元素 (B, 3) loss
- `.sum(dim=1)` → 沿坐标维度求和，得到逐样本的"两点总误差" (B,)
- `torch.minimum(l_id, l_sw)` → **element-wise** 取小者（不是 `min(scalar, scalar)`）
- 最后 `.mean()` → 才得到反向传播用的标量 loss

每个样本独立选 identity 或 swap：有些样本 (A, B) 顺序"恰好对齐"训练标签，另一些恰好反过来，
loss 自动选对那一边。模型不会因为预测出"反过来"的顺序而被错误惩罚。

#### 4.3.3 为什么不能用 mean 而要用 min？

- `mean(l_id, l_sw)` = "两种排列都得对" → 模型被强迫学一个不可学的目标 → 梯度自相抵消，loss 卡死
- `min(l_id, l_sw)` = "至少一种排列对就行" → 给模型自由度，让它自己选最自然的顺序

**和图像识别里 set prediction 的关系**：DETR、object detection 里"模型预测 N 个 box，
ground truth 也是 N 个 box，怎么对应？" 用的是匈牙利算法（Hungarian matching）。
这里 K=2 太小，直接枚举 2! = 2 种排列即可。K=3 时就是 6 种（参见 §8）。

#### 4.3.4 sanity check

训练开始前会跑一段：

```python
def sanity_check_perm_invariance(device):
    pred = torch.randn(64, 6, device=device)
    target = torch.randn(64, 6, device=device)
    target_swapped = torch.cat([target[:, 3:], target[:, :3]], dim=1)
    l1 = perm_invariant_loss(pred, target).item()
    l2 = perm_invariant_loss(pred, target_swapped).item()
    print(f"[sanity] perm-invariance: loss={l1:.6f}  swapped-target loss={l2:.6f}  "
          f"diff={abs(l1 - l2):.2e}  (应≈0)")
```

把 target 的 (A, B) 整体调换后 loss 必须**完全不变**（diff = 0），这是置换不变性的定义。
如果你以后想魔改 loss，把这段当成回归测试。

#### 4.3.5 其他细节（与 flow3 train2.py 几乎一致）

- **数据一次性 `.to(cuda)`** + 手写 batch 循环 + `index_select` → 没有 H2D 拷贝
- **AdamW** 而非 Adam（带正确的 L2 实现）
- **ReduceLROnPlateau** patience=8、factor=0.5
- **早停** patience=30
- 归一化、80/10/10 划分、`split_idx.npz` 保存等套路同 flow3，参考 `flow3_NN/GUIDE.md §4.3`

### 4.4 `evaluate_double.py` —— best-match 对齐评估

#### 4.4.1 numpy 版 best-match

```python
def best_match(pred, true):
    pa, pb = pred[:, :3], pred[:, 3:]
    ya, yb = true[:, :3], true[:, 3:]
    err_id = np.linalg.norm(pa - ya, axis=1) + np.linalg.norm(pb - yb, axis=1)
    err_sw = np.linalg.norm(pa - yb, axis=1) + np.linalg.norm(pb - ya, axis=1)
    swap = err_sw < err_id
    out = pred.copy()
    out[swap, :3] = pb[swap]
    out[swap, 3:] = pa[swap]
    return out
```

逻辑和训练 loss 完全对偶：每个样本独立判断"哪种匹配总误差更小"，把 pred 的 A/B 槽位
换成那一种。**不做这一步，所有指标都会被 A↔B 顺序错位污染**。

#### 4.4.2 算什么

- **per-axis MAE/RMSE** for A 和 B（共 6 个数）
- **3D-RMSE** for A 和 B（每点的 √(mean(d²))）
- **separation 误差**：`|A_pred − B_pred| − |A_true − B_true|`
  - 这个数衡量"两点距离判得对不对"，比单独看 A、B 的 RMSE 多一层信息
  - 如果 A、B 都准但 separation 偏大 / 偏小 → 两点同时被推到外侧或内侧

#### 4.4.3 hexbin 而非 scatter

数据量上去（>1k 样本）后，scatter 会堆成黑色一坨；hexbin 用六角格显示密度，对角线偏离
更直观。`gridsize=60`、`mincnt=1` 是经验值。**和 flow3 同款**。

### 4.5 `predict_double.py` —— 单 config 推理

整体流程和 evaluate 镜像，区别在于**输入是单个 DP_* 目录**而非测试集。

**取均值的合理性**：每个样本是 N 个光子聚合的"γ 事件等价物"，多 sample 平均 = 抑制统计涨落。
如果一个 config 跑了 35000 events、`--photons-per-sample 2000` → 17 个 sample → 平均后
等价于看了 34000 个光子，比单 sample 误差小 √17 ≈ 4 倍。

**best-match 对齐到 ground truth**：当 metadata.csv 存在时，会把 (pa, pb) 重排成与
(ya, yb) 总误差更小的顺序。如果 metadata 缺失（真实实验数据），跳过对齐，输出原始预测。

`nn_double_predicted.csv` 列：`sample_id, ax_pred, ay_pred, az_pred, bx_pred, by_pred, bz_pred`。
这套预测可以直接拿去和 ground truth 散点图对比、或喂给下游分析脚本。

---

## 5. 与 flow3_NN 的关键差异

| 维度 | flow3_NN（单点） | flow4_NN_double（双点） |
|---|---|---|
| **输出维度** | 3 (x, y, z) | **6** (ax, ay, az, bx, by, bz) |
| **Loss** | 普通 SmoothL1 | **置换不变** `min(identity, swap)` |
| **训练数据来源** | 单 root（一个 `Output/SiPM6_Output_*/`） | **多 root**（`--input-roots A B C ...` 合并多档 separation） |
| **单位** | cm | **mm** |
| **y_scale** | 1.25 cm | **12.5 mm**（同 = 12.5 mm，只是单位换） |
| **典型样本量** | 百万级（千万 events / 5000 photons） | **万级**（千 events × 几十 config / 2000 photons） |
| **过拟合风险** | 低 | **高**（数据少 100×，要更猛的 dropout / weight_decay） |
| **评估** | 直接比 (x, y, z) | **best-match 对齐后**才能比 |
| **推理输出 CSV** | 行级 (event → x, y, z) | **样本级** 6 维，外加 mean prediction |
| **对称镜像推理** | predict2 / predict3 利用立方体对称翻转到 +++ 卦限 | **没做** —— 训练数据已覆盖整个立方体（±12.5 mm 全空间随机），用不上 |
| **上游 split + Histo 步骤** | 必须 | **跳过** —— build 直接读 `_t*.csv` |
| **特殊 metadata 字段** | `source.fp_source_*` (cm) | 多 `source.fp_source_b_*_mm` + `source.fraction_a` |
| **额外辅助标签** | 无 | `fraction_a` / `separations` / `config_id` 都存进 npz，留给将来分析 |

**特别要记住的两条**：

1. **flow4 单位是 mm**。所有命令行参数、metadata 字段、图轴标签都是 mm。**不要混 cm**。
2. **每次改 build 参数都要同步改 predict 参数**。`--photons-per-sample` 和 `--normalize-counts`
   两边必须一致，否则输入尺度差几个量级，Tanh 直接饱和。

---

## 6. 关键概念词典

> 通用 NN 术语（epoch / batch / loss / dropout / weight_decay / 早停 / cosine schedule / ...
> 等等）见 `workflow/flow3_NN/GUIDE.md §5`，本节只列**双点 / set prediction 专属**词条。

| 词 | 意思 |
|---|---|
| **置换不变 (permutation-invariant)** | 函数 `f(A, B) = f(B, A)` —— 输入元素顺序不影响输出。本项目的 loss 就是这种 |
| **set prediction** | 模型输出"一组无序元素"而非"一个有序列表"。双点重建是 K=2 的 set prediction |
| **identity match / swap match** | K=2 时只有两种排列：原顺序 和 交换 (A, B)。loss 取这两个的 min |
| **best-match alignment** | 推理后，把模型输出的 (A, B) 顺序重排成"和真值总误差最小"的那种。仅用于评估和打印，不参与训练 |
| **Hungarian / 二分图匹配** | 匈牙利算法，K ≥ 3 时找最优排列的 O(K³) 算法。K=2 时退化为枚举 2! = 2 |
| **Chamfer distance** | 另一种 set-to-set 距离：每个 pred 找最近 truth + 每个 truth 找最近 pred。K=2 时和 best-match 等价 |
| **fraction_a** | 仿真里 A 点发射光子的概率（B 是 1 − fraction_a），范围通常 [1/3, 2/3] |
| **二项分布涨落** | 35000 events × p=0.5 → 实际 A 命中数 ≈ 17500 ± √(35000·0.5·0.5) ≈ 17500 ± 93。这是物理上不可避免的"哪个点亮多少"随机性 |
| **separation** | 两点 3D 距离 `‖A − B‖`，flow4 训练数据里的 sep 范围一般是 0.1 ~ 20 mm |
| **separation 误差** | 模型预测的 separation 和真实 separation 之差，是除 per-point 误差外的第二维评价指标 |
| **y_scale_mm** | 标签归一化常数（默认 12.5 mm，= 晶体半边长）。Tanh 输出 (-1, 1) 乘上它得到 mm |

---

## 7. 排错 checklist

### 通用问题

复用 flow3 的 checklist：CUDA 不可用、loss NaN、Windows DataLoader 卡死、推理饱和等等。
详见 `flow3_NN/GUIDE.md §6`。

### 双点专属

#### 训练 loss 卡在 ~y_scale²/2 不下降

**症状**：epoch 1 loss = 0.5，epoch 50 loss = 0.49，几乎不动。

可能原因：
- ☐ **模型只学到"输出 (0, 0, 0, 0, 0, 0)"**：`y_test.var()` ≈ 0.5，预测全 0 时 loss 就是这个值。
  检查归一化对不对（`y_train` 应在 [-1, 1] 内、X 标准化到方差 1 附近）
- ☐ **数据集里 fraction_a 极端**（绝大多数 ≈ 0.99 或 0.01）：模型只学到主导那个点，
  另一个点的梯度被淹没。重跑 Stage 0 把 `--ratio-min/--ratio-max` 收紧到 [0.4, 0.6]
- ☐ sanity check 没过（`diff` 不是 0）：loss 实现有 bug，回 §4.3.4 对一遍

#### predict 出来 A 和 B 总是相同

**症状**：mean prediction 输出 `A = B = (0.x, 0.x, 0.x)`。

可能原因：
- ☐ **模型缺少区分能力**：hidden 太小（试 `--hidden 512`），或输入归一化把所有差异抹平
- ☐ **训练数据 separation 普遍 < 1 mm**：模型从来没见过分得开的两点，自然学到"两点重合"
  的捷径。**多加几档大 separation 的 sweep 文件夹**
- ☐ Dropout=0.5 之类太大 → 模型容量被打残，回到 0.1

#### evaluate 时一半样本对得很准、另一半离谱

**症状**：散点图上 50% 点贴对角线、50% 点贴反对角线（即 ax_pred 和 bx_true 相关）。

**多半是没用 best-match alignment**！pred (A, B) 直接和 truth (A, B) 比，但模型输出的
A↔B 顺序对不上 → 一半样本恰好对齐、一半完全反过来。

确认 `evaluate_double.py` 里调了 `best_match(y_pred, y_test)` —— 当前脚本是有的。
如果你魔改了脚本忘了调，加回来。

#### 多 sweep 合并训练后某档表现特别差

**症状**：在 sep0p1mm 上 RMSE = 5 mm，在 sep5mm 上 RMSE = 0.8 mm。

正常现象：**近重合的两点物理上就难分**（光子分布几乎一样）。可以做的：
- 把不同 separation 分开评估：用 `dataset.npz` 里的 `separations` 数组分桶画 RMSE 曲线
- 训练时不平衡：远 sep 的 config 多、近 sep 的 config 少，让模型偏向"分得开"的场景
- 接受这是物理极限：< 1 mm 的两点本来就只能拟合到加权重心，无解

#### Tanh 饱和（推理预测都是 ±12.5）

完全等价于 flow3 的同名问题。99% 是 build 和 predict 的 `--photons-per-sample` /
`--normalize-counts` 不一致。再读一遍 §3 Step 1 末尾的警告。

---

## 8. 下一步可以怎么做

按"收益 / 工作量"排：

1. **多任务学习：把 `fraction_a` 当辅助 head**
   - 在 model 里加一个 `Linear(hidden, 1) + Sigmoid` 输出 fraction_a 估计
   - loss 加一项 `BCE(pred_frac, true_frac)`
   - 物理意义：fraction_a 是"谁亮谁暗"的指示，预测它能让模型把"主点"和"副点"显式区分
2. **K=3+ 的 set prediction → Hungarian matching**
   - 当前 min(2!) 暴力枚举只对 K=2 work。K=3 要 O(K³) 的匈牙利算法（`scipy.optimize.linear_sum_assignment`）
   - 用 detr 那种 set 损失：对每个 pred 找最近 truth，全集对应再求和
3. **对称增强（48 倍数据）**
   - 立方体的 48 种对称操作对输入和 (A, B) 标签做同样的变换
   - flow4 数据量小 100×，对称增强能直接让有效样本数 × 48
   - 实现要点：变换矩阵作用到 face/j/k 重排同时作用到 A、B 两点坐标
4. **CNN 重新组织 96 维输入**
   - 把 96 维 reshape 成 `(6, 4, 4)` 喂小 CNN（或 transformer attention）
   - 4×4 太小 CNN 收益有限，但 6 个 face 之间的"对面相关性"用 attention 可能更好捕捉
5. **预测 (中点, 间距向量) 6 维 vs (A, B) 6 维**
   - 中点 `m = (A + B) / 2` 和差向量 `d = A − B` 也是 6 维表示
   - 中点没有置换歧义、差向量有 `±d` 的 ±1 歧义
   - 损失变成 `SmoothL1(m_pred, m_true) + min(SmoothL1(d_pred, d_true), SmoothL1(d_pred, -d_true))`
   - 实验对比哪种参数化更稳
6. **加入 separation 范围过滤的 curriculum learning**
   - 先训远 separation 的（容易），逐步加进近 separation 的（难）
   - 比一次性混合训练收敛更稳
7. **不确定性估计**
   - 输出 `(mu_A, sigma_A, mu_B, sigma_B)` 共 12 维，loss 用置换不变 NLL
   - 物理意义：模型告诉你"我对这次定位有多大把握"，工程上很有用

每做一次改动：
- 固定 seed (`--seed 42`) 和 split (`split_idx.npz` 不重生成)
- 对比同一测试集 RMSE_A / RMSE_B / separation MAE
- 把 `test_metrics.csv` 另存到 `artifacts/exp_<日期>_<描述>/` 备查

---

## 附：最常见的命令

```powershell
# 0. Stage 0 上游仿真（一次双点扫描，按 separation 分档保存）
uv run python workflow/flow4_NN_double/run_batch_double_point.py `
    --random-count 100 `
    --min-separation 2 --max-separation 8 `
    --random-seed 1234 `
    --results-dir Results_flow4_trainingdata_sep5mm

# 1. 构建数据集（多 root 合并，推荐设置：聚合 + 归一化）
uv run python workflow/flow4_NN_double/build_dataset_double.py `
    --input-roots Results_flow4_trainingdata `
                  Results_flow4_trainingdata_sep5mm `
                  Results_flow4_trainingdata_sep3mm `
                  Results_flow4_trainingdata_sep1mm `
                  Results_flow4_trainingdata_sep0p1mm `
    --photons-per-sample 2000 --normalize-counts

# 2. 训练
uv run python workflow/flow4_NN_double/train_double.py --epochs 300

# 3. 评估（best-match 对齐后算 MAE/RMSE/separation）
uv run python workflow/flow4_NN_double/evaluate_double.py

# 4. 单 config 推理
uv run python workflow/flow4_NN_double/predict_double.py `
    --config-dir Results_flow4_trainingdata_sep5mm/DP_0001_S0p3_..._<ts> `
    --photons-per-sample 2000 --normalize-counts

# 5. 看模型结构 / 参数量
uv run python workflow/flow4_NN_double/model_double.py
```

**记住**：训练时传了 `--normalize-counts`、`--photons-per-sample 2000`，predict 命令必须
**完全一致**。否则 Tanh 直接饱和，结果全是 ±12.5 mm。

祝玩得开心 —— 看到 `[sanity] perm-invariance: ... diff=0.00e+00` 的那行打印出来时，
就知道你已经把"两点是个集合"这件事用 4 行 PyTorch 教给了模型。
