# 神经网络回归做位置重建 —— 从零上手笔记

**创建日期**：2026-04-18
**场景**：立方体 GAGG 晶体 + 6 面 × 4×4 SiPM（共 96 路）读出。已有 Geant4 仿真生成的扫描数据（真值位置 → 每面 4×4 光子计数）。目标：用神经网络直接从 96 维 SiPM 响应回归 3D 入射位置 `(x, y, z)`。

---

## 0. 为什么 NN 比公式法强

- **公式法（linear_scaled、half_side_ratio 等）** 把 96 个数压成 2–6 个求和再做比值，**信息损失极大**。
- 每面的 4×4 分布里藏着入射点的横向信息（类似图像重心），NN 能自动学会怎么用。
- 光学响应在边缘非线性很强，NN 是**万能函数逼近器**，能把 S 曲线、边缘爆炸这些 hand-crafted 公式搞不定的非线性自动拟合掉。
- 一旦训练完，推理就是一次矩阵乘法，**比最小二乘前向模型快得多**。

代价：需要一份"带真值标签"的数据集 —— 这个你已经有了（就是你扫描网格跑出来的仿真）。

---

## 1. 大局观（你要理解的 5 件事）

### 1.1 监督学习的基本流程

```
原始仿真数据 (CSV)
    ↓  解析 + 组织成 (X, y)
数据集 (X: 输入特征, y: 真值标签)
    ↓  划分
训练集 / 验证集 / 测试集
    ↓  喂给网络
模型 f_θ(X) ≈ y
    ↓  评估
测试集 MAE / RMSE / 分辨率
```

- **X**：输入特征，这里 shape = `(N_events, 96)` 或 `(N_events, 6, 4, 4)`。
- **y**：真值标签，shape = `(N_events, 3)`（x, y, z 三个坐标）。
- **θ**：网络内部的权重参数（几千到几万个数），训练就是找到让预测最接近真值的 θ。

### 1.2 张量（Tensor）

就是多维数组，和 numpy 的 `ndarray` 几乎一样，但能在 GPU 上跑、能自动求导。

- `tensor.shape` 像 numpy 的 `.shape`
- `tensor.to("cuda")` 搬到 GPU
- PyTorch 和 numpy 互转：`torch.from_numpy(arr)`、`t.numpy()`

### 1.3 前向传播 & 反向传播

- **前向**：输入 X 经过每一层，算出预测 `ŷ`。
- **损失 Loss**：`L(ŷ, y_true)`，数值越小越好。
- **反向**：PyTorch 自动算 `∂L/∂θ`（链式法则，叫 autograd）。
- **优化器**（Optimizer，如 Adam）：根据梯度把 θ 往变小的方向挪一小步。

一次这样的循环叫一个 **iteration**；把整份训练集过一遍叫一个 **epoch**。通常训练几十到几百个 epoch。

### 1.4 过拟合 vs 欠拟合

- **欠拟合**：训练集 loss 都降不下去 → 模型太小或训练不够。
- **过拟合**：训练集 loss 很低，但验证集 loss 不降反升 → 模型背下了训练集。
- **应对**：更多数据、更小模型、正则化（dropout、weight decay）、early stopping。

### 1.5 归一化（最容易忽略但最关键）

- 输入 X：每路 SiPM 计数量级差很大 → **对每个特征按训练集均值/方差标准化**：`X' = (X - μ) / σ`。
- 输出 y：真值坐标 ∈ [−1.25, 1.25] cm → 归一化到 [−1, 1]：`y' = y / 1.25`。
- **必须保存 μ、σ、scale**，推理阶段要用同一套。

---

## 2. 选 PyTorch（不用犹豫）

理由：
- 社区大、资料多、调试直观；
- 动态图，错误信息人类可读；
- Geant4/物理圈子基本都用 PyTorch；
- 你后期想迁到 TensorFlow 随时可以。

### 2.1 安装（在你的 conda/venv 里）

```bash
# CPU 版，完全够用（96→3 的小网络 CPU 秒级）
pip install torch torchvision

# 有 NVIDIA 显卡想用 GPU（看 CUDA 版本，例 12.1）
pip install torch --index-url https://download.pytorch.org/whl/cu121
```

验证：
```python
import torch
print(torch.__version__, torch.cuda.is_available())
```

---

## 3. 数据准备（第一大关）

### 3.1 现有数据长什么样

你的 `Output/<config>/merged_face_jk.csv` 长这样：

```
Face,j,k,Count
0,0,0,123
0,0,1,156
...
5,3,3,88
```

每个配置目录对应**一个真值位置**（从 `metadata.csv` 的 `source.position_*` 能读到）。

### 3.2 目标数据集格式

组织成两个 numpy 数组：

```python
# X.shape == (N_samples, 96)         —— flatten 后的 SiPM 响应
# y.shape == (N_samples, 3)          —— 对应 (x_true, y_true, z_true) in cm
```

若用 CNN，则 `X.shape == (N_samples, 6, 4, 4)`。MLP 起步够用，先走 96 维。

### 3.3 构造脚本（伪代码）

```python
from pathlib import Path
import numpy as np
import pandas as pd

def face_jk_to_vec(agg_csv: Path) -> np.ndarray:
    """返回 96 维向量，顺序 face 0..5 × j 0..3 × k 0..3"""
    df = pd.read_csv(agg_csv)
    vec = np.zeros((6, 4, 4), dtype=np.float32)
    for _, r in df.iterrows():
        vec[int(r.Face), int(r.j), int(r.k)] = r.Count
    return vec.reshape(-1)  # 96

def read_true_position(meta_csv: Path) -> np.ndarray:
    meta = dict(zip(*pd.read_csv(meta_csv).values.T.astype(str)))
    return np.array([
        float(meta["source.position_x_mm"]) / 10.0,  # mm → cm
        float(meta["source.position_y_mm"]) / 10.0,
        float(meta["source.position_z_mm"]) / 10.0,
    ], dtype=np.float32)

def build_dataset(output_root: Path):
    X, y = [], []
    for cfg in sorted(output_root.iterdir()):
        agg = cfg / "merged_face_jk.csv"
        meta = cfg / "metadata.csv"
        if not agg.exists() or not meta.exists():
            continue
        X.append(face_jk_to_vec(agg))
        y.append(read_true_position(meta))
    return np.stack(X), np.stack(y)
```

### 3.4 数据集要多大？

- MLP 96→3：**1 k–10 k 样本** 就能看到效果；
- CNN：建议 **10 k+** 才能稳；
- 训练样本就是"不同入射位置的配置目录数量" —— 所以**扫描网格密度越大越好**。若你现在只有几百个配置目录，先考虑加密扫描（或者用每个目录内的**单事件**作样本，而不是整个目录聚合 —— 这样每个目录能拆成几百上千个样本）。

**关键点**：单事件级数据信息量大得多，但噪声也大 —— 这才是 NN 真正能发挥的战场。公式法对单事件数据更难处理，NN 能学到事件级的分布模式。

### 3.5 划分训练/验证/测试

```python
from sklearn.model_selection import train_test_split
X_trainval, X_test, y_trainval, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
X_train,    X_val,  y_train,    y_val   = train_test_split(X_trainval, y_trainval, test_size=0.1, random_state=42)
```

比例常用 80/10/10 或 70/15/15。**测试集要一直锁起来到最后才看**。

### 3.6 标准化

```python
mean = X_train.mean(axis=0)
std  = X_train.std(axis=0) + 1e-6   # 防除零
X_train = (X_train - mean) / std
X_val   = (X_val   - mean) / std
X_test  = (X_test  - mean) / std

y_scale = 1.25   # 半长 cm
y_train_n = y_train / y_scale
y_val_n   = y_val   / y_scale
y_test_n  = y_test  / y_scale

# 保存 mean/std/y_scale 到 npz，推理必须复用
np.savez("norm.npz", mean=mean, std=std, y_scale=y_scale)
```

---

## 4. 写模型（MLP 版，60 行）

### 4.1 完整最小可运行代码

```python
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

# ---------- 1. 包成 PyTorch Dataset ----------
def to_loader(X, y, batch_size=128, shuffle=True):
    ds = TensorDataset(torch.from_numpy(X).float(), torch.from_numpy(y).float())
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle)

train_loader = to_loader(X_train, y_train_n, shuffle=True)
val_loader   = to_loader(X_val,   y_val_n,   shuffle=False)

# ---------- 2. 模型 ----------
class PositionMLP(nn.Module):
    def __init__(self, in_dim=96, hidden=256, out_dim=3, dropout=0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, out_dim),
            nn.Tanh(),  # 输出强制落到 [-1, 1]（因为 y 已归一化）
        )
    def forward(self, x):
        return self.net(x)

device = "cuda" if torch.cuda.is_available() else "cpu"
model = PositionMLP().to(device)

# ---------- 3. 损失 + 优化器 ----------
criterion = nn.SmoothL1Loss()            # Huber loss，对离群点比 MSE 稳
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.5, patience=10)

# ---------- 4. 训练循环 ----------
best_val = float("inf")
for epoch in range(200):
    # 训练
    model.train()
    train_loss = 0.0
    for xb, yb in train_loader:
        xb, yb = xb.to(device), yb.to(device)
        optimizer.zero_grad()
        pred = model(xb)
        loss = criterion(pred, yb)
        loss.backward()
        optimizer.step()
        train_loss += loss.item() * xb.size(0)
    train_loss /= len(train_loader.dataset)

    # 验证
    model.eval()
    val_loss = 0.0
    with torch.no_grad():
        for xb, yb in val_loader:
            xb, yb = xb.to(device), yb.to(device)
            val_loss += criterion(model(xb), yb).item() * xb.size(0)
    val_loss /= len(val_loader.dataset)
    scheduler.step(val_loss)

    if val_loss < best_val:
        best_val = val_loss
        torch.save(model.state_dict(), "best.pt")

    if epoch % 10 == 0:
        print(f"epoch {epoch:3d} | train {train_loss:.4f} | val {val_loss:.4f}")

print(f"best val loss: {best_val:.4f}")
```

### 4.2 每一行在干什么

| 代码 | 含义 |
|---|---|
| `nn.Linear(96, 256)` | 全连接层，96 维输入 → 256 维输出 |
| `nn.ReLU()` | 激活函数，让网络有非线性能力 |
| `nn.Dropout(0.1)` | 训练时随机丢 10% 神经元，防过拟合 |
| `nn.Tanh()` | 最后一层，把输出压到 [−1, 1] —— 配合归一化后的 y |
| `model.train()` / `model.eval()` | 切换训练/推理模式（影响 dropout、batchnorm） |
| `optimizer.zero_grad()` | 清空上一轮梯度，不清会累加 |
| `loss.backward()` | 自动微分，算出所有参数的梯度 |
| `optimizer.step()` | 用梯度更新参数 |
| `with torch.no_grad():` | 推理时关掉 autograd，省内存加速 |

### 4.3 常见超参数起点

| 名字 | 起点 | 说明 |
|---|---|---|
| hidden 层宽 | 256 | 太小欠拟合，太大过拟合 |
| 层数 | 3（两个隐藏层） | 96→3 足够 |
| batch_size | 128 | 大 batch 稳但收敛慢 |
| lr | 1e-3 | Adam 默认 |
| dropout | 0.1 | 验证集在抖就加 |
| weight_decay | 1e-5 | L2 正则 |
| epochs | 200 + early stop | 看 val loss 平台期 |

---

## 5. 评估

```python
model.load_state_dict(torch.load("best.pt"))
model.eval()
with torch.no_grad():
    X_test_t = torch.from_numpy(X_test).float().to(device)
    y_pred_n = model(X_test_t).cpu().numpy()
y_pred = y_pred_n * y_scale   # 反归一化回 cm

err = y_pred - y_test
print("MAE  per axis (cm):", np.mean(np.abs(err), axis=0))
print("RMSE per axis (cm):", np.sqrt(np.mean(err**2, axis=0)))
print("3D RMSE (cm):", np.sqrt(np.mean(np.sum(err**2, axis=1))))
```

与 linear_scaled / half_side_ratio 同一个测试集对比，画同样的 scatter 图直观感受。

---

## 6. 进阶：CNN 版（等 MLP 跑通再上）

6 个面每个是 4×4 图像，适合 CNN。保持 shape 为 `(N, 6, 4, 4)`：

```python
class PositionCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(6, 32, 3, padding=1), nn.ReLU(),
            nn.Conv2d(32, 64, 3, padding=1), nn.ReLU(),
            nn.Flatten(),
        )
        self.head = nn.Sequential(
            nn.Linear(64 * 4 * 4, 128), nn.ReLU(),
            nn.Linear(128, 3), nn.Tanh(),
        )
    def forward(self, x):
        return self.head(self.conv(x))
```

4×4 太小，实际上 MLP 和 CNN 差距不会很大。CNN 更有意义是在 6 面拼成单张"展开图"或 SiPM 阵列更大时。

---

## 7. 坑和 checklist

**训练诊断**

- [ ] Loss 第一个 epoch 就 NaN → 学习率太大 or 输入没归一化
- [ ] train loss 降、val loss 涨 → 过拟合，加 dropout / 减宽
- [ ] train loss 就是不降 → 模型太小 or lr 太小 or 数据没归一化
- [ ] val loss 震荡剧烈 → batch size 太小、数据噪声大
- [ ] 推理精度比训练时差很多 → 归一化统计量没存下来/推理用错了

**物理上的坑**

- 真值位置标签来源要**绝对可信**；读 metadata 的单位（mm vs cm）检查两次。
- 训练集必须**覆盖你想推理的位置范围**：NN 对外推（extrapolation）很烂，你扫描网格边界在哪里，NN 就只在那个范围里靠谱。
- 训练集里 `(x, y, z)` 的分布应尽量均匀，不均匀会让网络偏向高密度区。
- 对称性：立方体几何有对称性，数据增强可以镜像翻转（同时翻输入和 y 标签），等效于 8× 数据。

**复现性**

```python
import random
random.seed(42); np.random.seed(42); torch.manual_seed(42)
```

---

## 8. 你要系统学习的关键点（按优先级）

### 必学（做完上面代码就懂一半）
1. **张量操作**：shape、reshape、view、permute、广播
2. **Dataset / DataLoader**：为什么要 batch、shuffle
3. **损失函数**：MSE / SmoothL1 / 什么时候用哪个
4. **Adam 优化器**：学习率是什么，learning rate schedule
5. **训练/验证/测试的划分**：为什么要三分、为什么不能用测试集调参
6. **归一化**：为什么必须、训练和推理要一致

### 次学（训得不理想的时候一定要学）
7. **过拟合 & 正则化**：dropout、weight decay、early stopping
8. **激活函数**：ReLU / GELU / Tanh 的区别和选择
9. **BatchNorm / LayerNorm**：稳定训练
10. **梯度诊断**：grad norm、是否爆炸/消失

### 再进阶（想做到 SOTA 必学）
11. **CNN 基础**：卷积、池化、感受野
12. **数据增强**：物理对称性利用
13. **贝叶斯 / 不确定性估计**：MC Dropout、Ensemble —— 对物理测量很重要
14. **前向模型 + NN 混合**：把 NN 当"残差修正器"叠在公式法之上

### 推荐学习材料
- **Andrej Karpathy 的 "Neural Networks: Zero to Hero"**（YouTube） —— 从头搓 NN，理解最彻底
- **PyTorch 官方 60-minute blitz** —— 极速上手
- **《深度学习入门》（斋藤康毅）** —— 中文友好，讲清原理
- **fast.ai 课程** —— 实战导向

---

## 9. 推荐实施顺序

1. **Day 1**：跑通上面 60 行 MLP 代码，哪怕输入是随机数据 —— 先保证环境没问题
2. **Day 2**：写数据组装脚本，把你的仿真结果变成 `(X, y)`，保存成 `.npz`
3. **Day 3**：训 MLP，跟 linear_scaled 对比 MAE
4. **Day 4–5**：调网络结构、超参数，做消融
5. **Day 6**：引入事件级数据（不是聚合级），重做
6. **Day 7**：CNN / 对称性数据增强 / 集成

**卡住时的经验法则**：永远先问"数据对不对？归一化对不对？label 单位对不对？"—— 90% 的 NN 问题是数据问题，不是模型问题。

---

## 10. 最终产出清单

做完之后你的代码仓库应该多出这些文件：

```
workflow/nn_recon/
├── build_dataset.py        # 从 Output/ 生成 dataset.npz
├── train.py                # 训练主脚本，产出 best.pt + norm.npz
├── model.py                # 模型定义（MLP、CNN）
├── predict.py              # 推理：接受 merged_face_jk.csv → 输出 (x,y,z)
├── evaluate.py             # 在测试集上算 MAE/RMSE，画 scatter
└── dataset.npz             # 缓存的训练数据
```

然后在 `Histo10_Cubic.py` 的 `compute_reconstruction_rows` 里新增一个 `("nn_mlp", predict_nn(...))`，和其它 7 个算法并列比较。

---

## 11. 事后复盘（2026-04-18 全天迭代）

实际落地时踩到几个上面没提的坑，记录如下。最终脚本见 `workflow/flow3_NN/`，上手文档见 `workflow/flow3_NN/GUIDE.md`。

### 11.1 数据坑：1 个 EventID = 1 个光子，不能直接当样本

我们仿真用的是 `/source/mode optical`，每个 Geant4 event 只产 1 个光子。所以 `merged_event.csv` 里每行只有 `PhotonCount=1`，整个文件 EventID 稀疏、总行数巨大。

如果直接把每个 EventID 当样本：
- 样本数暴涨到 800 多万，每个样本只有 1 路 SiPM 非零
- 信息量极少，NN 学不到任何空间分布
- 第一版训练 3D RMSE 卡在 0.55 cm，loss 曲线几乎水平

**解决**：`build_dataset.py` 加 `--photons-per-sample N` 选项，用 `np.add.at` 把同一 config 内 N 个光子行向量化累加到一个样本。N=5000 对应单次 γ 事件典型量级。加上后 3D RMSE 立刻降到 0.028 cm。

**实现关键**（速度差几十倍）：
```python
sample_idx = np.repeat(np.arange(n_samples), photons_per_sample)
np.add.at(X, (sample_idx, flat_ch[use]), counts[use])
```

### 11.2 数据坑：光子总数在训练/推理之间会变

仿真时每个 config 光子数固定（比如 5000），但真实事件的光子数从几百到几万都有可能。如果模型输入是**原始计数**，推理时数量级一变 Tanh 立刻饱和。

**解决**：`--normalize-counts` 把每个样本除以总光子数，转成"光子**分数**"（每样本 sum=1）。模型输入只看分布形状，对总光子数免疫。加了这个选项的 artifacts 目录用后缀 `_norm` 标记。

### 11.3 推理坑 1：predict 漏传 `--normalize-counts`

训练和推理的预处理**必须完全一致**。我们遇到一次超级典型的 bug：
- 训练用的 artifacts 是 `artifacts_p5000_norm`（含归一化）
- 但推理命令没加 `--normalize-counts`
- 现象：任何 config 的预测都是 (±1.25, ±1.25, ±1.25)，无论真值

原因：原始计数 vs 归一化分数差 1000×，第一层激活爆炸 → 最终 Tanh 彻底饱和到 ±1 → 输出永远卡在 ±y_scale。**这种"整齐极值"的症状记下来，以后见到就知道问题在哪里**。

### 11.4 推理坑 2：训练只在一个卦限，其它卦限全失效

训练数据扫描只覆盖 +X/+Y/+Z 卦限（y ∈ [0, 10] mm）。第一版 `predict.py` 不做任何处理，负卦限全错。

**解决**：立方体几何关于三轴完全对称，利用这一点做"对称镜像推理"——
1. 判源所在卦限：比较对面 SiPM 总光子数，`sign_x = sign(N(+X face) - N(-X face))`
2. 把输入 96 维通道按几何镜像重排成 +++ 等价形式（`predict2.py._build_mirror_perm`）
3. 模型在 +++ 卦限推理
4. 输出坐标按 signs 反翻

通道置换的几何依据（重要，写脚本时反复 debug）：
- X 镜像：swap face 0↔1；face 2/3/4/5 上 **j 方向**翻转（因为 j 轴沿 +X）
- Y 镜像：swap face 2↔3；face 0/1 翻 j；face 4/5 翻 k
- Z 镜像：swap face 4↔5；face 0/1/2/3 翻 k

**自检方法**：对 +++ 卦限的 config，`predict.py` 和 `predict2.py` 结果必须完全一致（identity case）。

### 11.5 训练加速：train2.py 把每 epoch 从几十秒压到 1–2 秒

`train.py` 用标准 DataLoader，每 batch 都要 H2D 拷贝 + Python 调度，GPU 利用率低得可怜。

`train2.py` 做了三件事：
1. 整个 dataset 一次性 `.to(cuda)` 常驻显存
2. 不用 DataLoader，手写 `X.index_select(0, sel)` 切 batch
3. 默认 batch=8192（train.py 是 512），配合 lr=2e-3（大 batch + 大 lr 的经验法则）

显存足够时 train2 完胜，RTX 5060 上百万级样本 < 1 分钟训完。

### 11.6 批量推理 + 汇总：predict3.py

为了和 `flow2/analyze_position_accuracy.py` 的输出格式对齐，加了 `predict3.py`：输入一个总文件夹，自动遍历所有含 `merged_event.csv` 的子目录，跑对称镜像推理，最后写：
- `accuracy_summary_nn.csv`：逐 config 的真值/预测/误差 + 末尾按算法分组的 BIAS/STD
- `scatter_nn_mlp_sym.png`：XY/XZ/YZ 三视图 × 各深度切片（格式同 flow2）

单 config 用 `predict2.py`，整个实验批次用 `predict3.py`。

### 11.7 还没解决的问题：训练网格太稀疏

当前训练位置是 `{0, 2, 4, 6, 8, 10} mm` 的 6×6×6 = 216 个离散点，如果测试位置在两个 grid 点之间（比如 6.25 mm），模型行为不稳定。下一步可能的做法：
- Geant4 侧补更密的扫描
- 在 `build_dataset.py` 里加 `--mirror-augment`，用 8 卦限对称把数据 ×8
- 检查训练 y 分布是否真的触及 10 mm，必要时调 y_scale 贴合实际范围

留给下次。
