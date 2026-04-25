# uv：Python 项目管理上手笔记

**创建日期**：2026-04-18
**适用场景**：Windows + PowerShell，本项目（Scintillator-Detector-Continuous）
**目标读者**：从未用过 uv、目前靠裸 pip 管理 Python 的开发者

---

## 0. 为什么要切到 uv

### 起因：昨天踩的坑

- 机器上装了两个 Python：3.10（`AppData\Local\Programs\Python\python310\`）和 3.14（`AppData\Local\Python\pythoncore-3.14-64\`）。
- `pip install torch ...` 装到了 3.10，但运行脚本用的是 3.14 —— 装的东西白装。
- pip 看到"已有 torch"就跳过，`--index-url cu121` 根本没起作用，两个 Python 都是 `+cpu`。
- 没有 lockfile，换机器/换时间依赖版本对不齐。

这些问题**用 uv 一条命令解决**。

### uv 是什么

一句话：`uv = pyenv + venv + pip + pip-tools + poetry`，Rust 写的，速度是 pip 的 **10–100×**。

出自 Astral（ruff 的同公司），已是 2025 年 Python 工具链事实标准。

### uv 能解决的痛点

| 原来的问题 | uv 的解法 |
|---|---|
| pip 和 python 指向不同解释器 | 项目绑定一个 venv，`uv run` 自动用对的 Python |
| 已装了就不重装 | `uv add --reinstall`，状态可见 |
| 装包慢 | 并行 + 全局缓存，快 10–100× |
| 依赖不可复现 | 自动生成 `uv.lock` |
| 切 Python 版本麻烦 | `uv python install 3.11` 一条命令搞定 |

---

## 1. 安装 uv

### Windows PowerShell

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

装完**重开 PowerShell**（让 PATH 更新），然后验证：

```powershell
uv --version
```

应输出类似 `uv 0.5.x (xxxxxx)`。

### 如果 irm 被企业网络拦截

改用 winget：
```powershell
winget install --id=astral-sh.uv -e
```

或下载 zip：<https://github.com/astral-sh/uv/releases>

---

## 2. 核心心智模型（只需记住 3 点）

### 2.1 项目 = 一个目录 + `pyproject.toml` + `.venv`

- `pyproject.toml`：人类可读的依赖清单（取代 `requirements.txt`）
- `uv.lock`：机器生成的精确版本锁（取代手工维护）
- `.venv/`：专属虚拟环境，uv 自动管理，**不需要手动 `activate`**
- `.python-version`：锁定项目用哪个 Python 版本

### 2.2 `uv run <cmd>` 就是入口

原来你写 `python foo.py`，在 uv 项目里都换成：
```powershell
uv run python foo.py
```

uv 会自动：
1. 确认 Python 版本（没装就下载）
2. 确认依赖都装齐（没装就 sync）
3. 在 `.venv` 里执行命令

不用手动 activate，不用担心 PATH。

### 2.3 `uv add` 代替 `pip install`

```powershell
uv add numpy                 # 添加依赖到 pyproject.toml 并装好
uv add torch --index https://download.pytorch.org/whl/cu121
uv remove numpy              # 卸载并从 pyproject.toml 移除
```

**和 pip 的关键区别**：`uv add` 会同时更新 `pyproject.toml` 和 `uv.lock`，依赖变化对团队/未来的自己可见。

---

## 3. 本项目迁移步骤（实操）

### Step 1：在项目根目录初始化

```powershell
cd D:\Geant4\Projects\Scintillator-Detector-Continuous
uv init --python 3.11
```

效果：
- 生成 `pyproject.toml`、`.python-version`、`.gitignore`、`README.md`（若无）
- 如果系统没 3.11，uv 自动下一个到 `%USERPROFILE%\AppData\Roaming\uv\python\`（**不动现有 3.10/3.14**）

如果 `uv init` 说已有 `pyproject.toml`，跳过这步。

### Step 2：添加依赖

根据本项目目前 import 的包，逐步加：

```powershell
# 基础数据/绘图
uv add numpy pandas matplotlib scikit-learn

# PyTorch（CUDA 12.1 版；没 NVIDIA 卡就去掉 --index）
uv add torch torchvision --index https://download.pytorch.org/whl/cu121

# 其它你用到的
uv add scipy   # 如有需要
```

装完可以查看：
```powershell
uv pip list
```

### Step 3：跑你现有脚本

原来：
```powershell
python workflow\flow2\Histo10_Cubic.py Results
```

换成：
```powershell
uv run python workflow\flow2\Histo10_Cubic.py Results
```

### Step 4：让 VS Code 用这个环境

VS Code 右下角 Python interpreter → 选：
```
D:\Geant4\Projects\Scintillator-Detector-Continuous\.venv\Scripts\python.exe
```

之后点 ▶ 运行、调试、Jupyter notebook 全都走 uv 的 venv。

### Step 5：把 `.venv` 排除出 git

`.gitignore` 里加（uv init 默认已加）：
```
.venv/
__pycache__/
*.pyc
```

`pyproject.toml` 和 `uv.lock` **要提交**。

---

## 4. 常用命令速查

| 命令 | 作用 |
|---|---|
| `uv init` | 初始化项目 |
| `uv add <pkg>` | 加依赖（装包 + 更新 lock） |
| `uv add <pkg> --dev` | 加开发依赖（不进生产） |
| `uv remove <pkg>` | 卸载依赖 |
| `uv sync` | 按 `uv.lock` 装全部依赖（新机器克隆代码后就跑这个） |
| `uv run <cmd>` | 在项目环境里跑任何命令 |
| `uv run python script.py` | 跑 Python 脚本 |
| `uv pip list` | 看已装包 |
| `uv pip show <pkg>` | 看某包详细信息 |
| `uv lock` | 重新生成 lock（改了 pyproject 后） |
| `uv lock --upgrade` | 升级所有依赖到最新兼容版本 |
| `uv python list` | 看本机所有可用 Python |
| `uv python install 3.12` | 装新版 Python |
| `uv cache clean` | 清全局包缓存（节省磁盘） |
| `uv tree` | 看依赖树 |

---

## 5. 典型场景

### 5.1 克隆代码后如何快速启动

```powershell
git clone <repo>
cd <repo>
uv sync          # 自动读 pyproject.toml + uv.lock，建 .venv，装依赖
uv run python main.py
```

三条命令，新电脑 5 分钟出结果。

### 5.2 升级某个包

```powershell
uv add numpy --upgrade         # 升级到最新兼容版
uv add "numpy>=2.0"            # 升级到指定约束
```

### 5.3 加装 GPU 版 torch（本项目大概率会遇到）

```powershell
# 先卸掉可能已装的 cpu 版
uv remove torch torchvision

# 再装 CUDA 版
uv add torch torchvision --index https://download.pytorch.org/whl/cu121

# 验证
uv run python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```

### 5.4 给不同项目用不同 Python 版本

每个项目目录里 `.python-version` 文件内容不同即可，uv 自动切换。
无需全局修改任何东西。

### 5.5 临时跑个 Python 片段

```powershell
uv run python -c "import torch; print(torch.cuda.is_available())"
```

### 5.6 一次性工具（不装进项目）

```powershell
uvx ruff check .            # 临时跑 ruff，不写进依赖
uvx black workflow/         # 类似 pipx
```

---

## 6. pyproject.toml 长什么样

uv init 后大致这样：

```toml
[project]
name = "scintillator-detector-continuous"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "matplotlib>=3.9.0",
    "numpy>=2.0.0",
    "pandas>=2.2.0",
    "scikit-learn>=1.5.0",
    "torch>=2.5.0",
    "torchvision>=0.20.0",
]

[tool.uv]
# 可选：固定 pytorch 索引
[[tool.uv.index]]
name = "pytorch-cu121"
url = "https://download.pytorch.org/whl/cu121"
explicit = true

[tool.uv.sources]
torch = { index = "pytorch-cu121" }
torchvision = { index = "pytorch-cu121" }
```

这个文件是**项目的单一真相源**。直接编辑它，然后跑 `uv sync`，也能生效。

---

## 7. 常见问题 (FAQ)

### Q1：会不会影响我现有的 3.10、3.14？
**不会**。uv 只在项目目录下创建 `.venv`，系统 Python 完全独立，可以随时删掉 uv 再回到 pip。

### Q2：和 conda 能共存吗？
**能**，但别在同一个项目里混用。conda 环境里可以跑 `uv pip install`，当成快速 pip 用。新项目建议纯 uv。

### Q3：CUDA 版 torch 怎么确认装对了？
```powershell
uv run python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```
期望：`2.5.x+cu121 True`；如果是 `+cpu` 或 `False`，见下面 §5.3。

### Q4：有没有公司代理/防火墙问题？
设环境变量：
```powershell
$env:HTTPS_PROXY = "http://your.proxy:port"
uv sync
```
或写进 `pyproject.toml`：
```toml
[tool.uv]
index-url = "https://your.mirror/simple"
```

### Q5：删除项目环境重来？
```powershell
Remove-Item -Recurse -Force .venv
uv sync
```

### Q6：如何查某命令到底用的是哪个 Python？
```powershell
uv run python -c "import sys; print(sys.executable)"
```

### Q7：本项目已有 `requirements.txt` 怎么办？
```powershell
uv add -r requirements.txt
```
uv 会把里面每一行都加到 `pyproject.toml`。

---

## 8. 和 pip/conda 的对比一眼看清

| 动作 | pip | conda | uv |
|---|---|---|---|
| 装依赖 | `pip install X` | `conda install X` | `uv add X` |
| 建环境 | `python -m venv .venv` | `conda create -n foo` | 自动（uv init 时） |
| 激活 | `.venv\Scripts\activate` | `conda activate foo` | **不需要**（uv run） |
| 锁版本 | 手写 requirements.txt | 手写 environment.yml | **自动** uv.lock |
| 切 Python 版本 | 另装一个 Python | `conda install python=3.11` | `uv python install 3.11` |
| 速度 | 慢 | 慢 | 快 10–100× |
| 复现性 | 差 | 中 | 强 |

---

## 9. 本项目建议的 3 条命令收尾

```powershell
cd D:\Geant4\Projects\Scintillator-Detector-Continuous
uv init --python 3.11
uv add numpy pandas matplotlib scikit-learn torch torchvision --index https://download.pytorch.org/whl/cu121
uv run python workflow\flow2\Histo10_Cubic.py Results --no-heatmap
```

跑通这三条，你就从"两个 Python 打架的混沌"迁到"干净 uv 项目"了。

---

## 10. 学习材料

- **官方文档**：<https://docs.astral.sh/uv/>
- **快速上手（推荐）**：<https://docs.astral.sh/uv/guides/projects/>
- **PyTorch 官方 uv 指引**：<https://docs.astral.sh/uv/guides/integration/pytorch/>
- **视频教程**：YouTube 搜 "uv Astral Python" —— ArjanCodes 有一个 15 分钟的讲得很清楚

---

## 附录：和本笔记配套的另一份笔记

- `Note-2026-04-18-NN-Position-Reconstruction.md` —— 神经网络重建位置实施笔记。
  两份结合：用 uv 管项目环境 → 在环境里跑 NN 训练脚本，互补。
